"""Local web server for OligoVisualizer 2.1: linear and IDT structure views."""

import math
import socket
import threading
import webbrowser
from collections.abc import Mapping

from flask import Flask, jsonify, request

from idt import (
    DEFAULT_CONDITIONS,
    FoldConditions,
    IdtAuthError,
    IdtError,
    IdtSession,
    IdtSessionExpired,
    NoStructureError,
)
from notation import normalize, parse, summary_line
from render import (
    legend_items,
    render_homodimer_svg,
    render_structure_svg,
    render_svg,
    structure_legend_items,
)

app = Flask(__name__, static_folder="static", static_url_path="/static")

STRUCTURE_CACHE_LIMIT = 64
# IDT Analyzer's default mode is DNA, even when the sequence contains U.
STRUCTURE_NUCLEOTIDE_TYPE = "DNA"
_session = IdtSession()
_structure_cache: dict[tuple[str, str, FoldConditions], dict] = {}
_structure_lock = threading.Lock()

_CONDITION_SPECS = (
    {
        "query": "temp",
        "field": "temperature_c",
        "display": "Temperature",
        "minimum": 0,
        "maximum": 100,
        "kind": "float",
        "suffix": " °C",
    },
    {
        "query": "na",
        "field": "sodium_mm",
        "display": "Na⁺",
        "minimum": 0,
        "maximum": 1500,
        "kind": "float",
        "suffix": " mM",
    },
    {
        "query": "mg",
        "field": "magnesium_mm",
        "display": "Mg²⁺",
        "minimum": 0,
        "maximum": 600,
        "kind": "float",
        "suffix": " mM",
    },
    {
        "query": "suboptimality",
        "field": "suboptimality_percent",
        "display": "Suboptimality",
        "minimum": 0,
        "maximum": 100,
        "kind": "float",
        "suffix": "%",
    },
    {
        "query": "max_foldings",
        "field": "max_foldings",
        "display": "Maximum foldings",
        "minimum": 1,
        "maximum": 100,
        "kind": "int",
        "suffix": "",
    },
    {
        "query": "oligo",
        "field": "oligo_conc_um",
        "display": "Oligo concentration",
        "minimum": 0.0001,
        "maximum": 10000,
        "kind": "float",
        "suffix": " µM",
    },
    {
        "query": "dntps",
        "field": "dntps_mm",
        "display": "dNTP concentration",
        "minimum": 0,
        "maximum": 1000,
        "kind": "float",
        "suffix": " mM",
    },
)


def _condition_problem(spec: dict) -> dict:
    whole = "whole " if spec["kind"] == "int" else ""
    return {
        "message": (
            f"{spec['display']} must be a {whole}number between "
            f"{spec['minimum']:g} and {spec['maximum']:g}{spec['suffix']}."
        ),
        "position": 0,
        "token": "",
    }


def _parse_conditions(
    query: Mapping[str, str],
) -> tuple[FoldConditions | None, dict | None]:
    """Validate the folding and self-dimer inputs; missing values use defaults."""
    values: dict[str, float | int] = {}
    for spec in _CONDITION_SPECS:
        raw = query.get(spec["query"])
        if raw is None:
            values[spec["field"]] = getattr(DEFAULT_CONDITIONS, spec["field"])
            continue
        try:
            number = float(raw)
        except (TypeError, ValueError):
            return None, _condition_problem(spec)
        if not math.isfinite(number):
            return None, _condition_problem(spec)
        if number < spec["minimum"] or number > spec["maximum"]:
            return None, _condition_problem(spec)
        if spec["kind"] == "int":
            if not number.is_integer():
                return None, _condition_problem(spec)
            values[spec["field"]] = int(number)
        else:
            values[spec["field"]] = number
    return FoldConditions(**values), None


def _problem_json(problem) -> dict:
    return {
        "message": problem.message,
        "position": problem.position,
        "token": problem.token,
    }


def build_payload(raw: str) -> dict:
    """Turn a raw IDT string into the JSON payload the page consumes."""
    if raw == "":
        return {
            "ok": False,
            "svg": "",
            "summary": "",
            "legend": [],
            "warnings": [],
            "errors": [
                {"message": "no sequence supplied", "position": 0, "token": ""},
            ],
        }

    seq = parse(raw)
    warnings = [_problem_json(problem) for problem in seq.warnings]
    errors = [_problem_json(problem) for problem in seq.errors]

    if not seq.nucleotides:
        if not errors:
            errors = [
                {"message": "no nucleotide found", "position": 0, "token": ""},
            ]
        return {
            "ok": False,
            "svg": "",
            "summary": "",
            "legend": [],
            "warnings": warnings,
            "errors": errors,
        }

    return {
        "ok": True,
        "svg": render_svg(seq),
        "summary": summary_line(seq),
        "legend": legend_items(seq),
        "warnings": warnings,
        "errors": errors,
    }


def _structure_error(
    message: str,
    *,
    signed_in: bool = True,
    login_required: bool = False,
    warnings: list[dict] | None = None,
    errors: list[dict] | None = None,
) -> dict:
    return {
        "ok": False,
        "svg": "",
        "summary": "",
        "detail": "",
        "legend": [],
        "warnings": warnings or [],
        "errors": (
            [{"message": message, "position": 0, "token": ""}]
            if errors is None
            else errors
        ),
        "structure": None,
        "homodimer": None,
        "properties": None,
        "signed_in": signed_in,
        "login_required": login_required,
    }


def build_structure_payload(raw: str, query=None) -> dict:
    """Fold one raw IDT string with the signed-in Analyzer and draw it."""
    if raw == "":
        return _structure_error("no sequence supplied")

    seq = parse(raw)
    parse_warnings = [_problem_json(problem) for problem in seq.warnings]
    parse_errors = [_problem_json(problem) for problem in seq.errors]
    if not seq.nucleotides:
        if not parse_errors:
            parse_errors = [
                {"message": "no nucleotide found", "position": 0, "token": ""},
            ]
        return _structure_error(
            "no nucleotide found",
            warnings=parse_warnings,
            errors=parse_errors,
        )

    conditions, condition_problem = _parse_conditions(query or {})
    if condition_problem is not None:
        return _structure_error(
            condition_problem["message"],
            signed_in=_session.signed_in,
            login_required=False,
            warnings=parse_warnings,
            errors=parse_errors + [condition_problem],
        )
    assert conditions is not None

    requested_type = STRUCTURE_NUCLEOTIDE_TYPE
    if not _session.signed_in:
        return _structure_error(
            "Sign in to IDT to use the Analyzer.",
            signed_in=False,
            login_required=True,
            warnings=parse_warnings,
            errors=parse_errors
            + [
                {
                    "message": "Sign in to IDT to use the Analyzer.",
                    "position": 0,
                    "token": "",
                }
            ],
        )

    key = (normalize(raw), requested_type, conditions)
    with _structure_lock:
        cached = _structure_cache.get(key)
    if cached is not None:
        return cached

    sugar_warning: dict | None = None
    properties = None
    properties_warning: dict | None = None
    homodimer = None
    homodimer_svg = ""
    homodimer_warning: dict | None = None
    try:
        fold = _session.predict(raw, requested_type, conditions)
        sugars = tuple(nt.sugar for nt in seq.nucleotides)
        if len(sugars) != len(fold.bases):
            sugars = ()
            sugar_warning = {
                "message": (
                    "IDT returned a different number of bases than the parsed "
                    "sequence; 2' modification icons are hidden."
                ),
                "position": 0,
                "token": "",
            }
        svg, layout = render_structure_svg(fold.bases, fold.pairs, sugars)

        try:
            properties = _session.analyze_properties(
                raw, requested_type, conditions
            )
        except (IdtAuthError, IdtSessionExpired):
            raise
        except IdtError as exc:
            properties_warning = {
                "message": f"IDT property analysis failed: {exc}",
                "position": exc.position,
                "token": exc.token,
            }
        except Exception as exc:
            properties_warning = {
                "message": f"IDT property analysis failed: {exc}",
                "position": 0,
                "token": "",
            }

        try:
            homodimer = _session.predict_homodimer(
                raw, requested_type, conditions
            )
            dimer_sugars = sugars
            if len(dimer_sugars) != len(homodimer.sequence):
                dimer_sugars = ()
                homodimer_warning = {
                    "message": (
                        "IDT returned a different self-dimer sequence length "
                        "than the parsed sequence; 2' modification icons are "
                        "hidden on the homodimer."
                    ),
                    "position": 0,
                    "token": "",
                }
            homodimer_svg = render_homodimer_svg(
                homodimer.sequence, homodimer.dimer, dimer_sugars
            )
        except (IdtAuthError, IdtSessionExpired):
            raise
        except NoStructureError as exc:
            homodimer = None
            homodimer_warning = {
                "message": str(exc),
                "position": 0,
                "token": "",
            }
        except IdtError as exc:
            homodimer = None
            homodimer_warning = {
                "message": f"IDT self-dimer analysis failed: {exc}",
                "position": exc.position,
                "token": exc.token,
            }
        except Exception as exc:
            homodimer = None
            homodimer_warning = {
                "message": f"IDT self-dimer analysis failed: {exc}",
                "position": 0,
                "token": "",
            }
    except IdtAuthError as exc:
        return _structure_error(
            str(exc),
            signed_in=False,
            login_required=True,
            warnings=parse_warnings,
            errors=parse_errors
            + [{"message": str(exc), "position": exc.position, "token": exc.token}],
        )
    except IdtSessionExpired as exc:
        _session.sign_out()
        return _structure_error(
            str(exc),
            signed_in=False,
            login_required=True,
            warnings=parse_warnings,
            errors=parse_errors
            + [{"message": str(exc), "position": exc.position, "token": exc.token}],
        )
    except NoStructureError as exc:
        return _structure_error(
            str(exc),
            signed_in=True,
            warnings=parse_warnings
            + [{"message": str(exc), "position": 0, "token": ""}],
            errors=parse_errors,
        )
    except IdtError as exc:
        return _structure_error(
            str(exc),
            signed_in=True,
            warnings=parse_warnings,
            errors=parse_errors
            + [
                {
                    "message": str(exc),
                    "position": exc.position,
                    "token": exc.token,
                }
            ],
        )
    except Exception as exc:
        message = f"IDT prediction failed: {exc}"
        return _structure_error(
            message,
            signed_in=True,
            warnings=parse_warnings,
            errors=parse_errors
            + [{"message": message, "position": 0, "token": ""}],
        )

    if layout == "none":
        message = "IDT predicted no base pairs for this sequence."
        return _structure_error(
            message,
            signed_in=True,
            warnings=parse_warnings
            + [{"message": message, "position": 0, "token": ""}],
            errors=parse_errors,
        )

    summary = f"ΔG {fold.delta_g:+.2f} kcal/mol"
    if fold.tm_c is not None:
        summary += f" · Tm {fold.tm_c:.1f} °C"
    summary += (
        f" · structure {fold.structure_name} of {fold.structure_count} "
        "(lowest ΔG)"
    )

    detail = (
        f"IDT OligoAnalyzer · {fold.nucleotide_type} · {len(fold.bases)} nt · "
        f"{len(fold.pairs)} base pairs"
    )
    if fold.delta_h is not None:
        detail += f" · ΔH {fold.delta_h:.2f} kcal/mol"
    if fold.delta_s is not None:
        detail += f" · ΔS {fold.delta_s:.2f} cal·K⁻¹·mol⁻¹"
    detail += (
        f" · {fold.temperature_c:g} °C · {fold.sodium_mm:g} mM Na⁺ · "
        f"{fold.magnesium_mm:g} mM Mg²⁺ · {fold.suboptimality_percent:g}% "
        f"suboptimality · max {fold.max_foldings} foldings"
    )
    if layout == "arc":
        detail += " · complex fold drawn as base-pair arcs"
    elif layout == "tree":
        detail += " · branched fold drawn as connected stem-loops"

    properties_payload = None
    if properties is not None:
        properties_payload = {
            "melt_temp_c": properties.melt_temp_c,
            "molecular_weight_g_mol": properties.molecular_weight_g_mol,
            "extinction_coefficient_l_mol_cm": (
                properties.extinction_coefficient_l_mol_cm
            ),
            "nmol_per_od260": properties.nmol_per_od260,
            "ug_per_od260": properties.ug_per_od260,
        }

    homodimer_payload = None
    if homodimer is not None:
        dimer_summary = (
            f"ΔG {homodimer.delta_g:+.2f} kcal/mol · "
            f"{homodimer.base_pairs} base pairs · most likely homodimer"
        )
        dimer_detail = (
            f"IDT OligoAnalyzer self-dimer · 2 × "
            f"{len(homodimer.sequence)} nt · {homodimer.result_count} "
            f"candidate{'s' if homodimer.result_count != 1 else ''} · "
            f"{homodimer.base_pairs} counted base pairs"
        )
        if homodimer.additional_pairs:
            dimer_detail += f" · {homodimer.additional_pairs} additional"
        if homodimer.max_delta_g is not None:
            dimer_detail += (
                f" · perfect-complement ΔG "
                f"{homodimer.max_delta_g:+.2f} kcal/mol"
            )
        dimer_detail += (
            f" · {conditions.oligo_conc_um:g} µM oligo · "
            f"{conditions.dntps_mm:g} mM dNTPs · "
            f"{conditions.sodium_mm:g} mM Na⁺ · "
            f"{conditions.magnesium_mm:g} mM Mg²⁺"
        )
        homodimer_payload = {
            "ok": True,
            "svg": homodimer_svg,
            "summary": dimer_summary,
            "detail": dimer_detail,
            "delta_g_kcal_mol": homodimer.delta_g,
            "base_pairs": homodimer.base_pairs,
            "additional_pairs": homodimer.additional_pairs,
            "result_count": homodimer.result_count,
            "perfect_delta_g_kcal_mol": homodimer.max_delta_g,
            "sequence": homodimer.sequence,
            "layout": "homodimer",
        }

    legend = structure_legend_items(fold.bases, sugars)
    if homodimer is not None and homodimer.additional_pairs:
        legend.append(
            {
                "kind": "pair",
                "label": "Additional, not counted",
                "style": "additional",
            }
        )

    warnings = list(parse_warnings)
    if fold.delta_g >= 0:
        warnings.append(
            {
                "message": (
                    "Lowest predicted ΔG is not negative; the oligo is not "
                    "predicted to fold stably under these conditions."
                ),
                "position": 0,
                "token": "",
            }
        )
    if sugar_warning is not None:
        warnings.append(sugar_warning)
    if properties_warning is not None:
        warnings.append(properties_warning)
    if homodimer_warning is not None:
        warnings.append(homodimer_warning)

    payload = {
        "ok": True,
        "svg": svg,
        "summary": summary,
        "detail": detail,
        "legend": legend,
        "warnings": warnings,
        "errors": parse_errors,
        "structure": {
            "name": fold.structure_name,
            "delta_g_kcal_mol": fold.delta_g,
            "tm_c": fold.tm_c,
            "delta_h_kcal_mol": fold.delta_h,
            "delta_s_cal_k_mol": fold.delta_s,
            "structure_count": fold.structure_count,
            "nucleotide_type": fold.nucleotide_type,
            "temperature_c": fold.temperature_c,
            "sodium_mm": fold.sodium_mm,
            "magnesium_mm": fold.magnesium_mm,
            "suboptimality_percent": fold.suboptimality_percent,
            "max_foldings": fold.max_foldings,
            "base_pairs": [[i, j] for i, j in fold.pairs],
            "layout": layout,
        },
        "homodimer": homodimer_payload,
        "properties": properties_payload,
        "signed_in": True,
        "login_required": False,
    }

    with _structure_lock:
        _structure_cache[key] = payload
        while len(_structure_cache) > STRUCTURE_CACHE_LIMIT:
            del _structure_cache[next(iter(_structure_cache))]
    return payload


def _no_store(payload: dict, status: int = 200):
    response = jsonify(payload)
    response.status_code = status
    response.headers["Cache-Control"] = "no-store"
    return response


@app.get("/")
def index():
    """Serve the viewer page verbatim; it is static and never templated."""
    return app.send_static_file("index.html")


@app.get("/api/render")
def api_render():
    """Return the drawing payload for one pasted sequence; always HTTP 200."""
    return jsonify(build_payload(request.args.get("seq", "")))


@app.get("/api/idt/status")
def api_idt_status():
    return _no_store({"signed_in": _session.signed_in})


@app.post("/api/idt/login")
def api_idt_login():
    origin = request.headers.get("Origin")
    if origin and origin.rstrip("/") != request.host_url.rstrip("/"):
        return _no_store(
            {"ok": False, "error": "cross-origin request rejected"}, status=403
        )

    body = request.get_json(silent=True)
    if not isinstance(body, dict):
        body = {}
    username = body.get("username")
    password = body.get("password")
    username = "" if username is None else str(username)
    password = "" if password is None else str(password)
    if not username or not password:
        return _no_store(
            {"ok": False, "error": "Enter your IDT username and password."}
        )

    try:
        _session.sign_in(username, password)
    except IdtAuthError as exc:
        return _no_store({"ok": False, "error": str(exc)})
    except Exception:
        return _no_store({"ok": False, "error": "IDT sign-in failed."})
    return _no_store({"ok": True})


@app.post("/api/idt/logout")
def api_idt_logout():
    _session.sign_out()
    return _no_store({"ok": True})


@app.get("/api/structure")
def api_structure():
    """Fold one pasted sequence with IDT; always HTTP 200."""
    return _no_store(
        build_structure_payload(request.args.get("seq", ""), request.args)
    )


def main() -> None:
    """Take the first free port in 8765-8785, announce it, open the browser, serve."""
    port = None
    for candidate in range(8765, 8786):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
                probe.bind(("127.0.0.1", candidate))
            port = candidate
            break
        except OSError:
            continue

    if port is None:
        print("No free port found between 8765 and 8785.")
        return

    url = f"http://127.0.0.1:{port}/"
    print(f"Oligo viewer running at {url}")
    print("Close this window to stop the viewer.")
    webbrowser.open(url)
    app.run(host="127.0.0.1", port=port, debug=False, use_reloader=False)


if __name__ == "__main__":
    main()
