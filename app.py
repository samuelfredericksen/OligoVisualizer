"""Local web server for OligoVisualizer 2.0: linear and IDT hairpin views."""

import socket
import threading
import webbrowser

from flask import Flask, jsonify, request

from idt import (
    IdtAuthError,
    IdtError,
    IdtSession,
    IdtSessionExpired,
    NoStructureError,
)
from notation import normalize, parse, summary_line
from render import (
    legend_items,
    render_structure_svg,
    render_svg,
    structure_legend_items,
)

app = Flask(__name__, static_folder="static", static_url_path="/static")

STRUCTURE_CACHE_LIMIT = 64
# IDT Analyzer's default mode is DNA, even when the sequence contains U.
STRUCTURE_NUCLEOTIDE_TYPE = "DNA"
_session = IdtSession()
_structure_cache: dict[tuple[str, str], dict] = {}
_structure_lock = threading.Lock()


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
        "signed_in": signed_in,
        "login_required": login_required,
    }


def build_structure_payload(raw: str) -> dict:
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

    key = (normalize(raw), requested_type)
    with _structure_lock:
        cached = _structure_cache.get(key)
    if cached is not None:
        return cached

    try:
        fold = _session.predict(raw, requested_type)
        svg, layout = render_structure_svg(fold.bases, fold.pairs)
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

    payload = {
        "ok": True,
        "svg": svg,
        "summary": summary,
        "detail": detail,
        "legend": structure_legend_items(fold.bases),
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
    return _no_store(build_structure_payload(request.args.get("seq", "")))


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
