"""IDT OligoAnalyzer client: sign in, predict a hairpin, fetch its CT map.

Standard library only.  The password and username are never stored; the IDT
cookies and anti-forgery tokens live in this process's memory.
"""

from __future__ import annotations

import base64
import binascii
import html
import json
import re
import socket
import threading
from dataclasses import dataclass
from html.parser import HTMLParser
from http.cookiejar import CookieJar
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import unquote, urlencode, urljoin, urlsplit
from urllib.request import (
    HTTPCookieProcessor,
    HTTPRedirectHandler,
    OpenerDirector,
    Request,
    build_opener,
)
from xml.etree import ElementTree

from notation import (
    MOD_CODE_CHARS,
    Sequence,
    _looks_like_sequence,
    normalize,
)

LOGIN_PAGE = "https://www.idtdna.com/site/account/account?returnurl=/calc/analyzer"
LOGIN_POST = "https://www.idtdna.com/site/account/login/gatekeepernew"
HAIRPIN_PATH = "/calc/analyzer/home/hairpin"
SELFDIMER_PATH = "/calc/analyzer/home/selfdimer"
ANALYZE_PATH = "/api/OligoAnalyzer/Analyze"
DATAOUT_PATH = "/api/scitools/dataout"
PUBLIC_CT_PATH = "/AnalyzerService/analyzerservice.asmx/UNAFoldRetrieveImage"
BASE_URL = "https://www.idtdna.com"

USER_AGENT = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) OligoVisualizer/2"


class IdtError(Exception):
    """An IDT or transport error that can be shown to the user."""

    def __init__(self, message: str, position: int = 0, token: str = "") -> None:
        super().__init__(message)
        self.position = position
        self.token = token


class IdtAuthError(IdtError):
    """The user must sign in, or the supplied credentials were rejected."""


class IdtSessionExpired(IdtError):
    """IDT redirected a signed-in request back to the login page."""


class NoStructureError(IdtError):
    """IDT completed the request but returned no usable base pairs."""


@dataclass(frozen=True)
class UnafoldStructure:
    order: int
    name: str
    delta_g: float
    tm_c: float | None
    delta_h: float | None
    delta_s: float | None
    connectivity_url: str


@dataclass(frozen=True)
class Fold:
    bases: str
    pairs: tuple[tuple[int, int], ...]
    structure_name: str
    delta_g: float
    tm_c: float | None
    delta_h: float | None
    delta_s: float | None
    structure_count: int
    nucleotide_type: str
    temperature_c: float
    sodium_mm: float
    magnesium_mm: float
    suboptimality_percent: float
    max_foldings: int


@dataclass(frozen=True)
class OligoProperties:
    melt_temp_c: float
    molecular_weight_g_mol: float
    extinction_coefficient_l_mol_cm: float
    nmol_per_od260: float
    ug_per_od260: float


@dataclass(frozen=True)
class Homodimer:
    sequence: str
    delta_g: float
    base_pairs: int
    additional_pairs: int
    max_delta_g: float | None
    dimer: str
    result_count: int


@dataclass(frozen=True)
class FoldConditions:
    temperature_c: float = 25
    sodium_mm: float = 50
    magnesium_mm: float = 3
    suboptimality_percent: float = 50
    max_foldings: int = 20
    # The Analyzer's self-dimer endpoint also uses these two values.
    oligo_conc_um: float = 0.2
    dntps_mm: float = 0.8


DEFAULT_CONDITIONS = FoldConditions()

DEFAULTS = {
    "SequenceType": "Linear",
    "StartPos": 0,
    "StopPos": 0,
}


# --- pure helpers ---------------------------------------------------------


def _protect_code_chains(text: str) -> tuple[str, list[str]]:
    """Canonicalise slash-delimited code chains for the Analyzer payload.

    IDT accepts two adjacent codes as ``/i2FG//i2MOErG/``.  The local parser
    also accepts the shorthand ``/i2FG/i2MOErG/``, where one slash serves as
    both closing and opening delimiter; the Analyzer does not.  Rebuild a
    chain with two slashes between codes, and keep the whole chain out of the
    base-uppercasing pass so an ``i`` is never changed to inosine.
    """
    protected: list[str] = []
    pieces: list[str] = []
    i = 0
    while i < len(text):
        if text[i] != "/":
            pieces.append(text[i])
            i += 1
            continue
        closing = text.find("/", i + 1)
        if closing == -1:
            pieces.append(text[i:])
            break
        codes = [text[i + 1:closing]]
        end = closing + 1
        while True:
            next_slash = text.find("/", end)
            if next_slash == -1:
                break
            following = text[end:next_slash]
            if (
                not following
                or not all(ch in MOD_CODE_CHARS for ch in following)
                or _looks_like_sequence(following)
            ):
                break
            codes.append(following)
            end = next_slash + 1
        protected.append("/" + "//".join(codes) + "/")
        pieces.append(f"\x00{len(protected) - 1}\x00")
        i = end
    return "".join(pieces), protected


def clean_sequence(sequence: str) -> str:
    """Clean IDT notation without stripping any modification code."""
    text, protected = _protect_code_chains(normalize(sequence))
    uppercase = set("atcgiykswhbvdn")
    text = "".join(ch.upper() if ch in uppercase else ch for ch in text)
    for index, token in enumerate(protected):
        text = text.replace(f"\x00{index}\x00", token)
    return text


def nucleotide_type(seq: Sequence) -> str:
    """Return 'RNA' for ribose/uridine chains, otherwise 'DNA'."""
    if any(nt.sugar == "r" or nt.display == "U" for nt in seq.nucleotides):
        return "RNA"
    return "DNA"


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _number(value: Any) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value.strip())
        except ValueError:
            return None
    return None


def _int_or(value: Any, default: int) -> int:
    number = _number(value)
    if number is None:
        return default
    return int(round(number))


def _float_or(value: Any, default: float) -> float:
    number = _number(value)
    return default if number is None else number


def _struct_error(item: Any) -> IdtError:
    if isinstance(item, dict):
        message = (
            item.get("ErrorText")
            or item.get("ErrorMessage")
            or item.get("Message")
            or item.get("ErrorSymbol")
            or "IDT rejected the sequence."
        )
        position = _int_or(item.get("Position"), 0)
        symbol = str(item.get("ErrorSymbol") or "")
        return IdtError(_plain_text(str(message)), position, symbol)
    return IdtError(_plain_text(str(item)))


def select_structure(
    data: dict,
    requested_type: str,
    conditions: FoldConditions = DEFAULT_CONDITIONS,
) -> tuple[UnafoldStructure, int, dict]:
    """Validate the hairpin response and choose the lowest-G structure."""
    errors = _as_list(data.get("Errors"))
    if data.get("HasErrors"):
        if errors:
            error = _struct_error(errors[0])
            if "no structure" in str(error).lower():
                raise NoStructureError(str(error))
            raise error
        raise IdtError("IDT rejected the sequence.")

    validation = data.get("ValidationError")
    if validation:
        text = _plain_text(str(validation))
        if "no structure" in text.lower():
            raise NoStructureError(text)
        raise IdtError(text)

    model_errors = _as_list(data.get("ModelErrors"))
    if data.get("HasModelErrors") and model_errors:
        error = _struct_error(model_errors[0])
        if "no structure" in str(error).lower():
            raise NoStructureError(str(error))
        raise error

    output = data.get("OutputObj")
    if not isinstance(output, dict):
        raise NoStructureError("IDT found no hairpin structures for this sequence.")
    items = _as_list(output.get("Structures"))
    if not items:
        raise NoStructureError("IDT found no hairpin structures for this sequence.")

    structures: list[UnafoldStructure] = []
    for order, item in enumerate(items):
        if not isinstance(item, dict):
            continue
        delta_g = _number(item.get("GNode"))
        if delta_g is None:
            continue
        structures.append(
            UnafoldStructure(
                order=order,
                name=str(item.get("StructureName") or (order + 1)),
                delta_g=delta_g,
                tm_c=_number(item.get("TMNode")),
                delta_h=_number(item.get("HNode")),
                delta_s=_number(item.get("SNode")),
                connectivity_url=str(item.get("ConnectivityURL") or ""),
            )
        )
    if not structures:
        raise NoStructureError("IDT found no hairpin structures for this sequence.")

    chosen = min(structures, key=lambda item: (item.delta_g, item.order))

    def echoed(key: str, default):
        value = data.get(key)
        if value is None:
            value = output.get(key)
        return default if value is None else value

    meta = {
        "nucleotide_type": str(echoed("NucleotideType", requested_type)),
        "temperature_c": _float_or(
            echoed("Temperature", conditions.temperature_c),
            conditions.temperature_c,
        ),
        "sodium_mm": _float_or(
            echoed("SodiumConc", conditions.sodium_mm),
            conditions.sodium_mm,
        ),
        "magnesium_mm": _float_or(
            echoed("MagConc", conditions.magnesium_mm),
            conditions.magnesium_mm,
        ),
        "suboptimality_percent": _float_or(
            echoed("PercentSuboptimal", conditions.suboptimality_percent),
            conditions.suboptimality_percent,
        ),
        "max_foldings": _int_or(
            echoed("MaxFoldings", conditions.max_foldings),
            conditions.max_foldings,
        ),
    }
    return chosen, len(structures), meta


def _dimer_bar_counts(dimer: str) -> tuple[int, int]:
    """Return (counted, additional) base-pair counts from IDT's text map."""
    lines = [line for line in dimer.splitlines() if line.strip()]
    if len(lines) < 3:
        return 0, 0
    bars = lines[1]
    return bars.count("|"), bars.count(":")


def select_homodimer(data: dict) -> Homodimer:
    """Choose the self-dimer with the lowest predicted Gibbs free energy."""
    errors = _as_list(data.get("Errors"))
    if data.get("HasErrors"):
        if errors:
            error = _struct_error(errors[0])
            if "no dimer" in str(error).lower():
                raise NoStructureError(str(error))
            raise error
        raise IdtError("IDT rejected the sequence.")

    validation = data.get("ValidationError")
    if validation:
        text = _plain_text(str(validation))
        if "no dimer" in text.lower():
            raise NoStructureError(text)
        raise IdtError(text)

    model_errors = _as_list(data.get("ModelErrors"))
    if data.get("HasModelErrors") and model_errors:
        error = _struct_error(model_errors[0])
        if "no dimer" in str(error).lower():
            raise NoStructureError(str(error))
        raise error

    items = _as_list(data.get("Results"))
    matches: list[Homodimer] = []
    for item in items:
        if not isinstance(item, dict):
            continue
        delta_g = _number(item.get("DeltaG"))
        dimer = str(item.get("Dimer") or "")
        if delta_g is None or not dimer.strip():
            continue
        counted, additional = _dimer_bar_counts(dimer)
        base_pairs = _int_or(item.get("BasePairs"), counted)
        matches.append(
            Homodimer(
                sequence=str(
                    data.get("PrimarySequence")
                    or data.get("Sequence")
                    or ""
                ),
                delta_g=delta_g,
                base_pairs=base_pairs,
                additional_pairs=additional,
                max_delta_g=_number(data.get("MaxDeltaG")),
                dimer=dimer,
                result_count=0,
            )
        )
    if not matches:
        raise NoStructureError(
            "IDT found no self-dimer for this sequence."
        )

    chosen = min(enumerate(matches), key=lambda item: (item[1].delta_g, item[0]))
    result = chosen[1]
    return Homodimer(
        sequence=result.sequence,
        delta_g=result.delta_g,
        base_pairs=result.base_pairs,
        additional_pairs=result.additional_pairs,
        max_delta_g=result.max_delta_g,
        dimer=result.dimer,
        result_count=len(items),
    )


def select_properties(data: dict) -> OligoProperties:
    """Validate the Analyzer property response and return its five values."""
    errors = _as_list(data.get("Errors"))
    if data.get("HasErrors"):
        if errors:
            raise _struct_error(errors[0])
        raise IdtError("IDT rejected the sequence.")

    validation = data.get("ValidationError")
    if validation:
        raise IdtError(_plain_text(str(validation)))

    model_errors = _as_list(data.get("ModelErrors"))
    if data.get("HasModelErrors") and model_errors:
        raise _struct_error(model_errors[0])

    values = {
        "melt_temp_c": _number(data.get("MeltTemp")),
        "molecular_weight_g_mol": _number(data.get("MolecularWeight")),
        "extinction_coefficient_l_mol_cm": _number(data.get("ExtCoefficient")),
        "nmol_per_od260": _number(data.get("NmoleOD")),
        "ug_per_od260": _number(data.get("UgOD")),
    }
    if any(value is None for value in values.values()):
        raise IdtError("IDT returned incomplete Analyzer properties.")
    return OligoProperties(
        melt_temp_c=values["melt_temp_c"],
        molecular_weight_g_mol=values["molecular_weight_g_mol"],
        extinction_coefficient_l_mol_cm=values[
            "extinction_coefficient_l_mol_cm"
        ],
        nmol_per_od260=values["nmol_per_od260"],
        ug_per_od260=values["ug_per_od260"],
    )


def _checked_ct(text: str) -> str:
    """Return ``text`` only if it is a valid UNAFold CT map."""
    parse_ct(text)
    return text


def _decode_ct_text(text: str) -> str:
    """Accept CT text, a data: URL, URL-encoding, or base64-wrapped CT."""
    candidates = [text]
    stripped = text.strip()
    if "base64," in stripped:
        candidates.append(stripped.split("base64,", 1)[1])
    if "%" in text:
        candidates.append(unquote(text))
    if "<" in text or "&lt;" in text or "&#" in text:
        unwrapped = re.sub(r"<br\s*/?>", "\n", html.unescape(text), flags=re.I)
        unwrapped = re.sub(r"<[^>]*>", "", unwrapped)
        candidates.extend([unwrapped, unwrapped.strip()])
    for candidate in candidates:
        try:
            return _checked_ct(candidate)
        except IdtError:
            pass

    for candidate in candidates:
        compact = re.sub(r"\s+", "", candidate)
        if len(compact) < 16 or not re.fullmatch(r"[A-Za-z0-9+/=_-]+", compact):
            continue
        padded = compact + "=" * (-len(compact) % 4)
        for decoder in (base64.b64decode, base64.urlsafe_b64decode):
            try:
                decoded = decoder(padded).decode("utf-8-sig")
                return _checked_ct(decoded)
            except (binascii.Error, ValueError, UnicodeDecodeError, IdtError):
                continue
    raise IdtError("IDT returned an unreadable base-pair map.")


def parse_ct(text: str) -> tuple[str, tuple[tuple[int, int], ...]]:
    """Parse UNAFold CT text into bases and 0-based pairs."""
    lines = [line.strip() for line in text.lstrip("\ufeff").splitlines() if line.strip()]
    if not lines:
        raise IdtError("IDT returned an empty base-pair map.")
    header = lines[0].split()
    try:
        declared = int(header[0])
    except (IndexError, ValueError):
        raise IdtError(
            "IDT returned an unreadable base-pair map (bad header)."
        ) from None
    if declared < 1:
        raise IdtError("IDT returned an incomplete base-pair map (empty sequence).")
    if len(lines) - 1 != declared:
        raise IdtError(
            "IDT returned an incomplete base-pair map "
            f"(expected {declared} lines, got {len(lines) - 1})."
        )

    bases: list[str] = []
    partners: dict[int, int] = {}
    for expected, line in enumerate(lines[1:]):
        line_number = expected + 2
        fields = line.split()
        if len(fields) < 5:
            raise IdtError(
                "IDT returned an unreadable base-pair map "
                f"(line {line_number}: too few columns)."
            )
        try:
            raw_index = int(fields[0])
            raw_partner = int(fields[4])
        except ValueError:
            raise IdtError(
                "IDT returned an unreadable base-pair map "
                f"(line {line_number}: non-numeric index)."
            ) from None
        base = fields[1]
        if raw_index != expected + 1 or len(base) != 1:
            raise IdtError(
                "IDT returned an unreadable base-pair map "
                f"(line {line_number}: bad base index)."
            )
        index = raw_index - 1
        partner = raw_partner - 1 if raw_partner > 0 else -1
        if not -1 <= partner < declared or partner == index:
            raise IdtError(
                "IDT returned an unreadable base-pair map "
                f"(line {line_number}: partner out of range)."
            )
        bases.append(base)
        partners[index] = partner

    pairs: list[tuple[int, int]] = []
    for index, partner in partners.items():
        if partner < 0:
            continue
        if partners.get(partner) != index:
            raise IdtError(
                "IDT returned an unreadable base-pair map "
                f"(base {index + 1}: non-reciprocal pair)."
            )
        if partner > index:
            pairs.append((index, partner))
    return "".join(bases), tuple(sorted(pairs))


# --- transport helpers ----------------------------------------------------


class _SafeRedirectHandler(HTTPRedirectHandler):
    """Follow login-page refreshes but never follow an expired-session redirect."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        target = urljoin(req.full_url, newurl)
        lowered = target.lower()
        if "/500.aspx" in lowered or "/404.aspx" in lowered:
            return None
        if "/site/account" in lowered and "/site/account" not in req.full_url.lower():
            return None
        return super().redirect_request(req, fp, code, msg, headers, newurl)


class _InputParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.inputs: list[dict[str, str]] = []

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag.lower() != "input":
            return
        self.inputs.append({key.lower(): (value or "") for key, value in attrs})


def _hidden_input(page: str, name: str) -> str:
    parser = _InputParser()
    try:
        parser.feed(page)
    except Exception:
        return ""
    for attrs in parser.inputs:
        if attrs.get("name") == name:
            return attrs.get("value", "")
    return ""


def _plain_text(value: str) -> str:
    text = re.sub(r"<[^>]*>", " ", value)
    return " ".join(html.unescape(text).split())


class IdtSession:
    """One signed-in Analyzer session, serialised by a re-entrant lock."""

    def __init__(self) -> None:
        self._lock = threading.RLock()
        self._jar: CookieJar | None = None
        self._opener: OpenerDirector | None = None
        self._token = ""
        self._stateless_token = ""
        self._access_token = ""
        self._access_account = ""
        self._signed_in = False
        self._base = BASE_URL

    @property
    def signed_in(self) -> bool:
        with self._lock:
            return self._signed_in

    def _clear_locked(self) -> None:
        self._jar = None
        self._opener = None
        self._token = ""
        self._stateless_token = ""
        self._access_token = ""
        self._access_account = ""
        self._signed_in = False
        self._base = BASE_URL

    def sign_out(self) -> None:
        with self._lock:
            self._clear_locked()

    def sign_in(self, username: str, password: str) -> None:
        """Log in to IDT, then read the Analyzer's anti-forgery tokens."""
        if not username or not password:
            raise IdtAuthError("Enter your IDT username and password.")

        with self._lock:
            self._clear_locked()
            try:
                self._jar = CookieJar()
                self._opener = build_opener(
                    HTTPCookieProcessor(self._jar), _SafeRedirectHandler()
                )
                self._get(LOGIN_PAGE)
                body = self._post_form(
                    LOGIN_POST,
                    {
                        "UserName": username,
                        "Password": password,
                        "RememberMe": "false",
                        "ReturnURL": "/calc/analyzer",
                    },
                    headers={
                        "Accept": "application/json",
                        "X-Requested-With": "XMLHttpRequest",
                        "Referer": LOGIN_PAGE,
                        "Origin": BASE_URL,
                    },
                )
                try:
                    data = json.loads(body.decode("utf-8-sig"))
                except (UnicodeDecodeError, json.JSONDecodeError):
                    raise IdtAuthError(
                        "IDT sign-in returned an unexpected response."
                    ) from None
                if not isinstance(data, dict):
                    raise IdtAuthError("IDT sign-in returned an unexpected response.")
                if not data.get("userAuthenticated"):
                    message = _plain_text(str(data.get("JsError") or ""))
                    raise IdtAuthError(message or "Invalid IDT credentials.")
                if data.get("accessToken") is not None:
                    self._access_token = str(data["accessToken"])
                if data.get("accessAccount") is not None:
                    self._access_account = str(data["accessAccount"])

                app_url = urljoin(
                    self._base + "/", str(data.get("url") or "/calc/analyzer")
                )
                page = self._get(app_url).decode("utf-8-sig", errors="replace")
                token = _hidden_input(page, "__RequestVerificationToken")
                stateless = _hidden_input(
                    page, "__RequestVerificationStatelessToken"
                )
                if not token:
                    raise IdtAuthError(
                        "IDT sign-in did not return the Analyzer tokens."
                    )
                parts = urlsplit(app_url)
                if parts.scheme and parts.netloc:
                    self._base = f"{parts.scheme}://{parts.netloc}"
                self._token = token
                self._stateless_token = stateless
                self._signed_in = True
            except IdtAuthError:
                self._clear_locked()
                raise
            except IdtSessionExpired:
                self._clear_locked()
                raise IdtAuthError(
                    "IDT sign-in did not reach the Analyzer."
                ) from None
            except IdtError as exc:
                self._clear_locked()
                raise IdtAuthError(str(exc)) from None
            except Exception:
                self._clear_locked()
                raise IdtAuthError("IDT sign-in failed.") from None

    def _request(self, req: Request) -> tuple[bytes, str, str]:
        if self._opener is None:
            raise IdtError("Sign in to IDT first.")
        try:
            with self._opener.open(req, timeout=60) as response:
                body = response.read()
                final = response.geturl()
                ctype = response.headers.get_content_type()
        except HTTPError as exc:
            location = exc.headers.get("Location", "") if exc.headers else ""
            target = urljoin(req.full_url, location).lower() if location else ""
            if "/500.aspx" in target:
                raise IdtError("IDT returned HTTP 500.") from None
            if "/404.aspx" in target:
                raise IdtError("IDT returned HTTP 404.") from None
            if "/site/account" in target and "/site/account" not in req.full_url.lower():
                raise IdtSessionExpired(
                    "IDT session expired; sign in again."
                ) from None
            if exc.code in (401, 403):
                raise IdtSessionExpired(
                    "IDT session expired; sign in again."
                ) from None
            raise IdtError(f"IDT returned HTTP {exc.code}.") from None
        except (TimeoutError, socket.timeout):
            raise IdtError("IDT did not respond in time.") from None
        except URLError as exc:
            reason = str(getattr(exc, "reason", "")).lower()
            if isinstance(getattr(exc, "reason", None), (TimeoutError, socket.timeout)):
                raise IdtError("IDT did not respond in time.") from None
            if "timed out" in reason:
                raise IdtError("IDT did not respond in time.") from None
            raise IdtError(
                "Could not reach IDT; check the internet connection."
            ) from None
        except OSError:
            raise IdtError(
                "Could not reach IDT; check the internet connection."
            ) from None

        if "/site/account" in final.lower() and "/site/account" not in req.full_url.lower():
            raise IdtSessionExpired("IDT session expired; sign in again.")
        return body, final, ctype

    def _get(self, url: str, accept: str = "text/html") -> bytes:
        req = Request(
            url,
            headers={"User-Agent": USER_AGENT, "Accept": accept},
            method="GET",
        )
        body, _final, _ctype = self._request(req)
        return body

    def _post_form(
        self,
        url: str,
        fields: dict[str, str],
        headers: dict[str, str] | None = None,
    ) -> bytes:
        request_headers = {
            "User-Agent": USER_AGENT,
            "Content-Type": "application/x-www-form-urlencoded; charset=utf-8",
        }
        if headers:
            request_headers.update(headers)
        req = Request(
            url,
            data=urlencode(fields).encode("utf-8"),
            headers=request_headers,
            method="POST",
        )
        body, _final, _ctype = self._request(req)
        return body

    def _auth_headers(self) -> dict[str, str]:
        headers: dict[str, str] = {}
        if self._access_token:
            headers["Authorization"] = f"Bearer {self._access_token}"
            if self._access_account:
                headers["AccessAccount"] = self._access_account
        return headers

    def _post_json(
        self,
        url: str,
        payload: dict,
        headers: dict[str, str] | None = None,
    ) -> dict:
        request_headers = {
            "User-Agent": USER_AGENT,
            "Content-Type": "application/json; charset=utf-8",
            "Accept": "application/json",
        }
        request_headers.update(self._auth_headers())
        if headers:
            request_headers.update(headers)
        req = Request(
            url,
            data=json.dumps(payload, separators=(",", ":")).encode("utf-8"),
            headers=request_headers,
            method="POST",
        )
        body, final, ctype = self._request(req)
        text = body.decode("utf-8-sig", errors="replace")
        if (
            ctype == "text/html"
            or text.lstrip().startswith("<")
            or "/site/account" in final.lower()
        ):
            raise IdtSessionExpired("IDT session expired; sign in again.")
        if ctype != "application/json" and not ctype.endswith("+json"):
            raise IdtError("IDT returned an unexpected response.")
        try:
            data = json.loads(text)
        except json.JSONDecodeError:
            raise IdtError("IDT returned an unexpected response.") from None
        if not isinstance(data, dict):
            raise IdtError("IDT returned an unexpected response.")
        return data

    def _post_raw(
        self,
        url: str,
        body: bytes,
        headers: dict[str, str] | None = None,
    ) -> bytes:
        request_headers = {
            "User-Agent": USER_AGENT,
            "Content-Type": "application/json; charset=utf-8",
            "Accept": "text/plain, application/json, */*",
        }
        request_headers.update(self._auth_headers())
        if headers:
            request_headers.update(headers)
        req = Request(url, data=body, headers=request_headers, method="POST")
        data, _final, _ctype = self._request(req)
        return data

    # --- prediction -------------------------------------------------------

    def predict(
        self,
        sequence: str,
        nucleotide_type: str,
        conditions: FoldConditions = DEFAULT_CONDITIONS,
    ) -> Fold:
        """Run the Analyzer hairpin prediction and fetch the chosen CT map."""
        with self._lock:
            if not self._signed_in:
                raise IdtAuthError("Sign in to IDT first.")
            cleaned = clean_sequence(sequence)
            payload = {
                "settings": {
                    "Sequence": cleaned,
                    "NaConc": conditions.sodium_mm,
                    "MgConc": conditions.magnesium_mm,
                    "NucleotideType": nucleotide_type,
                },
                "hairpinSettings": {
                    "SequenceType": DEFAULTS["SequenceType"],
                    "MaxFoldings": conditions.max_foldings,
                    "StartPos": DEFAULTS["StartPos"],
                    "StopPos": DEFAULTS["StopPos"],
                    "Temp": conditions.temperature_c,
                    "Suboptimality": conditions.suboptimality_percent,
                },
            }
            data = self._post_json(
                self._base + HAIRPIN_PATH,
                payload,
                headers={
                    "__RequestVerificationToken": self._token,
                    "Referer": self._base + "/calc/analyzer",
                    "Origin": self._base,
                    "X-Requested-With": "XMLHttpRequest",
                },
            )
            chosen, count, meta = select_structure(data, nucleotide_type, conditions)
            if not chosen.connectivity_url:
                raise IdtError(
                    "IDT returned a structure without a base-pair map."
                )
            ct_text = self._retrieve_ct(chosen.connectivity_url)
            bases, pairs = parse_ct(ct_text)
            if not pairs:
                raise NoStructureError(
                    "IDT predicted no base pairs for this sequence."
                )
            return Fold(
                bases=bases,
                pairs=pairs,
                structure_name=chosen.name,
                delta_g=chosen.delta_g,
                tm_c=chosen.tm_c,
                delta_h=chosen.delta_h,
                delta_s=chosen.delta_s,
                structure_count=count,
                nucleotide_type=meta["nucleotide_type"],
                temperature_c=meta["temperature_c"],
                sodium_mm=meta["sodium_mm"],
                magnesium_mm=meta["magnesium_mm"],
                suboptimality_percent=meta["suboptimality_percent"],
                max_foldings=meta["max_foldings"],
            )

    def predict_homodimer(
        self,
        sequence: str,
        nucleotide_type: str,
        conditions: FoldConditions = DEFAULT_CONDITIONS,
    ) -> Homodimer:
        """Run the Analyzer self-dimer prediction and choose the best match."""
        with self._lock:
            if not self._signed_in:
                raise IdtAuthError("Sign in to IDT first.")
            cleaned = clean_sequence(sequence)
            payload = {
                "settings": {
                    "Sequence": cleaned,
                    "NaConc": conditions.sodium_mm,
                    "MgConc": conditions.magnesium_mm,
                    "DNTPsConc": conditions.dntps_mm,
                    "OligoConc": conditions.oligo_conc_um,
                    "NucleotideType": nucleotide_type,
                },
            }
            data = self._post_json(
                self._base + SELFDIMER_PATH,
                payload,
                headers={
                    "Referer": self._base + "/calc/analyzer",
                    "Origin": self._base,
                    "X-Requested-With": "XMLHttpRequest",
                },
            )
            return select_homodimer(data)

    def analyze_properties(
        self,
        sequence: str,
        nucleotide_type: str,
        conditions: FoldConditions = DEFAULT_CONDITIONS,
    ) -> OligoProperties:
        """Run the Analyzer's basic property calculation for one sequence."""
        with self._lock:
            if not self._signed_in:
                raise IdtAuthError("Sign in to IDT first.")
            payload = {
                "Sequence": clean_sequence(sequence),
                "NaConc": conditions.sodium_mm,
                "MgConc": conditions.magnesium_mm,
                "DNTPsConc": conditions.dntps_mm,
                "OligoConc": conditions.oligo_conc_um,
                "NucleotideType": nucleotide_type,
            }
            data = self._post_json(self._base + ANALYZE_PATH, payload)
            return select_properties(data)

    def _retrieve_ct(self, location: str) -> str:
        if not location:
            raise IdtError("IDT returned a structure without a base-pair map.")

        if self._stateless_token:
            try:
                body = self._post_raw(
                    self._base + DATAOUT_PATH,
                    json.dumps(
                        {"type": "hairpin", "data": location},
                        separators=(",", ":"),
                    ).encode("utf-8"),
                    headers={
                        "__RequestVerificationToken": self._stateless_token,
                        "Referer": self._base + "/calc/analyzer",
                        "Origin": self._base,
                        "X-Requested-With": "XMLHttpRequest",
                    },
                )
                text = body.decode("utf-8-sig", errors="replace")
                stripped = text.strip()
                if not stripped:
                    raise IdtError("IDT returned an unreadable base-pair map.")
                if stripped.startswith("<"):
                    try:
                        return _decode_ct_text(stripped)
                    except IdtError:
                        raise IdtSessionExpired(
                            "IDT session expired; sign in again."
                        ) from None
                if stripped.startswith(("{", '"', "[")):
                    try:
                        parsed = json.loads(stripped)
                    except json.JSONDecodeError:
                        raise IdtError(
                            "IDT returned an unreadable base-pair map."
                        ) from None
                    if isinstance(parsed, str):
                        return _decode_ct_text(parsed)
                    if isinstance(parsed, list):
                        parts = [item for item in parsed if isinstance(item, str)]
                        if parts:
                            return _decode_ct_text("\n".join(parts))
                        raise IdtError(
                            "IDT returned an unreadable base-pair map."
                        )
                    if isinstance(parsed, dict):
                        for key in (
                            "data",
                            "Data",
                            "content",
                            "Content",
                            "result",
                            "Result",
                            "d",
                            "value",
                            "ct",
                            "connectivityData",
                        ):
                            value = parsed.get(key)
                            if isinstance(value, str) and value.strip():
                                return _decode_ct_text(value)
                            if isinstance(value, list):
                                parts = [
                                    item for item in value if isinstance(item, str)
                                ]
                                if parts:
                                    return _decode_ct_text("\n".join(parts))
                        for value in parsed.values():
                            if isinstance(value, str) and value.strip():
                                try:
                                    return _decode_ct_text(value)
                                except IdtError:
                                    continue
                        for key in ("error", "Error", "message", "Message"):
                            value = parsed.get(key)
                            if isinstance(value, str) and value.strip():
                                raise IdtError(_plain_text(value))
                        raise IdtError(
                            "IDT returned an unreadable base-pair map."
                        )
                    raise IdtError("IDT returned an unreadable base-pair map.")
                return _decode_ct_text(text)
            except IdtSessionExpired:
                raise
            except IdtError:
                pass

        path = urlsplit(location).path if "://" in location else location
        path = path.split("?", 1)[0].lstrip("/")
        match = re.search(
            r"(?:[^/\\]+[/\\])?[^/\\]+\.ct$", path, re.IGNORECASE
        )
        if not match:
            raise IdtError(
                "IDT returned a structure but its base-pair map could not be retrieved."
            )
        try:
            body = self._post_form(
                self._base + PUBLIC_CT_PATH,
                {"location": match.group(0)},
                headers={"Accept": "text/xml, application/xml"},
            )
            root = ElementTree.fromstring(body.decode("utf-8-sig", errors="replace"))
        except IdtSessionExpired:
            raise
        except Exception:
            raise IdtError(
                "IDT returned a structure but its base-pair map could not be retrieved."
            ) from None

        encoded = ""
        for element in root.iter():
            if element.tag.rsplit("}", 1)[-1] == "base64Binary" and element.text:
                encoded = element.text
                break
        if not encoded:
            raise IdtError(
                "IDT returned a structure but its base-pair map could not be retrieved."
            )
        try:
            raw = base64.b64decode(re.sub(r"\s+", "", encoded), validate=True)
        except (binascii.Error, ValueError):
            raise IdtError(
                "IDT returned a structure but its base-pair map could not be retrieved."
            ) from None
        return _decode_ct_text(raw.decode("utf-8-sig", errors="replace"))
