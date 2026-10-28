"""IDT oligo notation -> Sequence model.

Standard library only.  Every error and warning string here is part of the
frozen contract with render.py, app.py and the tests.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Problem:
    """One issue found while parsing, located by its offset in the input."""

    message: str            # plain English, names the token
    position: int           # 0-based offset into the normalized string
    token: str              # exact offending text, "" if not applicable


@dataclass(frozen=True)
class Nucleotide:
    """One position in the chain: base, sugar, and the linkage that follows it."""

    position: int                 # 0-based, 5' end first
    display: str                  # letter drawn on the tile: "A".."U","I","N","R",…,"?" for a mix
    alternatives: tuple[str, ...]  # ("A",) normally; ("A","G") for R; () for a custom mix
    sugar: str                    # "d" | "r" | "m" | "LNA" | "F" | "MOE"
    ps_after: bool
    raw: str                      # exact substring consumed, e.g. "mA*"


@dataclass(frozen=True)
class Modification:
    """A catalog code as typed, plus the text between its slashes."""

    code: str    # text between the slashes, e.g. "5Biosg"
    raw: str     # exact text as typed, e.g. "/5Biosg/"


@dataclass(frozen=True)
class CustomMix:
    """A custom-mix position: optional name and A/C/G/T percentages summing to 100."""

    name: str                              # "" when unnamed
    percentages: tuple[int, int, int, int]  # A, C, G, T


@dataclass
class Sequence:
    """A parsed oligo, plus every problem noticed while parsing it."""

    nucleotides: list[Nucleotide]
    five_prime: list[Modification]
    three_prime: list[Modification]
    internal: dict[int, list[Modification]]  # key = index of the nucleotide the mod follows
    mixes: dict[int, CustomMix]              # key = nucleotide index; that nucleotide's
                                             # display is "?" and alternatives is ()
    warnings: list[Problem]
    errors: list[Problem]


# IUB mixed-base codes: uppercase only.  A, C, G, T, U and I are plain bases.
IUB_ALTERNATIVES = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "C", "G", "T"),
}

PLAIN_BASES = ("A", "C", "G", "T", "U", "I")
MOD_CODE_CHARS = set(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"
)


def _sugar_base_codes() -> dict[str, tuple[str, str]]:
    codes: dict[str, tuple[str, str]] = {}
    for base in "ACGUT":
        codes[f"i2F{base}"] = (base, "F")
    for base in "ACGUT":
        codes[f"i2MOEr{base}"] = (base, "MOE")
        codes[f"52MOEr{base}"] = (base, "MOE")
        codes[f"32MOEr{base}"] = (base, "MOE")
    return codes


SUGAR_BASE_CODES = _sugar_base_codes()

# Characters that can appear in a run of nucleotides (bases, IUB codes, the
# r/m prefixes and the LNA prefix).  Used only to decide whether text sitting
# between two slashes is a sequence rather than an adjacent catalog code.
_NUCLEOTIDE_RUN_CHARS = set("ACGTUIacgtuiRYSWKMBDHVNrm+")


def _looks_like_sequence(text: str) -> bool:
    return all(ch in _NUCLEOTIDE_RUN_CHARS for ch in text)


def normalize(raw: str) -> str:
    """Strip whitespace and the 5'/3' quoting IDT users paste along with the chain."""
    s = "".join(raw.split())
    if s.startswith("5'"):
        s = s[2:]
        if s.startswith("-"):
            s = s[1:]
    if s.endswith("3'"):
        s = s[:-2]
        if s.endswith("-"):
            s = s[:-1]
    return s


def parse(raw: str) -> Sequence:
    """Turn an IDT string into a Sequence, collecting problems instead of raising."""
    s = normalize(raw)
    n = len(s)

    nucleotides: list[Nucleotide] = []
    five_prime: list[Modification] = []
    three_prime: list[Modification] = []
    internal: dict[int, list[Modification]] = {}
    mixes: dict[int, CustomMix] = {}
    warnings: list[Problem] = []
    errors: list[Problem] = []
    pending: list[Modification] = []

    i = 0

    def emit(display: str, alternatives: tuple[str, ...], sugar: str, start: int, end: int) -> None:
        nonlocal i
        ps_after = False
        if end < n and s[end] == "*":
            ps_after = True
            end += 1
        index = len(nucleotides)
        nucleotides.append(
            Nucleotide(index, display, alternatives, sugar, ps_after, s[start:end])
        )
        if pending:
            if index == 0:
                five_prime.extend(pending)
            else:
                internal.setdefault(index - 1, []).extend(pending)
            pending.clear()
        i = end

    while i < n:
        c = s[i]

        if c == "/":
            j = s.find("/", i + 1)
            if j == -1:
                errors.append(Problem("unterminated modification code", i, s[i:]))
                break
            code = s[i + 1:j]
            if not code or any(ch not in MOD_CODE_CHARS for ch in code):
                errors.append(
                    Problem(
                        "modification code must contain only letters, digits, - and _: "
                        f"'{code}'",
                        i,
                        s[i:j + 1],
                    )
                )
                i = j + 1
                continue
            # IDT writes two adjacent codes as /5Phos/iSp18/, so a closing
            # slash can also open the next code.  Collect the whole chain,
            # then read base codes (/i2FA/, /i2MOErC/) as nucleotides rather
            # than as tags attached to a later base.
            codes = [(code, i, j + 1)]
            while True:
                k = s.find("/", j + 1)
                if k == -1:
                    break
                following = s[j + 1:k]
                if (
                    not following
                    or not all(ch in MOD_CODE_CHARS for ch in following)
                    or _looks_like_sequence(following)
                ):
                    break
                codes.append((following, j, k + 1))
                j = k
            last_is_base = codes[-1][0] in SUGAR_BASE_CODES
            for code_text, start, end in codes:
                sugar_base = SUGAR_BASE_CODES.get(code_text)
                if sugar_base is not None:
                    base, sugar = sugar_base
                    emit(base, (base,), sugar, start, end)
                else:
                    pending.append(Modification(code_text, s[start:end]))
            if not last_is_base:
                i = codes[-1][2]
            continue

        if c == "(":
            j = s.find(")", i + 1)
            if j == -1:
                errors.append(Problem("unterminated custom mix", i, s[i:]))
                break
            token = s[i:j + 1]
            inner = s[i + 1:j]
            name = ""
            digits = inner
            if ":" in inner:
                name, _, digits = inner.partition(":")
            if len(digits) != 8 or not digits.isdigit():
                errors.append(
                    Problem(
                        "custom mix needs 8 digits, two per base in A C G T order",
                        i,
                        token,
                    )
                )
                i = j + 1
                continue
            percentages = (
                int(digits[0:2]),
                int(digits[2:4]),
                int(digits[4:6]),
                int(digits[6:8]),
            )
            if sum(percentages) != 100:
                errors.append(
                    Problem("custom mix percentages must sum to 100", i, token)
                )
                i = j + 1
                continue
            emit("?", (), "d", i, j + 1)
            mixes[len(nucleotides) - 1] = CustomMix(name, percentages)
            continue

        if c == "+":
            nxt = s[i + 1] if i + 1 < n else ""
            if nxt.upper() in ("A", "C", "G", "T"):
                base = nxt.upper()
                emit(base, (base,), "LNA", i, i + 2)
            else:
                errors.append(
                    Problem("'+' must be followed by A, C, G or T", i, "+")
                )
                i += 1
            continue

        if c in "rm" and i + 1 < n and s[i + 1].upper() in ("A", "C", "G", "U", "T"):
            base = s[i + 1].upper()
            if base == "T":
                errors.append(
                    Problem(
                        f"{c}T is not valid IDT notation; write {c}U", i, s[i:i + 2]
                    )
                )
                emit("U", ("U",), c, i, i + 2)
            else:
                emit(base, (base,), c, i, i + 2)
            continue

        if c in "rm":
            errors.append(
                Problem(f"{c} must be followed by A, C, G or U", i, c)
            )
            i += 1
            continue

        upper = c.upper()
        if upper in PLAIN_BASES:
            emit(upper, (upper,), "d", i, i + 1)
            continue

        if c in IUB_ALTERNATIVES:
            emit(c, IUB_ALTERNATIVES[c], "d", i, i + 1)
            continue

        errors.append(Problem(f"unexpected character '{c}'", i, c))
        i += 1

    if pending:
        if nucleotides:
            three_prime.extend(pending)
        else:
            five_prime.extend(pending)
        pending.clear()

    if not nucleotides and s:
        warnings.append(Problem("sequence contains no nucleotides", 0, ""))

    return Sequence(
        nucleotides=nucleotides,
        five_prime=five_prime,
        three_prime=three_prime,
        internal=internal,
        mixes=mixes,
        warnings=warnings,
        errors=errors,
    )


def summary_line(seq: Sequence) -> str:
    """One-line description: length, PS linkage count, modification count."""
    if not seq.nucleotides:
        return ""
    ps = sum(1 for nt in seq.nucleotides if nt.ps_after)
    modifications = (
        len(seq.five_prime)
        + len(seq.three_prime)
        + sum(len(mods) for mods in seq.internal.values())
    )
    return (
        f"{len(seq.nucleotides)} nt · {ps} PS linkage"
        + ("s" if ps != 1 else "")
        + f" · {modifications} modification"
        + ("s" if modifications != 1 else "")
    )
