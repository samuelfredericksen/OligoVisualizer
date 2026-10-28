"""Render a parsed IDT oligonucleotide as an SVG schematic.

Standard library only: this module must stay importable with bare ``py``.
"""

from __future__ import annotations

import math

from notation import Sequence

# --- geometry -------------------------------------------------------------

TILE_W = 64
TILE_H = 64
LINK_W = 30
TAG_H = 26
TAG_GAP = 4
LANE_GAP = 8

# --- palette (Okabe-Ito, colourblind-safe) --------------------------------

BASE_COLOURS = {
    "A": ("#009E73", "#FFFFFF"),
    "C": ("#0072B2", "#FFFFFF"),
    "G": ("#E69F00", "#1A1A1A"),
    "T": ("#D55E00", "#FFFFFF"),
    "U": ("#CC79A7", "#1A1A1A"),
    "I": ("#56B4E9", "#1A1A1A"),
}
DEGENERATE_COLOURS = ("#7F7F7F", "#FFFFFF")
MIX_COLOURS = ("#E0E0E0", "#1A1A1A")

LINK_COLOUR = "#6B6B6B"
PS_COLOUR = "#B00020"
TAG_FILL = "#FFF4D6"
TAG_STROKE = "#B58900"
TAG_TEXT = "#1A1A1A"
INK = "#1A1A1A"

SUGAR_GLYPH = {"r": "r", "m": "m", "LNA": "+", "MOE": "MOE"}
SUGAR_MEANING = {
    "r": "ribose (2'-OH)",
    "m": "2'-O-methyl",
    "+": "locked nucleic acid (LNA)",
    "MOE": "2'-MOE",
}


def tag_width(raw: str) -> float:
    """Width in px of the tag box for a modification written as ``raw``."""
    return 10 + 7.2 * len(raw)


# --- small helpers --------------------------------------------------------


def _fmt(value: float) -> str:
    """Deterministic number formatting: no trailing ``.0``, at most 2 decimals."""
    if float(value).is_integer():
        return str(int(value))
    return ("%.2f" % value).rstrip("0").rstrip(".")


def _esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _tile_x(index: int, pad_x: int) -> int:
    return pad_x + index * (TILE_W + LINK_W)


def _pad_x(raws: list[str]) -> int:
    if not raws:
        return 40
    return max(40, math.ceil(max(tag_width(raw) for raw in raws) / 2) + 12)


def _tile_colours(nucleotide) -> tuple[str, str]:
    display = nucleotide.display
    if display == "?":
        return MIX_COLOURS
    if display in BASE_COLOURS:
        return BASE_COLOURS[display]
    return DEGENERATE_COLOURS


def _is_degenerate(nucleotide) -> bool:
    return len(nucleotide.alternatives) > 1


def _all_raws(seq: Sequence) -> list[str]:
    raws = [mod.raw for mod in seq.five_prime]
    for index in sorted(seq.internal):
        raws.extend(mod.raw for mod in seq.internal[index])
    raws.extend(mod.raw for mod in seq.three_prime)
    return raws


def _tag_entries(seq: Sequence, pad_x: int) -> list[dict]:
    """Every tag as ``{"raw", "x", "kind"}`` in sequence order."""
    n = len(seq.nucleotides)
    entries: list[dict] = []
    for mod in seq.five_prime:
        entries.append(
            {"raw": mod.raw, "x": _tile_x(0, pad_x) + TILE_W / 2, "kind": "5'"}
        )
    for index in sorted(seq.internal):
        for mod in seq.internal[index]:
            if index + 1 < n:
                x = _tile_x(index, pad_x) + TILE_W + LINK_W / 2
            else:
                # A modification trailing the last nucleotide rides the stub.
                x = _tile_x(index, pad_x) + TILE_W + 7
            entries.append({"raw": mod.raw, "x": x, "kind": "internal"})
    for mod in seq.three_prime:
        entries.append(
            {
                "raw": mod.raw,
                "x": _tile_x(n - 1, pad_x) + TILE_W / 2,
                "kind": "3'",
            }
        )
    return entries


def _pack_tags(entries: list[dict]) -> int:
    """Assign each tag the lowest free lane; returns the number of lanes used."""
    lane_right_edge: list[float] = []
    for entry in sorted(entries, key=lambda item: item["x"]):
        half = tag_width(entry["raw"]) / 2
        left = entry["x"] - half
        right = entry["x"] + half
        for lane, last_right in enumerate(lane_right_edge):
            if last_right + LANE_GAP < left:
                lane_right_edge[lane] = right
                entry["lane"] = lane
                break
        else:
            lane_right_edge.append(right)
            entry["lane"] = len(lane_right_edge) - 1
    return len(lane_right_edge)


# --- markup builders ------------------------------------------------------


def _tile_markup(nucleotide, x: int, y: int, hatch: bool, mix) -> list[str]:
    fill, text_colour = _tile_colours(nucleotide)
    parts = ['  <g class="tile">']
    parts.append(
        f'    <rect x="{_fmt(x)}" y="{_fmt(y)}" width="{TILE_W}" '
        f'height="{TILE_H}" rx="6" fill="{fill}"/>'
    )
    if hatch:
        parts.append(
            f'    <rect x="{_fmt(x)}" y="{_fmt(y)}" width="{TILE_W}" '
            f'height="{TILE_H}" rx="6" fill="url(#hatch)"/>'
        )
    parts.append(
        f'    <text x="{_fmt(x + TILE_W / 2)}" y="{_fmt(y + TILE_H / 2)}" '
        f'text-anchor="middle" dominant-baseline="central" font-family="sans-serif" '
        f'font-size="30" font-weight="bold" fill="{text_colour}">'
        f"{_esc(nucleotide.display)}</text>"
    )
    glyph = SUGAR_GLYPH.get(nucleotide.sugar)
    if glyph:
        parts.append(
            f'    <text x="{_fmt(x + 56)}" y="{_fmt(y + 16)}" '
            f'text-anchor="middle" font-family="sans-serif" font-size="12" '
            f'font-weight="bold" fill="{text_colour}">{_esc(glyph)}</text>'
        )
    if mix is not None:
        parts.extend(_mix_strips(mix.percentages, x, y))
    parts.append("  </g>")
    return parts


def _mix_strips(percentages: tuple[int, int, int, int], x: int, y: int) -> list[str]:
    """Four proportional strips under the '?' glyph, A C G T left to right."""
    parts = []
    cursor = x + 4
    strip_y = y + TILE_H / 2 + 20
    for glyph, percent in zip("ACGT", percentages):
        width = max(2, percent * 46 // 100)
        parts.append(
            f'    <rect x="{_fmt(cursor)}" y="{_fmt(strip_y)}" '
            f'width="{width}" height="4" fill="{BASE_COLOURS[glyph][0]}"/>'
        )
        cursor += width
    return parts


def _link_markup(start: int, end: int, y: int, ps: bool) -> list[str]:
    if not ps:
        return [
            '  <g class="link">',
            f'    <line x1="{_fmt(start)}" y1="{_fmt(y)}" x2="{_fmt(end)}" '
            f'y2="{_fmt(y)}" stroke="{LINK_COLOUR}" stroke-width="2"/>',
            "  </g>",
        ]
    mid = (start + end) / 2
    return [
        '  <g class="link ps">',
        f'    <line x1="{_fmt(start)}" y1="{_fmt(y)}" x2="{_fmt(end)}" '
        f'y2="{_fmt(y)}" stroke="{PS_COLOUR}" stroke-width="4"/>',
        f'    <circle cx="{_fmt(mid)}" cy="{_fmt(y)}" r="9" fill="#FFFFFF" '
        f'stroke="{PS_COLOUR}" stroke-width="2"/>',
        f'    <text x="{_fmt(mid)}" y="{_fmt(y)}" text-anchor="middle" '
        f'dominant-baseline="central" font-family="sans-serif" font-size="10" '
        f'font-weight="bold" fill="{PS_COLOUR}">S</text>',
        "  </g>",
    ]


def _tag_markup(entry: dict, bottom: float, anchor_y: float) -> list[str]:
    width = tag_width(entry["raw"])
    left = entry["x"] - width / 2
    top = bottom - TAG_H
    centre = bottom - TAG_H / 2
    return [
        '  <g class="tag">',
        f'    <line x1="{_fmt(entry["x"])}" y1="{_fmt(bottom)}" '
        f'x2="{_fmt(entry["x"])}" y2="{_fmt(anchor_y)}" '
        f'stroke="{TAG_STROKE}" stroke-width="1"/>',
        f'    <rect x="{_fmt(left)}" y="{_fmt(top)}" width="{_fmt(width)}" '
        f'height="{TAG_H}" rx="4" fill="{TAG_FILL}" stroke="{TAG_STROKE}" '
        f'stroke-width="1"/>',
        f'    <text x="{_fmt(entry["x"])}" y="{_fmt(centre)}" text-anchor="middle" '
        f'dominant-baseline="central" font-family="monospace" font-size="12" '
        f'fill="{TAG_TEXT}">{_esc(entry["raw"])}</text>',
        "  </g>",
    ]


# --- public surface -------------------------------------------------------


def render_svg(seq: Sequence) -> str:
    """Return the SVG for ``seq``; raises for an empty sequence."""
    n = len(seq.nucleotides)
    if n == 0:
        raise ValueError("render_svg needs at least one nucleotide")

    pad_x = _pad_x(_all_raws(seq))
    entries = _tag_entries(seq, pad_x)
    lanes = _pack_tags(entries) if entries else 0

    top = 24 + lanes * (TAG_H + TAG_GAP) + 16
    tile_y = top
    width = pad_x + n * (TILE_W + LINK_W) - LINK_W + pad_x
    height = top + TILE_H + 56

    has_hatch = any(_is_degenerate(nt) for nt in seq.nucleotides)

    lines: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" class="oligo-svg" role="img" '
        f'aria-label="Oligonucleotide schematic">'
    ]
    if has_hatch:
        lines.extend(
            [
                "  <defs>",
                '    <pattern id="hatch" width="8" height="8" '
                'patternUnits="userSpaceOnUse" patternTransform="rotate(45)">',
                '      <line x1="0" y1="0" x2="0" y2="8" stroke="#FFFFFF" '
                'stroke-width="2" opacity="0.55"/>',
                "    </pattern>",
                "  </defs>",
            ]
        )

    # connectors first, so tiles and their labels sit on top
    connector_y = tile_y + TILE_H / 2
    for index in range(n - 1):
        lines.extend(
            _link_markup(
                _tile_x(index, pad_x) + TILE_W,
                _tile_x(index + 1, pad_x),
                connector_y,
                bool(seq.nucleotides[index].ps_after),
            )
        )
    if seq.nucleotides[-1].ps_after:
        stub_start = _tile_x(n - 1, pad_x) + TILE_W
        lines.extend(_link_markup(stub_start, stub_start + 14, connector_y, True))

    for entry in entries:
        anchor_y = tile_y if entry["kind"] != "internal" else connector_y
        bottom = tile_y - 12 - entry["lane"] * (TAG_H + TAG_GAP)
        lines.extend(_tag_markup(entry, bottom, anchor_y))

    for index, nucleotide in enumerate(seq.nucleotides):
        x = _tile_x(index, pad_x)
        lines.extend(
            _tile_markup(
                nucleotide,
                x,
                tile_y,
                _is_degenerate(nucleotide),
                seq.mixes.get(index),
            )
        )

    lines.extend(
        [
            f'  <text x="{pad_x}" y="{_fmt(tile_y + TILE_H + 26)}" '
            f'text-anchor="start" font-family="sans-serif" font-size="14" '
            f'fill="{INK}">5\'</text>',
            f'  <text x="{_fmt(_tile_x(n - 1, pad_x) + TILE_W)}" '
            f'y="{_fmt(tile_y + TILE_H + 26)}" text-anchor="end" '
            f'font-family="sans-serif" font-size="14" fill="{INK}">3\'</text>',
        ]
    )
    lines.append("</svg>")
    return "\n".join(lines)


def legend_items(seq: Sequence) -> list[dict]:
    """Legend entries for exactly what this sequence contains."""
    if not seq.nucleotides:
        return []

    items: list[dict] = []

    seen_glyphs: set[str] = set()
    for nucleotide in seq.nucleotides:
        glyph = nucleotide.display
        if glyph in seen_glyphs:
            continue
        seen_glyphs.add(glyph)
        fill, text = _tile_colours(nucleotide)
        if glyph == "?":
            label = "custom mix"
        elif _is_degenerate(nucleotide):
            label = f"{glyph} ({'/'.join(nucleotide.alternatives)})"
        else:
            label = glyph
        items.append({"kind": "base", "label": label, "color": fill, "text": text})

    seen_sugars: set[str] = set()
    for nucleotide in seq.nucleotides:
        if nucleotide.sugar == "d":
            continue
        glyph = SUGAR_GLYPH.get(nucleotide.sugar, nucleotide.sugar)
        if glyph in seen_sugars:
            continue
        seen_sugars.add(glyph)
        items.append(
            {
                "kind": "sugar",
                "label": glyph,
                "meaning": SUGAR_MEANING.get(glyph, nucleotide.sugar),
            }
        )

    if any(nucleotide.ps_after for nucleotide in seq.nucleotides):
        items.append({"kind": "link", "label": "phosphate"})
        items.append({"kind": "link", "label": "phosphorothioate"})

    seen_codes: set[str] = set()
    for mod in _ordered_modifications(seq):
        if mod.code in seen_codes:
            continue
        seen_codes.add(mod.code)
        items.append({"kind": "tag", "label": mod.raw})

    return items


def _ordered_modifications(seq: Sequence):
    for mod in seq.five_prime:
        yield mod
    for index in sorted(seq.internal):
        for mod in seq.internal[index]:
            yield mod
    for mod in seq.three_prime:
        yield mod
