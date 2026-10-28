"""Render a parsed IDT oligonucleotide as an SVG schematic.

Standard library only: this module must stay importable with bare ``py``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

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

SUGAR_ICON_TEXT = {"r": "OH", "m": "OMe", "LNA": "LNA", "F": "F", "MOE": "MOE"}
SUGAR_LABEL = {
    "r": "2'-OH",
    "m": "2'-OMe",
    "LNA": "LNA",
    "F": "2'-F",
    "MOE": "2'-MOE",
}
SUGAR_MEANING = {
    "r": "ribose (2'-OH)",
    "m": "2'-O-methyl",
    "LNA": "locked nucleic acid (LNA)",
    "F": "2'-fluoro",
    "MOE": "2'-O-methoxyethyl",
}
SUGAR_ICON_GAP = 4
SUGAR_ICON_W = 44
SUGAR_ICON_H = 30
SUGAR_BOTTOM_EXTRA = 20

# Secondary-structure geometry.  Coordinates are SVG user units.
STRUCT_TILE_R = 12
STRUCT_ROW_H = 28
STRUCT_ARM_GAP = 64
STRUCT_BULGE_GAP = 22
STRUCT_LOOP_SPACING = 32
STRUCT_PAD = 28
STRUCT_BASE_STEP = 28
STRUCT_ITEM_GAP = 24
STRUCT_LOOP_GAP = 48
STRUCT_LOOP_SIDE = 22
STRUCT_STEM_WIDTH = STRUCT_ARM_GAP + 2 * STRUCT_TILE_R


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
    cx = x + TILE_W / 2
    cy = y + TILE_H / 2
    parts = ['  <g class="tile">']
    parts.append(
        f'    <circle cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="32" fill="{fill}"/>'
    )
    if hatch:
        parts.append(
            f'    <circle cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="32" fill="url(#hatch)"/>'
        )
    parts.append(
        f'    <text x="{_fmt(cx)}" y="{_fmt(cy)}" '
        f'text-anchor="middle" dominant-baseline="central" font-family="sans-serif" '
        f'font-size="30" font-weight="bold" fill="{text_colour}">'
        f"{_esc(nucleotide.display)}</text>"
    )
    if nucleotide.sugar != "d":
        parts.append(
            _sugar_icon_markup(
                nucleotide.sugar, cx, y + TILE_H + SUGAR_ICON_GAP
            )
        )
    if mix is not None:
        parts.extend(_mix_strips(mix.percentages, x, y))
    parts.append("  </g>")
    return parts


def _mix_strips(percentages: tuple[int, int, int, int], x: int, y: int) -> list[str]:
    """Four proportional strips inside the '?' circle, A C G T left to right."""
    parts = []
    cursor = x + 12
    strip_y = y + TILE_H / 2 + 20
    for glyph, percent in zip("ACGT", percentages):
        width = max(2, percent * 40 // 100)
        parts.append(
            f'    <rect x="{_fmt(cursor)}" y="{_fmt(strip_y)}" '
            f'width="{width}" height="4" fill="{BASE_COLOURS[glyph][0]}"/>'
        )
        cursor += width
    return parts


def _sugar_icon_body(sugar: str) -> str:
    """Mini pentose ring; origin at C1', the base attachment, y grows downward."""
    parts = [
        '<line x1="0" y1="-4" x2="0" y2="0" stroke="#1A1A1A" stroke-width="1.2"/>',
        '<polygon points="0,0 9.51,6.91 5.88,18.09 -5.88,18.09 -9.51,6.91" '
        'fill="#FFFFFF" stroke="#1A1A1A" stroke-width="1.2" '
        'stroke-linejoin="round"/>',
        '<line x1="5.88" y1="18.09" x2="8.2" y2="19.3" '
        'stroke="#1A1A1A" stroke-width="1.2"/>',
    ]
    if sugar == "LNA":
        parts.append(
            '<line x1="5.88" y1="18.09" x2="-9.51" y2="6.91" stroke="#0072B2" '
            'stroke-width="1.8" stroke-linecap="round"/>'
        )
    text = _esc(SUGAR_ICON_TEXT.get(sugar, sugar))
    parts.append(
        f'<text x="9.5" y="19.3" dominant-baseline="central" '
        f'font-family="sans-serif" font-size="9" font-weight="600" '
        f'fill="#1A1A1A">{text}</text>'
    )
    return "\n      ".join(parts)


def _sugar_icon_markup(sugar: str, cx: float, top: float) -> str:
    return (
        f'    <g class="sugar-icon" '
        f'transform="translate({_fmt(cx)}, {_fmt(top)})">\n'
        f"      {_sugar_icon_body(sugar)}\n"
        f"    </g>"
    )


def sugar_icon_svg(sugar: str) -> str:
    """Standalone icon for the page legend."""
    meaning = _esc(SUGAR_MEANING.get(sugar, sugar))
    return (
        '<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{SUGAR_ICON_W}" height="{SUGAR_ICON_H}" '
        'viewBox="-12 -4 44 30" class="sugar-icon" role="img" '
        f'aria-label="{meaning} sugar icon">'
        + _sugar_icon_body(sugar)
        + "</svg>"
    )


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
    has_sugar = any(nucleotide.sugar != "d" for nucleotide in seq.nucleotides)
    bottom_extra = SUGAR_BOTTOM_EXTRA if has_sugar else 0
    height = top + TILE_H + 56 + bottom_extra

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
            f'  <text x="{pad_x}" y="{_fmt(tile_y + TILE_H + 26 + bottom_extra)}" '
            f'text-anchor="start" font-family="sans-serif" font-size="14" '
            f'fill="{INK}">5\'</text>',
            f'  <text x="{_fmt(_tile_x(n - 1, pad_x) + TILE_W)}" '
            f'y="{_fmt(tile_y + TILE_H + 26 + bottom_extra)}" text-anchor="end" '
            f'font-family="sans-serif" font-size="14" fill="{INK}">3\'</text>',
        ]
    )
    lines.append("</svg>")
    return "\n".join(lines)


def _structure_colours(base: str) -> tuple[str, str]:
    return BASE_COLOURS.get(base.upper(), DEGENERATE_COLOURS)


def _structure_tile(base: str, cx: float, cy: float, extra: str) -> list[str]:
    fill, text_colour = _structure_colours(base)
    classes = "base" + (f" {extra}" if extra else "")
    return [
        f'  <g class="{classes}">',
        f'    <circle cx="{_fmt(cx)}" cy="{_fmt(cy)}" r="{STRUCT_TILE_R}" '
        f'fill="{fill}"/>',
        f'    <text x="{_fmt(cx)}" y="{_fmt(cy)}" text-anchor="middle" '
        f'dominant-baseline="central" font-family="sans-serif" font-size="12" '
        f'font-weight="bold" fill="{text_colour}">{_esc(base.upper())}</text>',
        "  </g>",
    ]


def _validated_pairs(
    n: int, pairs: tuple[tuple[int, int], ...]
) -> tuple[tuple[int, int], ...]:
    cleaned: list[tuple[int, int]] = []
    seen: set[int] = set()
    for pair in pairs:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError("structure pairs must be (i, j) tuples")
        try:
            i, j = int(pair[0]), int(pair[1])
        except (TypeError, ValueError):
            raise ValueError("structure pairs must be integer indices") from None
        if not (0 <= i < j < n):
            raise ValueError("structure pair index out of range")
        if i in seen or j in seen:
            raise ValueError("a base cannot appear in two base pairs")
        seen.add(i)
        seen.add(j)
        cleaned.append((i, j))

    cleaned.sort()
    for index, (i, j) in enumerate(cleaned):
        for k, ell in cleaned[index + 1:]:
            if i < k < j < ell:
                raise ValueError("structure pairs must not cross")
    return tuple(cleaned)


@dataclass
class _TreeHelix:
    pairs: list[tuple[int, int]]
    outer_i: int
    outer_j: int
    inner_i: int
    inner_j: int
    parent: _TreeHelix | None = None
    children: list[_TreeHelix] = field(default_factory=list)
    segments: list[list[int]] = field(default_factory=list)
    width: float = 0.0


def _hairpin_radius(loop_count: int) -> float:
    return max(30.0, (loop_count + 1) * STRUCT_LOOP_SPACING / math.pi)


def _loop_segments(
    start: int, end: int, children: list[_TreeHelix]
) -> list[list[int]]:
    segments: list[list[int]] = []
    cursor = start
    for child in children:
        if not (start <= child.outer_i and child.outer_j < end):
            raise ValueError("child helix is outside its parent loop")
        segments.append(list(range(cursor, child.outer_i)))
        cursor = child.outer_j + 1
    segments.append(list(range(cursor, end)))
    return segments


def _build_helix_tree(
    n: int, pairs: tuple[tuple[int, int], ...]
) -> list[_TreeHelix]:
    helices: list[_TreeHelix] = []
    current: list[tuple[int, int]] = []
    for pair in pairs:
        if current and pair == (current[-1][0] + 1, current[-1][1] - 1):
            current.append(pair)
        else:
            if current:
                helices.append(_new_helix(current))
            current = [pair]
    if current:
        helices.append(_new_helix(current))

    for helix in helices:
        enclosing = [
            other
            for other in helices
            if other is not helix
            and other.outer_i < helix.outer_i
            and helix.outer_j < other.outer_j
            and other.inner_i < helix.outer_i
            and helix.outer_j < other.inner_j
        ]
        if enclosing:
            helix.parent = max(enclosing, key=lambda other: other.outer_i)
            helix.parent.children.append(helix)
    for helix in helices:
        helix.children.sort(key=lambda child: child.outer_i)
        helix.segments = _loop_segments(
            helix.inner_i + 1, helix.inner_j, helix.children
        )
    top_level = [helix for helix in helices if helix.parent is None]
    top_level.sort(key=lambda helix: helix.outer_i)
    return top_level


def _new_helix(pairs: list[tuple[int, int]]) -> _TreeHelix:
    outer_i, outer_j = pairs[0]
    inner_i, inner_j = pairs[-1]
    return _TreeHelix(
        pairs=list(pairs),
        outer_i=outer_i,
        outer_j=outer_j,
        inner_i=inner_i,
        inner_j=inner_j,
    )


def _loop_items(
    segments: list[list[int]], children: list[_TreeHelix]
) -> list[tuple[str, object, float]]:
    items: list[tuple[str, object, float]] = []
    for index, segment in enumerate(segments):
        for base_index in segment:
            items.append(("base", base_index, STRUCT_BASE_STEP))
        if index < len(children):
            child = children[index]
            items.append(("helix", child, child.width))
    return items


def _content_width(items: list[tuple[str, object, float]]) -> float:
    if not items:
        return 0.0
    return sum(item[2] for item in items) + STRUCT_ITEM_GAP * (len(items) - 1)


def _compute_helix_widths(helices: list[_TreeHelix]) -> None:
    for child in helices:
        _compute_helix_widths(child.children)
        if child.children:
            content = _content_width(_loop_items(child.segments, child.children))
            loop_width = max(
                STRUCT_STEM_WIDTH, content + 2 * STRUCT_LOOP_SIDE
            )
        else:
            loop_count = child.inner_j - child.inner_i - 1
            if loop_count:
                radius = _hairpin_radius(loop_count)
                loop_width = max(
                    STRUCT_STEM_WIDTH, 2 * radius + 2 * STRUCT_TILE_R
                )
            else:
                loop_width = STRUCT_STEM_WIDTH
        child.width = max(STRUCT_STEM_WIDTH, loop_width)


def _layout_helix(
    helix: _TreeHelix,
    center_x: float,
    base_y: float,
    positions: dict[int, tuple[float, float]],
    classes: dict[int, str],
) -> None:
    left_x = center_x - STRUCT_ARM_GAP / 2
    right_x = center_x + STRUCT_ARM_GAP / 2
    for row, (left_base, right_base) in enumerate(helix.pairs):
        y = base_y - row * STRUCT_ROW_H
        positions[left_base] = (left_x, y)
        positions[right_base] = (right_x, y)
        classes[left_base] = "paired"
        classes[right_base] = "paired"

    top_y = base_y - (len(helix.pairs) - 1) * STRUCT_ROW_H
    if not helix.children:
        loop_count = helix.inner_j - helix.inner_i - 1
        if loop_count:
            radius = _hairpin_radius(loop_count)
            for offset in range(loop_count):
                angle = math.pi - (offset + 1) * math.pi / (loop_count + 1)
                base_index = helix.inner_i + 1 + offset
                positions[base_index] = (
                    center_x + radius * math.cos(angle),
                    top_y - radius * math.sin(angle),
                )
                classes[base_index] = "loop"
        return

    items = _loop_items(helix.segments, helix.children)
    content = _content_width(items)
    loop_y = top_y - STRUCT_LOOP_GAP
    cursor = center_x - content / 2
    for kind, item, _item_width in items:
        if kind == "base":
            base_index = int(item)
            positions[base_index] = (
                cursor + STRUCT_BASE_STEP / 2,
                loop_y,
            )
            classes[base_index] = "loop"
            cursor += STRUCT_BASE_STEP
        else:
            child = item
            assert isinstance(child, _TreeHelix)
            child_center = cursor + child.width / 2
            _layout_helix(child, child_center, loop_y, positions, classes)
            cursor += child.width
        cursor += STRUCT_ITEM_GAP


def _layout_exterior(
    n: int,
    top_level: list[_TreeHelix],
    positions: dict[int, tuple[float, float]],
    classes: dict[int, str],
) -> None:
    segments = _loop_segments(0, n, top_level)
    items = _loop_items(segments, top_level)
    cursor = 0.0
    for kind, item, _item_width in items:
        if kind == "base":
            base_index = int(item)
            positions[base_index] = (
                cursor + STRUCT_BASE_STEP / 2,
                0.0,
            )
            classes[base_index] = ""
            cursor += STRUCT_BASE_STEP
        else:
            child = item
            assert isinstance(child, _TreeHelix)
            child_center = cursor + child.width / 2
            _layout_helix(child, child_center, 0.0, positions, classes)
            cursor += child.width
        cursor += STRUCT_ITEM_GAP


def _render_structure_tree(
    bases: str, pairs: tuple[tuple[int, int], ...]
) -> str:
    top_level = _build_helix_tree(len(bases), pairs)
    if not top_level:
        raise ValueError("no helices to draw")
    _compute_helix_widths(top_level)

    positions: dict[int, tuple[float, float]] = {}
    classes: dict[int, str] = {}
    _layout_exterior(len(bases), top_level, positions, classes)
    if len(positions) != len(bases):
        raise ValueError("structure tree does not cover every base")

    xs = [position[0] for position in positions.values()]
    ys = [position[1] for position in positions.values()]
    min_cx, max_cx = min(xs), max(xs)
    min_cy, max_cy = min(ys), max(ys)
    left = min_cx - STRUCT_TILE_R - STRUCT_PAD - 22
    right = max_cx + STRUCT_TILE_R + STRUCT_PAD + 22
    top = min_cy - STRUCT_TILE_R - STRUCT_PAD
    bottom = max_cy + STRUCT_TILE_R + STRUCT_PAD
    width = right - left
    height = bottom - top

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_fmt(width)}" '
        f'height="{_fmt(height)}" viewBox="{_fmt(left)} {_fmt(top)} '
        f'{_fmt(width)} {_fmt(height)}" class="structure-svg" '
        f'data-layout="tree" role="img" '
        f'aria-label="Predicted secondary structure">'
    ]

    for base_index in range(len(bases) - 1):
        x1, y1 = positions[base_index]
        x2, y2 = positions[base_index + 1]
        lines.append(
            f'  <line x1="{_fmt(x1)}" y1="{_fmt(y1)}" '
            f'x2="{_fmt(x2)}" y2="{_fmt(y2)}" '
            f'stroke="{LINK_COLOUR}" stroke-width="2" class="backbone"/>'
        )

    for left_base, right_base in pairs:
        x1, y1 = positions[left_base]
        x2, _y2 = positions[right_base]
        lines.append(
            f'  <line x1="{_fmt(x1 + STRUCT_TILE_R)}" y1="{_fmt(y1)}" '
            f'x2="{_fmt(x2 - STRUCT_TILE_R)}" y2="{_fmt(y1)}" '
            f'stroke="#0072B2" stroke-width="2" class="base-pair"/>'
        )

    for base_index in range(len(bases)):
        cx, cy = positions[base_index]
        lines.extend(
            _structure_tile(bases[base_index], cx, cy, classes[base_index])
        )

    first_x, _first_y = positions[0]
    last_x, _last_y = positions[len(bases) - 1]
    lines.extend(
        [
            f'  <text x="{_fmt(first_x - 20)}" y="4" '
            f'text-anchor="end" font-family="sans-serif" font-size="14" '
            f'fill="{INK}" class="terminus">5\'</text>',
            f'  <text x="{_fmt(last_x + 20)}" y="4" '
            f'text-anchor="start" font-family="sans-serif" font-size="14" '
            f'fill="{INK}" class="terminus">3\'</text>',
        ]
    )
    lines.append("</svg>")
    return "\n".join(lines)


def render_structure_svg(
    bases: str, pairs: tuple[tuple[int, int], ...]
) -> tuple[str, str]:
    """Return (svg, layout); layout is ``stem-loop``, ``tree``, ``arc``, or ``none``."""
    n = len(bases)
    if n == 0:
        raise ValueError("structure needs bases and valid pairs")
    ordered_pairs = _validated_pairs(n, tuple(pairs))
    if not ordered_pairs:
        return "", "none"

    hairpin_loops = [
        pair
        for pair in ordered_pairs
        if not any(
            pair[0] < other[0] and other[1] < pair[1]
            for other in ordered_pairs
        )
    ]
    stem_loop = (
        len(hairpin_loops) == 1
        and all(
            ordered_pairs[k + 1][1] < ordered_pairs[k][1]
            for k in range(len(ordered_pairs) - 1)
        )
        and ordered_pairs[-1] == hairpin_loops[0]
    )
    if not stem_loop:
        try:
            return _render_structure_tree(bases, ordered_pairs), "tree"
        except Exception:
            return _render_structure_arc(bases, ordered_pairs), "arc"

    left_x = 0
    right_x = STRUCT_ARM_GAP
    positions: dict[int, tuple[float, float]] = {}
    classes: dict[int, str] = {}
    pair_rows: list[float] = []

    for index, (left_base, right_base) in enumerate(ordered_pairs):
        if index == 0:
            y = 0
        else:
            prev_left, prev_right = ordered_pairs[index - 1]
            gap5 = left_base - prev_left - 1
            gap3 = prev_right - right_base - 1
            rows = max(gap5, gap3)
            for row in range(rows):
                row_y = pair_rows[-1] - (row + 1) * STRUCT_ROW_H
                if row < gap5:
                    base_index = prev_left + 1 + row
                    positions[base_index] = (left_x - STRUCT_BULGE_GAP, row_y)
                    classes[base_index] = "bulge"
                if row < gap3:
                    base_index = prev_right - 1 - row
                    positions[base_index] = (right_x + STRUCT_BULGE_GAP, row_y)
                    classes[base_index] = "bulge"
            y = pair_rows[-1] - (rows + 1) * STRUCT_ROW_H
        positions[left_base] = (left_x, y)
        positions[right_base] = (right_x, y)
        classes[left_base] = "paired"
        classes[right_base] = "paired"
        pair_rows.append(y)

    inner_left, inner_right = ordered_pairs[-1]
    inner_y = pair_rows[-1]
    loop_count = inner_right - inner_left - 1
    if loop_count:
        radius = max(
            30.0,
            (loop_count + 1) * STRUCT_LOOP_SPACING / math.pi,
        )
        mid_x = STRUCT_ARM_GAP / 2
        for offset in range(loop_count):
            angle = math.pi - (offset + 1) * math.pi / (loop_count + 1)
            base_index = inner_left + 1 + offset
            positions[base_index] = (
                mid_x + radius * math.cos(angle),
                inner_y - radius * math.sin(angle),
            )
            classes[base_index] = "loop"

    outer_left, outer_right = ordered_pairs[0]
    for base_index in range(outer_left):
        distance = outer_left - 1 - base_index
        positions[base_index] = (
            left_x - 32 - distance * 28,
            0,
        )
        classes[base_index] = ""
    for base_index in range(outer_right + 1, n):
        distance = base_index - outer_right - 1
        positions[base_index] = (
            right_x + 32 + distance * 28,
            0,
        )
        classes[base_index] = ""

    if len(positions) != n:
        raise ValueError("structure pairs do not cover every base")

    xs = [position[0] for position in positions.values()]
    ys = [position[1] for position in positions.values()]
    min_cx, max_cx = min(xs), max(xs)
    min_cy, max_cy = min(ys), max(ys)
    left = min_cx - STRUCT_TILE_R - STRUCT_PAD - 22
    right = max_cx + STRUCT_TILE_R + STRUCT_PAD + 22
    top = min_cy - STRUCT_TILE_R - STRUCT_PAD
    bottom = max_cy + STRUCT_TILE_R + STRUCT_PAD
    width = right - left
    height = bottom - top

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_fmt(width)}" '
        f'height="{_fmt(height)}" viewBox="{_fmt(left)} {_fmt(top)} '
        f'{_fmt(width)} {_fmt(height)}" class="structure-svg" '
        f'data-layout="stem-loop" role="img" '
        f'aria-label="Predicted secondary structure">'
    ]

    for base_index in range(n - 1):
        x1, y1 = positions[base_index]
        x2, y2 = positions[base_index + 1]
        lines.append(
            f'  <line x1="{_fmt(x1)}" y1="{_fmt(y1)}" '
            f'x2="{_fmt(x2)}" y2="{_fmt(y2)}" '
            f'stroke="{LINK_COLOUR}" stroke-width="2" class="backbone"/>'
        )

    for left_base, right_base in ordered_pairs:
        x1, y = positions[left_base]
        x2, _ = positions[right_base]
        lines.append(
            f'  <line x1="{_fmt(x1 + STRUCT_TILE_R)}" y1="{_fmt(y)}" '
            f'x2="{_fmt(x2 - STRUCT_TILE_R)}" y2="{_fmt(y)}" '
            f'stroke="#0072B2" stroke-width="2" class="base-pair"/>'
        )

    for base_index in range(n):
        cx, cy = positions[base_index]
        lines.extend(
            _structure_tile(bases[base_index], cx, cy, classes[base_index])
        )

    lines.extend(
        [
            f'  <text x="{_fmt(min_cx - STRUCT_TILE_R - 8)}" y="{_fmt(4)}" '
            f'text-anchor="end" font-family="sans-serif" font-size="14" '
            f'fill="{INK}" class="terminus">5\'</text>',
            f'  <text x="{_fmt(max_cx + STRUCT_TILE_R + 8)}" y="{_fmt(4)}" '
            f'text-anchor="start" font-family="sans-serif" font-size="14" '
            f'fill="{INK}" class="terminus">3\'</text>',
        ]
    )
    lines.append("</svg>")
    return "\n".join(lines), "stem-loop"


def _render_structure_arc(
    bases: str, pairs: tuple[tuple[int, int], ...]
) -> str:
    n = len(bases)
    heights = [min(260, 20 + (j - i) * 8) for i, j in pairs]
    max_height = max(heights)
    baseline = max_height + 20
    width = 2 * STRUCT_PAD + (n - 1) * 26
    height = baseline + 40

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{_fmt(width)}" '
        f'height="{_fmt(height)}" viewBox="0 0 {_fmt(width)} {_fmt(height)}" '
        f'class="structure-svg" data-layout="arc" role="img" '
        f'aria-label="Predicted secondary structure">'
    ]

    for i, j in sorted(pairs, key=lambda pair: (pair[1] - pair[0], pair), reverse=True):
        x_i = STRUCT_PAD + i * 26
        x_j = STRUCT_PAD + j * 26
        curve_height = min(260, 20 + (j - i) * 8)
        lines.append(
            f'  <path d="M {_fmt(x_i)} {_fmt(baseline)} C {_fmt(x_i)} '
            f'{_fmt(baseline - curve_height)}, {_fmt(x_j)} '
            f'{_fmt(baseline - curve_height)}, {_fmt(x_j)} {_fmt(baseline)}" '
            f'fill="none" stroke="#0072B2" stroke-width="1.5" '
            f'opacity="0.7" class="pair-arc"/>'
        )

    first_x = STRUCT_PAD
    last_x = STRUCT_PAD + (n - 1) * 26
    lines.append(
        f'  <line x1="{_fmt(first_x)}" y1="{_fmt(baseline)}" '
        f'x2="{_fmt(last_x)}" y2="{_fmt(baseline)}" '
        f'stroke="{LINK_COLOUR}" stroke-width="2" class="backbone"/>'
    )

    for index, base in enumerate(bases):
        cx = STRUCT_PAD + index * 26
        fill, text_colour = _structure_colours(base)
        lines.extend(
            [
                '  <g class="base">',
                f'    <circle cx="{_fmt(cx)}" cy="{_fmt(baseline)}" r="10" '
                f'fill="{fill}"/>',
                f'    <text x="{_fmt(cx)}" y="{_fmt(baseline)}" '
                f'text-anchor="middle" dominant-baseline="central" '
                f'font-family="sans-serif" font-size="10" font-weight="bold" '
                f'fill="{text_colour}">{_esc(base.upper())}</text>',
                "  </g>",
            ]
        )

    lines.extend(
        [
            f'  <text x="2" y="{_fmt(baseline + 28)}" text-anchor="start" '
            f'font-family="sans-serif" font-size="14" fill="{INK}" '
            f'class="terminus">5\'</text>',
            f'  <text x="{_fmt(width - 2)}" y="{_fmt(baseline + 28)}" '
            f'text-anchor="end" font-family="sans-serif" font-size="14" '
            f'fill="{INK}" class="terminus">3\'</text>',
        ]
    )
    lines.append("</svg>")
    return "\n".join(lines)


def structure_legend_items(bases: str) -> list[dict]:
    """One base-legend entry per distinct base, in first-appearance order."""
    items: list[dict] = []
    seen: set[str] = set()
    for base in bases:
        glyph = base.upper()
        if glyph in seen:
            continue
        seen.add(glyph)
        fill, text = _structure_colours(glyph)
        items.append(
            {"kind": "base", "label": glyph, "color": fill, "text": text}
        )
    return items


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
        sugar = nucleotide.sugar
        if sugar == "d" or sugar in seen_sugars:
            continue
        seen_sugars.add(sugar)
        items.append(
            {
                "kind": "sugar",
                "label": SUGAR_LABEL.get(sugar, sugar),
                "meaning": SUGAR_MEANING.get(sugar, sugar),
                "icon": sugar_icon_svg(sugar),
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
