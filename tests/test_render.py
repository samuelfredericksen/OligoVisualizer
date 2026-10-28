"""Tests for render.py (stdlib unittest; importable with bare `py`)."""

import re
import unittest
from unittest.mock import patch

from notation import parse
from render import (
    base_pair_kind,
    legend_items,
    pair_legend_items,
    render_homodimer_svg,
    render_structure_svg,
    render_svg,
    structure_legend_items,
    sugar_icon_svg,
)


class RenderSvgTest(unittest.TestCase):
    def test_plain_sequence_uses_base_palette_and_exact_width(self):
        svg = render_svg(parse("ACGT"))
        for colour in ("#009E73", "#0072B2", "#E69F00", "#D55E00"):
            self.assertIn(colour, svg)
        self.assertNotIn('class="link ps"', svg)
        self.assertIn('width="426"', svg)
        self.assertEqual(svg.count("<circle"), 4)

    def test_phosphorothioate_link_gets_badge(self):
        svg = render_svg(parse("A*C"))
        self.assertEqual(svg.count('class="link ps"'), 1)
        self.assertEqual(svg.count(">S<"), 1)

    def test_plain_link_has_no_ps_class(self):
        self.assertNotIn('class="link ps"', render_svg(parse("AC")))

    def test_modification_codes_drawn_as_tags(self):
        svg = render_svg(parse("/5Biosg/TTTT/3Bio/"))
        self.assertIn("/5Biosg/", svg)
        self.assertIn("/3Bio/", svg)

    def test_hatch_only_for_degenerate_tiles(self):
        self.assertIn("url(#hatch)", render_svg(parse("ACGTN")))
        self.assertNotIn("url(#hatch)", render_svg(parse("ACGT")))

    def test_sugar_icons(self):
        svg = render_svg(parse("mA+CrG"))
        self.assertEqual(svg.count('class="sugar-icon"'), 3)
        self.assertIn(">OMe<", svg)
        self.assertIn(">LNA<", svg)
        self.assertIn(">OH<", svg)
        self.assertNotIn(">m<", svg)
        self.assertNotIn(">+<", svg)

    def test_sugar_glyph_has_ring_oxygen_and_base(self):
        svg = sugar_icon_svg("m")
        self.assertIn(">O<", svg)
        self.assertIn(">N<", svg)
        self.assertIn(">OMe<", svg)
        self.assertIn(">3'<", svg)
        self.assertIn(">5'<", svg)

    def test_bases_are_circles(self):
        svg = render_svg(parse("A"))
        self.assertIn("<circle", svg)
        self.assertNotIn("<rect", svg)

    def test_modified_base_codes_draw_sugar_icons_not_tags(self):
        svg = render_svg(parse("/i2FA/i2FC/"))
        self.assertEqual(svg.count('class="sugar-icon"'), 2)
        self.assertIn(">F<", svg)
        self.assertNotIn("/i2FA/", svg)
        self.assertNotIn("/i2FC/", svg)

    def test_fluoro_thymidine_draws_f_icon_not_tag(self):
        svg = render_svg(parse("/i2FT/"))
        self.assertEqual(svg.count('<g class="tile">'), 1)
        self.assertEqual(svg.count('class="sugar-icon"'), 1)
        self.assertIn(">F<", svg)
        self.assertNotIn("/i2FT/", svg)

    def test_moe_codes_draw_moe_icons_and_ps(self):
        svg = render_svg(parse("/52MOErA/*/i2MOErC/*/32MOErT/"))
        self.assertEqual(svg.count('class="sugar-icon"'), 3)
        self.assertEqual(svg.count(">MOE<"), 3)
        self.assertEqual(svg.count('class="link ps"'), 2)

    def test_output_is_deterministic(self):
        self.assertEqual(
            render_svg(parse("mA*+C/3Bio/")), render_svg(parse("mA*+C/3Bio/"))
        )

    def test_empty_sequence_raises(self):
        with self.assertRaises(ValueError):
            render_svg(parse(""))


class PairClassTest(unittest.TestCase):
    def test_watson_crick_pairs_include_t_u_equivalence(self):
        for left, right in (
            ("A", "T"), ("T", "A"), ("A", "U"), ("U", "A"),
            ("G", "C"), ("C", "G"),
        ):
            self.assertEqual(base_pair_kind(left, right), "wc")

    def test_g_u_and_g_t_are_wobble(self):
        for left, right in (("G", "U"), ("U", "G"), ("G", "T"), ("T", "G")):
            self.assertEqual(base_pair_kind(left, right), "wobble")

    def test_other_known_bases_are_noncanonical(self):
        for left, right in (("G", "A"), ("A", "C"), ("C", "C"), ("U", "U")):
            self.assertEqual(base_pair_kind(left, right), "noncanonical")

    def test_degenerate_base_is_unknown(self):
        self.assertEqual(base_pair_kind("N", "A"), "unknown")
        self.assertEqual(base_pair_kind("A", "?"), "unknown")


class HomodimerRenderTest(unittest.TestCase):
    SEQUENCE = "GCGAAAAAGCGCTTTTGCGCAAAAGCGC"
    DIMER = (
        "5' GCGAAAAAGCGCTTTTGCGCAAAAGCGC\n"
        "           ||||||||||||||||||||        \n"
        "3'         CGCGAAAACGCGTTTTCGCGAAAAAGCG"
    )

    def test_two_antiparallel_rows_and_solid_pairs(self):
        svg = render_homodimer_svg(self.SEQUENCE, self.DIMER)
        self.assertIn('data-layout="homodimer"', svg)
        self.assertEqual(svg.count("<circle"), 56)
        self.assertEqual(svg.count('class="base-pair '), 20)
        self.assertEqual(svg.count('class="backbone"'), 54)
        self.assertIn(">5'<", svg)
        self.assertIn(">3'<", svg)
        self.assertNotIn("NaN", svg)

    def test_additional_pairs_are_faded(self):
        dimer = "5' ACGT\n   ||::\n3' TGCA"
        svg = render_homodimer_svg("ACGT", dimer)
        self.assertEqual(svg.count('class="base-pair '), 4)
        self.assertEqual(
            svg.count('class="base-pair wc additional"'), 2
        )
        self.assertEqual(svg.count('opacity="0.45"'), 2)
        self.assertEqual(svg.count('stroke-width="1.4"'), 2)

    def test_wobble_and_noncanonical_rungs_are_distinct(self):
        dimer = "5' GA\n   ||\n3' UC"
        svg = render_homodimer_svg("GA", dimer)
        self.assertIn('class="base-pair wobble"', svg)
        self.assertIn('stroke-dasharray="7 3"', svg)
        self.assertIn('class="base-pair noncanonical"', svg)
        self.assertIn('stroke-dasharray="2 4"', svg)

    def test_sugars_appear_on_both_copies(self):
        sugars = ["d"] * len(self.SEQUENCE)
        sugars[10] = "m"
        svg = render_homodimer_svg(
            self.SEQUENCE, self.DIMER, tuple(sugars)
        )
        self.assertEqual(svg.count('class="structure-sugar"'), 2)
        self.assertEqual(svg.count(">OMe<"), 2)
        origins = re.findall(
            r'transform="translate\(([^,]+), ([^)]+)\) scale\(0\.55\)"',
            svg,
        )
        self.assertEqual(len(origins), 2)
        self.assertLess(float(origins[0][1]), 0)
        self.assertGreater(float(origins[1][1]), 72)

    def test_sugar_length_mismatch_raises(self):
        with self.assertRaisesRegex(
            ValueError, "homodimer sugars must match bases"
        ):
            render_homodimer_svg(self.SEQUENCE, self.DIMER, ("m",))

    def test_malformed_alignment_raises(self):
        with self.assertRaisesRegex(ValueError, "unreadable self-dimer"):
            render_homodimer_svg("ACGT", "not a dimer")


class StructureRenderTest(unittest.TestCase):
    def test_stacked_stem_loop(self):
        pairs = tuple((index, 19 - index) for index in range(8))
        svg, layout = render_structure_svg("ACGTACGTACGTACGTACGT", pairs)
        self.assertEqual(layout, "stem-loop")
        self.assertEqual(svg.count('class="base-pair '), 8)
        self.assertEqual(svg.count('class="base loop"'), 4)
        self.assertEqual(svg, render_structure_svg("ACGTACGTACGTACGTACGT", pairs)[0])
        self.assertIn('width="', svg)
        self.assertIn('height="', svg)
        self.assertNotIn("NaN", svg)

    def test_stem_loop_wobble_line_is_dashed_teal(self):
        bases = "G" + "A" * 18 + "U"
        pairs = tuple((index, 19 - index) for index in range(8))
        svg, layout = render_structure_svg(bases, pairs)
        self.assertEqual(layout, "stem-loop")
        self.assertIn('class="base-pair wobble"', svg)
        self.assertIn('stroke-dasharray="7 3"', svg)

    def test_arc_wobble_line_is_dashed_teal(self):
        bases = "G" + "A" * 6 + "U" + "C" * 4 + "A" * 4
        pairs = ((0, 7), (1, 6), (8, 15), (9, 14))
        with patch("render._render_structure_tree", side_effect=ValueError):
            svg, layout = render_structure_svg(bases, pairs)
        self.assertEqual(layout, "arc")
        self.assertIn('class="pair-arc wobble"', svg)
        self.assertIn('stroke-dasharray="7 3"', svg)

    def test_structure_sugar_icons_and_legend(self):
        bases = "ACGTACGTACGTACGTACGT"
        pairs = tuple((index, 19 - index) for index in range(8))
        sugars = ["d"] * 20
        sugars[0] = "m"
        sugars[7] = "F"
        sugars[12] = "r"
        svg, layout = render_structure_svg(bases, pairs, tuple(sugars))
        self.assertEqual(layout, "stem-loop")
        self.assertEqual(svg.count('class="structure-sugar"'), 3)
        self.assertEqual(svg.count('class="base-pair '), 8)
        self.assertIn(">OMe<", svg)
        self.assertIn(">F<", svg)
        self.assertIn(">OH<", svg)

        items = structure_legend_items(bases, tuple(sugars))
        self.assertEqual(
            [item["label"] for item in items],
            [
                "A", "C", "G", "T", "2'-OMe", "2'-F", "2'-OH",
                "Watson–Crick", "G·U / G·T wobble",
                "Other noncanonical", "Unknown pair",
            ],
        )
        self.assertEqual([item["kind"] for item in items[:4]], ["base"] * 4)
        self.assertEqual([item["kind"] for item in items[4:7]], ["sugar"] * 3)
        self.assertEqual([item["kind"] for item in items[7:]], ["pair"] * 4)

    def test_paired_sugars_sit_outside_the_stem(self):
        bases = "ACGTACGTACGTACGTACGT"
        pairs = tuple((index, 19 - index) for index in range(8))
        sugars = ("m",) + ("d",) * 18 + ("F",)
        svg, layout = render_structure_svg(bases, pairs, sugars)
        self.assertEqual(layout, "stem-loop")
        origins = [
            float(value)
            for value in re.findall(
                r'transform="translate\(([^,]+), [^)]+\) scale\(0\.55\)"',
                svg,
            )
        ]
        self.assertEqual(len(origins), 2)
        self.assertLess(origins[0], 0)
        self.assertGreater(origins[1], 64)

    def test_tree_structure_carries_sugar_icons(self):
        bases = "AAAACCCCAAAAGGGG"
        pairs = ((0, 7), (1, 6), (8, 15), (9, 14))
        sugars = ["d"] * 16
        sugars[0] = "m"
        svg, layout = render_structure_svg(bases, pairs, tuple(sugars))
        self.assertEqual(layout, "tree")
        self.assertEqual(svg.count('class="structure-sugar"'), 1)
        self.assertIn(">OMe<", svg)

    def test_arc_fallback_carries_sugar_icons(self):
        bases = "AAAACCCCAAAAGGGG"
        pairs = ((0, 7), (1, 6), (8, 15), (9, 14))
        sugars = ["d"] * 16
        sugars[0] = "m"
        with patch("render._render_structure_tree", side_effect=ValueError):
            svg, layout = render_structure_svg(bases, pairs, tuple(sugars))
        self.assertEqual(layout, "arc")
        self.assertEqual(svg.count('class="structure-sugar"'), 1)
        self.assertIn(">OMe<", svg)

    def test_structure_sugars_must_match_bases(self):
        with self.assertRaisesRegex(
            ValueError, "structure sugars must match bases"
        ):
            render_structure_svg("ACGT", ((0, 3),), ("m",))

    def test_two_disparate_hairpins_use_tree_layout(self):
        pairs = ((0, 7), (1, 6), (8, 15), (9, 14))
        svg, layout = render_structure_svg("AAAACCCCAAAAGGGG", pairs)
        self.assertEqual(layout, "tree")
        self.assertEqual(svg.count('class="base-pair '), 4)
        self.assertEqual(svg.count('class="base loop"'), 8)
        self.assertIn('data-layout="tree"', svg)

    def test_two_hairpins_with_exterior_linker_use_tree(self):
        pairs = (
            (0, 24), (1, 23), (2, 22), (3, 21), (4, 20), (5, 19),
            (6, 18), (7, 17), (8, 16), (9, 15),
            (30, 42), (31, 41), (32, 40), (33, 39),
        )
        svg, layout = render_structure_svg("A" * 48, pairs)
        self.assertEqual(layout, "tree")
        self.assertEqual(svg.count('class="base-pair '), 14)
        self.assertEqual(svg.count('class="base loop"'), 10)

    def test_multiloop_uses_tree(self):
        pairs = (
            (0, 29), (1, 28),
            (4, 12), (5, 11),
            (16, 24), (17, 23),
        )
        svg, layout = render_structure_svg("A" * 30, pairs)
        self.assertEqual(layout, "tree")
        self.assertEqual(svg.count('class="base-pair '), 6)
        self.assertEqual(svg.count('class="base loop"'), 18)

    def test_bulge_bases_are_labelled(self):
        svg, layout = render_structure_svg("ACGTACGTAA", ((0, 9), (2, 7)))
        self.assertEqual(layout, "stem-loop")
        self.assertEqual(svg.count('class="base bulge"'), 2)

    def test_no_pairs_returns_empty(self):
        self.assertEqual(render_structure_svg("ACGT", ()), ("", "none"))

    def test_invalid_pairs_raise(self):
        with self.assertRaises(ValueError):
            render_structure_svg("ACGT", ((0, 9),))
        with self.assertRaises(ValueError):
            render_structure_svg("ACGT", ((0, 3), (1, 2), (0, 1)))
        with self.assertRaises(ValueError):
            render_structure_svg("ACGT", ((0, 2), (1, 3)))

    def test_structure_legend_first_appearance(self):
        items = structure_legend_items("ACGTAC")
        self.assertEqual(
            [item["label"] for item in items[:4]], ["A", "C", "G", "T"]
        )
        self.assertEqual(items[0]["kind"], "base")
        self.assertEqual(items[0]["color"], "#009E73")
        self.assertEqual(
            [item["style"] for item in pair_legend_items()],
            ["wc", "wobble", "noncanonical", "unknown"],
        )


class LegendItemsTest(unittest.TestCase):
    def test_legend_for_mixed_sequence(self):
        items = legend_items(parse("mA*+C/3Bio/"))
        bases = [item["label"] for item in items if item["kind"] == "base"]
        sugars = [item["label"] for item in items if item["kind"] == "sugar"]
        links = [item["label"] for item in items if item["kind"] == "link"]
        tags = [item["label"] for item in items if item["kind"] == "tag"]

        self.assertEqual(bases, ["A", "C"])
        self.assertEqual(sugars, ["2'-OMe", "LNA"])
        self.assertEqual(links, ["phosphate", "phosphorothioate"])
        self.assertEqual(tags, ["/3Bio/"])
        for item in items:
            if item["kind"] == "sugar":
                self.assertTrue(item["icon"])

    def test_legend_for_fluoro_and_moe(self):
        items = legend_items(parse("/i2FA/i2MOErC/"))
        sugars = [item for item in items if item["kind"] == "sugar"]
        self.assertEqual([item["label"] for item in sugars], ["2'-F", "2'-MOE"])
        self.assertIn(">F<", sugars[0]["icon"])
        self.assertIn(">MOE<", sugars[1]["icon"])

    def test_base_entry_carries_tile_colours(self):
        entry = legend_items(parse("A"))[0]
        self.assertEqual(entry["kind"], "base")
        self.assertEqual(entry["label"], "A")
        self.assertEqual(entry["color"], "#009E73")
        self.assertEqual(entry["text"], "#FFFFFF")

    def test_degenerate_label_lists_alternatives(self):
        items = legend_items(parse("ACGTN"))
        self.assertIn("N (A/C/G/T)", [item["label"] for item in items])

    def test_empty_sequence_has_empty_legend(self):
        self.assertEqual(legend_items(parse("")), [])


if __name__ == "__main__":
    unittest.main()
