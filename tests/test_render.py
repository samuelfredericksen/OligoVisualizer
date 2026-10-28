"""Tests for render.py (stdlib unittest; importable with bare `py`)."""

import unittest

from notation import parse
from render import legend_items, render_svg


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
        self.assertEqual(svg.count("<circle"), 1)
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
