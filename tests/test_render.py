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

    def test_sugar_superscripts(self):
        svg = render_svg(parse("mA+CrG"))
        self.assertIn(">m<", svg)
        self.assertIn(">+<", svg)
        self.assertIn(">r<", svg)

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
        self.assertEqual(sugars, ["m", "+"])
        self.assertEqual(links, ["phosphate", "phosphorothioate"])
        self.assertEqual(tags, ["/3Bio/"])

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
