import unittest

from notation import (
    Modification,
    parse,
    normalize,
    summary_line,
)


class TestNotationTable(unittest.TestCase):
    def test_plain_bases(self):
        seq = parse("ACGT")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertTrue(all(nt.sugar == "d" for nt in seq.nucleotides))
        self.assertEqual(seq.errors, [])
        self.assertEqual(seq.warnings, [])

    def test_lowercase_bases(self):
        seq = parse("acgt")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertEqual([nt.display for nt in seq.nucleotides], ["A", "C", "G", "T"])

    def test_rna(self):
        seq = parse("rArCrGrU")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertTrue(all(nt.sugar == "r" for nt in seq.nucleotides))

    def test_rT_is_rejected_but_drawn_as_u(self):
        seq = parse("rT")
        self.assertTrue(any("rT is not valid" in p.message for p in seq.errors))
        self.assertEqual(len(seq.nucleotides), 1)
        self.assertEqual(seq.nucleotides[0].display, "U")

    def test_2_o_methyl_with_ps(self):
        seq = parse("mA*")
        self.assertEqual(len(seq.nucleotides), 1)
        self.assertEqual(seq.nucleotides[0].sugar, "m")
        self.assertIs(seq.nucleotides[0].ps_after, True)

    def test_lna_with_ps(self):
        seq = parse("+A*+C")
        self.assertEqual(len(seq.nucleotides), 2)
        self.assertTrue(all(nt.sugar == "LNA" for nt in seq.nucleotides))
        self.assertIs(seq.nucleotides[0].ps_after, True)
        self.assertIs(seq.nucleotides[1].ps_after, False)

    def test_ps_suffix_mid_chain(self):
        seq = parse("TTTTTT*GGGGGG")
        self.assertEqual(len(seq.nucleotides), 12)
        ps_positions = [nt.position for nt in seq.nucleotides if nt.ps_after]
        self.assertEqual(ps_positions, [5])

    def test_terminal_modifications(self):
        seq = parse("/5Biosg/TTTT/3Bio/")
        self.assertEqual([m.raw for m in seq.five_prime], ["/5Biosg/"])
        self.assertEqual([m.raw for m in seq.three_prime], ["/3Bio/"])
        self.assertEqual(len(seq.nucleotides), 4)

    def test_two_five_prime_mods_keep_order(self):
        seq = parse("/5Phos/iSp18/ACGT")
        self.assertEqual([m.raw for m in seq.five_prime], ["/5Phos/", "/iSp18/"])

    def test_internal_modification(self):
        seq = parse("AC/5Me-dC/GT")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertEqual(
            seq.internal, {1: [Modification("5Me-dC", "/5Me-dC/")]}
        )

    def test_modification_position_comes_from_placement(self):
        seq = parse("/ZEN/ACGT")
        self.assertEqual([m.code for m in seq.five_prime], ["ZEN"])

    def test_degenerate_n(self):
        seq = parse("ACGTN")
        last = seq.nucleotides[-1]
        self.assertEqual(last.display, "N")
        self.assertEqual(last.alternatives, ("A", "C", "G", "T"))

    def test_degenerate_r(self):
        seq = parse("R")
        self.assertEqual(seq.nucleotides[0].alternatives, ("A", "G"))

    def test_custom_mix(self):
        seq = parse("A(10203040)T")
        self.assertEqual(len(seq.nucleotides), 3)
        self.assertEqual(seq.nucleotides[1].display, "?")
        self.assertEqual(seq.mixes[1].percentages, (10, 20, 30, 40))

    def test_custom_mix_must_sum_to_100(self):
        seq = parse("A(10903040)T")
        self.assertTrue(
            any(
                p.message == "custom mix percentages must sum to 100"
                for p in seq.errors
            )
        )

    def test_unexpected_character(self):
        seq = parse("ACGT!")
        matches = [p for p in seq.errors if p.message == "unexpected character '!'"]
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].position, 4)
        self.assertEqual(matches[0].token, "!")

    def test_unterminated_modification(self):
        seq = parse("/Unclosed")
        self.assertTrue(
            any(p.message == "unterminated modification code" for p in seq.errors)
        )

    def test_quoted_and_dashed_sequence(self):
        seq = parse("5'-ACGT-3'")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertEqual(seq.errors, [])
        self.assertEqual(seq.warnings, [])

    def test_whitespace_is_stripped(self):
        seq = parse(" AC \n GT ")
        self.assertEqual(len(seq.nucleotides), 4)
        self.assertEqual(normalize(" AC \n GT "), "ACGT")

    def test_empty_string(self):
        seq = parse("")
        self.assertEqual(len(seq.nucleotides), 0)
        self.assertEqual(seq.errors, [])
        self.assertEqual(seq.warnings, [])

    def test_two_ps_links_single_star_each(self):
        seq = parse("mA*mA*")
        self.assertEqual(len(seq.nucleotides), 2)
        self.assertTrue(all(nt.ps_after for nt in seq.nucleotides))
        self.assertEqual(seq.nucleotides[0].raw, "mA*")
        self.assertEqual(seq.nucleotides[1].raw, "mA*")


class TestSummaryLine(unittest.TestCase):
    def test_one_ps_no_mods(self):
        self.assertEqual(
            summary_line(parse("TTTTTT*GGGGGG")),
            "12 nt · 1 PS linkage · 0 modifications",
        )

    def test_two_ps_no_mods(self):
        self.assertEqual(
            summary_line(parse("mA*mA*")),
            "2 nt · 2 PS linkages · 0 modifications",
        )

    def test_no_ps_two_mods(self):
        self.assertEqual(
            summary_line(parse("/5Biosg/ACGT/3Bio/")),
            "4 nt · 0 PS linkages · 2 modifications",
        )


if __name__ == "__main__":
    unittest.main()
