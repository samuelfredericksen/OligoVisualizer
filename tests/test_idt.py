"""Tests for idt.py; every network method is replaced with a fake."""

import base64
import json
import unittest

from idt import (
    Fold,
    FoldConditions,
    Homodimer,
    IdtAuthError,
    IdtError,
    IdtSession,
    IdtSessionExpired,
    NoStructureError,
    OligoProperties,
    clean_sequence,
    nucleotide_type,
    parse_ct,
    select_homodimer,
    select_properties,
    select_structure,
)
from notation import parse


def ct_text(bases: str, pairs: dict[int, int]) -> str:
    partners = dict(pairs)
    for left, right in pairs.items():
        partners[right] = left
    lines = [f"  {len(bases)}"]
    for index, base in enumerate(bases, start=1):
        partner = partners.get(index - 1, -1)
        lines.append(
            f"{index} {base} 0 {index + 1} "
            f"{partner + 1 if partner >= 0 else 0} {index}"
        )
    return "\n".join(lines) + "\n"


class CleanSequenceTest(unittest.TestCase):
    def test_strips_quotes_and_keeps_modification_codes(self):
        self.assertEqual(
            clean_sequence("5'-/5Biosg/acgt/3Bio/-3'"),
            "/5Biosg/ACGT/3Bio/",
        )

    def test_keeps_ps_sugar_and_custom_mix_notation(self):
        self.assertEqual(
            clean_sequence("a*c+g*mr(10203040)t"),
            "A*C+G*mr(10203040)T",
        )

    def test_keeps_modified_base_tokens_verbatim(self):
        self.assertEqual(
            clean_sequence("/i2FA/*/i2FC/"),
            "/i2FA/*/i2FC/",
        )

    def test_separates_adjacent_single_slash_codes(self):
        self.assertEqual(
            clean_sequence("aaa/i2FG/i2MOErG/ggg"),
            "AAA/i2FG//i2MOErG/GGG",
        )
        self.assertEqual(
            clean_sequence("/5Phos/iSp18/ttt"),
            "/5Phos//iSp18/TTT",
        )

    def test_reported_adjacent_fluoro_moe_sequence_is_not_mangled(self):
        raw = (
            "CCmC CmCC C*C*C* CTrU TrUT G*G*G* G*G*mG GGG GAC TGA CTG "
            "ACT GTT TTG G*G*A AAA/i2FG/i2MOErG/GG"
        )
        cleaned = clean_sequence(raw)
        self.assertIn("/i2FG//i2MOErG/", cleaned)
        self.assertNotIn("/i2FG/I2MOErG/", cleaned)

    def test_whitespace_is_removed(self):
        self.assertEqual(clean_sequence(" AC \n GT "), "ACGT")


class NucleotideTypeTest(unittest.TestCase):
    def test_dna(self):
        self.assertEqual(nucleotide_type(parse("ACGT")), "DNA")

    def test_rna_ribose(self):
        self.assertEqual(nucleotide_type(parse("rArCrGrU")), "RNA")

    def test_rna_uridine(self):
        self.assertEqual(nucleotide_type(parse("ACGU")), "RNA")


class ParseCtTest(unittest.TestCase):
    def test_twelve_base_four_pair_map(self):
        text = ct_text(
            "GCGCAAAAGCGC",
            {0: 11, 1: 10, 2: 9, 3: 8},
        )
        bases, pairs = parse_ct(text)
        self.assertEqual(bases, "GCGCAAAAGCGC")
        self.assertEqual(pairs, ((0, 11), (1, 10), (2, 9), (3, 8)))

    def test_real_unafold_header_line_is_accepted(self):
        text = (
            "  4\tdG = -1.234\tUNAFOLDRDM_example\n"
            "1 A 0 2 4 1 0 2\n"
            "2 C 1 3 3 2 1 3\n"
            "3 G 2 4 2 3 2 4\n"
            "4 T 3 5 1 4 3 5\n"
        )
        self.assertEqual(parse_ct(text), ("ACGT", ((0, 3), (1, 2))))

    def test_truncated_map_raises(self):
        text = "  4\n1 A 0 2 4 1\n2 C 0 3 3 2\n3 G 0 2 2 3\n"
        with self.assertRaises(IdtError):
            parse_ct(text)

    def test_duplicate_index_raises(self):
        text = "  2\n1 A 0 2 2 1\n1 C 0 2 1 2\n"
        with self.assertRaises(IdtError):
            parse_ct(text)

    def test_non_reciprocal_pair_raises(self):
        text = "  4\n1 A 0 2 4 1\n2 C 0 3 3 2\n3 G 0 2 2 3\n4 T 0 1 2 4\n"
        with self.assertRaises(IdtError):
            parse_ct(text)


def hairpin_payload() -> dict:
    return {
        "HasErrors": False,
        "HasModelErrors": False,
        "ValidationError": None,
        "NucleotideType": "DNA",
        "Temperature": 25,
        "SodiumConc": 50,
        "MagConc": 3,
        "PercentSuboptimal": 50,
        "MaxFoldings": 20,
        "OutputObj": {
            "Structures": [
                {
                    "StructureName": "1",
                    "GNode": -8.96,
                    "TMNode": 68.1,
                    "HNode": -70.0,
                    "SNode": -210.0,
                    "ConnectivityURL": "/AnalyzerService/x/UNAFOLDRDM_A/one.ct",
                },
                {
                    "StructureName": "2",
                    "GNode": -9.47,
                    "TMNode": 68.5,
                    "HNode": -74.3,
                    "SNode": -217.45,
                    "ConnectivityURL": "/AnalyzerService/x/UNAFOLDRDM_A/two.ct",
                },
            ],
        },
    }


class SelectStructureTest(unittest.TestCase):
    def test_chooses_lowest_delta_g(self):
        chosen, count, meta = select_structure(hairpin_payload(), "DNA")
        self.assertEqual(chosen.delta_g, -9.47)
        self.assertEqual(chosen.name, "2")
        self.assertEqual(count, 2)
        self.assertEqual(meta["temperature_c"], 25)
        self.assertEqual(meta["sodium_mm"], 50)
        self.assertEqual(meta["magnesium_mm"], 3)
        self.assertEqual(meta["nucleotide_type"], "DNA")

    def test_has_errors_carries_position_and_symbol(self):
        data = {
            "HasErrors": True,
            "Errors": {
                "Position": 3,
                "ErrorSymbol": "!",
                "ErrorText": "Invalid bases were detected.",
            },
        }
        with self.assertRaises(IdtError) as caught:
            select_structure(data, "DNA")
        self.assertEqual(caught.exception.position, 3)
        self.assertEqual(caught.exception.token, "!")

    def test_no_structure_error_is_a_no_structure_result(self):
        data = {
            "HasErrors": True,
            "Errors": {
                "Position": 0,
                "ErrorSymbol": "",
                "ErrorText": "No structure was found for this sequence!",
            },
        }
        with self.assertRaises(NoStructureError):
            select_structure(data, "DNA")

    def test_validation_error_raises(self):
        with self.assertRaisesRegex(IdtError, "rejected"):
            select_structure({"ValidationError": "<b>rejected</b>"}, "DNA")

    def test_validation_no_structure_is_a_no_structure_result(self):
        with self.assertRaises(NoStructureError):
            select_structure(
                {"ValidationError": "No structure was found for this sequence!"},
                "DNA",
            )

    def test_missing_output_raises_no_structure(self):
        with self.assertRaises(NoStructureError):
            select_structure({}, "DNA")

    def test_single_structure_dict_is_accepted(self):
        data = {
            "OutputObj": {
                "Structures": {
                    "StructureName": "solo",
                    "GNode": -1.5,
                    "ConnectivityURL": "/x.ct",
                }
            }
        }
        chosen, count, _meta = select_structure(data, "DNA")
        self.assertEqual(chosen.name, "solo")
        self.assertEqual(count, 1)

    def test_null_g_node_is_skipped(self):
        data = {
            "OutputObj": {
                "Structures": [
                    {"GNode": None},
                    {"GNode": -2.0, "StructureName": "ok"},
                ]
            }
        }
        chosen, count, _meta = select_structure(data, "DNA")
        self.assertEqual(chosen.name, "ok")
        self.assertEqual(count, 1)

    def test_model_error_raises(self):
        with self.assertRaisesRegex(IdtError, "model"):
            select_structure(
                {"HasModelErrors": True, "ModelErrors": ["model failed"]},
                "DNA",
            )

    def test_conditions_fill_missing_idt_metadata(self):
        conditions = FoldConditions(37, 100, 1, 25, 10)
        data = {
            "OutputObj": {
                "Structures": [
                    {"StructureName": "solo", "GNode": -1.5},
                ]
            }
        }
        _chosen, _count, meta = select_structure(data, "DNA", conditions)
        self.assertEqual(meta["temperature_c"], 37)
        self.assertEqual(meta["sodium_mm"], 100)
        self.assertEqual(meta["magnesium_mm"], 1)
        self.assertEqual(meta["suboptimality_percent"], 25)
        self.assertEqual(meta["max_foldings"], 10)

    def test_echoed_idt_metadata_wins_over_conditions(self):
        data = hairpin_payload()
        data.update(
            Temperature=30,
            SodiumConc=80,
            MagConc=2,
            PercentSuboptimal=10,
            MaxFoldings=5,
        )
        _chosen, _count, meta = select_structure(
            data, "DNA", FoldConditions(37, 100, 1, 25, 10)
        )
        self.assertEqual(meta["temperature_c"], 30)
        self.assertEqual(meta["sodium_mm"], 80)
        self.assertEqual(meta["magnesium_mm"], 2)
        self.assertEqual(meta["suboptimality_percent"], 10)
        self.assertEqual(meta["max_foldings"], 5)


class SignInTest(unittest.TestCase):
    def test_success_stores_tokens_but_not_password(self):
        session = IdtSession()
        app_page = (
            '<html><input type="hidden" name="__RequestVerificationToken" '
            'value="verify">'
            '<input type="hidden" name="__RequestVerificationStatelessToken" '
            'value="stateless"></html>'
        )
        captured = {}

        def fake_get(url, accept="text/html"):
            if "returnurl" in url:
                return b"<html>login</html>"
            return app_page.encode()

        def fake_post_form(url, fields, headers=None):
            captured.update(fields)
            return (
                b'{"userAuthenticated":true,"url":"/calc/analyzer",'
                b'"accessToken":"access","accessAccount":"account",'
                b'"JsError":""}'
            )

        session._get = fake_get
        session._post_form = fake_post_form
        session.sign_in("me@example.com", "secret")

        self.assertTrue(session.signed_in)
        self.assertEqual(session._token, "verify")
        self.assertEqual(session._stateless_token, "stateless")
        self.assertEqual(session._access_token, "access")
        self.assertEqual(session._access_account, "account")
        self.assertNotIn("secret", repr(session.__dict__))
        self.assertEqual(captured["Password"], "secret")

    def test_invalid_credentials_raise_and_clear(self):
        session = IdtSession()

        def fake_post_form(url, fields, headers=None):
            return (
                b'{"userAuthenticated":false,'
                b'"JsError":"<div>Invalid Credentials. Try again.</div>"}'
            )

        session._get = lambda url, accept="text/html": b"<html></html>"
        session._post_form = fake_post_form
        with self.assertRaisesRegex(IdtAuthError, "Invalid Credentials"):
            session.sign_in("me@example.com", "wrong")
        self.assertFalse(session.signed_in)

    def test_missing_token_raises_and_clears(self):
        session = IdtSession()

        def fake_post_form(url, fields, headers=None):
            return b'{"userAuthenticated":true,"url":"/calc/analyzer"}'

        session._get = lambda url, accept="text/html": b"<html></html>"
        session._post_form = fake_post_form
        with self.assertRaises(IdtAuthError):
            session.sign_in("me@example.com", "secret")
        self.assertFalse(session.signed_in)


def properties_payload() -> dict:
    return {
        "MeltTemp": 73.2,
        "MolecularWeight": 8649.6,
        "ExtCoefficient": 269700.0,
        "NmoleOD": 3.71,
        "UgOD": 32.07,
        "HasErrors": False,
        "HasModelErrors": False,
    }


class SelectPropertiesTest(unittest.TestCase):
    def test_returns_all_five_analyzer_values(self):
        properties = select_properties(properties_payload())
        self.assertIsInstance(properties, OligoProperties)
        self.assertEqual(properties.melt_temp_c, 73.2)
        self.assertEqual(properties.molecular_weight_g_mol, 8649.6)
        self.assertEqual(
            properties.extinction_coefficient_l_mol_cm, 269700.0
        )
        self.assertEqual(properties.nmol_per_od260, 3.71)
        self.assertEqual(properties.ug_per_od260, 32.07)

    def test_missing_property_is_an_error(self):
        data = properties_payload()
        data["UgOD"] = None
        with self.assertRaisesRegex(IdtError, "incomplete"):
            select_properties(data)

    def test_model_error_is_reported(self):
        data = {
            "HasErrors": False,
            "HasModelErrors": True,
            "ModelErrors": ["Oligo Sequence is required."],
        }
        with self.assertRaisesRegex(IdtError, "Oligo Sequence"):
            select_properties(data)


def homodimer_payload() -> dict:
    return {
        "PrimarySequence": "ACGTACGT",
        "MaxDeltaG": -20.0,
        "Results": [
            {
                "DeltaG": -4.0,
                "BasePairs": 3,
                "Dimer": "5' ACGTACGT\n   |||     \n3' TGCATGCA",
            },
            {
                "DeltaG": -8.0,
                "BasePairs": 5,
                "Dimer": "5' ACGTACGT\n   ||||:   \n3' TGCATGCA",
            },
        ],
    }


class SelectHomodimerTest(unittest.TestCase):
    def test_chooses_lowest_delta_g(self):
        chosen = select_homodimer(homodimer_payload())
        self.assertEqual(chosen.delta_g, -8.0)
        self.assertEqual(chosen.base_pairs, 5)
        self.assertEqual(chosen.additional_pairs, 1)
        self.assertEqual(chosen.result_count, 2)
        self.assertEqual(chosen.sequence, "ACGTACGT")
        self.assertEqual(chosen.max_delta_g, -20.0)

    def test_no_results_is_no_structure(self):
        with self.assertRaises(NoStructureError):
            select_homodimer({"PrimarySequence": "ACGT", "Results": []})

    def test_errors_are_reported(self):
        data = {
            "HasErrors": True,
            "Errors": {
                "Position": 2,
                "ErrorSymbol": "!",
                "ErrorText": "Invalid bases were detected.",
            },
        }
        with self.assertRaises(IdtError) as caught:
            select_homodimer(data)
        self.assertEqual(caught.exception.position, 2)
        self.assertEqual(caught.exception.token, "!")


class PredictTest(unittest.TestCase):
    def test_payload_headers_and_chosen_fold(self):
        session = IdtSession()
        session._signed_in = True
        session._token = "verify"
        session._stateless_token = "stateless"
        session._access_token = "access"
        session._access_account = "account"
        session._base = "https://www.idtdna.com"
        captured = {}

        def fake_post_json(url, payload, headers=None):
            captured["url"] = url
            captured["payload"] = payload
            captured["headers"] = headers
            return hairpin_payload()

        def fake_post_raw(url, body, headers=None):
            captured["raw_url"] = url
            captured["raw_body"] = body
            captured["raw_headers"] = headers
            return ct_text(
                "GCGCAAAAGCGC",
                {0: 11, 1: 10, 2: 9, 3: 8},
            ).encode()

        session._post_json = fake_post_json
        session._post_raw = fake_post_raw
        fold = session.predict("/5Biosg/GCGCAAAAGCGC/3Bio/", "DNA")

        self.assertIsInstance(fold, Fold)
        self.assertEqual(fold.delta_g, -9.47)
        self.assertEqual(fold.pairs, ((0, 11), (1, 10), (2, 9), (3, 8)))
        self.assertEqual(fold.structure_count, 2)
        self.assertEqual(captured["url"], "https://www.idtdna.com/calc/analyzer/home/hairpin")
        self.assertEqual(
            captured["payload"]["settings"]["Sequence"],
            "/5Biosg/GCGCAAAAGCGC/3Bio/",
        )
        self.assertEqual(captured["payload"]["settings"]["NaConc"], 50)
        self.assertEqual(captured["payload"]["settings"]["MgConc"], 3)
        self.assertEqual(captured["payload"]["settings"]["NucleotideType"], "DNA")
        self.assertEqual(captured["payload"]["hairpinSettings"]["Temp"], 25)
        self.assertEqual(
            captured["payload"]["hairpinSettings"]["Suboptimality"], 50
        )
        self.assertEqual(
            captured["payload"]["hairpinSettings"]["MaxFoldings"], 20
        )
        self.assertEqual(captured["headers"]["__RequestVerificationToken"], "verify")
        self.assertEqual(
            session._auth_headers(),
            {"Authorization": "Bearer access", "AccessAccount": "account"},
        )
        self.assertIn(b'"type":"hairpin"', captured["raw_body"])
        self.assertEqual(
            captured["raw_headers"]["__RequestVerificationToken"], "stateless"
        )

    def test_properties_payload_and_values(self):
        session = IdtSession()
        session._signed_in = True
        session._base = "https://www.idtdna.com"
        captured = {}

        def fake_post_json(url, payload, headers=None):
            captured["url"] = url
            captured["payload"] = payload
            return properties_payload()

        session._post_json = fake_post_json
        conditions = FoldConditions(37, 100, 1, 25, 10, 0.3, 0.0)
        properties = session.analyze_properties("ACGTACGT", "DNA", conditions)

        self.assertEqual(properties.melt_temp_c, 73.2)
        self.assertEqual(
            captured["url"],
            "https://www.idtdna.com/api/OligoAnalyzer/Analyze",
        )
        self.assertEqual(captured["payload"]["Sequence"], "ACGTACGT")
        self.assertEqual(captured["payload"]["NaConc"], 100)
        self.assertEqual(captured["payload"]["MgConc"], 1)
        self.assertEqual(captured["payload"]["OligoConc"], 0.3)
        self.assertEqual(captured["payload"]["DNTPsConc"], 0.0)
        self.assertEqual(captured["payload"]["NucleotideType"], "DNA")

    def test_homodimer_payload_and_selection(self):
        session = IdtSession()
        session._signed_in = True
        session._base = "https://www.idtdna.com"
        captured = {}

        def fake_post_json(url, payload, headers=None):
            captured["url"] = url
            captured["payload"] = payload
            captured["headers"] = headers
            return homodimer_payload()

        session._post_json = fake_post_json
        conditions = FoldConditions(37, 100, 1, 25, 10, 0.25, 0.0)
        dimer = session.predict_homodimer("ACGTACGT", "DNA", conditions)

        self.assertIsInstance(dimer, Homodimer)
        self.assertEqual(dimer.delta_g, -8.0)
        self.assertEqual(
            captured["url"],
            "https://www.idtdna.com/calc/analyzer/home/selfdimer",
        )
        settings = captured["payload"]["settings"]
        self.assertEqual(settings["Sequence"], "ACGTACGT")
        self.assertEqual(settings["NaConc"], 100)
        self.assertEqual(settings["MgConc"], 1)
        self.assertEqual(settings["OligoConc"], 0.25)
        self.assertEqual(settings["DNTPsConc"], 0.0)
        self.assertEqual(settings["NucleotideType"], "DNA")

    def test_custom_conditions_reach_payload_and_fold(self):
        session = IdtSession()
        session._signed_in = True
        session._token = "verify"
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        data = hairpin_payload()
        for key in (
            "Temperature",
            "SodiumConc",
            "MagConc",
            "PercentSuboptimal",
            "MaxFoldings",
        ):
            data.pop(key)
        captured = {}

        def fake_post_json(url, payload, headers=None):
            captured["payload"] = payload
            return data

        def fake_post_raw(url, body, headers=None):
            return ct_text(
                "GCGCAAAAGCGC",
                {0: 11, 1: 10, 2: 9, 3: 8},
            ).encode()

        session._post_json = fake_post_json
        session._post_raw = fake_post_raw
        conditions = FoldConditions(37, 100, 1, 25, 10)
        fold = session.predict("GCGCAAAAGCGC", "DNA", conditions)

        self.assertEqual(captured["payload"]["settings"]["NaConc"], 100)
        self.assertEqual(captured["payload"]["settings"]["MgConc"], 1)
        self.assertEqual(captured["payload"]["hairpinSettings"]["Temp"], 37)
        self.assertEqual(
            captured["payload"]["hairpinSettings"]["Suboptimality"], 25
        )
        self.assertEqual(
            captured["payload"]["hairpinSettings"]["MaxFoldings"], 10
        )
        self.assertEqual(fold.temperature_c, 37)
        self.assertEqual(fold.sodium_mm, 100)
        self.assertEqual(fold.magnesium_mm, 1)
        self.assertEqual(fold.suboptimality_percent, 25)
        self.assertEqual(fold.max_foldings, 10)

    def test_auth_headers_reach_json_and_raw_requests(self):
        session = IdtSession()
        session._access_token = "access"
        session._access_account = "account"
        captured = []

        def fake_request(req):
            captured.append(req)
            if req.full_url.endswith("/json"):
                return b'{"ok": true}', req.full_url, "application/json"
            return b"raw", req.full_url, "text/plain"

        session._request = fake_request
        session._post_json("https://www.idtdna.com/json", {"x": 1})
        session._post_raw("https://www.idtdna.com/raw", b"x")
        self.assertEqual(
            [req.get_header("Authorization") for req in captured],
            ["Bearer access", "Bearer access"],
        )
        self.assertEqual(
            [req.get_header("Accessaccount") for req in captured],
            ["account", "account"],
        )

    def test_not_signed_in_raises(self):
        session = IdtSession()
        with self.assertRaises(IdtAuthError):
            session.predict("ACGT", "DNA")

    def test_no_pairs_raises_no_structure(self):
        session = IdtSession()
        session._signed_in = True
        session._token = "verify"
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        session._post_json = lambda url, payload, headers=None: hairpin_payload()
        session._post_raw = lambda url, body, headers=None: ct_text(
            "ACGT", {}
        ).encode()
        with self.assertRaises(NoStructureError):
            session.predict("ACGT", "DNA")


class RetrieveCtTest(unittest.TestCase):
    def test_public_fallback_accepts_modern_guid_path(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})
        captured = {}

        def fake_form(url, fields, headers=None):
            captured.update(fields)
            encoded = base64.b64encode(expected.encode()).decode()
            return (
                '<Envelope xmlns="http://services.idtdna.com/">'
                f"<base64Binary>{encoded}</base64Binary></Envelope>"
            ).encode()

        session._post_raw = lambda url, body, headers=None: (_ for _ in ()).throw(
            IdtError("dataout failed")
        )
        session._post_form = fake_form
        location = (
            "62106641-8c36-4e3c-8401-b886460570c1/"
            "62106641-8c36-4e3c-8401-b886460570c1_1.ct"
        )
        self.assertEqual(session._retrieve_ct(location), expected)
        self.assertEqual(captured["location"], location)

    def test_public_fallback_decodes_base64_xml(self):
        session = IdtSession()
        session._signed_in = True
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})
        captured = {}

        def fake_raw(url, body, headers=None):
            raise IdtError("dataout failed")

        def fake_form(url, fields, headers=None):
            captured.update(fields)
            encoded = base64.b64encode(expected.encode()).decode()
            return (
                '<Envelope xmlns="http://services.idtdna.com/">'
                f"<base64Binary>{encoded}</base64Binary></Envelope>"
            ).encode()

        session._post_raw = fake_raw
        session._post_form = fake_form
        text = session._retrieve_ct(
            "/AnalyzerService/x/UNAFOLDRDM_A/one.ct"
        )
        self.assertEqual(text, expected)
        self.assertEqual(captured["location"], "UNAFOLDRDM_A/one.ct")

    def test_raw_pre_wrapped_ct_is_unwrapped(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})
        session._post_raw = lambda url, body, headers=None: (
            "<pre>" + expected + "</pre>"
        ).encode()
        session._post_form = lambda url, fields, headers=None: self.fail(
            "public fallback must not run"
        )
        self.assertEqual(
            session._retrieve_ct("/AnalyzerService/x/one.ct"),
            expected,
        )

    def test_json_wrapped_html_ct_is_unwrapped(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})
        wrapped = (
            "<pre style='white-space: pre-wrap;'>" + expected + "</pre>"
        )
        session._post_raw = lambda url, body, headers=None: json.dumps(
            wrapped
        ).encode()
        session._post_form = lambda url, fields, headers=None: self.fail(
            "public fallback must not run"
        )
        self.assertEqual(
            session._retrieve_ct("/AnalyzerService/x/one.ct"),
            expected,
        )

    def test_base64_dataout_is_decoded_without_fallback(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})
        session._post_raw = lambda url, body, headers=None: base64.b64encode(
            expected.encode()
        )
        session._post_form = lambda url, fields, headers=None: self.fail(
            "public fallback must not run"
        )
        self.assertEqual(
            session._retrieve_ct("/AnalyzerService/x/UNAFOLDRDM_A/one.ct"),
            expected,
        )

    def test_malformed_dataout_falls_back_to_public_ct(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        expected = ct_text("ACGT", {0: 3, 1: 2})

        session._post_raw = lambda url, body, headers=None: b"not a ct map"
        session._post_form = lambda url, fields, headers=None: (
            '<Envelope xmlns="http://services.idtdna.com/">'
            f"<base64Binary>{base64.b64encode(expected.encode()).decode()}"
            "</base64Binary></Envelope>"
        ).encode()
        self.assertEqual(
            session._retrieve_ct("/AnalyzerService/x/UNAFOLDRDM_A/one.ct"),
            expected,
        )

    def test_html_dataout_is_session_expired_and_does_not_fall_back(self):
        session = IdtSession()
        session._stateless_token = "stateless"
        session._base = "https://www.idtdna.com"
        session._post_raw = lambda url, body, headers=None: b"<html>login</html>"
        session._post_form = lambda url, fields, headers=None: self.fail(
            "public fallback must not run"
        )
        with self.assertRaises(IdtSessionExpired):
            session._retrieve_ct("/AnalyzerService/x/UNAFOLDRDM_A/one.ct")


if __name__ == "__main__":
    unittest.main()
