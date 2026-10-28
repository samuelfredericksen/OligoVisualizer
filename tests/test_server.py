"""Tests for the Flask layer, using the test client (no real port)."""

import unittest

import app
from idt import (
    Fold,
    FoldConditions,
    Homodimer,
    IdtAuthError,
    IdtError,
    IdtSessionExpired,
    NoStructureError,
    OligoProperties,
)


class FakeSession:
    def __init__(self, signed_in=True):
        self._signed_in = signed_in
        self.predict_calls = 0
        self.predict_args = []
        self.sign_in_calls = []
        self.sign_out_calls = 0
        self.error = None
        self.login_error = None
        self.fold = make_fold()
        self.homodimer = None
        self.homodimer_error = None
        self.homodimer_calls = 0
        self.homodimer_args = []
        self.properties = make_properties()
        self.properties_error = None
        self.properties_calls = 0
        self.properties_args = []

    @property
    def signed_in(self):
        return self._signed_in

    def sign_in(self, username, password):
        self.sign_in_calls.append((username, password))
        if self.login_error is not None:
            raise self.login_error
        self._signed_in = True

    def sign_out(self):
        self.sign_out_calls += 1
        self._signed_in = False

    def predict(self, sequence, nucleotide_type, conditions):
        self.predict_calls += 1
        self.predict_args.append((sequence, nucleotide_type, conditions))
        if self.error is not None:
            raise self.error
        return self.fold

    def predict_homodimer(self, sequence, nucleotide_type, conditions):
        self.homodimer_calls += 1
        self.homodimer_args.append((sequence, nucleotide_type, conditions))
        if self.homodimer_error is not None:
            raise self.homodimer_error
        return self.homodimer or make_homodimer(self.fold.bases)

    def analyze_properties(self, sequence, nucleotide_type, conditions):
        self.properties_calls += 1
        self.properties_args.append((sequence, nucleotide_type, conditions))
        if self.properties_error is not None:
            raise self.properties_error
        return self.properties


def make_properties():
    return OligoProperties(
        melt_temp_c=73.2,
        molecular_weight_g_mol=8649.6,
        extinction_coefficient_l_mol_cm=269700.0,
        nmol_per_od260=3.71,
        ug_per_od260=32.07,
    )


def make_homodimer(bases):
    top = bases
    bottom = bases[::-1]
    dimer = "5' " + top + "\n   " + "|" * len(top) + "\n3' " + bottom
    return Homodimer(
        sequence=bases,
        delta_g=-12.0,
        base_pairs=len(bases),
        additional_pairs=0,
        max_delta_g=-20.0,
        dimer=dimer,
        result_count=1,
    )


def make_fold():
    return Fold(
        bases="ACGTACGTACGTACGTACGT",
        pairs=tuple((index, 19 - index) for index in range(8)),
        structure_name="2",
        delta_g=-9.467,
        tm_c=68.5,
        delta_h=-74.30,
        delta_s=-217.45,
        structure_count=2,
        nucleotide_type="DNA",
        temperature_c=25,
        sodium_mm=50,
        magnesium_mm=3,
        suboptimality_percent=50,
        max_foldings=20,
    )


class ServerTests(unittest.TestCase):
    def setUp(self):
        self.client = app.app.test_client()
        self.saved_session = app._session
        self.saved_cache = app._structure_cache
        self.fake = FakeSession()
        app._session = self.fake
        app._structure_cache = {}

    def tearDown(self):
        app._session = self.saved_session
        app._structure_cache = self.saved_cache

    def test_index_serves_page(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        body = response.get_data(as_text=True)
        response.close()
        self.assertIn('id="seq"', body)
        self.assertIn('id="idtform"', body)
        self.assertIn('id="predict"', body)
        self.assertIn('id="homodimerpanel"', body)
        self.assertIn('id="cond-oligo"', body)
        self.assertIn('id="cond-dntps"', body)
        self.assertIn('id="properties"', body)
        self.assertIn('id="prop-tm"', body)

    def test_empty_sequence(self):
        response = self.client.get("/api/render?seq=")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], False)
        self.assertEqual(payload["errors"][0]["message"], "no sequence supplied")

    def test_renders_sequence(self):
        response = self.client.get("/api/render?seq=/5Biosg/TTTTTT*GGGGGG/3Bio/")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertEqual(
            payload["summary"], "12 nt · 1 PS linkage · 2 modifications"
        )
        self.assertEqual(payload["svg"].count('class="tile"'), 12)
        self.assertEqual(payload["svg"].count('class="link ps"'), 1)
        self.assertEqual(payload["errors"], [])

    def test_errors_do_not_suppress_drawing(self):
        response = self.client.get("/api/render?seq=ACGT!")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertEqual(len(payload["errors"]), 1)
        self.assertEqual(payload["svg"].count('class="tile"'), 4)

    def test_modified_base_codes_render(self):
        response = self.client.get(
            "/api/render", query_string={"seq": "/i2FA/*/i2FC/"}
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertEqual(
            payload["summary"], "2 nt · 1 PS linkage · 0 modifications"
        )
        self.assertIn('class="sugar-icon"', payload["svg"])
        self.assertNotIn("/i2FA/", payload["svg"])

    def test_structure_not_signed_in(self):
        self.fake._signed_in = False
        response = self.client.get("/api/structure?seq=ACGT")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], False)
        self.assertIs(payload["login_required"], True)
        self.assertIs(payload["signed_in"], False)
        self.assertIn("Sign in to IDT", payload["errors"][0]["message"])
        self.assertEqual(self.fake.predict_calls, 0)

    def test_structure_success(self):
        response = self.client.get("/api/structure?seq=ACGTACGTACGTACGTACGT")
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertIs(payload["signed_in"], True)
        self.assertIn("ΔG -9.47", payload["summary"])
        self.assertIn("25 °C", payload["detail"])
        self.assertIn("50 mM", payload["detail"])
        self.assertIn("3 mM", payload["detail"])
        self.assertIn('class="base-pair wc"', payload["svg"])
        self.assertEqual(response.headers["Cache-Control"], "no-store")
        self.assertEqual(payload["structure"]["name"], "2")
        self.assertEqual(payload["structure"]["delta_g_kcal_mol"], -9.467)
        self.assertEqual(payload["structure"]["layout"], "stem-loop")
        self.assertEqual(len(payload["structure"]["base_pairs"]), 8)
        dimer = payload["homodimer"]
        self.assertIsNotNone(dimer)
        self.assertIs(dimer["ok"], True)
        self.assertEqual(dimer["layout"], "homodimer")
        self.assertIn('data-layout="homodimer"', dimer["svg"])
        self.assertEqual(dimer["base_pairs"], 20)
        self.assertIn("most likely homodimer", dimer["summary"])
        properties = payload["properties"]
        self.assertEqual(properties["melt_temp_c"], 73.2)
        self.assertEqual(properties["molecular_weight_g_mol"], 8649.6)
        self.assertEqual(
            properties["extinction_coefficient_l_mol_cm"], 269700.0
        )
        self.assertEqual(properties["nmol_per_od260"], 3.71)
        self.assertEqual(properties["ug_per_od260"], 32.07)

    def test_structure_uses_idt_dna_mode_even_with_u(self):
        response = self.client.get(
            "/api/structure?seq=CCCCCCCCCCTUTUTGGGGGGGGGGACTGACTGACTGTTTTGGGAAAA"
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(self.fake.predict_args[-1][1], "DNA")

    def test_structure_tree_layout_reaches_payload(self):
        self.fake.fold = Fold(
            bases="A" * 16,
            pairs=((0, 7), (1, 6), (8, 15), (9, 14)),
            structure_name="1",
            delta_g=-9.0,
            tm_c=60.0,
            delta_h=-50.0,
            delta_s=-150.0,
            structure_count=1,
            nucleotide_type="DNA",
            temperature_c=25,
            sodium_mm=50,
            magnesium_mm=3,
            suboptimality_percent=50,
            max_foldings=20,
        )
        response = self.client.get("/api/structure?seq=AAAAAAAAAAAAAAAA")
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertEqual(payload["structure"]["layout"], "tree")
        self.assertIn("branched fold", payload["detail"])
        self.assertIn('data-layout="tree"', payload["svg"])

    def test_structure_cache_avoids_second_prediction(self):
        url = "/api/structure?seq=ACGTACGTACGTACGTACGT"
        self.client.get(url)
        self.client.get(url)
        self.assertEqual(self.fake.predict_calls, 1)
        self.assertEqual(self.fake.properties_calls, 1)
        self.assertEqual(self.fake.homodimer_calls, 1)

    def test_structure_conditions_reach_idt(self):
        response = self.client.get(
            "/api/structure",
            query_string={
                "seq": "ACGTACGTACGTACGTACGT",
                "temp": "37",
                "na": "100",
                "mg": "1",
                "suboptimality": "25",
                "max_foldings": "10",
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            self.fake.predict_args[-1][2],
            FoldConditions(37, 100, 1, 25, 10),
        )
        self.assertEqual(self.fake.properties_args[-1][2], self.fake.predict_args[-1][2])
        self.assertEqual(self.fake.homodimer_args[-1][2], self.fake.predict_args[-1][2])

    def test_structure_conditions_separate_cache_entries(self):
        url = "/api/structure?seq=ACGTACGTACGTACGTACGT"
        changed = url + "&temp=37"
        self.client.get(url)
        self.client.get(url)
        self.client.get(changed)
        self.client.get(changed)
        self.assertEqual(self.fake.predict_calls, 2)

    def test_structure_invalid_conditions_do_not_call_idt(self):
        cases = (
            (
                {"temp": "101"},
                "Temperature must be a number between 0 and 100 °C.",
            ),
            (
                {"max_foldings": "2.5"},
                "Maximum foldings must be a whole number between 1 and 100.",
            ),
            (
                {"oligo": "0"},
                "Oligo concentration must be a number between 0.0001 and 10000 µM.",
            ),
            (
                {"dntps": "-1"},
                "dNTP concentration must be a number between 0 and 1000 mM.",
            ),
        )
        for query, message in cases:
            query["seq"] = "ACGT"
            payload = self.client.get(
                "/api/structure", query_string=query
            ).get_json()
            self.assertIs(payload["ok"], False)
            self.assertEqual(payload["errors"][0]["message"], message)
        self.assertEqual(self.fake.predict_calls, 0)
        self.assertEqual(self.fake.properties_calls, 0)
        self.assertEqual(self.fake.homodimer_calls, 0)

    def test_homodimer_conditions_reach_idt(self):
        response = self.client.get(
            "/api/structure",
            query_string={
                "seq": "ACGTACGTACGTACGTACGT",
                "oligo": "0.25",
                "dntps": "0",
            },
        )
        self.assertEqual(response.status_code, 200)
        conditions = self.fake.homodimer_args[-1][2]
        self.assertEqual(conditions.oligo_conc_um, 0.25)
        self.assertEqual(conditions.dntps_mm, 0.0)

    def test_property_failure_keeps_structures(self):
        self.fake.properties_error = IdtError("Property model failed.")
        response = self.client.get(
            "/api/structure?seq=ACGTACGTACGTACGTACGT"
        )
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertIsNone(payload["properties"])
        self.assertTrue(
            any(
                "Property model failed" in warning["message"]
                for warning in payload["warnings"]
            )
        )

    def test_homodimer_additional_pairs_get_legend_key(self):
        self.fake.homodimer = Homodimer(
            sequence="ACGT",
            delta_g=-8.0,
            base_pairs=2,
            additional_pairs=2,
            max_delta_g=-12.0,
            dimer="5' ACGT\n   ||::\n3' TGCA",
            result_count=1,
        )
        response = self.client.get(
            "/api/structure?seq=ACGTACGTACGTACGTACGT"
        )
        payload = response.get_json()
        self.assertIn(
            {"kind": "pair", "label": "Additional, not counted", "style": "additional"},
            payload["legend"],
        )

    def test_homodimer_failure_keeps_hairpin(self):
        self.fake.homodimer_error = NoStructureError(
            "IDT found no self-dimer for this sequence."
        )
        response = self.client.get(
            "/api/structure?seq=ACGTACGTACGTACGTACGT"
        )
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertIsNone(payload["homodimer"])
        self.assertTrue(
            any(
                "no self-dimer" in warning["message"]
                for warning in payload["warnings"]
            )
        )

    def test_structure_sugar_icons_from_local_sequence(self):
        self.fake.fold = Fold(
            bases="GCGCAAAAGCG",
            pairs=((0, 10), (1, 9), (2, 8)),
            structure_name="1",
            delta_g=-3.2,
            tm_c=40.0,
            delta_h=-30.0,
            delta_s=-90.0,
            structure_count=1,
            nucleotide_type="DNA",
            temperature_c=25,
            sodium_mm=50,
            magnesium_mm=3,
            suboptimality_percent=50,
            max_foldings=20,
        )
        response = self.client.get(
            "/api/structure",
            query_string={"seq": "rGrCrGrCrArArArArGrCrG"},
        )
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertEqual(payload["svg"].count('class="structure-sugar"'), 11)
        self.assertIn(">OH<", payload["svg"])
        sugars = [item for item in payload["legend"] if item["kind"] == "sugar"]
        self.assertEqual([item["label"] for item in sugars], ["2'-OH"])

    def test_structure_sugar_mismatch_hides_icons_and_warns(self):
        response = self.client.get(
            "/api/structure",
            query_string={"seq": "rGrCrGrCrArArArArGrCrG"},
        )
        payload = response.get_json()
        self.assertIs(payload["ok"], True)
        self.assertNotIn('class="structure-sugar"', payload["svg"])
        messages = [warning["message"] for warning in payload["warnings"]]
        self.assertIn(
            (
                "IDT returned a different number of bases than the parsed "
                "sequence; 2' modification icons are hidden."
            ),
            messages,
        )
        self.assertEqual(
            [item for item in payload["legend"] if item["kind"] == "sugar"],
            [],
        )

    def test_structure_no_base_pairs_is_warning_not_error(self):
        self.fake.error = NoStructureError(
            "IDT predicted no base pairs for this sequence."
        )
        response = self.client.get("/api/structure?seq=ACGT")
        payload = response.get_json()
        self.assertIs(payload["ok"], False)
        self.assertEqual(payload["errors"], [])
        self.assertIn("no base pairs", payload["warnings"][0]["message"])

    def test_structure_idt_error_is_red(self):
        self.fake.error = IdtError("IDT rejected the sequence.", 2, "!")
        response = self.client.get("/api/structure?seq=ACGT")
        payload = response.get_json()
        self.assertIs(payload["ok"], False)
        self.assertEqual(payload["errors"][-1]["message"], "IDT rejected the sequence.")
        self.assertEqual(payload["errors"][-1]["position"], 2)
        self.assertEqual(payload["errors"][-1]["token"], "!")

    def test_structure_expired_session_is_cleared(self):
        self.fake.error = IdtSessionExpired("IDT session expired; sign in again.")
        response = self.client.get("/api/structure?seq=ACGT")
        payload = response.get_json()
        self.assertIs(payload["ok"], False)
        self.assertIs(payload["signed_in"], False)
        self.assertIs(payload["login_required"], True)
        self.assertEqual(self.fake.sign_out_calls, 1)
        self.assertFalse(self.fake.signed_in)

    def test_login_status_and_logout(self):
        self.fake._signed_in = False
        response = self.client.post(
            "/api/idt/login",
            json={"username": "me@example.com", "password": "secret"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIs(response.get_json()["ok"], True)
        self.assertEqual(self.fake.sign_in_calls, [("me@example.com", "secret")])
        self.assertEqual(response.headers["Cache-Control"], "no-store")

        status = self.client.get("/api/idt/status").get_json()
        self.assertIs(status["signed_in"], True)

        self.client.post("/api/idt/logout")
        status = self.client.get("/api/idt/status").get_json()
        self.assertIs(status["signed_in"], False)

    def test_login_failure_reports_error(self):
        self.fake.login_error = IdtAuthError("Invalid IDT credentials.")
        response = self.client.post(
            "/api/idt/login",
            json={"username": "me@example.com", "password": "wrong"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIs(response.get_json()["ok"], False)
        self.assertEqual(response.get_json()["error"], "Invalid IDT credentials.")

    def test_login_empty_values_do_not_call_sign_in(self):
        response = self.client.post("/api/idt/login", json={})
        self.assertEqual(response.status_code, 200)
        self.assertIs(response.get_json()["ok"], False)
        self.assertEqual(self.fake.sign_in_calls, [])

    def test_login_cross_origin_is_rejected(self):
        response = self.client.post(
            "/api/idt/login",
            json={"username": "me@example.com", "password": "secret"},
            headers={"Origin": "http://evil.example"},
        )
        self.assertEqual(response.status_code, 403)
        self.assertEqual(self.fake.sign_in_calls, [])


if __name__ == "__main__":
    unittest.main()
