"""Tests for the Flask layer, using the test client (no real port)."""

import unittest

import app
from idt import Fold, IdtAuthError, IdtError, IdtSessionExpired, NoStructureError


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

    def predict(self, sequence, nucleotide_type):
        self.predict_calls += 1
        self.predict_args.append((sequence, nucleotide_type))
        if self.error is not None:
            raise self.error
        return self.fold


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
        self.assertIn('class="base-pair"', payload["svg"])
        self.assertEqual(response.headers["Cache-Control"], "no-store")
        self.assertEqual(payload["structure"]["name"], "2")
        self.assertEqual(payload["structure"]["delta_g_kcal_mol"], -9.467)
        self.assertEqual(payload["structure"]["layout"], "stem-loop")
        self.assertEqual(len(payload["structure"]["base_pairs"]), 8)

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
