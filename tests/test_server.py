"""Tests for the Flask layer, using the test client (no real port)."""

import unittest

import app


class ServerTests(unittest.TestCase):
    def setUp(self):
        self.client = app.app.test_client()

    def test_index_serves_page(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        body = response.get_data(as_text=True)
        response.close()
        self.assertIn('id="seq"', body)

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


if __name__ == "__main__":
    unittest.main()
