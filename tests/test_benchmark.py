"""The benchmark's gate-only verdicts are pinned, so a change to a threshold or a gate shows up here."""

import importlib.util
import json
import unittest

from helpers import ROOT

from llomics import validate


def load_score():
    spec = importlib.util.spec_from_file_location("score", ROOT / "benchmark" / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Benchmark(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.score = load_score()
        cls.rows = {r["id"]: r for r in cls.score.evaluate()}

    def test_every_scenario_has_a_valid_payload(self):
        for sid in self.rows:
            path = ROOT / "benchmark" / "payloads" / f"{sid}.llomics.json"
            with self.subTest(scenario=sid):
                self.assertEqual(validate(json.loads(path.read_text(encoding="utf-8")))[0], [])

    def test_no_miss_and_no_false_alarm_from_the_gates(self):
        wrong = {sid: r["without_model"] for sid, r in self.rows.items()
                 if r["without_model"] in ("MISSED", "FALSE ALARM")}
        self.assertEqual(wrong, {})

    def test_clean_runs_reach_the_model(self):
        for sid in ("clean_interferon", "clean_hypoxia"):
            self.assertEqual(self.rows[sid]["reasons"], [])

    def test_planted_problems_are_caught_by_the_expected_gate(self):
        expected = {
            "null": ["too_few_hits"],
            "batch_aliased": ["trust:batch_confounding:fail"],
            "no_batch_info": ["trust:batch_confounding:not_run"],
            "unmodelled_batch": ["trust:pvalue_distribution:warn"],
            "two_replicates": ["trust:replicates:warn"],
        }
        for sid, reasons in expected.items():
            with self.subTest(scenario=sid):
                self.assertEqual(self.rows[sid]["reasons"], reasons)

    def test_known_blind_spot_is_still_reported_as_one(self):
        self.assertTrue(self.rows["technical_replicates"]["without_model"].startswith("not caught"))

    def test_old_gate_settings_would_miss_two_replicates(self):
        old = {"blocking_checks": ["replicates", "pvalue_distribution", "batch_confounding"],
               "blocking_statuses": ["fail", "not_run"], "facts": [{"fact": "n_hits", "op": "<", "value": 10, "reason": "too_few_hits"}]}
        rows = {r["id"]: r for r in self.score.evaluate(old)}
        self.assertEqual(rows["two_replicates"]["without_model"], "MISSED")
        self.assertEqual(rows["unmodelled_batch"]["without_model"], "MISSED")


if __name__ == "__main__":
    unittest.main()
