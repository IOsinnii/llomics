"""Repository-level checks: synced copies, committed example payloads, and the ask step's logic."""

import copy
import importlib.util
import json
import subprocess
import sys
import threading
import unittest
from types import SimpleNamespace as NS

from helpers import ROOT, build_dge

from llomics import validate


def load_judge():
    path = ROOT / "skills" / "llomics-judge" / "scripts" / "judge_bedrock.py"
    spec = importlib.util.spec_from_file_location("judge_bedrock", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Repo(unittest.TestCase):
    def test_copies_are_in_sync(self):
        result = subprocess.run([sys.executable, str(ROOT / "tools" / "sync.py"), "--check"],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_committed_example_payloads_are_valid(self):
        for path in sorted((ROOT / "examples").rglob("*.json")):
            payload = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(payload, dict) and "llomics" in payload:
                with self.subTest(path=path.name):
                    self.assertEqual(validate(payload)[0], [])

    def test_skill_validator_runs_standalone(self):
        script = ROOT / "skills" / "llomics-core" / "scripts" / "validate_payload.py"
        example = ROOT / "examples" / "pbmc10k" / "payload.json"
        result = subprocess.run([sys.executable, str(script), str(example)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout)
        self.assertIn("VALID", result.stdout)


class Judge(unittest.TestCase):
    """The ask step with a stand-in for the model: parsing, vote shares and the three kinds of gate."""

    GOOD = {"dominant_signal": "response to the comparison", "up_is_one_program": True, "up_program_label": "x",
            "down_is_one_program": False, "down_program_label": "none", "evidence": "F0001"}

    @classmethod
    def setUpClass(cls):
        cls.jb = load_judge()
        cls.cfg = json.loads((ROOT / "skills" / "llomics-limma" / "references" / "questions.contrast.json")
                             .read_text(encoding="utf-8"))

    def setUp(self):
        payload = build_dge()
        payload["units"] = payload["units"][:1]
        for chk in payload["trust"] + payload["units"][0]["trust"]:
            chk["status"] = "pass"  # start from a payload that passes every trust gate
        self.payload = payload

    def response(self, text=None, stop="end_turn"):
        content = [NS(type="thinking", thinking="")] + ([NS(type="text", text=text)] if text is not None else [])
        details = NS(category="bio") if stop == "refusal" else None
        return NS(stop_reason=stop, content=content, model="stub", stop_details=details)

    def sequence(self, results):
        it, lock = iter(results), threading.Lock()

        def call(_request):
            with lock:
                return next(it)
        return call

    def test_parse_response(self):
        jb = self.jb
        self.assertTrue(jb.parse_response(self.response(json.dumps(self.GOOD)))["ok"])
        self.assertEqual(jb.parse_response(self.response(stop="refusal")), {"ok": False, "error": "refusal (bio)"})
        self.assertFalse(jb.parse_response(self.response("{", stop="max_tokens"))["ok"])
        self.assertFalse(jb.parse_response(self.response("not json"))["ok"])
        self.assertFalse(jb.parse_response(self.response(None))["ok"])

    def test_unanimous_answers_pass(self):
        call = lambda _r: {"ok": True, "answers": self.GOOD}
        result = self.jb.judge_units(self.payload, self.cfg, call, 5, 3)[0]
        self.assertEqual((result["status"], result["reasons"]), ("auto", []))
        self.assertEqual(result["summary"]["dominant_signal"]["share"], 1.0)
        self.assertEqual(len(result["runs"]), 5)

    def test_split_vote_goes_to_review(self):
        bad = dict(self.GOOD, dominant_signal="no recognisable pattern")
        runs = [{"ok": True, "answers": a} for a in (self.GOOD, self.GOOD, self.GOOD, bad, bad)]
        result = self.jb.judge_units(self.payload, self.cfg, self.sequence(runs), 5, 1)[0]
        self.assertEqual(result["status"], "review")
        self.assertEqual(result["reasons"], ["signal_is_not_the_comparison"])
        self.assertEqual(result["summary"]["dominant_signal"]["share"], 0.6)

    def test_failed_request_is_never_an_answer(self):
        runs = [{"ok": True, "answers": self.GOOD}, {"ok": False, "error": "http_403"}, {"ok": True, "answers": self.GOOD}]
        result = self.jb.judge_units(self.payload, self.cfg, self.sequence(runs), 3, 1)[0]
        self.assertEqual(result["status"], "review")
        self.assertTrue(result["reasons"][0].startswith("no_answer"))

    def test_trust_and_fact_gates_override_a_confident_model(self):
        payload = copy.deepcopy(self.payload)
        for chk in payload["trust"]:
            if chk["check"] == "batch_confounding":
                chk["status"] = "not_run"
        payload["units"][0]["facts"]["n_hits"] = 4
        call = lambda _r: {"ok": True, "answers": self.GOOD}
        result = self.jb.judge_units(payload, self.cfg, call, 3, 2)[0]
        self.assertEqual(result["reasons"], ["trust:batch_confounding:not_run", "too_few_hits"])

    def test_saved_runs_gate_as_when_they_were_asked(self):
        bad = dict(self.GOOD, dominant_signal="no recognisable pattern")
        runs = [{"ok": True, "answers": a} for a in (self.GOOD, self.GOOD, bad)]
        asked = self.jb.judge_units(self.payload, self.cfg, self.sequence(runs), 3, 1)[0]
        regated = self.jb.gate_units(self.payload, self.cfg, [runs])[0]
        self.assertEqual(regated, asked)
        with self.assertRaises(ValueError):
            self.jb.gate_units(self.payload, self.cfg, [runs, runs])

    def test_requests_need_no_model(self):
        requests = self.jb.build_requests(self.payload, self.cfg)
        self.assertEqual(requests, [self.jb.judge_units(self.payload, self.cfg, None, 0, 1)[0]["request"]])

    def test_dry_run_sends_nothing(self):
        result = self.jb.judge_units(self.payload, self.cfg, None, 0, 1)[0]
        self.assertEqual((result["status"], result["runs"]), ("not_asked", []))
        self.assertIn('"unit"', result["request"]["user"])

    def test_schema(self):
        schema = self.jb.build_schema(self.cfg["questions"])
        self.assertFalse(schema["additionalProperties"])
        self.assertEqual(set(schema["required"]), set(self.GOOD))


if __name__ == "__main__":
    unittest.main()
