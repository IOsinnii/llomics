"""The browser demo's Python side (web/bridge.py) and the site build (tools/build_site.py).

The page runs bridge.py under Pyodide; here it runs under CPython with the same files. Pyodide
itself and the JavaScript are not exercised by these tests.
"""

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from helpers import ROOT


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


bridge = load("bridge", ROOT / "web" / "bridge.py")


def call(fn, **args):
    return json.loads(bridge.call(fn, json.dumps(args)))


class Gallery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = call("gallery")["result"]
        cls.score = load("score", ROOT / "benchmark" / "score.py")

    def test_summary_matches_the_scorer(self):
        rows = self.score.evaluate()
        s = self.result["summary"]
        self.assertEqual(s["n_scenarios"], len(rows))
        self.assertEqual(s["caught_by_gates"], sum(1 for r in rows if r["expected"] == "review" and r["reasons"]))
        self.assertEqual((s["false_alarms"], s["pending_model"], s["blind_spots"]), (0, 1, 1))

    def test_each_scenario_carries_the_scorer_reasons(self):
        rows = {r["id"]: r for r in self.score.evaluate()}
        for item in self.result["scenarios"]:
            with self.subTest(scenario=item["id"]):
                self.assertEqual(item["reasons"], rows[item["id"]]["reasons"])
                self.assertTrue((ROOT / item["payload_path"]).exists())

    def test_trust_badges_name_each_check_once_with_the_worst_status(self):
        for item in self.result["scenarios"]:
            names = [c["check"] for c in item["trust"]]
            self.assertEqual(len(names), len(set(names)), item["id"])
        unmodelled = next(i for i in self.result["scenarios"] if i["id"] == "unmodelled_batch")
        status = {c["check"]: c["status"] for c in unmodelled["trust"]}
        self.assertEqual(status["pvalue_distribution"], "warn")

    def test_saved_answers_are_gated_again(self):
        cfg = bridge._questions()
        rows = {r["id"]: r for r in self.score.evaluate()}
        good = {"dominant_signal": "response to the comparison", "up_is_one_program": True, "up_program_label": "x",
                "down_is_one_program": False, "down_program_label": "none", "evidence": "x"}
        composition = dict(good, dominant_signal="cell-type composition")
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / "stub-model"
            folder.mkdir()
            for sid, answers in (("clean_interferon", good), ("composition", composition)):
                runs = [{"ok": True, "answers": answers, "model": "stub"}] * 3
                # status and reasons in the file are deliberately wrong: the page must recompute them
                record = {"units": [{"unit_id": "x", "status": "review", "reasons": ["stale"], "runs": runs}]}
                (folder / f"{sid}.answers.json").write_text(json.dumps(record), encoding="utf-8")
            for sid, status in (("clean_interferon", "auto"), ("composition", "review")):
                payload = json.loads((ROOT / "benchmark" / "payloads" / f"{sid}.llomics.json").read_text(encoding="utf-8"))
                replay = bridge._replay(payload, cfg, rows[sid], folder)
                with self.subTest(scenario=sid):
                    self.assertEqual((replay["status"], replay["verdict"], replay["n_runs"]), (status, "correct", 3))
                    self.assertNotIn("stale", replay["reasons"])
            missing = bridge._replay(payload, cfg, rows["null"], folder)
            self.assertIsNone(missing)


class Build(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.example = call("example")["result"]
        cls.built = call("build", form=cls.example["form"], table_text=cls.example["text"],
                         filename=cls.example["filename"])

    def test_example_is_deterministic(self):
        self.assertEqual(call("example")["result"]["text"], self.example["text"])

    def test_example_builds_a_valid_payload_that_passes_the_gates(self):
        self.assertTrue(self.built["ok"], self.built.get("error"))
        r = self.built["result"]
        self.assertEqual(r["inspect"]["errors"], [])
        self.assertEqual(r["table"]["tool_guess"], "limma")
        self.assertEqual({c["status"] for c in r["inspect"]["trust"]}, {"pass"})
        unit = r["inspect"]["units"][0]
        self.assertEqual(unit["gate"], {"status": "not_asked", "reasons": []})
        self.assertGreaterEqual(unit["unit"]["facts"]["n_up"], 40)
        self.assertIsNotNone(unit["unit"]["panel"]["ISG15"])
        self.assertIn('"unit"', unit["request"]["user"])
        self.assertLess(r["inspect"]["chars"], r["input_chars"] / 10)

    def test_python_snippet_compiles(self):
        compile(self.built["result"]["python"], "snippet", "exec")

    def test_built_payload_matches_a_direct_call(self):
        from llomics import core
        parsed = bridge.read_table(self.example["text"], self.example["filename"])
        form = self.example["form"]
        payload = core.dge_payload(
            {form["unit_id"]: parsed["table"]}, {form["unit_id"]: form["definition"]},
            analysis_id=form["analysis_id"], question=form["question"], tool="limma",
            steps=form["steps"].split("\n"), n_observations=12, group_sizes={"ifn": 6, "control": 6},
            batch_crosstab={"ifn": {"b1": 3, "b2": 3}, "control": {"b1": 3, "b2": 3}}, batch_in_design=True,
            mean_variance="trend", panel=["ISG15", "MX1", "ACTB", "GAPDH"],
            sizes={form["unit_id"]: {"n_per_group": {"ifn": 6, "control": 6}}},
            source_tables={form["unit_id"]: self.example["filename"]})
        built = json.loads(self.built["result"]["payload_text"])
        self.assertEqual(built, json.loads(json.dumps(core.jsonable(payload))))

    def test_form_errors_are_messages(self):
        form = dict(self.example["form"], definition="")
        r = call("build", form=form, table_text=self.example["text"])
        self.assertFalse(r["ok"])
        self.assertIn("definition", r["error"])
        r = call("build", form=dict(self.example["form"], groups="treated six"), table_text=self.example["text"])
        self.assertIn("name=number", r["error"])
        r = call("build", form=dict(self.example["form"], exclude="(["), table_text=self.example["text"])
        self.assertIn("exclude", r["error"])


class ReadTable(unittest.TestCase):
    def test_tool_formats(self):
        cases = {
            "DESeq2": ('"","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"\n'
                       '"G1",10.5,2.1,0.3,7,1e-10,1e-8\n"G2",8,0.1,0.3,0.3,0.7,NA\n', "res.csv"),
            "Seurat": ('"","p_val","avg_log2FC","pct.1","pct.2","p_val_adj"\n"CD3D",0,2.5,0.9,0.1,0\n', "markers.csv"),
            "scanpy": ("names\tscores\tlogfoldchanges\tpvals\tpvals_adj\nCD3D\t30\t2.5\t1e-50\t1e-48\n", "df.tsv"),
            "edgeR": ("logFC\tlogCPM\tF\tPValue\tFDR\nG1\t2\t5\t40\t1e-5\t1e-3\n", "toptags.txt"),
        }
        for tool, (text, name) in cases.items():
            with self.subTest(tool=tool):
                r = bridge.read_table(text, name)
                self.assertEqual(r["tool_guess"], tool)
                self.assertEqual(set(r["table"]) >= {"feature", "effect", "p", "p_adj"}, True)
        deseq = bridge.read_table(*cases["DESeq2"])
        self.assertEqual(deseq["table"]["feature"], ["G1", "G2"])
        self.assertEqual(deseq["table"]["p_adj"], [1e-8, None])
        self.assertIn("baseMean", deseq["unmapped"])
        # write.table: the header is one column short of the rows
        self.assertEqual(bridge.read_table(*cases["edgeR"])["table"]["feature"], ["G1"])

    def test_missing_adjusted_p_is_computed(self):
        r = bridge.read_table("gene,logFC,P.Value\nA,1,0.01\nB,-1,0.04\nC,0.1,NA\n", "x.csv")
        self.assertEqual(r["table"]["p_adj"], [0.02, 0.04, None])
        self.assertTrue(any("Benjamini-Hochberg" in n for n in r["notes"]))

    def test_unreadable_tables(self):
        for text, fragment in (("gene,score\nA,1\n", "'effect'"), ("gene,logFC,p\nA,big,0.1\n", "line 2"),
                               ("logFC,p\n1,0.1\n", "feature names"), ("", "empty")):
            with self.subTest(fragment=fragment):
                r = call("build", form={"definition": "a - b"}, table_text=text, filename="t.csv")
                self.assertFalse(r["ok"])
                self.assertIn(fragment, r["error"])

    def test_batch_table_forms(self):
        expected = {"a": {"b1": 2, "b2": 1}, "c": {"b1": 0, "b2": 3}}
        self.assertEqual(bridge._crosstab("b1 b2\na 2 1\nc 0 3"), expected)
        self.assertEqual(bridge._crosstab("group,b1,b2\na,2,1\nc,0,3"), expected)
        self.assertIsNone(bridge._crosstab("  "))
        for bad in ("b1 b2\na 2 1\nc 3", "b1 b2\na 1 2 3", "b1 b2\na x y", "b1 b2"):
            with self.subTest(table=bad), self.assertRaises(ValueError):
                bridge._crosstab(bad)


class Inspect(unittest.TestCase):
    def test_bad_inputs(self):
        self.assertIn("parse_error", call("inspect", text="{")["result"])
        r = call("inspect", text=json.dumps({"llomics": "0.1"}))["result"]
        self.assertTrue(r["errors"])

    def test_cluster_payload_is_shown_without_a_request(self):
        r = call("inspect", text=(ROOT / "examples" / "pbmc10k" / "payload.json").read_text(encoding="utf-8"))["result"]
        self.assertEqual(r["errors"], [])
        self.assertFalse(r["askable"])
        self.assertIn("cluster units", r["questions_note"])
        self.assertEqual(len(r["units"]), 19)
        self.assertIsNone(r["units"][0]["request"])

    def test_unknown_function(self):
        self.assertFalse(call("nope")["ok"])


class Site(unittest.TestCase):
    def test_built_site_runs_on_its_own_files(self):
        build_site = load("build_site", ROOT / "tools" / "build_site.py")
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "site"
            files = build_site.build(out, with_pyodide=False)
            for name in build_site.PAGE:
                self.assertTrue((out / name).exists(), name)
            manifest = json.loads((out / "repo" / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["files"], files)
            for rel in files:
                self.assertTrue((out / "repo" / rel).is_file(), rel)
            # -I: no PYTHONPATH, no user site, so only the copied files can be imported
            script = ("import sys, json; sys.path.insert(0, sys.argv[1]); import bridge; "
                      "print(json.loads(bridge.call('gallery'))['ok'], json.loads(bridge.call('example'))['ok'])")
            result = subprocess.run([sys.executable, "-I", "-c", script, str(out / "repo" / "web")],
                                    capture_output=True, text=True, cwd=tmp)
            self.assertEqual(result.stdout.split(), ["True", "True"], result.stderr)

    def test_pyodide_version_matches_the_worker(self):
        build_site = load("build_site", ROOT / "tools" / "build_site.py")
        worker = (ROOT / "web" / "worker.js").read_text(encoding="utf-8")
        self.assertIn("./pyodide/pyodide.mjs", worker)
        self.assertIn("pyodide.mjs", build_site.PYODIDE_FILES)
        self.assertTrue(all(len(h) == 64 for h in build_site.PYODIDE_FILES.values()))


if __name__ == "__main__":
    unittest.main()
