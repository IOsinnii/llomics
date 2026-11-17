"""Core behaviour: golden payloads from fixtures, and the rules of the contract one at a time."""

import copy
import unittest

from helpers import build_clusters, build_dge, differences, golden

from llomics import core, validate


class GoldenPayloads(unittest.TestCase):
    def test_dge_matches_expected(self):
        payload = build_dge()
        self.assertEqual(differences(payload, golden("expected_dge.json", payload), rel=1e-12), [])

    def test_clusters_matches_expected(self):
        payload = build_clusters()
        self.assertEqual(differences(payload, golden("expected_clusters.json", payload), rel=1e-12), [])

    def test_fixture_payloads_are_valid(self):
        for payload in (build_dge(), build_clusters()):
            errors, _ = validate(payload)
            self.assertEqual(errors, [])


class ContractRules(unittest.TestCase):
    def setUp(self):
        self.dge = build_dge()
        self.unit = self.dge["units"][0]
        self.clusters = build_clusters()

    def test_never_padded_and_tail_counted(self):
        small = self.dge["units"][1]
        shown = len(small["top_features"]["up"]) + len(small["top_features"]["down"])
        self.assertEqual(small["tail"]["n_shown"], shown)
        self.assertLessEqual(shown, small["tail"]["n_eligible"])
        table = {"feature": ["a", "b", "c", "d"], "effect": [2.0, 1.5, 0.1, -0.1],
                 "p": [1e-6, 1e-5, 0.6, 0.7], "p_adj": [1e-5, 1e-4, 0.7, 0.7]}
        unit = core.contrast_unit(table, unit_id="x", definition="x - y", n_features=5)
        self.assertEqual(len(unit["top_features"]["up"]), 2)  # two eligible, five asked for: not padded
        self.assertEqual(unit["top_features"]["down"], [])
        self.assertEqual(unit["tail"], {"n_shown": 2, "n_eligible": 2, "rule": unit["tail"]["rule"]})

    def test_null_is_not_zero(self):
        self.assertIsNone(self.unit["panel"]["NOT_MEASURED"])
        self.assertIsNone(self.unit["panel"]["F0400"])  # present in the table but not tested (p missing)
        self.assertIsNotNone(self.unit["panel"]["F0002"])
        self.assertIsNone(self.clusters["units"][0]["panel"]["NOT_MEASURED"])
        self.assertEqual(self.clusters["units"][0]["panel"]["MS4A1"], 0.0)

    def test_panel_matches_symbols_too(self):
        self.assertIsNotNone(self.unit["panel"]["SYM0050"])

    def test_no_zero_p_values(self):
        top = self.unit["top_features"]["up"][0]
        self.assertEqual(top["feature"], "F0001")  # the underflowed row has the largest statistic
        self.assertEqual(top["p_adj"], self.dge["run"]["p_floor"])
        self.assertGreaterEqual(self.dge["run"]["n_p_values_floored"], 1)

    def test_missing_adjusted_p_is_never_eligible(self):
        self.assertEqual(self.unit["facts"]["n_p_adj_missing"], 1)
        shown = [r["feature"] for d in ("up", "down") for r in self.unit["top_features"][d]]
        self.assertNotIn("F0040", shown)

    def test_excluded_names_are_dropped_and_declared(self):
        symbols = [r["symbol"] for d in ("up", "down") for r in self.unit["top_features"][d]]
        self.assertFalse(any(s and s.startswith("MT-") for s in symbols))
        self.assertTrue(any("exclude" in f["rule"] for f in self.dge["run"]["filters"]))

    def test_ties_keep_table_order(self):
        up = [r["feature"] for r in self.unit["top_features"]["up"]]
        self.assertLess(up.index("F0002"), up.index("F0003"))

    def test_uninformative_features_leave_the_histogram(self):
        facts = self.unit["facts"]
        self.assertEqual(facts["n_uninformative"], 6)
        self.assertEqual(sum(facts["p_histogram"]), self.unit["size"]["n_tested"] - 6)

    def test_trust_layer_is_complete(self):
        names = {c["check"] for c in self.dge["trust"]}
        self.assertEqual(names, set(core.CONTRACT["methods"]["dge"]["required_checks"]))
        by = {c["check"]: c for c in self.dge["trust"]}
        self.assertEqual(by["multiple_testing"]["status"], "warn")  # two contrasts
        self.assertIn("declared by the analyst", by["mean_variance"]["note"])

    def test_batch_dominated_cluster_is_flagged(self):
        by = {c["check"]: c for c in self.clusters["trust"]}
        self.assertEqual(by["batch_integration"]["status"], "warn")
        self.assertIn("c2", by["batch_integration"]["note"])
        c2 = [u for u in self.clusters["units"] if u["unit_id"] == "c2"][0]
        self.assertGreater(c2["facts"]["doublet_fraction_vs_reference"], 1)
        self.assertGreater(c2["facts"]["n_genes_vs_reference"], 1.5)


class TrustChecks(unittest.TestCase):
    def test_p_shape(self):
        self.assertEqual(core.p_shape([800, 480, 470, 440, 470, 480, 500, 460, 460, 450])[0], "pass")
        self.assertEqual(core.p_shape([100, 100, 100, 100, 100, 100, 100, 100, 100, 100])[0], "pass")
        self.assertEqual(core.p_shape([100, 120, 150, 200, 260, 300, 340, 400, 520, 900])[0], "warn")  # peak near 1
        self.assertEqual(core.p_shape([300, 310, 450, 520, 480, 400, 330, 300, 290, 280])[0], "warn")  # hill
        self.assertEqual(core.p_shape([900, 300, 250, 240, 230, 240, 235, 230, 240, 520])[0], "warn")  # U shape
        self.assertEqual(core.p_shape([0] * 10)[0], "warn")

    def test_replicates(self):
        self.assertEqual(core.replicates_check({"a": 6, "b": 6})["status"], "pass")
        self.assertEqual(core.replicates_check({"a": 2, "b": 6})["status"], "warn")
        self.assertEqual(core.replicates_check({"a": 1, "b": 6})["status"], "fail")
        self.assertEqual(core.replicates_check(None)["status"], "not_run")
        self.assertEqual(core.replicates_check({"a": 6, "b": 6}, units_are_replicates=False)["status"], "warn")
        self.assertEqual(core.replicates_check({"a": 1, "b": 1}, units_are_replicates=False)["status"], "fail")

    def test_batch_confounding(self):
        balanced = {"a": {"b1": 3, "b2": 3}, "b": {"b1": 3, "b2": 3}}
        aliased = {"a": {"b1": 6}, "b": {"b2": 6}}
        partial = {"a": {"b1": 3, "b2": 3}, "b": {"b1": 6}}
        status = lambda *a, **k: core.batch_confounding_check(*a, **k)["status"]
        self.assertEqual(status(balanced, True), "pass")
        self.assertEqual(status(balanced, False), "warn")
        self.assertEqual(status(balanced, None), "warn")
        self.assertEqual(status(aliased, True), "fail")
        self.assertEqual(status(partial, True), "warn")
        self.assertEqual(status(balanced, True, adjustable=False), "warn")
        self.assertEqual(status({"a": {"b1": 6}, "b": {"b1": 6}}, True), "not_applicable")
        self.assertEqual(status(None), "not_run")
        self.assertEqual(status(None, single_batch=True), "not_applicable")

    def test_p_shape_skip_records_not_run(self):
        table = {"feature": ["a", "b"], "effect": [1.0, -1.0], "p": [0.001, 0.5], "p_adj": [0.002, 0.5]}
        payload = core.dge_payload({"x": table}, {"x": "x - y"}, analysis_id="t", question="q", tool="t",
                                   steps=["s"], n_observations=4, p_shape_skip="features were pre-filtered")
        by = {c["check"]: c for c in payload["trust"]}
        self.assertEqual(by["pvalue_distribution"]["status"], "not_run")
        self.assertEqual(payload["units"][0]["trust"][0]["status"], "not_run")

    def test_ranks_without_a_statistic(self):
        table = {"feature": ["a", "b", "c"], "effect": [1.0, 3.0, 2.0], "p": [1e-4, 1e-4, 1e-6], "p_adj": [1e-3, 1e-3, 1e-5]}
        unit = core.contrast_unit(table, unit_id="x", definition="x - y")
        self.assertEqual([r["feature"] for r in unit["top_features"]["up"]], ["c", "b", "a"])
        self.assertIn("p then |effect|", unit["tail"]["rule"])

    def test_empty_table_gives_an_empty_unit(self):
        empty = {"feature": [], "effect": [], "p": [], "p_adj": []}
        unit = core.cluster_unit(empty, unit_id="c9")
        self.assertEqual(unit["top_features"], [])
        self.assertEqual(unit["tail"]["n_eligible"], 0)

    def test_missing_definition_is_an_error(self):
        table = {"feature": ["a"], "effect": [1.0], "p": [0.01], "p_adj": [0.02]}
        with self.assertRaises(ValueError):
            core.dge_payload({"x": table}, {}, analysis_id="t", question="q", tool="t", steps=[], n_observations=2)

    def test_missing_column_is_an_error(self):
        with self.assertRaises(KeyError):
            core.contrast_unit({"feature": ["a"], "effect": [1.0], "p": [0.1]}, unit_id="x", definition="d")


class TechnicalPanelsAndEnrichment(unittest.TestCase):
    def setUp(self):
        self.dge = build_dge()
        self.a, self.c = self.dge["units"]

    def test_sex_signature_is_flagged_where_planted(self):
        self.assertEqual(self.c["facts"]["technical_flags"], ["sex"])
        self.assertEqual(self.c["facts"]["n_technical_flags"], 1)
        sex = self.c["facts"]["technical_signatures"]["sex"]
        self.assertEqual((sex["n_measured"], sex["n_hits"]), (7, 6))  # includes the mouse-case symbol Ddx3y
        self.assertEqual(self.a["facts"]["technical_flags"], [])
        self.assertEqual(self.a["facts"]["n_technical_flags"], 0)

    def test_a_few_panel_hits_are_reported_but_not_flagged(self):
        mito = self.a["facts"]["technical_signatures"]["mitochondrial"]
        self.assertEqual(mito["n_hits"], 2)
        self.assertNotIn("stress", self.a["facts"]["technical_signatures"])  # measured, but no hits
        self.assertEqual(self.a["facts"]["n_panel_features_measured"], 6)

    def test_unmeasurable_is_null_not_zero(self):
        table = {"feature": ["m1", "m2", "m3"], "effect": [2.0, 1.0, 0.1], "p": [1e-6, 1e-4, 0.5], "p_adj": [1e-5, 1e-3, 0.5]}
        facts = core.contrast_unit(table, unit_id="x", definition="x - y")["facts"]
        self.assertEqual(facts["n_panel_features_measured"], 0)
        self.assertIsNone(facts["n_technical_flags"])

    def test_panels_can_be_switched_off_or_replaced(self):
        n = 60  # five hits among sixty features, all five in the custom panel
        table = {"feature": [f"g{i}" for i in range(n)], "effect": [2.0] * 5 + [0.1] * (n - 5),
                 "p": [1e-6] * 5 + [0.5] * (n - 5), "p_adj": [1e-5] * 5 + [0.5] * (n - 5)}
        off = core.contrast_unit(table, unit_id="x", definition="d", technical_panels=False)["facts"]
        self.assertNotIn("n_technical_flags", off)
        mine = {"mine": {"genes": ["G0", "G1", "G2", "G3", "G4"]}}
        custom = core.contrast_unit(table, unit_id="x", definition="d", technical_panels=mine)["facts"]
        self.assertEqual(custom["technical_flags"], ["mine"])
        # the same three-of-three overlap in a four-feature table is not evidence: p = 0.25
        tiny = {k: v[:3] + v[-1:] for k, v in table.items()}
        small = core.contrast_unit(tiny, unit_id="x", definition="d", technical_panels={"mine": {"genes": ["G0", "G1", "G2"]}})["facts"]
        self.assertEqual(small["technical_flags"], [])

    def test_enrichment_finds_planted_sets_only(self):
        e = self.a["enrichment"]
        self.assertEqual(e["n_sets_tested"], 3)  # too-small, unmeasured and other-universe sets are not tested
        self.assertEqual([s["set"] for s in e["up"]], ["SET_UP"])
        self.assertEqual([s["set"] for s in e["down"]], ["SET_DOWN"])  # the lower-case set matched
        up = e["up"][0]
        self.assertGreater(up["n_overlap"], 5 * up["expected"])
        self.assertLessEqual(len(up["genes"]), core.CONTRACT["enrichment"]["max_genes"])
        self.assertEqual(self.c["enrichment"]["up"], [])
        self.assertEqual(self.dge["run"]["enrichment"]["gene_sets"], "sets.gmt")

    def test_cluster_enrichment(self):
        clusters = build_clusters()
        self.assertEqual([s["set"] for s in clusters["units"][0]["enrichment"]["up"]], ["SET_CLUSTER"])
        self.assertIn("n_technical_flags", clusters["units"][0]["facts"])

    def test_statistics(self):
        from llomics import enrich

        self.assertAlmostEqual(enrich.hypergeom_sf(1, 10, 5, 2), 1 - 10 / 45, places=12)
        self.assertEqual(enrich.hypergeom_sf(0, 10, 5, 2), 1.0)
        self.assertEqual(enrich.hypergeom_sf(3, 10, 5, 2), 0.0)
        self.assertEqual(enrich.bh([0.01, 0.04, 0.03, 0.5]), [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])


class Validator(unittest.TestCase):
    def setUp(self):
        self.payload = build_dge()

    def broken(self, mutate):
        payload = copy.deepcopy(self.payload)
        mutate(payload)
        return validate(payload)[0]

    def test_missing_required_check(self):
        errors = self.broken(lambda p: p["trust"].pop(0))
        self.assertTrue(any("required check" in e for e in errors))

    def test_zero_p_value(self):
        errors = self.broken(lambda p: p["units"][0]["top_features"]["up"][0].update(p_adj=0.0))
        self.assertTrue(any("exactly 0" in e for e in errors))

    def test_wrong_tail_count(self):
        errors = self.broken(lambda p: p["units"][0]["tail"].update(n_shown=1))
        self.assertTrue(any("n_shown" in e for e in errors))

    def test_nan(self):
        errors = self.broken(lambda p: p["units"][0]["facts"].update(pi0=float("nan")))
        self.assertTrue(any("NaN" in e for e in errors))

    def test_unknown_method_and_wrong_unit_type(self):
        self.assertTrue(any("not in the contract" in e for e in self.broken(lambda p: p.update(method="nope"))))
        errors = self.broken(lambda p: p["units"][0].update(unit_type="cluster"))
        self.assertTrue(any("unit_type" in e for e in errors))

    def test_status_needs_note(self):
        errors = self.broken(lambda p: p["trust"][0].update(status="fail", note=""))
        self.assertTrue(any("needs a 'note'" in e for e in errors))

    def test_duplicate_unit_id(self):
        errors = self.broken(lambda p: p["units"][1].update(unit_id=p["units"][0]["unit_id"]))
        self.assertTrue(any("duplicate" in e for e in errors))


if __name__ == "__main__":
    unittest.main()
