"""Scanpy adapter on a small synthetic AnnData. Skipped when scanpy is not installed."""

import importlib.util
import unittest

import helpers  # noqa: F401  (puts python/ on the path)

HAVE_SCANPY = importlib.util.find_spec("scanpy") is not None


@unittest.skipUnless(HAVE_SCANPY, "scanpy is not installed")
class ScanpyAdapter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import anndata as ad
        import numpy as np
        import pandas as pd
        import scanpy as sc

        rng = np.random.default_rng(0)
        n_cells, n_genes = 300, 120
        groups = np.repeat(["a", "b", "c"], n_cells // 3)
        counts = rng.poisson(1.0, size=(n_cells, n_genes)).astype(float)
        for k, g in enumerate(["a", "b", "c"]):
            counts[np.ix_(groups == g, range(k * 10, k * 10 + 10))] += rng.poisson(6.0, size=(100, 10))
        counts[:, -5:] = 0  # features with no counts at all
        obs = pd.DataFrame({
            "group": pd.Categorical(groups),
            "donor": np.tile(["d1", "d2"], n_cells // 2),
            "run": np.where(groups == "c", "r2", "r1"),
            "predicted_doublet": rng.random(n_cells) < 0.05,
            "n_genes_by_counts": (counts > 0).sum(axis=1),
        }, index=[f"cell{i}" for i in range(n_cells)])
        adata = ad.AnnData(counts, obs=obs, var=pd.DataFrame(index=[f"G{i:03d}" for i in range(n_genes)]))
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        cls.adata = adata

    def test_cluster_mode(self):
        import llomics

        payload = llomics.rank_genes_groups(self.adata, "group", method="wilcoxon", n_features=6,
                                            panel=["G000", "G015", "NOT_A_GENE"], sample_key="donor",
                                            batch_key="run", question="Which group is which?")
        self.assertEqual(llomics.validate(payload)[0], [])
        self.assertEqual(payload["method"], "clusters")
        self.assertEqual([u["unit_id"] for u in payload["units"]], ["a", "b", "c"])
        unit_a = payload["units"][0]
        self.assertTrue({r["feature"] for r in unit_a["top_features"]} <= {f"G{i:03d}" for i in range(10)})
        self.assertIn("pct_in", unit_a["top_features"][0])
        self.assertIsNone(unit_a["panel"]["NOT_A_GENE"])
        self.assertGreater(unit_a["panel"]["G000"], unit_a["panel"]["G015"])
        by = {c["check"]: c for c in payload["trust"]}
        self.assertEqual(by["batch_integration"]["status"], "warn")  # group c is entirely run r2
        self.assertEqual(by["doublet_detection"]["status"], "pass")
        self.assertEqual(payload["dataset"]["n_biological_samples"], 2)

    def test_contrast_mode(self):
        import llomics

        payload = llomics.rank_genes_groups(self.adata, "group", groups=["a"], reference="b", method="wilcoxon",
                                            n_features=5, panel=["G000", "NOT_A_GENE"], sample_key="donor",
                                            key="a_vs_b")
        self.assertEqual(llomics.validate(payload)[0], [])
        self.assertEqual(payload["method"], "dge")
        unit = payload["units"][0]
        self.assertEqual(unit["unit_id"], "a_vs_b")
        self.assertIn("higher in a", unit["definition"])
        self.assertTrue(all(r["effect"] > 0 for r in unit["top_features"]["up"]))
        self.assertTrue(all(r["effect"] < 0 for r in unit["top_features"]["down"]))
        self.assertGreaterEqual(unit["facts"]["n_uninformative"], 5)  # the all-zero features
        self.assertIsNotNone(unit["panel"]["G000"])
        self.assertIsNone(unit["panel"]["NOT_A_GENE"])
        by = {c["check"]: c for c in payload["trust"]}
        self.assertEqual(by["replicates"]["status"], "warn")  # two donors, but cells were the test units
        self.assertEqual(by["mean_variance"]["status"], "not_applicable")

    def test_reuse_requires_matching_groupby(self):
        import llomics

        llomics.rank_genes_groups(self.adata, "group", method="wilcoxon", key="stored")
        reused = llomics.rank_genes_groups(self.adata, "group", rerun=False, key="stored", n_features=3)
        self.assertEqual(len(reused["units"]), 3)
        with self.assertRaises(ValueError):
            llomics.rank_genes_groups(self.adata, "donor", rerun=False, key="stored")


if __name__ == "__main__":
    unittest.main()
