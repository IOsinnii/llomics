"""R and Python must build the same payload from the same tables.

Needs Rscript with the llomics R package installed (R CMD INSTALL r/llomics). Skipped otherwise.
Set LLOMICS_RSCRIPT to the Rscript executable if it is not on the path.
"""

import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from helpers import FIXTURES, ROOT, build_clusters, build_dge, differences

RSCRIPT = os.environ.get("LLOMICS_RSCRIPT") or shutil.which("Rscript")


def r_has_llomics() -> bool:
    if not RSCRIPT:
        return False
    probe = subprocess.run([RSCRIPT, "-e", "cat(requireNamespace('llomics', quietly = TRUE))"],
                           capture_output=True, text=True)
    return probe.stdout.strip().endswith("TRUE")


@unittest.skipUnless(r_has_llomics(), "Rscript with the llomics package is not available")
class Parity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="llomics-parity-")
        result = subprocess.run(
            [RSCRIPT, str(ROOT / "tests" / "parity" / "build_from_fixtures.R"), str(FIXTURES), cls.tmp],
            capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError("R build failed:\n" + result.stderr[-3000:])

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def from_r(self, name: str):
        return json.loads((Path(self.tmp) / name).read_text(encoding="utf-8"))

    def test_dge_payloads_match(self):
        self.assertEqual(differences(build_dge(), self.from_r("dge.json")), [])

    def test_clusters_payloads_match(self):
        self.assertEqual(differences(build_clusters(), self.from_r("clusters.json")), [])


if __name__ == "__main__":
    unittest.main()
