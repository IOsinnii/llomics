"""llomics: turn omics pipeline results into small, self-describing payloads for LLMs.

    import llomics
    payload = llomics.rank_genes_groups(adata, "leiden", n_features=15, panel=["CD3D", "MS4A1"],
                                        question="Which cell type is each cluster?", out="payload.json")
"""

from ._core import CONTRACT_VERSION, write_payload
from ._scanpy import rank_genes_groups

__version__ = "0.1.0"
__all__ = ["CONTRACT_VERSION", "rank_genes_groups", "write_payload"]
