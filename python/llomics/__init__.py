"""llomics: turn omics pipeline results into small, self-describing payloads for LLMs.

    import llomics
    payload = llomics.rank_genes_groups(adata, "leiden", n_features=15, panel=["CD3D", "MS4A1"],
                                        question="Which cell type is each cluster?", out="payload.json")

Layout: ``core`` turns plain tables into payload units (standard library only), ``adapters`` map
scanpy and other tools onto those tables, ``validate`` checks a payload against the contract.
"""

from . import core
from .contract import CONTRACT, VERSION
from .core import clusters_payload, dge_payload, write_payload
from .validate import validate

__version__ = "0.1.0"
__all__ = ["CONTRACT", "VERSION", "clusters_payload", "core", "dge_payload", "rank_genes_groups", "validate",
           "write_payload"]


def rank_genes_groups(*args, **kwargs):
    """Run ``sc.tl.rank_genes_groups`` and build the payload. See ``llomics.adapters.scanpy``."""
    from .adapters.scanpy import rank_genes_groups as _impl

    return _impl(*args, **kwargs)
