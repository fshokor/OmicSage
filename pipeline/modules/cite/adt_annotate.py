"""
adt_annotate.py — ADT Leiden clustering and cell type annotation for OmicSage.

Runs Leiden clustering on the post-Harmony neighbor graph, computes
rank_genes_groups and a dendrogram, and optionally maps cluster IDs to
cell type labels via an annotation_map.

Reference:
  https://www.sc-best-practices.org/surface_protein/annotation.html

Input
-----
mdata["adt"].obsp["connectivities"]   — post-Harmony kNN graph
mdata["adt"].obsp["distances"]        — post-Harmony kNN distances
                                        (both written by adt_harmony.py)

Outputs added to mdata["adt"]
------------------------------
obs["leiden"]               str  Leiden cluster ID for every cell (always)
obs["adt_celltype"]         str  Cell type label (only when annotation_map is given)
uns["rank_genes_groups"]         Marker proteins per cluster
uns["dendrogram_leiden"]         Cluster dendrogram (for dotplot ordering)
uns["omicsage_adt_annotate"]     Provenance

Naming convention
-----------------
Labels go into obs["adt_celltype"] — never obs["celltype"] or obs["cell_type"].
This prevents collision with the RNA pipeline's obs["cell_type_vote"] when both
AnnData objects are joined into a MuData for integration.

API
---
annotate_adt(
    mdata,
    annotation_map=None,
    resolution=0.1,
    n_iterations=2,
    random_state=0,
    inplace=False,
)
→ (AnnData, dict)
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from mudata import MuData


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def annotate_adt(
    mdata: MuData,
    annotation_map: Optional[dict] = None,
    resolution: float = 0.1,
    n_iterations: int = 2,
    random_state: int = 0,
    inplace: bool = False,
    preset: Optional[str] = None,   # reserved for future use; ignored
) -> tuple[AnnData, dict]:
    """
    Run Leiden clustering and optional cell type annotation on ADT data.

    Parameters
    ----------
    mdata : MuData
        Must contain ``mdata["adt"]`` with ``obsp["connectivities"]`` present
        (output of ``adt_harmony.py`` or ``adt_reduce.py``).
    annotation_map : dict[str, str] or None
        Maps Leiden cluster ID strings to cell type labels.
        Integer keys are coerced to str automatically.
        Unknown keys (not in the actual Leiden output) emit a UserWarning
        but do not raise.
        When ``None`` (default), clustering is performed but no labels are
        written and ``metrics["annotated"]`` is ``False``.
        When ``{}`` (empty dict), ``obs["adt_celltype"]`` is written as a
        copy of ``obs["leiden"]`` and ``metrics["annotated"]`` is ``True``.
    resolution : float
        Leiden resolution.  Lower = fewer, broader clusters.  Default: 0.1.
    n_iterations : int
        Leiden iterations (igraph flavor).  Default: 2.
    random_state : int
        Leiden random seed.  Default: 0.
    inplace : bool
        If True, modify ``mdata["adt"]`` in place and return the same object.
        If False (default), operate on a copy.

    Returns
    -------
    adata : AnnData
        ``mdata["adt"]`` (copy or in-place) with clustering applied:
          ``obs["leiden"]``                — Leiden cluster IDs (always)
          ``obs["adt_celltype"]``          — labels (if annotation_map given)
          ``uns["rank_genes_groups"]``     — marker proteins
          ``uns["dendrogram_leiden"]``     — cluster dendrogram
          ``uns["omicsage_adt_annotate"]`` — provenance
    metrics : dict
        n_cells, n_clusters, cluster_sizes, annotation_map (applied, as dict),
        annotated, leiden_key, celltype_key, resolution, random_state.

    Raises
    ------
    KeyError
        If ``mdata["adt"]`` does not exist or ``obsp["connectivities"]``
        is missing.
    """
    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    if "adt" not in mdata.mod:
        raise KeyError(
            "'adt' modality not found in mdata. "
            f"Available modalities: {list(mdata.mod.keys())}"
        )

    adt_src = mdata["adt"]

    if "connectivities" not in adt_src.obsp:
        raise KeyError(
            "mdata['adt'].obsp['connectivities'] not found. "
            "Run adt_harmony.py (or adt_reduce.py) first to build the "
            "neighbor graph."
        )

    # ------------------------------------------------------------------
    # 2. Copy or in-place
    # ------------------------------------------------------------------
    adata = adt_src if inplace else adt_src.copy()

    # ------------------------------------------------------------------
    # 3. Normalise annotation_map keys to str
    # ------------------------------------------------------------------
    if annotation_map is not None:
        annotation_map = {str(k): str(v) for k, v in annotation_map.items()}

    # ------------------------------------------------------------------
    # 4. Leiden clustering
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            directed=False,
            n_iterations=n_iterations,
            random_state=random_state,
        )

    # Leiden writes as Categorical — cast to plain str for contract consistency
    adata.obs["leiden"] = adata.obs["leiden"].astype(str)

    leiden_clusters = sorted(adata.obs["leiden"].unique().tolist())
    n_clusters = len(leiden_clusters)

    # ------------------------------------------------------------------
    # 5. rank_genes_groups (marker proteins per cluster)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(
            adata,
            groupby="leiden",
            method="wilcoxon",
        )

    # ------------------------------------------------------------------
    # 6. Dendrogram (requires >= 2 groups)
    # ------------------------------------------------------------------
    if n_clusters >= 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.dendrogram(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 7. Apply annotation_map → obs["adt_celltype"]
    # ------------------------------------------------------------------
    annotated = annotation_map is not None
    celltype_key: Optional[str] = None

    if annotated:
        # Warn about unknown keys
        unknown_keys = set(annotation_map.keys()) - set(leiden_clusters)
        if unknown_keys:
            warnings.warn(
                f"annotation_map keys not found in Leiden clusters and will be "
                f"ignored: {sorted(unknown_keys)}. "
                f"Valid cluster IDs are: {leiden_clusters}",
                UserWarning,
                stacklevel=2,
            )

        # Map clusters → labels; unmapped clusters keep their numeric string ID
        adata.obs["adt_celltype"] = (
            adata.obs["leiden"]
            .map(lambda cid: annotation_map.get(cid, cid))
            .astype(str)
        )
        celltype_key = "adt_celltype"

    # ------------------------------------------------------------------
    # 8. cluster_sizes
    # ------------------------------------------------------------------
    cluster_sizes: dict[str, int] = {
        cid: int((adata.obs["leiden"] == cid).sum())
        for cid in leiden_clusters
    }

    # ------------------------------------------------------------------
    # 9. Provenance
    # ------------------------------------------------------------------
    provenance = {
        "module": "adt_annotate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "resolution": resolution,
            "flavor": "igraph",
            "directed": False,
            "n_iterations": n_iterations,
            "random_state": random_state,
            "annotation_map_provided": annotation_map is not None,
        },
        "outputs": {
            "leiden_key": "leiden",
            "n_clusters": n_clusters,
            "celltype_key": celltype_key,
        },
    }
    adata.uns["omicsage_adt_annotate"] = provenance

    # ------------------------------------------------------------------
    # 10. Metrics
    # ------------------------------------------------------------------
    metrics: dict = {
        "n_cells": int(adata.n_obs),
        "n_clusters": n_clusters,
        "cluster_sizes": cluster_sizes,
        "annotation_map": annotation_map if annotation_map is not None else {},
        "annotated": annotated,
        "leiden_key": "leiden",
        "celltype_key": celltype_key,
        "resolution": resolution,
        "random_state": random_state,
    }

    return adata, metrics
