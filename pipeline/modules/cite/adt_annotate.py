"""
adt_annotate.py — ADT Leiden clustering and marker-based annotation
for the OmicSage CITE-seq pipeline.

Takes mdata["adt"] (output of adt_harmony.py) and performs:
  1. Leiden clustering at low resolution on the harmony-corrected neighbor
     graph (obsp["connectivities"] / obsp["distances"])
  2. rank_genes_groups on Leiden clusters (for marker dotplot support)
  3. dendrogram on Leiden clusters (for ordered dotplot)
  4. Optional: maps cluster IDs → cell type labels stored in
     obs["adt_celltype"] when annotation_map is provided

OmicSage naming conventions enforced:
  obs["leiden"]        — Leiden cluster IDs (always written)
  obs["adt_celltype"]  — cell type labels (written only when annotation_map
                         is supplied; avoids collision with RNA pipeline's
                         obs["cell_type_vote"])

Reference:
  https://www.sc-best-practices.org/surface_protein/annotation.html

  Key calls from the reference notebook (ch.39):
    sc.tl.leiden(mdata["prot"], resolution=0.1, flavor="igraph",
                 n_iterations=2, directed=False, random_state=0)
    sc.tl.rank_genes_groups(mdata["prot"], groupby="leiden")
    sc.tl.dendrogram(mdata["prot"], groupby="leiden")
    mdata["prot"].obs["celltype"] = mdata["prot"].obs.leiden.copy()
    mdata["prot"].obs.celltype.replace({...}, inplace=True)

  OmicSage change: obs["celltype"] → obs["adt_celltype"] to avoid
  collision with RNA annotation column obs["cell_type_vote"].

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

Prerequisite pipeline order:
  adt_normalize.py → adt_doublets.py → adt_reduce.py → adt_harmony.py
  → adt_annotate.py (this module) → wnn.py
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Optional

import scanpy as sc
from anndata import AnnData
from mudata import MuData


# ---------------------------------------------------------------------------
# Constants (sc-best-practices ch.39 defaults)
# ---------------------------------------------------------------------------
_DEFAULT_RESOLUTION   = 0.1   # low resolution → broad immune cell types
_DEFAULT_N_ITERATIONS = 2     # igraph flavor uses n_iterations
_DEFAULT_RANDOM_STATE = 0


def annotate_adt(
    mdata: MuData,
    annotation_map: Optional[dict[str, str]] = None,
    resolution: float = _DEFAULT_RESOLUTION,
    n_iterations: int = _DEFAULT_N_ITERATIONS,
    random_state: int = _DEFAULT_RANDOM_STATE,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Leiden clustering and optional marker-based annotation on ADT data.

    Runs Leiden clustering on the harmony-corrected neighbor graph already
    present in ``mdata["adt"].obsp``, computes ``rank_genes_groups`` and a
    dendrogram for dotplot support, and optionally maps cluster IDs to cell
    type labels.

    Parameters
    ----------
    mdata : MuData
        MuData object containing ``mdata["adt"]`` with:
          - ``obsp["connectivities"]`` and ``obsp["distances"]`` — neighbor
            graph produced by ``adt_harmony.py``
          - ``obsm["X_umap_adt"]`` — harmony-corrected UMAP (used for
            provenance; not required for Leiden itself)
    annotation_map : dict[str, str] or None
        Mapping from Leiden cluster ID strings to cell type label strings,
        e.g. ``{"0": "CD4 T", "1": "B", "2": "CD8 T"}``.
        If ``None`` (default), clustering is performed but no cell type
        labels are assigned; ``obs["adt_celltype"]`` is not written.
        Cluster IDs not present in the map are left as their numeric string.
    resolution : float
        Leiden resolution parameter.  Lower values → fewer, broader clusters.
        Default: 0.1 (sc-best-practices ch.39).
    n_iterations : int
        Number of Leiden iterations (used with ``flavor="igraph"``).
        Default: 2 (sc-best-practices ch.39).
    random_state : int
        Random seed for Leiden.  Default: 0.
    inplace : bool
        If ``True``, modify ``mdata["adt"]`` in place.
        If ``False`` (default), operate on a copy.

    Returns
    -------
    adata : AnnData
        ``mdata["adt"]`` (copy or in-place) with:
          ``obs["leiden"]``                — Leiden cluster IDs (always)
          ``obs["adt_celltype"]``          — cell type labels (annotation_map
                                            provided only)
          ``uns["rank_genes_groups"]``     — marker gene results
          ``uns["dendrogram_leiden"]``     — dendrogram for dotplot
          ``uns["omicsage_adt_annotate"]`` — provenance dict
    metrics : dict
        Summary metrics — see Notes.

    Raises
    ------
    KeyError
        If ``mdata["adt"]`` is missing, or if the neighbor graph
        (``obsp["connectivities"]``) is absent.
    ValueError
        If ``annotation_map`` contains keys that are not valid cluster ID
        strings (i.e. keys that do not appear in ``obs["leiden"]`` after
        clustering).  A warning is issued instead of an error so that
        partially-matching maps still work.

    Notes
    -----
    Metrics dict keys:
      n_cells           int    number of cells processed
      n_clusters        int    number of Leiden clusters found
      cluster_sizes     dict   {cluster_id_str: cell_count} for every cluster
      annotation_map    dict   the map applied (copy), or {} if None
      annotated         bool   True when annotation_map was applied
      leiden_key        str    always "leiden"
      celltype_key      str    "adt_celltype" if annotated, else None
      resolution        float  Leiden resolution used
      random_state      int    seed used
    """
    # ------------------------------------------------------------------
    # 1. Validate inputs
    # ------------------------------------------------------------------
    _validate_inputs(mdata)

    # ------------------------------------------------------------------
    # 2. Copy or in-place
    # ------------------------------------------------------------------
    adata = mdata["adt"] if inplace else mdata["adt"].copy()

    # ------------------------------------------------------------------
    # 3. Leiden clustering
    #    flavor="igraph" + directed=False matches sc-best-practices ch.39
    #    The neighbor graph from adt_harmony.py is already in obsp.
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            n_iterations=n_iterations,
            directed=False,
            random_state=random_state,
        )

    # ------------------------------------------------------------------
    # 4. rank_genes_groups — marker genes per cluster (for dotplot)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 5. Dendrogram — ordered cluster relationships (for dotplot)
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.dendrogram(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 6. Optional annotation: cluster ID → cell type label
    # ------------------------------------------------------------------
    annotated = False
    applied_map: dict[str, str] = {}

    if annotation_map is not None:
        # Start from the leiden column
        adata.obs["adt_celltype"] = adata.obs["leiden"].copy()

        # Identify any keys in annotation_map that do not match any cluster
        found_clusters = set(adata.obs["leiden"].unique().tolist())
        unknown_keys = set(str(k) for k in annotation_map.keys()) - found_clusters
        if unknown_keys:
            warnings.warn(
                f"annotation_map contains keys not found in Leiden clusters: "
                f"{sorted(unknown_keys)}. These keys will be ignored.",
                UserWarning,
                stacklevel=2,
            )

        # Apply only keys that exist in the cluster set
        applied_map = {
            str(k): str(v)
            for k, v in annotation_map.items()
            if str(k) in found_clusters
        }
        if applied_map:
            adata.obs["adt_celltype"] = adata.obs["adt_celltype"].astype(str).replace(applied_map)

        annotated = True

    # ------------------------------------------------------------------
    # 7. Metrics
    # ------------------------------------------------------------------
    cluster_sizes: dict[str, int] = (
        adata.obs["leiden"].value_counts().to_dict()
    )
    # Ensure keys are plain strings (pandas may return category dtype)
    cluster_sizes = {str(k): int(v) for k, v in cluster_sizes.items()}

    n_clusters = len(cluster_sizes)

    metrics: dict = {
        "n_cells": int(adata.n_obs),
        "n_clusters": n_clusters,
        "cluster_sizes": cluster_sizes,
        "annotation_map": applied_map,
        "annotated": annotated,
        "leiden_key": "leiden",
        "celltype_key": "adt_celltype" if annotated else None,
        "resolution": resolution,
        "random_state": random_state,
    }

    # ------------------------------------------------------------------
    # 8. Provenance
    # ------------------------------------------------------------------
    adata.uns["omicsage_adt_annotate"] = {
        "module": "adt_annotate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "resolution": resolution,
            "flavor": "igraph",
            "n_iterations": n_iterations,
            "directed": False,
            "random_state": random_state,
            "annotation_map_provided": annotation_map is not None,
        },
        "outputs": {
            "leiden_key": "leiden",
            "celltype_key": "adt_celltype" if annotated else None,
            "n_clusters": n_clusters,
        },
    }

    return adata, metrics


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_inputs(mdata: MuData) -> None:
    """Raise informative errors for common input mistakes."""
    if "adt" not in mdata.mod:
        raise KeyError(
            "mdata must contain an 'adt' modality. "
            "Found: " + str(list(mdata.mod.keys()))
        )

    adata = mdata["adt"]

    if "connectivities" not in adata.obsp:
        raise KeyError(
            "mdata['adt'].obsp['connectivities'] not found. "
            "Run run_harmony_adt() (or reduce_adt()) before annotate_adt() "
            "to build the neighbor graph."
        )
