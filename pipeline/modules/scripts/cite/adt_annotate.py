"""
adt_annotate.py — ADT Leiden clustering and cell type annotation for OmicSage.

Runs Leiden clustering on the post-Harmony neighbor graph, computes
rank_genes_groups and a dendrogram, then produces two independent annotation
columns:

  obs["adt_celltype_score"]   — auto-scored labels from marker panel
                                 fold-change scoring (always written when a
                                 marker_panel or preset is provided)
  obs["adt_celltype_manual"]  — manual labels from an explicit annotation_map
                                 (written only when annotation_map is given)

Having two separate columns lets you compare the scoring result against your
manual map in the report side by side, and iterate on the map without losing
the auto-scored baseline.

Reference:
  https://www.sc-best-practices.org/surface_protein/annotation.html

Input
-----
mdata["adt"].obsp["connectivities"]   — post-Harmony kNN graph
mdata["adt"].obsp["distances"]        — post-Harmony kNN distances
                                        (both written by adt_harmony.py)
mdata["adt"].layers["adt_clr"]        — CLR-normalised ADT (for scoring)

Outputs added to mdata["adt"]
------------------------------
obs["leiden"]                str  Leiden cluster ID for every cell (always)
obs["adt_celltype_score"]    str  Auto-scored cell type per cluster (when
                                   marker_panel or preset is given)
obs["adt_celltype_manual"]   str  Manual cell type from annotation_map (when
                                   annotation_map is given)
uns["rank_genes_groups"]          Marker proteins per cluster
uns["dendrogram_leiden"]          Cluster dendrogram (for dotplot ordering)
uns["omicsage_adt_marker_panel"]  Marker panel used for scoring (if any)
uns["omicsage_adt_annotate"]      Provenance

Naming convention
-----------------
Labels go into obs["adt_celltype_score"] / obs["adt_celltype_manual"] —
never obs["celltype"] or obs["cell_type"] or obs["adt_celltype"].
This prevents collision with the RNA pipeline's obs["cell_type_vote"] when
both AnnData objects are joined into a MuData for integration.

API
---
annotate_adt(
    mdata,
    annotation_map=None,
    marker_panel=None,
    preset=None,
    resolution=0.1,
    n_iterations=2,
    random_state=0,
    inplace=False,
)
-> (AnnData, dict)
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
# Built-in marker panels
# ---------------------------------------------------------------------------

_PRESET_PANELS: dict[str, dict[str, list[str]]] = {
    "bmmc": {
        # Based on BioLegend TotalSeq-A Universal Cocktail V1.0 (GSE194122).
        # Keys = cell type labels; values = positive-marker protein prefixes
        # (case-insensitive prefix match against var_names).
        "CD4 T":     ["CD4",   "CD5",   "CD3"],
        "CD8 T":     ["CD8a",  "CD5",   "CD3"],
        "NK":        ["CD56",  "CD16",  "CD45RA"],
        "B":         ["CD19",  "CD20",  "CD72"],
        "Mono CD14": ["CD14",  "CD11c", "CD172a"],
        "Mono CD16": ["CD16",  "CD11c", "CD172a"],
        "DC":        ["CD123", "CD304", "CD1c"],
        "Plasma":    ["CD38",  "CD138"],
        "Erythroid": ["CD71",  "CD36"],
        "Platelet":  ["CD41",  "CD62P"],
    },
}


# ---------------------------------------------------------------------------
# Scoring helper
# ---------------------------------------------------------------------------

def _score_clusters(
    adata: AnnData,
    marker_panel: dict[str, list[str]],
    leiden_clusters: list[str],
) -> np.ndarray:
    """
    Assign a cell type label to each cell by scoring every Leiden cluster
    against the marker panel.

    Algorithm (per cluster):
      1. For each cell type in the panel, compute the mean CLR expression of
         its positive markers across all cells in the cluster.
      2. Normalise each protein by subtracting the panel-wide median of that
         protein (fold-change above background, not absolute expression).
      3. Sum the normalised means across all markers for that cell type.
      4. Assign the cell type with the highest total score to the cluster.
         Ties and "no markers resolved" fall back to "Unknown".

    Returns
    -------
    np.ndarray of str, shape (n_cells,)
        One label per cell; all cells in a cluster share the same label.
    """
    # Prefer CLR layer; fall back to .X
    if "adt_clr" in adata.layers:
        X = adata.layers["adt_clr"]
    else:
        X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X, dtype=float)

    var_names = list(adata.var_names)

    # Panel-wide protein medians for background normalisation
    protein_median: dict[str, float] = {
        vn: float(np.median(X[:, i])) for i, vn in enumerate(var_names)
    }

    def _resolve(name: str) -> Optional[int]:
        """First var_name that starts with `name` (case-insensitive)."""
        name_up = name.upper()
        for idx, vn in enumerate(var_names):
            if vn.upper().startswith(name_up):
                return idx
        return None

    # Resolve panel marker names to column indices once
    resolved_panel: dict[str, list[int]] = {}
    for ct, markers in marker_panel.items():
        idxs = [_resolve(m) for m in markers]
        idxs = [i for i in idxs if i is not None]
        if idxs:
            resolved_panel[ct] = idxs

    if not resolved_panel:
        return np.full(adata.n_obs, "Unknown", dtype=object)

    leiden_col = adata.obs["leiden"].astype(str).values
    labels = np.full(adata.n_obs, "Unknown", dtype=object)

    for cluster_id in leiden_clusters:
        cell_mask = leiden_col == cluster_id
        if not cell_mask.any():
            continue

        cluster_X = X[cell_mask, :]

        best_ct    = "Unknown"
        best_score = -np.inf

        for ct, idxs in resolved_panel.items():
            score = sum(
                float(cluster_X[:, idx].mean()) - protein_median[var_names[idx]]
                for idx in idxs
            )
            if score > best_score:
                best_score = score
                best_ct    = ct

        labels[cell_mask] = best_ct

    return labels


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def annotate_adt(
    mdata: MuData,
    annotation_map: Optional[dict] = None,
    marker_panel: Optional[dict[str, list[str]]] = None,
    preset: Optional[str] = None,
    resolution: float = 0.1,
    n_iterations: int = 2,
    random_state: int = 0,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Run Leiden clustering, auto-scoring, and optional manual annotation on
    ADT data.

    Parameters
    ----------
    mdata : MuData
        Must contain mdata["adt"] with obsp["connectivities"] present.
    annotation_map : dict[str, str] or None
        Maps Leiden cluster ID strings -> cell type labels.
        Integer keys are coerced to str automatically.
        Unknown keys emit UserWarning but do not raise.
        Writes obs["adt_celltype_manual"] and sets metrics["annotated"]=True.
        When None: no manual column is written.
        When {}: obs["adt_celltype_manual"] is written as a copy of leiden.
    marker_panel : dict[str, list[str]] or None
        Custom panel for auto-scoring. Keys = cell type labels; values = lists
        of protein name prefixes (case-insensitive prefix match).
        Writes obs["adt_celltype_score"]. Takes precedence over preset.
    preset : str or None
        Built-in panel for auto-scoring. Currently: "bmmc".
        Ignored when marker_panel is also provided.
        Writes obs["adt_celltype_score"].
    resolution : float
        Leiden resolution. Default: 0.1.
    n_iterations : int
        Leiden iterations (igraph flavor). Default: 2.
    random_state : int
        Leiden random seed. Default: 0.
    inplace : bool
        Modify mdata["adt"] in place if True; default makes a copy.

    Returns
    -------
    adata : AnnData
        mdata["adt"] with:
          obs["leiden"]               -- Leiden cluster IDs (always)
          obs["adt_celltype_score"]   -- auto-scored labels (if panel/preset)
          obs["adt_celltype_manual"]  -- manual labels (if annotation_map)
          uns["rank_genes_groups"]    -- marker proteins per cluster
          uns["dendrogram_leiden"]    -- cluster dendrogram
          uns["omicsage_adt_marker_panel"] -- panel used for scoring
          uns["omicsage_adt_annotate"]     -- provenance
    metrics : dict
        n_cells, n_clusters, cluster_sizes, annotation_map, annotated,
        leiden_key, celltype_key, score_key, marker_panel,
        resolution, random_state.

    Raises
    ------
    KeyError
        If mdata["adt"] or obsp["connectivities"] are missing.
    """
    # ------------------------------------------------------------------
    # 1. Validate
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
    # 3. Resolve marker panel  (custom > preset > None)
    # ------------------------------------------------------------------
    active_panel: Optional[dict[str, list[str]]] = None
    if marker_panel is not None:
        active_panel = marker_panel
    elif preset is not None:
        key = preset.lower().strip()
        if key not in _PRESET_PANELS:
            warnings.warn(
                f"Unknown preset '{preset}'. "
                f"Available presets: {list(_PRESET_PANELS.keys())}. "
                "Scoring will be skipped.",
                UserWarning, stacklevel=2,
            )
        else:
            active_panel = _PRESET_PANELS[key]

    # ------------------------------------------------------------------
    # 4. Normalise annotation_map keys to str
    # ------------------------------------------------------------------
    if annotation_map is not None:
        annotation_map = {str(k): str(v) for k, v in annotation_map.items()}

    # ------------------------------------------------------------------
    # 5. Leiden clustering
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

    adata.obs["leiden"] = adata.obs["leiden"].astype(str)
    leiden_clusters = sorted(adata.obs["leiden"].unique().tolist())
    n_clusters = len(leiden_clusters)

    # ------------------------------------------------------------------
    # 6. rank_genes_groups
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # ------------------------------------------------------------------
    # 7. Dendrogram (requires >= 2 groups)
    # ------------------------------------------------------------------
    if n_clusters >= 2:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.dendrogram(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 8. Auto-scoring -> obs["adt_celltype_score"]
    # ------------------------------------------------------------------
    score_key: Optional[str] = None
    if active_panel is not None:
        score_labels = _score_clusters(adata, active_panel, leiden_clusters)
        adata.obs["adt_celltype_score"] = score_labels.astype(str)
        score_key = "adt_celltype_score"
        adata.uns["omicsage_adt_marker_panel"] = active_panel

    # ------------------------------------------------------------------
    # 9. Manual annotation -> obs["adt_celltype_manual"]
    # ------------------------------------------------------------------
    annotated = annotation_map is not None
    celltype_key: Optional[str] = None

    if annotated:
        unknown_keys = set(annotation_map.keys()) - set(leiden_clusters)
        if unknown_keys:
            warnings.warn(
                f"annotation_map keys not found in Leiden clusters and will be "
                f"ignored: {sorted(unknown_keys)}. "
                f"Valid cluster IDs are: {leiden_clusters}",
                UserWarning, stacklevel=2,
            )

        adata.obs["adt_celltype_manual"] = (
            adata.obs["leiden"]
            .map(lambda cid: annotation_map.get(cid, cid))
            .astype(str)
        )
        celltype_key = "adt_celltype_manual"

    # ------------------------------------------------------------------
    # 10. cluster_sizes
    # ------------------------------------------------------------------
    cluster_sizes: dict[str, int] = {
        cid: int((adata.obs["leiden"] == cid).sum())
        for cid in leiden_clusters
    }

    # ------------------------------------------------------------------
    # 11. Provenance
    # ------------------------------------------------------------------
    adata.uns["omicsage_adt_annotate"] = {
        "module": "adt_annotate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "resolution": resolution,
            "flavor": "igraph",
            "directed": False,
            "n_iterations": n_iterations,
            "random_state": random_state,
            "annotation_map_provided": annotation_map is not None,
            "preset": preset,
            "marker_panel_provided": marker_panel is not None,
        },
        "outputs": {
            "leiden_key": "leiden",
            "n_clusters": n_clusters,
            "celltype_key": celltype_key,
            "score_key": score_key,
        },
    }

    # ------------------------------------------------------------------
    # 12. Metrics
    # ------------------------------------------------------------------
    metrics: dict = {
        "n_cells": int(adata.n_obs),
        "n_clusters": n_clusters,
        "cluster_sizes": cluster_sizes,
        "annotation_map": annotation_map if annotation_map is not None else {},
        "annotated": annotated,
        "leiden_key": "leiden",
        "celltype_key": celltype_key,
        "score_key": score_key,
        "marker_panel": active_panel,
        "resolution": resolution,
        "random_state": random_state,
    }

    return adata, metrics
