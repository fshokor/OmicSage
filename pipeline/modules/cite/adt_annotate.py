"""
adt_annotate.py — ADT Leiden clustering and marker-based annotation
for the OmicSage CITE-seq pipeline.

Takes mdata["adt"] (output of adt_harmony.py) and performs:
  1. Leiden clustering at low resolution on the harmony-corrected neighbor
     graph (obsp["connectivities"] / obsp["distances"])
  2. rank_genes_groups on Leiden clusters (for marker dotplot support)
  3. dendrogram on Leiden clusters (for ordered dotplot)
  4. Auto-generates obs["adt_celltype"] by scoring each cluster against a
     marker panel: the cell type whose markers show the highest mean
     expression in a cluster wins (plurality vote).

OmicSage naming conventions enforced:
  obs["leiden"]        — Leiden cluster IDs (always written)
  obs["adt_celltype"]  — cell type labels (written when a panel is active)

Built-in preset
---------------
``preset="bmmc"``
    Uses ``BMMC_MARKER_PANEL`` — lineage markers for the NeurIPS 2021 BMMC
    CITE-seq benchmark (GSE194122), BioLegend TotalSeq-A Human Universal
    Cocktail V1.0 (~134 proteins after QC).
    Cell type names follow sc-best-practices ch.39 naming.

    The annotation map is **derived from the data at runtime**: for each
    Leiden cluster, the cell type whose panel markers have the highest
    mean DSB/CLR expression in that cluster is assigned.  Works for any
    number of clusters.

Custom workflow
---------------
Supply your own ``marker_panel``::

    marker_panel = {
        "CD4 T":   ["CD3", "CD4"],
        "B":       ["CD19", "CD20"],
        "Monocyte":["CD14", "HLA-DR"],
    }
    annotate_adt(mdata, marker_panel=marker_panel)

Or supply an explicit ``annotation_map`` to bypass scoring entirely::

    annotation_map = {"0": "CD4 T", "1": "B", "2": "Monocyte"}
    annotate_adt(mdata, annotation_map=annotation_map)

Reference:
  https://www.sc-best-practices.org/surface_protein/annotation.html

API
---
annotate_adt(
    mdata,
    preset=None,
    annotation_map=None,
    marker_panel=None,
    resolution=0.1,
    n_iterations=2,
    random_state=0,
    inplace=False,
)
→ (AnnData, dict)

Prerequisite pipeline order:
  adt_normalize.py → adt_doublets.py → adt_reduce.py → adt_harmony.py
  → adt_annotate.py (this module) → cite_integration.py
"""

from __future__ import annotations

import warnings
from datetime import datetime, timezone
from typing import Literal, Optional

import numpy as np
import scanpy as sc
from anndata import AnnData
from mudata import MuData


# ---------------------------------------------------------------------------
# Built-in preset: NeurIPS 2021 BMMC CITE-seq (GSE194122)
# Protein names match var_names in the processed AnnData (no prefix).
# Only marker proteins are listed — isotype controls excluded.
# ---------------------------------------------------------------------------

BMMC_MARKER_PANEL: dict[str, list[str]] = {
    "CD4 T":     ["CD3", "CD4", "CD45RA", "CD45RO", "CD27", "CD28",
                  "CD127", "CD197"],
    "CD8 T":     ["CD3", "CD8a", "CD45RA", "CD45RO", "CD27", "CD57"],
    "Treg":      ["CD3", "CD4", "CD25", "CD127"],
    "NK":        ["CD56", "CD16", "CD3", "CD94", "CD161"],
    "B":         ["CD19", "CD20", "CD21", "CD22", "CD24", "IgD"],
    "Plasma":    ["CD19", "CD38", "CD27", "IgM"],
    "CD14 Mono": ["CD14", "CD16", "HLA-DR", "CD11b", "CD33"],
    "CD16 Mono": ["CD14", "CD16", "HLA-DR", "CD11c"],
    "DC":        ["HLA-DR", "CD11c", "CD141", "CD1c", "CD303"],
    "Progenitor":["CD34", "CD90", "CD117", "CD135", "CD133"],
    "Platelet":  ["CD41", "CD61", "CD42b"],
    # Added — clusters 0+1 in GSE194122 BMMC are erythroid progenitors
    # (~30% of cells); absent from original panel caused mislabelling.
    # CD235a may appear as "CD235a" or "GYPA" depending on feature reference.
    "Erythroid": ["CD71", "CD36", "CD235a", "CD88"],
    # HSPC — haematopoietic stem and progenitor cells
    "HSPC":      ["CD34", "CD38", "CD90", "CD117", "CD133", "CD135"],
}

# Registry — extend with new presets here
_PRESETS: dict[str, dict[str, list[str]]] = {
    "bmmc": BMMC_MARKER_PANEL,
}

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_DEFAULT_RESOLUTION   = 0.1
_DEFAULT_N_ITERATIONS = 2
_DEFAULT_RANDOM_STATE = 0


def annotate_adt(
    mdata: MuData,
    preset: Optional[Literal["bmmc"]] = None,
    annotation_map: Optional[dict[str, str]] = None,
    marker_panel: Optional[dict[str, list[str]]] = None,
    resolution: float = _DEFAULT_RESOLUTION,
    n_iterations: int = _DEFAULT_N_ITERATIONS,
    random_state: int = _DEFAULT_RANDOM_STATE,
    inplace: bool = False,
) -> tuple[AnnData, dict]:
    """
    Leiden clustering and marker-based annotation on ADT data.

    Parameters
    ----------
    mdata : MuData
        MuData object containing ``mdata["adt"]`` with:
          - ``obsp["connectivities"]`` / ``obsp["distances"]`` — from
            ``adt_harmony.py``
          - ``obsm["X_umap_adt"]`` — for provenance
    preset : {"bmmc"} or None
        Load a built-in marker panel and derive the annotation map from
        the data.  ``None`` means no preset is used.
        Explicit ``marker_panel`` overrides the preset panel.
        Explicit ``annotation_map`` bypasses scoring entirely.
    annotation_map : dict[str, str] or None
        Hard override: cluster ID → cell type label.
        When provided, scoring is skipped and this map is applied directly.
        Overrides both preset and marker_panel-derived maps.
    marker_panel : dict[str, list[str]] or None
        Cell type → list of marker protein names.
        Used to score clusters and derive the annotation map automatically.
        Overrides the preset panel when both are supplied.
        If neither preset nor marker_panel is given, clustering runs but
        ``obs["adt_celltype"]`` is not written.
    resolution : float
        Leiden resolution.  Default: 0.1.
    n_iterations : int
        Leiden iterations (igraph flavor).  Default: 2.
    random_state : int
        Random seed.  Default: 0.
    inplace : bool
        Modify ``mdata["adt"]`` in place.  Default: False.

    Returns
    -------
    adata : AnnData
        ``mdata["adt"]`` with:
          ``obs["leiden"]``                    — Leiden cluster IDs
          ``obs["adt_celltype_score"]``        — auto-scored labels (always, when panel active)
          ``obs["adt_celltype_manual"]``       — manual map labels (only when annotation_map provided)
          ``uns["rank_genes_groups"]``         — marker results
          ``uns["dendrogram_leiden"]``         — dendrogram
          ``uns["omicsage_adt_marker_panel"]`` — resolved marker panel (if any)
          ``uns["omicsage_adt_annotate"]``     — provenance dict
    metrics : dict
        Summary metrics — see Notes.

    Raises
    ------
    KeyError
        If ``mdata["adt"]`` is missing or neighbor graph is absent.
    ValueError
        If ``preset`` is not a recognised preset name.

    Notes
    -----
    Auto-annotation scoring
    ~~~~~~~~~~~~~~~~~~~~~~~
    For each Leiden cluster *c* and each cell type *t* in the marker panel:

    1. Collect the mean expression (DSB → CLR → .X, in priority order)
       of each marker protein listed for *t* that exists in ``var_names``.
    2. Compute the mean of those per-protein means → a scalar score S(c, t).
    3. Assign ``adt_celltype[c] = argmax_t S(c, t)``.

    Ties are broken by alphabetical order of cell type name.  Clusters
    where no marker protein is found in ``var_names`` are labelled
    ``"Unknown"``.

    Metrics dict keys
    ~~~~~~~~~~~~~~~~~
      n_cells              int    cells processed
      n_clusters           int    Leiden clusters found
      cluster_sizes        dict   {cluster_id: cell_count}
      annotation_map       dict   the map applied, or {} if no annotation
      annotation_source    str    "explicit" | "scoring" | "none"
      annotated            bool   True when adt_celltype was written
      leiden_key           str    "leiden"
      celltype_key         str    "adt_celltype" if annotated, else None
      resolution           float  Leiden resolution used
      random_state         int    seed used
      marker_panel         dict   resolved panel, or {}
      preset               str    preset name, or None
    """
    # ------------------------------------------------------------------
    # 1. Validate
    # ------------------------------------------------------------------
    _validate_inputs(mdata)

    # ------------------------------------------------------------------
    # 2. Resolve marker panel: explicit arg > preset > None
    # ------------------------------------------------------------------
    if preset is not None and preset not in _PRESETS:
        raise ValueError(
            f"Unknown preset '{preset}'. "
            f"Available: {sorted(_PRESETS.keys())}"
        )

    resolved_panel: Optional[dict[str, list[str]]] = marker_panel
    if resolved_panel is None and preset is not None:
        resolved_panel = _PRESETS[preset]

    # ------------------------------------------------------------------
    # 3. Copy or in-place
    # ------------------------------------------------------------------
    adata = mdata["adt"] if inplace else mdata["adt"].copy()

    # ------------------------------------------------------------------
    # 4. Leiden clustering
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
    # 5. rank_genes_groups on leiden — top marker proteins per cluster.
    #    Computed on leiden (not adt_celltype) so the dotplot is available
    #    as a diagnostic to guide annotation decisions.
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.rank_genes_groups(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 6. Dendrogram on leiden — cluster ordering for the dotplot.
    # ------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.tl.dendrogram(adata, groupby="leiden")

    # ------------------------------------------------------------------
    # 7. Store marker panel in uns
    # ------------------------------------------------------------------
    if resolved_panel is not None:
        adata.uns["omicsage_adt_marker_panel"] = dict(resolved_panel)

    # ------------------------------------------------------------------
    # 8. Annotation
    #    Priority: explicit annotation_map > score from panel > none
    # ------------------------------------------------------------------
    annotated         = False
    annotation_source = "none"
    applied_map: dict[str, str] = {}
    score_map:   dict[str, str] = {}

    # --- 8a. Always score when a panel is available ----------------------
    #     Score-derived labels are stored in obs["adt_celltype_score"]
    #     regardless of whether an explicit annotation_map was also given.
    #     This lets the report always show the scoring result for comparison.
    if resolved_panel is not None:
        score_map = _score_clusters(adata, resolved_panel)
        adata.obs["adt_celltype_score"] = (
            adata.obs["leiden"].astype(str).replace(score_map)
        )

    if annotation_map is not None:
        # --- 8b. Explicit map: apply directly to obs["adt_celltype_manual"]
        adata.obs["adt_celltype_manual"] = adata.obs["leiden"].copy()
        found_clusters = set(adata.obs["leiden"].unique().tolist())
        unknown_keys   = set(str(k) for k in annotation_map.keys()) - found_clusters
        if unknown_keys:
            warnings.warn(
                f"annotation_map keys not found in clusters: "
                f"{sorted(unknown_keys)} — ignored.",
                UserWarning, stacklevel=2,
            )
        applied_map = {
            str(k): str(v)
            for k, v in annotation_map.items()
            if str(k) in found_clusters
        }
        if applied_map:
            adata.obs["adt_celltype_manual"] = (
                adata.obs["adt_celltype_manual"].astype(str).replace(applied_map)
            )
        annotated         = True
        annotation_source = "explicit"

    elif resolved_panel is not None:
        # --- 8c. No explicit map — score is the only annotation ----------
        applied_map       = score_map
        annotated         = True
        annotation_source = "scoring"

    # ------------------------------------------------------------------
    # 9. Metrics
    # ------------------------------------------------------------------
    cluster_sizes: dict[str, int] = {
        str(k): int(v)
        for k, v in adata.obs["leiden"].value_counts().items()
    }
    n_clusters = len(cluster_sizes)

    metrics: dict = {
        "n_cells":              int(adata.n_obs),
        "n_clusters":           n_clusters,
        "cluster_sizes":        cluster_sizes,
        "annotation_map":       applied_map,
        "score_map":            score_map,
        "annotation_source":    annotation_source,
        "annotated":            annotated,
        "score_celltype_key":   "adt_celltype_score" if score_map else None,
        "manual_celltype_key":  "adt_celltype_manual" if annotation_map is not None else None,
        "leiden_key":           "leiden",
        "resolution":           resolution,
        "n_iterations":         n_iterations,
        "random_state":         random_state,
        "marker_panel":         dict(resolved_panel) if resolved_panel else {},
        "preset":               preset,
    }

    # ------------------------------------------------------------------
    # 10. Provenance
    # ------------------------------------------------------------------
    adata.uns["omicsage_adt_annotate"] = {
        "module":    "adt_annotate",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "params": {
            "preset":                  preset,
            "resolution":              resolution,
            "flavor":                  "igraph",
            "n_iterations":            n_iterations,
            "directed":                False,
            "random_state":            random_state,
            "annotation_map_provided": annotation_map is not None,
            "marker_panel_provided":   marker_panel is not None,
        },
        "outputs": {
            "leiden_key":          "leiden",
            "score_celltype_key":  "adt_celltype_score" if score_map else None,
            "manual_celltype_key": "adt_celltype_manual" if annotation_map is not None else None,
            "n_clusters":          n_clusters,
            "annotation_source":   annotation_source,
        },
    }

    return adata, metrics


# ---------------------------------------------------------------------------
# Auto-annotation scoring
# ---------------------------------------------------------------------------

def _score_clusters(
    adata: AnnData,
    marker_panel: dict[str, list[str]],
) -> dict[str, str]:
    """
    Score each Leiden cluster against marker_panel and return a
    cluster_id → cell_type mapping.

    Scoring uses fold-change: for each (cluster, cell_type) pair, compute
    the mean expression of that cell type's markers *in the cluster* minus
    the mean expression of those same markers *in all other cells*.
    The cell type with the highest fold-change score wins.

    This avoids the raw-mean trap where ubiquitously expressed proteins
    (e.g. CD47, CD52) dominate the score regardless of cluster identity.

    Clusters where no panel marker is found in var_names → "Unknown".
    Ties broken alphabetically by cell type name.
    """
    X, var_names = _get_expression_matrix(adata)
    var_index = {v: i for i, v in enumerate(var_names)}

    cluster_ids = sorted(
        adata.obs["leiden"].unique().tolist(),
        key=lambda x: int(x) if x.isdigit() else x,
    )

    result: dict[str, str] = {}
    for cid in cluster_ids:
        mask = (adata.obs["leiden"] == cid).values
        other_mask = ~mask
        cluster_expr = X[mask, :]
        other_expr   = X[other_mask, :]

        best_celltype = "Unknown"
        best_score    = -np.inf

        for celltype, markers in marker_panel.items():
            indices = [var_index[m] for m in markers if m in var_index]
            if not indices:
                continue

            mean_in    = float(cluster_expr[:, indices].mean())
            mean_out   = float(other_expr[:, indices].mean())
            # Fold-change in expression space; add small epsilon to avoid /0
            score = mean_in - mean_out

            if score > best_score or (
                np.isclose(score, best_score) and celltype < best_celltype
            ):
                best_score    = score
                best_celltype = celltype

        result[str(cid)] = best_celltype

    return result


def _get_expression_matrix(adata: AnnData) -> tuple[np.ndarray, list[str]]:
    """Return (cells × proteins ndarray, protein names).
    Prefers adt_dsb layer → adt_clr layer → .X.
    """
    var_names = list(adata.var_names)
    for layer in ("adt_dsb", "adt_clr"):
        if layer in adata.layers:
            X = adata.layers[layer]
            if hasattr(X, "toarray"):
                X = X.toarray()
            return np.array(X, dtype=float), var_names
    X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.array(X, dtype=float), var_names


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _validate_inputs(mdata: MuData) -> None:
    if "adt" not in mdata.mod:
        raise KeyError(
            "mdata must contain an 'adt' modality. "
            "Found: " + str(list(mdata.mod.keys()))
        )
    if "connectivities" not in mdata["adt"].obsp:
        raise KeyError(
            "mdata['adt'].obsp['connectivities'] not found. "
            "Run adt_harmony.py before adt_annotate.py."
        )
