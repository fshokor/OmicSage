"""
annotate.py
===========
OmicSage cell-type annotation module.

Supported methods (select via `methods` parameter)
---------------------------------------------------
celltypist   — CellTypist majority-voting at cluster level
               Models: Immune_All_High.pkl (coarse) and/or Immune_All_Low.pkl (fine)
               Models are downloaded on first use and cached in
               <celltypist_models_dir> (default: data/references/celltypist/)

markers      — Mean log-normalised expression of canonical marker gene sets
               per Leiden cluster; best-scoring type per cluster stored in
               obs['cell_type_markers']

vote         — Majority vote across all methods that were run.
               Today: 2-way (CellTypist fine + markers).
               ScType and SingleR slots are reserved but empty until Session B.
               Requires both "celltypist" and "markers" to also be in methods.

Methods planned for Session B (rpy2/R)
---------------------------------------
sctype       — ScType with tissue-specific DB (default tissue: "Immune")
singler      — SingleR with HumanPrimaryCellAtlasData() reference

Public API
----------
    from pipeline.modules.qc.annotate import annotate, MARKER_SETS

    adata_ann, ann_dict = annotate(
        adata_clustered,
        methods=["celltypist", "markers", "vote"],
        leiden_col="leiden",
        celltypist_models=["Immune_All_High.pkl", "Immune_All_Low.pkl"],
        celltypist_models_dir=None,   # → data/references/celltypist/
        marker_sets=None,             # → built-in MARKER_SETS
        inplace=False,
    )

obs columns written
-------------------
  celltypist_coarse      — High model majority-vote label per cluster
  celltypist_fine        — Low  model majority-vote label per cluster
  cell_type_markers      — best marker-score label per cluster
  cell_type_groundtruth  — copy of obs['cell_type'] if it existed before annotation
                           (preserves publication ground-truth from being overwritten)
  cell_type_vote         — final consensus label  (written when "vote" in methods)
  cell_type_confidence   — fraction of methods agreeing (0.0–1.0)

uns provenance
--------------
  adata.uns['omicsage_annotate']
    methods_requested, methods_run, leiden_col,
    celltypist_models, celltypist_models_dir,
    marker_sets_keys, n_clusters, n_cells,
    celltypist_version, scanpy_version,
    omicsage_module, omicsage_version, timestamp
"""

from __future__ import annotations

import logging
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

# ── Module metadata ────────────────────────────────────────────────────────────
OMICSAGE_VERSION = "0.1.0"
OMICSAGE_MODULE  = "annotate"

# ── Built-in marker sets ───────────────────────────────────────────────────────
# Covers both liver (HCC dataset) and bone-marrow/blood (BMMC CITE-seq dataset).
# Users can pass their own dict via the `marker_sets` parameter.
MARKER_SETS: Dict[str, List[str]] = {
    # ── Myeloid ───────────────────────────────────────────────────────────────
    "Macrophage":   ["CD68", "MARCO", "CSF1R", "MRC1", "VSIG4", "GPNMB",
                     "SPP1", "C1QA", "C1QB", "TIMD4"],
    "Monocyte":     ["CD14", "LYZ", "S100A8", "S100A9", "FCN1", "VCAN", "CXCL8"],
    "DC":           ["FCER1A", "CLEC10A", "CST3", "CLEC9A", "CD1C"],
    "pDC":          ["LILRA4", "CLEC4C", "IL3RA", "TCF4", "IRF7", "PLAC8"],
    "Mast_cell":    ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC"],
    # ── Lymphoid ──────────────────────────────────────────────────────────────
    "T_cell":       ["CD3D", "CD3E", "TRAC", "TRBC2", "IL7R", "CD2"],
    "CD8_T_cell":   ["CD8A", "CD8B", "GZMK", "GZMA", "GZMB", "PRF1", "CCL5"],
    "NK_ILC":       ["NKG7", "GNLY", "NCAM1", "KLRB1", "KLRD1", "TYROBP"],
    "B_cell":       ["MS4A1", "CD79A", "CD79B", "IGHM", "IGHD"],
    "Plasma_cell":  ["MZB1", "JCHAIN", "IGHG1", "XBP1", "SDC1", "PRDM1"],
    # ── Progenitors / erythroid (BMMC-specific) ───────────────────────────────
    "Progenitor":   ["CD34", "SPINK2", "CYTL1", "MLLT3", "MECOM", "AVP"],
    "Erythroid":    ["HBA1", "HBA2", "HBB", "GYPA", "TFRC", "KLF1", "GATA1"],
    "Platelet":     ["PPBP", "PF4", "GP1BA", "ITGA2B", "TUBB1"],
    # ── Non-immune / parenchymal ──────────────────────────────────────────────
    "Hepatocyte":   ["ALB", "APOC3", "TTR", "FGB", "FGG", "CYP3A4",
                     "GPC3", "APOE", "FABP1"],
    "Fibroblast":   ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
                     "ACTA2", "PDGFRB", "FAP"],
    "Endothelial":  ["PECAM1", "VWF", "CDH5", "CLDN5", "LYVE1", "ENG"],
}

# Cell types that are tissue-specific (get double weight in vote if ScType agrees)
# Reserved for Session B — not used in the 2-way vote today.
_PARENCHYMAL: frozenset = frozenset({"Hepatocyte", "Fibroblast", "Endothelial"})

# ── Default paths ──────────────────────────────────────────────────────────────
_DEFAULT_MODELS_DIR = Path("data/references/celltypist")
_VALID_METHODS = {"celltypist", "markers", "vote", "sctype", "singler"}


# ─────────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────────

def annotate(
    adata: sc.AnnData,
    methods: List[str] = None,
    leiden_col: str = "leiden",
    celltypist_models: List[str] = None,
    celltypist_models_dir: Optional[Path] = None,
    marker_sets: Optional[Dict[str, List[str]]] = None,
    inplace: bool = False,
) -> Tuple[sc.AnnData, dict]:
    """
    Annotate cell clusters using one or more automated methods.

    Parameters
    ----------
    adata : AnnData
        Clustered AnnData — must have adata.obs[leiden_col] and
        adata.layers["logcounts"] (log-normalised counts).
    methods : list of str, optional
        Subset of {"celltypist", "markers", "vote"}.
        "vote" requires both "celltypist" and "markers" also in the list.
        Defaults to ["celltypist", "markers", "vote"].
    leiden_col : str
        obs column with cluster labels (default: "leiden").
    celltypist_models : list of str, optional
        CellTypist model filenames to use.
        Default: ["Immune_All_High.pkl", "Immune_All_Low.pkl"].
    celltypist_models_dir : Path or None
        Directory to cache downloaded CellTypist models.
        Default: data/references/celltypist/ (relative to cwd).
    marker_sets : dict or None
        {cell_type: [gene_list]}. Defaults to built-in MARKER_SETS.
    inplace : bool
        If False (default), work on a copy; caller's object is unchanged.

    Returns
    -------
    adata_ann : AnnData
        Annotated AnnData with new obs columns and uns provenance.
    ann_dict : dict
        Per-cluster annotation summary (same data as uns provenance + vote tables).
    """
    # ── Defaults ───────────────────────────────────────────────────────────────
    if methods is None:
        methods = ["celltypist", "markers", "vote"]
    if celltypist_models is None:
        celltypist_models = ["Immune_All_High.pkl", "Immune_All_Low.pkl"]
    if celltypist_models_dir is None:
        celltypist_models_dir = _DEFAULT_MODELS_DIR
    if marker_sets is None:
        marker_sets = MARKER_SETS

    celltypist_models_dir = Path(celltypist_models_dir)

    # ── Validate ───────────────────────────────────────────────────────────────
    unknown = set(methods) - _VALID_METHODS
    if unknown:
        raise ValueError(f"Unknown methods: {unknown}. Valid: {_VALID_METHODS}")

    if leiden_col not in adata.obs.columns:
        raise KeyError(
            f"leiden_col='{leiden_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    if "vote" in methods:
        _check_vote_prerequisites(methods)

    # Warn about Session B methods requested but not yet implemented
    for m in ("sctype", "singler"):
        if m in methods:
            warnings.warn(
                f"Method '{m}' is planned for Session B (requires rpy2/R) "
                f"and is not yet implemented. It will be skipped.",
                UserWarning, stacklevel=2,
            )
    methods = [m for m in methods if m not in ("sctype", "singler")]

    # ── Copy guard ─────────────────────────────────────────────────────────────
    adata_ann = adata if inplace else adata.copy()

    # ── Preserve ground-truth before it can be overwritten ─────────────────────
    # Many public datasets (incl. GSE194122) ship with obs['cell_type'] as the
    # publication ground-truth.  Copy it now so the vote output never destroys it.
    if "cell_type" in adata_ann.obs.columns:
        adata_ann.obs["cell_type_groundtruth"] = adata_ann.obs["cell_type"].copy()
        logger.info("obs['cell_type'] preserved → obs['cell_type_groundtruth'].")

    # ── Run methods ────────────────────────────────────────────────────────────
    methods_run: List[str] = []
    ann_dict: dict = {}

    if "celltypist" in methods:
        try:
            _run_celltypist(
                adata_ann, leiden_col, celltypist_models, celltypist_models_dir
            )
            methods_run.append("celltypist")
            logger.info("CellTypist annotation complete.")
        except ImportError:
            warnings.warn(
                "celltypist is not installed — skipping CellTypist annotation. "
                "Install with: pip install celltypist",
                UserWarning, stacklevel=2,
            )

    if "markers" in methods:
        score_df = _run_marker_scoring(adata_ann, leiden_col, marker_sets)
        ann_dict["marker_score_df"] = score_df
        methods_run.append("markers")
        logger.info("Marker scoring complete.")

    if "vote" in methods and "celltypist" in methods_run and "markers" in methods_run:
        vote_df = _run_majority_vote(adata_ann, ann_dict["marker_score_df"], leiden_col)
        ann_dict["vote_df"] = vote_df
        methods_run.append("vote")
        logger.info("Majority vote complete.")
    elif "vote" in methods:
        warnings.warn(
            "Majority vote requested but not all prerequisite methods ran "
            "successfully — skipping vote.",
            UserWarning, stacklevel=2,
        )

    # ── Provenance ─────────────────────────────────────────────────────────────
    ct_version = _get_celltypist_version()
    provenance = {
        "methods_requested":      methods,
        "methods_run":            methods_run,
        "leiden_col":             leiden_col,
        "celltypist_models":      celltypist_models,
        "celltypist_models_dir":  str(celltypist_models_dir),
        "marker_sets_keys":       list(marker_sets.keys()),
        "n_clusters":             int(adata_ann.obs[leiden_col].nunique()),
        "n_cells":                int(adata_ann.n_obs),
        "celltypist_version":     ct_version,
        "scanpy_version":         sc.__version__,
        "omicsage_module":        OMICSAGE_MODULE,
        "omicsage_version":       OMICSAGE_VERSION,
        "timestamp":              datetime.now(timezone.utc).isoformat(),
    }
    adata_ann.uns["omicsage_annotate"] = provenance
    ann_dict["provenance"] = provenance

    return adata_ann, ann_dict


# ─────────────────────────────────────────────────────────────────────────────
# CellTypist
# ─────────────────────────────────────────────────────────────────────────────

def _run_celltypist(
    adata: sc.AnnData,
    leiden_col: str,
    model_names: List[str],
    models_dir: Path,
) -> None:
    """
    Run CellTypist majority-vote annotation and write results into adata.obs.

    Models are downloaded on first use into `models_dir` (inside the project
    repo) so that runs are reproducible without relying on system-wide caches.

    Writes
    ------
    obs['celltypist_coarse']  — if Immune_All_High.pkl in model_names
    obs['celltypist_fine']    — if Immune_All_Low.pkl  in model_names
    """
    import celltypist
    from celltypist import models as ct_models

    # ── Build a CP10K log-normalised dense copy for CellTypist ────────────────
    adata_ct = adata.copy()

    # Prefer raw counts layer for re-normalisation; fall back to .X
    if "counts" in adata_ct.layers:
        adata_ct.X = adata_ct.layers["counts"].copy()
    elif "logcounts" in adata_ct.layers:
        # Already log-normalised — use as-is (less ideal but workable)
        adata_ct.X = adata_ct.layers["logcounts"].copy()
        logger.warning(
            "No 'counts' layer found — using 'logcounts' directly. "
            "CellTypist prefers CP10K-normalised counts."
        )
    # else: use .X as-is and let CellTypist warn if needed

    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)
    if hasattr(adata_ct.X, "toarray"):
        adata_ct.X = adata_ct.X.toarray()

    # ── Ensure models directory exists inside the project repo ────────────────
    models_dir.mkdir(parents=True, exist_ok=True)

    # Override CellTypist's default cache to use our project-local directory.
    # This keeps downloaded models inside the repo (data/references/celltypist/)
    # rather than ~/.celltypist/ — making the project fully self-contained.
    import os
    _original_cache = os.environ.get("CELLTYPIST_FOLDER")
    os.environ["CELLTYPIST_FOLDER"] = str(models_dir.resolve())

    try:
        # Download only the models we need (force_update=False = skip if cached)
        ct_models.download_models(force_update=False, model=model_names)

        _MODEL_TO_OBS = {
            "Immune_All_High.pkl": "celltypist_coarse",
            "Immune_All_Low.pkl":  "celltypist_fine",
        }

        for model_name in model_names:
            obs_col = _MODEL_TO_OBS.get(model_name, f"celltypist_{model_name.replace('.pkl','')}")
            model   = ct_models.Model.load(model=model_name)
            pred    = celltypist.annotate(
                adata_ct,
                model=model,
                majority_voting=True,
                over_clustering=leiden_col,   # cluster-level majority voting
            )
            pred_adata = pred.to_adata()
            adata.obs[obs_col] = (
                pred_adata.obs.loc[adata.obs.index, "majority_voting"]
                .astype("category")
            )
            n_types = adata.obs[obs_col].nunique()
            logger.info(f"  {obs_col}: {n_types} types ({model_name})")
            print(f"  CellTypist {obs_col}: {n_types} cell types")

    finally:
        # Restore original env var so we don't pollute the caller's environment
        if _original_cache is None:
            os.environ.pop("CELLTYPIST_FOLDER", None)
        else:
            os.environ["CELLTYPIST_FOLDER"] = _original_cache

    # Guarantee both standard columns always exist
    if "celltypist_coarse" not in adata.obs.columns:
        adata.obs["celltypist_coarse"] = "not_run"
    if "celltypist_fine" not in adata.obs.columns:
        adata.obs["celltypist_fine"] = "not_run"


# ─────────────────────────────────────────────────────────────────────────────
# Marker gene scoring
# ─────────────────────────────────────────────────────────────────────────────

def _run_marker_scoring(
    adata: sc.AnnData,
    leiden_col: str,
    marker_sets: Dict[str, List[str]],
) -> pd.DataFrame:
    """
    Compute mean log-normalised expression of each marker set per cluster.

    Prefers adata.layers['logcounts'] if present, else adata.X.

    Returns
    -------
    score_df : DataFrame (clusters × cell-types + 'best_by_score' column)
    Writes obs['cell_type_markers'] into adata.
    """
    # Choose expression matrix
    if "logcounts" in adata.layers:
        X = adata.layers["logcounts"]
    else:
        X = adata.X

    clusters = sorted(adata.obs[leiden_col].unique(),
                      key=lambda x: (int(x) if str(x).isdigit() else x))
    rows = []
    for cl in clusters:
        mask = (adata.obs[leiden_col] == cl).values
        row  = {"cluster": str(cl)}
        for ct, markers in marker_sets.items():
            present = [g for g in markers if g in adata.var_names]
            if present:
                gene_idx = [adata.var_names.get_loc(g) for g in present]
                expr = X[np.where(mask)[0], :][:, gene_idx]
                if hasattr(expr, "toarray"):
                    expr = expr.toarray()
                row[ct] = float(np.asarray(expr).mean())
            else:
                row[ct] = 0.0
        rows.append(row)

    score_df = pd.DataFrame(rows).set_index("cluster")
    score_df["best_by_score"] = score_df.idxmax(axis=1)

    # Map best label back onto every cell
    cluster_map = score_df["best_by_score"].to_dict()
    adata.obs["cell_type_markers"] = (
        adata.obs[leiden_col]
        .astype(str)
        .map(cluster_map)
        .astype("category")
    )

    n_types = adata.obs["cell_type_markers"].nunique()
    print(f"  Marker scoring: {n_types} cell types across {len(clusters)} clusters")
    return score_df


# ─────────────────────────────────────────────────────────────────────────────
# Majority vote
# ─────────────────────────────────────────────────────────────────────────────

def _run_majority_vote(
    adata: sc.AnnData,
    score_df: pd.DataFrame,
    leiden_col: str,
) -> pd.DataFrame:
    """
    Combine available method labels into a consensus via majority vote.

    Evidence sources (today — 2-way):
      1. celltypist_fine  (CellTypist Low model, cluster-level majority)
      2. cell_type_markers (marker gene scoring)

    Reserved for Session B (slots present but skipped when column absent):
      3. sctype_cell_type  (ScType)
      4. SingleR_HPCA      (SingleR)

    ScType gets double weight for parenchymal types when present (Session B).

    Writes
    ------
    obs['cell_type']           — final consensus label
    obs['cell_type_confidence']— fraction of active methods agreeing (0.0–1.0)
    """
    def _cluster_majority(series: pd.Series) -> str:
        vc = series.dropna().value_counts()
        return vc.index[0] if len(vc) > 0 else "Unknown"

    # ── Build per-cluster evidence table ──────────────────────────────────────
    vote_rows: List[dict] = []

    # All potential evidence columns, in priority order
    # slot_name → (obs_col, weight, available)
    evidence_slots = [
        ("CellTypist",   "celltypist_fine",    1),
        ("Markers",      "cell_type_markers",  1),
        # Session B slots — will be non-None when those columns exist
        ("ScType",       "sctype_cell_type",   1),   # +1 extra for parenchymal
        ("SingleR_HPCA", "SingleR_HPCA",       1),
    ]

    clusters = sorted(adata.obs[leiden_col].unique(),
                      key=lambda x: (int(x) if str(x).isdigit() else x))

    for cl in clusters:
        mask   = adata.obs[leiden_col] == str(cl)
        row    = {"cluster": str(cl), "n_cells": int(mask.sum())}
        votes: List[str] = []

        for slot_name, obs_col, weight in evidence_slots:
            if obs_col not in adata.obs.columns:
                row[slot_name] = None
                continue
            label = _cluster_majority(adata.obs.loc[mask, obs_col])
            row[slot_name] = label

            # Add to vote pool (repeated weight times)
            if label and label not in ("not_run", "Unknown", "nan"):
                votes.extend([label] * weight)

                # Double weight for ScType parenchymal (Session B)
                if slot_name == "ScType" and label in _PARENCHYMAL:
                    votes.append(label)

        # Compute winner and confidence
        if votes:
            counts  = Counter(votes)
            winner, top_count = counts.most_common(1)[0]
            confidence = round(top_count / len(votes), 4)
        else:
            winner, confidence = "Unknown", 0.0

        row["final_label"]  = winner
        row["confidence"]   = confidence
        vote_rows.append(row)

    vote_df = pd.DataFrame(vote_rows).set_index("cluster")

    # ── Map labels back onto cells ─────────────────────────────────────────────
    label_map      = vote_df["final_label"].to_dict()
    confidence_map = vote_df["confidence"].to_dict()

    adata.obs["cell_type_vote"] = (
        adata.obs[leiden_col]
        .astype(str)
        .map(label_map)
        .astype("category")
    )
    adata.obs["cell_type_confidence"] = (
        adata.obs[leiden_col]
        .astype(str)
        .map(confidence_map)
        .astype(float)
    )

    n_types = adata.obs["cell_type_vote"].nunique()
    active_methods = [s for s, c, _ in evidence_slots if c in adata.obs.columns
                      and adata.obs[c].iloc[0] != "not_run"]
    print(f"  Majority vote ({len(active_methods)} methods): {n_types} consensus types")
    print(vote_df[["n_cells", "final_label", "confidence"]].to_string())

    return vote_df


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _check_vote_prerequisites(methods: List[str]) -> None:
    """Raise if vote is requested but its dependencies are missing."""
    required = {"celltypist", "markers"}
    missing  = required - set(methods)
    if missing:
        raise ValueError(
            f"'vote' requires {required} to also be in methods. "
            f"Missing: {missing}"
        )


def _get_celltypist_version() -> str:
    try:
        import celltypist
        return celltypist.__version__
    except ImportError:
        return "not_installed"
