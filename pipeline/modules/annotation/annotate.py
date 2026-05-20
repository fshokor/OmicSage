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
               Up to 3-way today (CellTypist + markers + ScType).
               SingleR slot reserved for a future session.
               Requires both "celltypist" and "markers" to also be in methods.

scanvi       — scANVI transfer learning via scvi-hub (optional, GPU-accelerated).
               Requires a pre-trained scANVI model directory (scanvi_model param).
               Posterior probability per cell is used as fractional vote weight.
               All scANVI tests are skipped in CI (OMICSAGE_CI=true).

Public API
----------
    from pipeline.modules.annotation.annotate import annotate, MARKER_SETS

    adata_ann, ann_dict = annotate(
        adata_clustered,
        methods=["celltypist", "markers", "vote"],
        leiden_col="leiden",
        celltypist_models=["Immune_All_High.pkl", "Immune_All_Low.pkl"],
        celltypist_models_dir=None,   # → data/references/celltypist/
        marker_sets=None,             # → built-in MARKER_SETS
        scanvi_model=None,            # → path to pre-trained scANVI model dir
        inplace=False,
    )

obs columns written
-------------------
  celltypist_coarse      — High model majority-vote label per cluster
  celltypist_fine        — Low  model majority-vote label per cluster
  cell_type_markers      — best marker-score label per cluster
  cell_type_sctype       — ScType best label per cluster (when "sctype" in methods)
  cell_type_scanvi       — scANVI transfer label per cell (when "scanvi" in methods)
  cell_type_groundtruth  — copy of obs['cell_type'] if it existed before annotation
                           (preserves publication ground-truth from being overwritten)
  cell_type_vote         — final consensus label  (written when "vote" in methods)
  cell_type_confidence   — fraction of methods agreeing (0.0–1.0)

uns provenance
--------------
  adata.uns['omicsage_annotate']
    methods_requested, methods_run, leiden_col,
    celltypist_models, celltypist_models_dir,
    marker_sets_keys, sctype_tissue,
    n_clusters, n_cells,
    celltypist_version, scanpy_version,
    omicsage_module, omicsage_version, timestamp
"""

from __future__ import annotations

import logging
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

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

# Cell types that are tissue-specific — get double weight in vote when ScType agrees
_PARENCHYMAL: frozenset = frozenset({"Hepatocyte", "Fibroblast", "Endothelial"})

# ── Default paths ──────────────────────────────────────────────────────────────
_DEFAULT_MODELS_DIR = Path("data/references/celltypist")
_VALID_METHODS = {"celltypist", "markers", "vote", "sctype", "singler", "scanvi"}

# ScType DB — fetched fresh at runtime, never cached to disk
_SCTYPE_DB_URL = (
    "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/"
    "master/ScTypeDB_full.xlsx"
)


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
    tissue: str = "Immune system",
    scanvi_model: Optional[str] = None,
    singler_ref: Optional[Union[str, Path]] = None,
    singler_ref_label_col: str = "cell_type",
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
        Subset of {"celltypist", "markers", "sctype", "singler", "vote"}.
        "vote" requires at least "celltypist" and "markers" also in the list.
        Defaults to ["celltypist", "markers", "sctype", "singler", "vote"].
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
    tissue : str
        Tissue filter for ScTypeDB (default: "Immune system").
        Valid values: "Immune system", "Liver", "Brain", "Lung", "Kidney",
        "Pancreas", "Heart", "Intestine", "Skin", "Adrenal", "Bladder",
        "Breast", "Cervix", "Eye", "Muscle", "Ovary", "Prostate",
        "Stomach", "Thyroid".
        Only used when "sctype" is in methods.
    scanvi_model : str or None
        Path to a pre-trained scANVI model directory (output of
        SCVI_Model.save() or SCANVI_Model.save()).  Required when "scanvi"
        is in methods; ignored otherwise.
        The model must have been trained with the same gene set as adata.
    singler_ref : str, Path, or None
        Reference source for SingleR annotation:
          None (default) → Human Primary Cell Atlas (HPCA) pseudobulk reference
                           (Mabbott et al. 2013, processed by Aran et al. 2019).
                           Downloaded via the ``celldex`` Python package on first
                           call and cached at
                           ``~/.cache/omicsage/singler/hpca_ref.h5ad``.
                           Requires: ``pip install celldex``
          "hca"          → HCA Census API (streamed at runtime, requires
                           cellxgene-census package).
          str/Path        → Path to a user-supplied H5AD pseudobulk reference
                           (obs = cell types, var = genes, X = mean log-
                           normalised expression).
        Only used when "singler" is in methods.
    singler_ref_label_col : str
        obs column in user-supplied H5AD that contains cell type labels.
        Only used when singler_ref is a file path. Default: "cell_type".
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
        methods = ["celltypist", "markers", "sctype", "singler", "vote"]
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

    # Validate scanvi_model is provided when scanvi is requested
    if "scanvi" in methods and scanvi_model is None:
        raise ValueError(
            "'scanvi' is in methods but scanvi_model=None. "
            "Provide the path to a pre-trained scANVI model directory."
        )

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

    if "sctype" in methods:
        try:
            _run_sctype(adata_ann, leiden_col, tissue)
            methods_run.append("sctype")
            logger.info("ScType annotation complete.")
        except Exception as e:
            warnings.warn(
                f"ScType annotation failed ({e}) — skipping. "
                f"Check your internet connection or tissue parameter.",
                UserWarning, stacklevel=2,
            )

    if "singler" in methods:
        try:
            _run_singler(
                adata_ann, leiden_col,
                singler_ref=singler_ref,
                singler_ref_label_col=singler_ref_label_col,
            )
            methods_run.append("singler")
            logger.info("SingleR annotation complete.")
        except Exception as e:
            warnings.warn(
                f"SingleR annotation failed ({e}) — skipping.",
                UserWarning, stacklevel=2,
            )

    if "scanvi" in methods:
        try:
            cluster_posteriors = _run_scanvi(adata_ann, leiden_col, scanvi_model)
            ann_dict["scanvi_posteriors"] = cluster_posteriors
            methods_run.append("scanvi")
            logger.info("scANVI annotation complete.")
        except ImportError:
            warnings.warn(
                "scvi-tools is not installed — skipping scANVI annotation. "
                "Install with: pip install scvi-tools",
                UserWarning, stacklevel=2,
            )
        except Exception as e:
            warnings.warn(
                f"scANVI annotation failed ({e}) — skipping.",
                UserWarning, stacklevel=2,
            )

    if "vote" in methods and "celltypist" in methods_run and "markers" in methods_run:
        # Collect per-cluster mean posterior probabilities from scANVI if run
        scanvi_posteriors = (
            ann_dict.get("scanvi_posteriors")  # dict: cluster → mean_prob
            if "scanvi" in methods_run else None
        )
        vote_df = _run_majority_vote(
            adata_ann, ann_dict["marker_score_df"], leiden_col,
            scanvi_posteriors=scanvi_posteriors,
        )
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
        "sctype_tissue":          tissue,
        "scanvi_model_path":      str(scanvi_model) if scanvi_model else None,
        "singler_ref":            str(singler_ref) if singler_ref is not None else "hpca",
        "singler_ref_label_col":  singler_ref_label_col,
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
    scanvi_posteriors: Optional[Dict[str, float]] = None,
) -> pd.DataFrame:
    """
    Combine available method labels into a consensus via N-way weighted vote.

    Evidence sources (up to 4-way):
      1. celltypist_fine    (CellTypist Low model, cluster-level majority)  weight=1
      2. cell_type_markers  (marker gene scoring)                           weight=1
      3. cell_type_sctype   (ScType)                                        weight=1
                            +1 extra for parenchymal types (double weight)
      4. cell_type_scanvi   (scANVI transfer)
                            weight = mean posterior probability × _SCANVI_MAX_VOTES
                            (fractional: high-confidence calls get more influence)
      5. cell_type_singler  (reserved; skipped when column absent)

    Writes
    ------
    obs['cell_type_vote']       — final consensus label
    obs['cell_type_confidence'] — fraction of active vote weight agreeing (0.0–1.0)
    """
    _SCANVI_MAX_VOTES = 2  # max fractional votes scANVI can contribute per cluster

    def _cluster_majority(series: pd.Series) -> str:
        vc = series.dropna().value_counts()
        return vc.index[0] if len(vc) > 0 else "Unknown"

    # ── Build per-cluster evidence table ──────────────────────────────────────
    vote_rows: List[dict] = []

    # Integer-weight slots (name, obs_col, weight)
    int_slots = [
        ("CellTypist",   "celltypist_fine",   1),
        ("Markers",      "cell_type_markers", 1),
        ("ScType",       "cell_type_sctype",  1),   # +1 extra for parenchymal
        ("SingleR",      "cell_type_singler", 1),   # future
    ]

    clusters = sorted(adata.obs[leiden_col].unique(),
                      key=lambda x: (int(x) if str(x).isdigit() else x))

    for cl in clusters:
        mask   = adata.obs[leiden_col] == str(cl)
        row    = {"cluster": str(cl), "n_cells": int(mask.sum())}
        vote_weight: Dict[str, float] = {}  # label → accumulated weight
        total_weight = 0.0

        # ── Integer-weight slots ─────────────────────────────────────────────
        for slot_name, obs_col, weight in int_slots:
            if obs_col not in adata.obs.columns:
                row[slot_name] = None
                continue
            label = _cluster_majority(adata.obs.loc[mask, obs_col])
            row[slot_name] = label

            if label and label not in ("not_run", "Unknown", "nan"):
                w = float(weight)
                # Double weight for parenchymal types when ScType labels them
                if slot_name == "ScType" and label in _PARENCHYMAL:
                    w += 1.0
                vote_weight[label] = vote_weight.get(label, 0.0) + w
                total_weight += w

        # ── Fractional scANVI slot ───────────────────────────────────────────
        if "cell_type_scanvi" in adata.obs.columns:
            scanvi_label = _cluster_majority(adata.obs.loc[mask, "cell_type_scanvi"])
            row["scANVI"] = scanvi_label

            if scanvi_label and scanvi_label not in ("not_run", "Unknown", "nan"):
                # Use per-cluster mean posterior as fractional weight
                mean_prob = (
                    float(scanvi_posteriors[str(cl)])
                    if scanvi_posteriors and str(cl) in scanvi_posteriors
                    else 0.5  # fallback if posteriors not available
                )
                w = mean_prob * _SCANVI_MAX_VOTES
                vote_weight[scanvi_label] = vote_weight.get(scanvi_label, 0.0) + w
                total_weight += w
        else:
            row["scANVI"] = None

        # ── Compute winner and confidence ────────────────────────────────────
        if vote_weight:
            winner = max(vote_weight, key=vote_weight.__getitem__)
            confidence = round(vote_weight[winner] / total_weight, 4)
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
    active_slots = [s for s, c, _ in int_slots if c in adata.obs.columns
                    and adata.obs[c].iloc[0] != "not_run"]
    if "cell_type_scanvi" in adata.obs.columns:
        active_slots.append("scANVI")
    print(f"  Majority vote ({len(active_slots)} methods): {n_types} consensus types")
    print(vote_df[["n_cells", "final_label", "confidence"]].to_string())

    return vote_df




# ─────────────────────────────────────────────────────────────────────────────
# ScType-py
# ─────────────────────────────────────────────────────────────────────────────

def _fetch_sctype_db() -> pd.DataFrame:
    """
    Fetch ScTypeDB_full.xlsx fresh from GitHub at runtime.
    No file is written to disk — data lives in memory only.
    Always returns the current version maintained by the ScType authors.
    """
    import io
    import requests

    logger.info(f"Fetching ScTypeDB from {_SCTYPE_DB_URL} ...")
    response = requests.get(_SCTYPE_DB_URL, timeout=30)
    response.raise_for_status()
    db = pd.read_excel(io.BytesIO(response.content), engine="openpyxl")
    logger.info(f"  ScTypeDB loaded: {len(db)} rows, tissues: {db['tissueType'].nunique()}")
    return db


def _parse_sctype_markers(
    db: pd.DataFrame,
    tissue: str,
) -> Dict[str, Dict[str, List[str]]]:
    """
    Parse ScTypeDB into {cell_type: {"pos": [...], "neg": [...]}} for a given tissue.

    Parameters
    ----------
    db : DataFrame
        Raw ScTypeDB loaded from Excel.
    tissue : str
        Tissue filter — must match a value in db['tissueType'] exactly.

    Returns
    -------
    marker_dict : dict
        {cell_type: {"pos": [gene, ...], "neg": [gene, ...]}}

    Raises
    ------
    ValueError
        If the requested tissue is not found in the database.
    """
    available = sorted(db["tissueType"].dropna().unique().tolist())
    if tissue not in available:
        raise ValueError(
            f"Tissue '{tissue}' not found in ScTypeDB. "
            f"Available tissues: {available}"
        )

    subset = db[db["tissueType"] == tissue].copy()
    marker_dict: Dict[str, Dict[str, List[str]]] = {}

    for _, row in subset.iterrows():
        cell_type = str(row["cellName"]).strip()

        # Positive markers — column: geneSymbolmore1
        pos_raw = str(row.get("geneSymbolmore1", ""))
        pos = [g.strip() for g in pos_raw.split(",") if g.strip() and g.strip() != "nan"]

        # Negative markers — column: geneSymbolmore2
        neg_raw = str(row.get("geneSymbolmore2", ""))
        neg = [g.strip() for g in neg_raw.split(",") if g.strip() and g.strip() != "nan"]

        if pos or neg:
            marker_dict[cell_type] = {"pos": pos, "neg": neg}

    logger.info(f"  ScType parsed {len(marker_dict)} cell types for tissue='{tissue}'")
    return marker_dict


def _run_sctype(
    adata: sc.AnnData,
    leiden_col: str,
    tissue: str = "Immune system",
) -> None:
    """
    Run ScType-py annotation and write results into adata.obs.

    Fetches ScTypeDB_full.xlsx fresh from GitHub at runtime (no local cache).
    Scores each cluster as: mean(positive marker expr) - mean(negative marker expr).
    Assigns the cell type with the highest score to each cluster.

    Parameters
    ----------
    adata : AnnData
        Must have adata.layers['logcounts'] or adata.X with log-normalised counts.
    leiden_col : str
        obs column with cluster labels.
    tissue : str
        Tissue filter for ScTypeDB (default: "Immune system").

    Writes
    ------
    obs['cell_type_sctype'] — ScType best label per cluster (category dtype)
    """
    # ── Fetch and parse DB ────────────────────────────────────────────────────
    db = _fetch_sctype_db()
    marker_dict = _parse_sctype_markers(db, tissue)

    # ── Choose expression matrix ──────────────────────────────────────────────
    if "logcounts" in adata.layers:
        X = adata.layers["logcounts"]
    else:
        X = adata.X

    # ── Score each cluster ────────────────────────────────────────────────────
    clusters = sorted(
        adata.obs[leiden_col].unique(),
        key=lambda x: (int(x) if str(x).isdigit() else x),
    )

    def _mean_expr(mask: np.ndarray, genes: List[str]) -> float:
        """Mean expression of a gene list across cells in mask."""
        present = [g for g in genes if g in adata.var_names]
        if not present:
            return 0.0
        gene_idx = [adata.var_names.get_loc(g) for g in present]
        expr = X[np.where(mask)[0], :][:, gene_idx]
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        return float(np.asarray(expr).mean())

    cluster_labels: Dict[str, str] = {}

    for cl in clusters:
        mask = (adata.obs[leiden_col] == cl).values
        best_type  = "Unknown"
        best_score = -np.inf

        for cell_type, markers in marker_dict.items():
            pos_score = _mean_expr(mask, markers["pos"])
            neg_score = _mean_expr(mask, markers["neg"])
            score     = pos_score - neg_score

            if score > best_score:
                best_score = score
                best_type  = cell_type

        cluster_labels[str(cl)] = best_type
        logger.debug(f"  Cluster {cl} → {best_type} (score={best_score:.4f})")

    # ── Map labels onto cells ─────────────────────────────────────────────────
    adata.obs["cell_type_sctype"] = (
        adata.obs[leiden_col]
        .astype(str)
        .map(cluster_labels)
        .astype("category")
    )

    n_types = adata.obs["cell_type_sctype"].nunique()
    print(f"  ScType ({tissue}): {n_types} cell types across {len(clusters)} clusters")


# ─────────────────────────────────────────────────────────────────────────────
# scANVI transfer learning annotation
# ─────────────────────────────────────────────────────────────────────────────

def _run_scanvi(
    adata: sc.AnnData,
    leiden_col: str,
    model_path: str,
) -> Dict[str, float]:
    """
    Run scANVI transfer label annotation on adata using a pre-trained model.

    The model is loaded with SCANVI.load() (scvi-tools ≥ 1.0).  If the model
    directory contains a SCVI base model rather than a SCANVI model, scvi-tools
    will raise; callers should ensure the model type matches.

    Posterior probabilities are extracted per cell via
    model.predict() which returns a DataFrame with columns:
      - 'prediction'    : predicted cell type label
      - 'probability'   : softmax probability for the winning class

    Parameters
    ----------
    adata : AnnData
        Must share the same gene set the model was trained on.
    leiden_col : str
        obs column with cluster labels (used for per-cluster posterior summary).
    model_path : str
        Path to a saved SCANVI model directory.

    Returns
    -------
    cluster_posteriors : dict
        {cluster_label: mean_posterior_probability} — used by vote as
        fractional weight for each cluster's scANVI call.

    Writes
    ------
    obs['cell_type_scanvi']  — per-cell predicted label (category dtype)
    """
    import scvi.model as _scvi_model  # type: ignore  # noqa: F401

    logger.info(f"Loading scANVI model from {model_path} ...")
    model = _scvi_model.SCANVI.load(model_path, adata=adata)

    # predict() returns a DataFrame: index = cell barcode, columns vary by version
    # scvi-tools ≥ 1.0: columns are 'prediction' and 'probability'
    pred_df = model.predict(adata, soft=False)

    if isinstance(pred_df, pd.DataFrame):
        if "prediction" in pred_df.columns and "probability" in pred_df.columns:
            labels      = pred_df["prediction"].values
            posteriors  = pred_df["probability"].values.astype(float)
        else:
            # Fallback: first column is label, second is probability
            labels     = pred_df.iloc[:, 0].values
            posteriors = pred_df.iloc[:, 1].values.astype(float) \
                         if pred_df.shape[1] > 1 else np.full(len(labels), 0.5)
    else:
        # Some versions return a plain array of labels
        labels     = np.asarray(pred_df)
        posteriors = np.full(len(labels), 0.5)

    adata.obs["cell_type_scanvi"] = pd.Categorical(labels)
    n_types = adata.obs["cell_type_scanvi"].nunique()
    print(f"  scANVI: {n_types} cell types predicted")
    logger.info(f"  scANVI mean posterior probability: {posteriors.mean():.3f}")

    # ── Compute per-cluster mean posterior ────────────────────────────────────
    cluster_posteriors: Dict[str, float] = {}
    for cl in adata.obs[leiden_col].unique():
        mask = (adata.obs[leiden_col] == cl).values
        cluster_posteriors[str(cl)] = float(posteriors[mask].mean())

    return cluster_posteriors


# ─────────────────────────────────────────────────────────────────────────────
# SingleR — Spearman correlation against a pseudobulk reference
# ─────────────────────────────────────────────────────────────────────────────

_SINGLER_CACHE_DIR  = Path.home() / ".cache" / "omicsage" / "singler"
_HPCA_CACHE         = _SINGLER_CACHE_DIR / "hpca_ref.h5ad"

# celldex dataset name and version for Human Primary Cell Atlas
_CELLDEX_HPCA_NAME    = "hpca"
_CELLDEX_HPCA_VERSION = "2024-02-26"

# Keep the old name as an alias so the public symbol _TABULA_SAPIENS_CACHE
# still resolves (used in the CI-skipped cache test).
_TABULA_SAPIENS_CACHE = _HPCA_CACHE


def _load_hpca_reference() -> sc.AnnData:
    """
    Download the Human Primary Cell Atlas (HPCA) pseudobulk reference via the
    Python ``celldex`` package on first call, then cache the result at
    ``~/.cache/omicsage/singler/hpca_ref.h5ad``.

    The HPCA reference (Mabbott et al. 2013, processed by Aran et al. 2019)
    contains 713 microarray samples covering 37 main cell types and is the
    standard SingleR reference for human data.  It is fetched from the same
    gypsum / ArtifactDB backend used by the R ``celldex`` / ``SingleR``
    packages, guaranteeing version-locked reproducibility.

    Requires
    --------
    celldex >= 0.3.0 (``pip install celldex``)

    Returns
    -------
    AnnData
        obs_names = cell type labels (``label.main``),
        var_names = gene symbols,
        X = mean log-normalised expression per cell type (dense float32).
    """
    if _HPCA_CACHE.exists():
        logger.info(f"Loading HPCA reference from cache: {_HPCA_CACHE}")
        return sc.read_h5ad(_HPCA_CACHE)

    try:
        from celldex import fetch_reference  # type: ignore
    except ImportError:
        raise ImportError(
            "celldex is required to download the HPCA reference. "
            "Install with: pip install celldex"
        )

    logger.info(
        f"Downloading HPCA reference via celldex "
        f"(name='{_CELLDEX_HPCA_NAME}', version='{_CELLDEX_HPCA_VERSION}') ..."
    )
    try:
        se = fetch_reference(_CELLDEX_HPCA_NAME, version=_CELLDEX_HPCA_VERSION)
    except Exception as e:
        raise RuntimeError(
            f"Failed to download HPCA reference via celldex: {e}. "
            f"Check your internet connection or run: pip install celldex"
        ) from e

    # se is a SummarizedExperiment — extract logcounts matrix and metadata
    # Gene names are in row_names; cell-type labels in col_data["label.main"]
    X_log = se.assay("logcounts")
    if hasattr(X_log, "toarray"):
        X_log = X_log.toarray()
    X_log = np.asarray(X_log, dtype=np.float32)  # shape: genes × samples

    col_data  = se.column_data
    row_names = list(se.row_names)

    # Each column is a sample; aggregate to pseudobulk per label.main
    labels_all = list(col_data["label.main"])
    unique_labels = list(dict.fromkeys(labels_all))  # preserve insertion order

    pb_rows = []
    for lbl in unique_labels:
        col_mask = np.array([l == lbl for l in labels_all])
        pb_rows.append(X_log[:, col_mask].mean(axis=1))  # mean over samples for this type

    X_pb = np.vstack(pb_rows).astype(np.float32)  # shape: n_types × n_genes

    ref = sc.AnnData(X=X_pb)
    ref.obs_names = unique_labels
    ref.var_names = row_names

    # Cache for future calls
    _SINGLER_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    ref.write_h5ad(_HPCA_CACHE)
    logger.info(
        f"  Cached {len(unique_labels)} cell types × {len(row_names)} genes "
        f"at: {_HPCA_CACHE}"
    )

    return ref


# Keep old name as internal alias (used in cache test via _TABULA_SAPIENS_CACHE symbol)
_load_tabula_sapiens = _load_hpca_reference


def _load_hca_census(tissue: str) -> sc.AnnData:
    """
    Stream mean log-normalised expression per cell type from the HCA Census API.

    Requires cellxgene-census (gated behind ImportError). Queries the census
    for cells matching the given tissue, computes pseudobulk mean per cell_type,
    and returns an AnnData in the standard pseudobulk format.

    Parameters
    ----------
    tissue : str
        Tissue name as it appears in the HCA Census obs metadata.

    Returns
    -------
    AnnData
        obs_names = cell type labels, var_names = genes,
        X = mean log-normalised expression per cell type (dense float32).
    """
    try:
        import cellxgene_census  # type: ignore
    except ImportError:
        raise ImportError(
            "cellxgene-census is required for singler_ref='hca'. "
            "Install with: pip install cellxgene-census"
        )

    logger.info(f"Streaming HCA Census reference for tissue='{tissue}' ...")

    with cellxgene_census.open_soma() as census:
        adata_tissue = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=f"tissue_general == '{tissue}'",
            column_names={"obs": ["cell_type"]},
        )

    if len(adata_tissue) == 0:
        raise ValueError(
            f"No cells found in HCA Census for tissue='{tissue}'. "
            f"Check the tissue name against available values."
        )

    # Build pseudobulk: mean expression per cell type
    sc.pp.normalize_total(adata_tissue, target_sum=1e4)
    sc.pp.log1p(adata_tissue)

    cell_types = adata_tissue.obs["cell_type"].unique().tolist()
    X_raw = adata_tissue.X
    if hasattr(X_raw, "toarray"):
        X_raw = X_raw.toarray()
    X_raw = np.asarray(X_raw, dtype=np.float32)

    pb_rows = []
    for ct in cell_types:
        mask = (adata_tissue.obs["cell_type"] == ct).values
        pb_rows.append(X_raw[mask, :].mean(axis=0))

    X_pb = np.vstack(pb_rows).astype(np.float32)
    ref  = sc.AnnData(X=X_pb)
    ref.obs_names  = cell_types
    ref.var_names  = adata_tissue.var_names.tolist()

    logger.info(f"  HCA Census reference: {len(cell_types)} cell types, {ref.n_vars} genes")
    return ref


def _load_singler_ref(
    singler_ref: Optional[Union[str, Path]],
    singler_ref_label_col: str,
    tissue: str = "",
) -> sc.AnnData:
    """
    Route to the correct reference loader based on singler_ref value.

    Returns
    -------
    AnnData
        obs_names = cell type labels, var_names = genes,
        X = mean log-normalised expression (pseudobulk).
    """
    if singler_ref is None:
        ref = _load_hpca_reference()

    elif isinstance(singler_ref, str) and singler_ref.lower() == "hca":
        ref = _load_hca_census(tissue)

    else:
        # User-supplied H5AD file
        path = Path(singler_ref)
        if not path.exists():
            raise FileNotFoundError(
                f"SingleR reference file not found: {path}. "
                f"Check the path passed to singler_ref."
            )
        logger.info(f"Loading user-supplied SingleR reference from {path} ...")
        ref = sc.read_h5ad(path)

        # If obs contains a cell_type column, use it as obs_names
        if singler_ref_label_col in ref.obs.columns:
            ref.obs_names = ref.obs[singler_ref_label_col].astype(str).values
        # else obs_names are already the cell type labels

    # Validate reference has enough cell types
    if ref.n_obs < 5:
        raise ValueError(
            f"SingleR reference has only {ref.n_obs} cell types — need at least 5. "
            f"Check the reference source."
        )

    return ref


def _run_singler(
    adata: sc.AnnData,
    leiden_col: str,
    singler_ref: Optional[Union[str, Path]] = None,
    singler_ref_label_col: str = "cell_type",
    num_threads: int = 1,
) -> None:
    """
    Annotate cells using the ``singler`` Python package (C++ SingleR bindings).

    Delegates all score computation, pairwise marker selection, iterative
    fine-tuning, and delta-based pruning to ``singler.annotate_single``.
    Cells whose delta score is ≤ 0 are labelled "Unassigned" (SingleR's own
    pruning convention); their ``singler_delta`` value is stored as NaN.

    Reference loading uses the 3-tier hierarchy defined in ``_load_singler_ref``:
      - None (default) → HPCA via ``celldex`` (``pip install celldex``)
      - "hca"          → HCA Census API (``pip install cellxgene-census``)
      - Path / str     → user-supplied H5AD pseudobulk reference

    Parameters
    ----------
    adata : AnnData
        Must have ``adata.obs[leiden_col]`` and ``adata.layers["logcounts"]``
        (or ``adata.X``) with log-normalised expression.
    leiden_col : str
        obs column with cluster labels.
    singler_ref : str, Path, or None
        Reference source — see ``_load_singler_ref()``.
    singler_ref_label_col : str
        obs column in user-supplied H5AD that contains cell type labels.
    num_threads : int
        Number of threads to pass to ``singler.annotate_single`` (default 1).

    Writes
    ------
    obs['cell_type_singler']  — per-cell label (str, category dtype).
                                "Unassigned" for low-confidence cells.
    obs['singler_delta']      — per-cell delta score (float).
                                NaN for "Unassigned" cells.
    """
    try:
        import singler as _singler  # type: ignore
    except ImportError:
        raise ImportError(
            "singler is required for SingleR annotation. "
            "Install with: pip install singler"
        )

    # ── Load reference (routes to HPCA / HCA Census / user H5AD) ─────────────
    ref = _load_singler_ref(singler_ref, singler_ref_label_col)

    # ── Build test matrix: genes × cells (singler convention) ─────────────────
    if "logcounts" in adata.layers:
        X_test = adata.layers["logcounts"]
    else:
        X_test = adata.X

    if hasattr(X_test, "toarray"):
        X_test = X_test.toarray()
    X_test = np.asarray(X_test, dtype=np.float32).T  # genes × cells

    test_features = list(adata.var_names)

    # ── Build ref matrix and labels ───────────────────────────────────────────
    # ref is an AnnData: obs = cell types (samples), var = genes.
    # singler expects ref as genes × samples.
    X_ref = ref.X
    if hasattr(X_ref, "toarray"):
        X_ref = X_ref.toarray()
    X_ref = np.asarray(X_ref, dtype=np.float32).T  # genes × cell-type-samples

    ref_features = list(ref.var_names)
    ref_labels   = list(ref.obs_names)  # one label per column of X_ref

    # ── Run SingleR ───────────────────────────────────────────────────────────
    logger.info(
        f"Running SingleR: {X_test.shape[1]} cells × "
        f"{len(ref_labels)} reference labels ..."
    )
    result = _singler.annotate_single(
        test_data=X_test,
        ref_data=X_ref,
        ref_labels=ref_labels,
        test_features=test_features,
        ref_features=ref_features,
        num_threads=num_threads,
    )

    # ── Extract per-cell results ───────────────────────────────────────────────
    best   = list(result.column("best"))    # per-cell label strings
    deltas = np.asarray(result.column("delta"), dtype=np.float64)

    # Pruning: singler signals low-confidence as delta ≤ 0 (or None)
    labels_final = []
    deltas_final = []
    for lbl, d in zip(best, deltas):
        if lbl is None or (np.isfinite(d) and d <= 0):
            labels_final.append("Unassigned")
            deltas_final.append(float("nan"))
        else:
            labels_final.append(str(lbl))
            deltas_final.append(float(d))

    # ── Write obs columns ─────────────────────────────────────────────────────
    adata.obs["cell_type_singler"] = pd.Categorical(labels_final)
    adata.obs["singler_delta"]     = np.array(deltas_final, dtype=np.float64)

    n_types    = adata.obs["cell_type_singler"].nunique()
    n_assigned = sum(l != "Unassigned" for l in labels_final)
    print(
        f"  SingleR (singler {_singler.__version__}): {n_types} types, "
        f"{n_assigned}/{len(adata)} cells assigned"
    )


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
