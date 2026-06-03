"""
OmicSage — CITE-seq Pipeline Runner
=====================================
Runs the Phase 4 CITE-seq pipeline from the terminal, with the same
step-control, checkpointing, and config-driven behaviour as run_pipeline.py.

Usage
-----
  # Run all enabled steps
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml

  # Run a single step
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --step normalize_adt

  # Stop before integration (run ADT steps only)
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --to-step annotate_adt

  # Resume from a specific step
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --from-step harmony_adt

  # Run a range of steps
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --from-step reduce_adt --to-step harmony_adt

  # Force re-run even if cached output exists
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --step normalize_adt --force

  # Force re-run from a step onward
  python run_cite_pipeline.py --config config/runs/GSE194122_cite.yaml --from-step harmony_adt --force

Step order
----------
  normalize_adt → doublets → reduce_adt → harmony_adt → annotate_adt → integration

Checkpointing
-------------
  Each step writes its output to processed_dir/cite_NN_<step>.h5ad (AnnData).
  The integration step writes a MuData (.h5mu).
  If the output already exists the step is skipped (cached).
  Use --force to override.

Input files (set in config under paths)
----------------------------------------
  adt_input   : output of QC step, e.g. data/processed/GSE194122/01_qc_adt.h5ad
  rna_input   : annotated RNA, e.g. data/processed/GSE194122/05_annotated.h5ad
                (only needed for the integration step)
"""

import sys
import io
import logging

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8",
                              errors="replace", line_buffering=True)

import argparse
import os
import warnings
from datetime import datetime
from pathlib import Path
from typing import Optional

import yaml

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

logging.getLogger("streamlit").setLevel(logging.ERROR)

# ── repo root ──────────────────────────────────────────────────────────────────
root = Path(__file__).resolve().parent
while not (root / "pipeline").exists() and root != root.parent:
    root = root.parent
os.chdir(root)
sys.path.insert(0, str(root))


# ── Step registry ──────────────────────────────────────────────────────────────

STEP_ORDER = [
    "normalize_adt",
    "doublets",
    "reduce_adt",
    "harmony_adt",
    "annotate_adt",
    "integration",
    "deg_cite",
    "gsea_cite",
]

# Output filename per step
STEP_OUTPUT = {
    "normalize_adt": "cite_01_normalized_adt.h5ad",
    "doublets":      "cite_02_doublets_adt.h5ad",
    "reduce_adt":    "cite_03_reduced_adt.h5ad",
    "harmony_adt":   "cite_04_harmony_adt.h5ad",
    "annotate_adt":  "cite_05_annotated_adt.h5ad",
    "integration":   "cite_06_integration.h5mu",   # MuData
    "deg_cite":      "cite_07_deg.h5mu",            # MuData (RNA + ADT, DEG results in uns)
    "gsea_cite":     "cite_08_gsea.h5mu",           # MuData (GSEA results in uns)
}

# Predecessor for each step (None = reads from adt_input path directly)
STEP_PREDECESSOR = {
    "normalize_adt": None,
    "doublets":      "normalize_adt",
    "reduce_adt":    "doublets",
    "harmony_adt":   "reduce_adt",
    "annotate_adt":  "harmony_adt",
    "integration":   "annotate_adt",
    "deg_cite":      "integration",   # reads cite_06 MuData (RNA + ADT + labels)
    "gsea_cite":     "deg_cite",      # reads cite_07 MuData (DEG results in uns)
}


# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f)


def get_step_cfg(cfg: dict, step: str) -> dict:
    step_cfg = cfg.get("steps", {}).get(step, {})
    if isinstance(step_cfg, bool):
        return {"enabled": step_cfg, "params": {}}
    return {
        "enabled": step_cfg.get("enabled", True),
        "params":  step_cfg.get("params", {}),
    }


# ══════════════════════════════════════════════════════════════════════════════
# INPUT RESOLUTION
# ══════════════════════════════════════════════════════════════════════════════

def resolve_input(step: str, cfg: dict, processed_dir: Path) -> Path:
    """
    Return the input path for a step.
      - normalize_adt reads from paths.adt_input (QC output)
      - integration reads BOTH paths.adt_input_for_integration (annotate_adt
        output, resolved automatically) AND paths.rna_input (annotated RNA)
      - all other steps read the predecessor's output
    Raises FileNotFoundError with a clear message if nothing is found.
    """
    pred = STEP_PREDECESSOR[step]

    # normalize_adt: read raw ADT input
    if pred is None:
        raw = cfg.get("paths", {}).get("adt_input")
        if not raw:
            raise ValueError(
                f"[{step}] paths.adt_input is not set in the config."
            )
        p = Path(raw)
        if not p.exists():
            raise FileNotFoundError(
                f"[{step}] adt_input not found: {p}\n"
                f"  Run the RNA pipeline QC step first to produce this file."
            )
        return p

    # all other steps: predecessor output
    pred_out = processed_dir / STEP_OUTPUT[pred]
    if pred_out.exists():
        return pred_out
    # Also check .h5mu extension (MuData checkpoints)
    pred_out_mu = processed_dir / STEP_OUTPUT[pred]
    if pred_out_mu.exists():
        return pred_out_mu
    raise FileNotFoundError(
        f"[{step}] requires output of '{pred}' at {pred_out}\n"
        f"  Run '{pred}' first, or use --from-step to include it in this run."
    )


def resolve_rna_input(cfg: dict, processed_dir: Path) -> Path:
    """Resolve the annotated RNA AnnData for the integration step."""
    rna_path = cfg.get("paths", {}).get("rna_input")
    if not rna_path:
        raise ValueError(
            "[integration] paths.rna_input is not set in the config.\n"
            "  Set it to the annotated RNA file, e.g.:\n"
            "    paths:\n"
            "      rna_input: data/processed/GSE194122/05_annotated.h5ad"
        )
    p = Path(rna_path)
    if not p.exists():
        raise FileNotFoundError(
            f"[integration] rna_input not found: {p}\n"
            f"  Run the RNA pipeline first (through the annotate step)."
        )
    return p


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def validate_plan(active_steps: list, cfg: dict, processed_dir: Path) -> None:
    errors = []
    for step in active_steps:
        pred = STEP_PREDECESSOR[step]
        if pred is not None and pred in active_steps:
            continue
        try:
            resolve_input(step, cfg, processed_dir)
        except (FileNotFoundError, ValueError) as e:
            errors.append(str(e))

    # Check rna_input separately if integration is active
    if "integration" in active_steps:
        try:
            resolve_rna_input(cfg, processed_dir)
        except (FileNotFoundError, ValueError) as e:
            errors.append(str(e))

    if errors:
        print("\n[OmicSage CITE] Validation failed — fix these before running:\n")
        for e in errors:
            print(f"  ✗  {e}\n")
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# STEP RUNNERS
# ══════════════════════════════════════════════════════════════════════════════

def _load_adt(path: Path):
    import anndata as ad
    return ad.read_h5ad(path)


def _wrap_mdata(adt):
    """Wrap a bare AnnData into a minimal MuData for modules that expect MuData."""
    import mudata as mu
    return mu.MuData({"adt": adt})


def _unwrap_adt(mdata):
    """Extract mdata['adt'] after a step that returns MuData."""
    return mdata["adt"]


# ── Step 1: normalize_adt ─────────────────────────────────────────────────────

def run_normalize_adt(input_path: Path, out: Path, reports_dir: Path,
                      params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[normalize_adt] cached → {out}")
        return out

    print("[normalize_adt] running …")
    from pipeline.modules.cite.adt_normalize import normalize_adt

    adt = _load_adt(input_path)
    adt_norm, metrics = normalize_adt(
        adt,
        clr_axis=params.get("clr_axis", 0),
        inplace=False,
    )

    adt_norm.write_h5ad(out)
    print(
        f"[normalize_adt] {metrics['n_cells']:,} cells"
        f"  {metrics['n_proteins']} proteins"
        f"  CLR axis={metrics['clr_axis']} → {out}"
    )
    try:
        from reports.templates.cite.cite_normalize_report import run_cite_normalize_report
        run_cite_normalize_report(
            adt=adt_norm,
            metrics=metrics,
            report_path=str(reports_dir / "cite_01_normalize_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[normalize_adt] WARNING: report failed: {e}")
    return out


# ── Step 2: doublets ──────────────────────────────────────────────────────────

def run_doublets(input_path: Path, out: Path, reports_dir: Path,
                 params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[doublets] cached → {out}")
        return out

    print("[doublets] running …")
    import mudata as mu
    from pipeline.modules.cite.adt_doublets import detect_adt_doublets

    adt = _load_adt(input_path)
    mdata = _wrap_mdata(adt)
    mdata, metrics = detect_adt_doublets(
        mdata,
        threshold=params.get("threshold", 2.5),
        filter_doublets=params.get("filter_doublets", False),
        inplace=True,
    )

    mdata["adt"].write_h5ad(out)
    print(
        f"[doublets] {metrics['n_doublets_detected']:,} doublets flagged"
        f"  ({metrics['pct_doublets']:.1f}%)"
        f"  filter={metrics['filter_doublets']} → {out}"
    )
    try:
        from reports.templates.cite.cite_doublets_report import run_cite_doublets_report
        run_cite_doublets_report(
            adt=mdata["adt"],
            metrics=metrics,
            report_path=str(reports_dir / "cite_02_doublets_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[doublets] WARNING: report failed: {e}")
    return out


# ── Step 3: reduce_adt ────────────────────────────────────────────────────────

def run_reduce_adt(input_path: Path, out: Path, reports_dir: Path,
                   params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[reduce_adt] cached → {out}")
        return out

    print("[reduce_adt] running …")
    import mudata as mu
    from pipeline.modules.cite.adt_reduce import reduce_adt

    adt = _load_adt(input_path)
    mdata = _wrap_mdata(adt)
    adt_reduced, metrics = reduce_adt(
        mdata,
        n_comps=params.get("n_comps", 50),
        n_pcs=params.get("n_pcs", 20),
        n_neighbors=params.get("n_neighbors", 15),
        isotype_controls=params.get("isotype_controls") or None,
        umap_color_keys=params.get("umap_color_keys") or None,
        inplace=False,
    )

    adt_reduced.write_h5ad(out)
    print(
        f"[reduce_adt] {metrics['n_comps_actual']} PCs"
        f"  UMAP computed → {out}"
    )
    try:
        from reports.templates.cite.cite_reduce_report import run_cite_reduce_report
        run_cite_reduce_report(
            adt=adt_reduced,
            metrics=metrics,
            report_path=str(reports_dir / "cite_03_reduce_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[reduce_adt] WARNING: report failed: {e}")
    return out


# ── Step 4: harmony_adt ───────────────────────────────────────────────────────

def run_harmony_adt(input_path: Path, out: Path, reports_dir: Path,
                    params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[harmony_adt] cached → {out}")
        return out

    print("[harmony_adt] running …")
    import mudata as mu
    from pipeline.modules.cite.adt_harmony import run_harmony_adt as _harmony

    adt = _load_adt(input_path)
    mdata = _wrap_mdata(adt)
    batch_key = params.get("batch_key",
                cfg.get("cite", {}).get("batch_key", "batch"))
    adt_harmonised, metrics = _harmony(
        mdata,
        batch_key=batch_key,
        n_pcs=params.get("n_pcs", 20),
        n_neighbors=params.get("n_neighbors", 15),
        random_state=params.get("random_state", 0),
        inplace=False,
    )

    adt_harmonised.write_h5ad(out)
    print(
        f"[harmony_adt] {metrics['n_batches']} batches"
        f"  batch_key={batch_key} → {out}"
    )
    try:
        from reports.templates.cite.cite_harmony_report import run_cite_harmony_report
        run_cite_harmony_report(
            adt=adt_harmonised,
            metrics=metrics,
            report_path=str(reports_dir / "cite_04_harmony_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[harmony_adt] WARNING: report failed: {e}")
    return out


# ── Step 5: annotate_adt ──────────────────────────────────────────────────────

def run_annotate_adt(input_path: Path, out: Path, reports_dir: Path,
                     params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[annotate_adt] cached → {out}")
        return out

    print("[annotate_adt] running …")
    import mudata as mu
    from pipeline.modules.cite.adt_annotate import annotate_adt

    adt = _load_adt(input_path)
    mdata = _wrap_mdata(adt)

    # annotation_map: convert string keys from YAML to the right type
    annotation_map = params.get("annotation_map")
    if annotation_map:
        annotation_map = {str(k): str(v) for k, v in annotation_map.items()}

    adt_annotated, metrics = annotate_adt(
        mdata,
        preset=params.get("preset"),  
        annotation_map=annotation_map,
        resolution=params.get("resolution", 0.1),
        n_iterations=params.get("n_iterations", 2),
        random_state=params.get("random_state", 0),
        inplace=False,
    )

    adt_annotated.write_h5ad(out)
    annotated_str = (
        f"  cell types annotated" if metrics["annotated"] else "  no annotation_map"
    )
    print(
        f"[annotate_adt] {metrics['n_clusters']} Leiden clusters"
        f"{annotated_str} → {out}"
    )
    try:
        from reports.templates.cite.cite_annotate_report import run_cite_annotate_report
        run_cite_annotate_report(
            adt=adt_annotated,
            metrics=metrics,
            report_path=str(reports_dir / "cite_05_annotate_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[annotate_adt] WARNING: report failed: {e}")
    return out


# ── Step 6: integration ───────────────────────────────────────────────────────

def run_integration(input_path: Path, out: Path, reports_dir: Path,
                    params: dict, cfg: dict, processed_dir: Path,
                    force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[integration] cached → {out}")
        return out

    print("[integration] running …")
    import anndata as ad
    import mudata as mu
    from pipeline.modules.cite.cite_integration import run_mofa, run_totalvi, run_both

    # Load ADT (annotated) and RNA (annotated with cell_type_vote)
    adt = _load_adt(input_path)
    rna_path = resolve_rna_input(cfg, processed_dir)
    rna = ad.read_h5ad(rna_path)
    print(f"[integration] RNA loaded: {rna.shape}  ADT: {adt.shape}")

    # Align barcodes — keep intersection
    common = adt.obs_names.intersection(rna.obs_names)
    if len(common) == 0:
        raise ValueError(
            "[integration] RNA and ADT share no cell barcodes. "
            "Check that both files come from the same dataset and QC run."
        )
    if len(common) < len(adt.obs_names):
        print(
            f"[integration] WARNING: {len(adt.obs_names) - len(common):,} ADT cells"
            f" not found in RNA — keeping {len(common):,} shared cells."
        )
    adt = adt[common].copy()
    rna = rna[common].copy()

    mdata = mu.MuData({"rna": rna, "adt": adt})

    method        = params.get("method", "mofa").lower()
    batch_key     = params.get("batch_key",
                    cfg.get("cite", {}).get("batch_key", "batch"))
    compute_scib  = params.get("compute_scib", False)
    cell_type_key = params.get("cell_type_key", None)

    if method == "mofa":
        print(f"[integration] MOFA+  batch_key={batch_key}  n_factors={params.get('n_factors', 15)}")
        mdata, metrics = run_mofa(
            mdata,
            batch_key=batch_key,
            n_factors=params.get("n_factors", 15),
            random_state=params.get("random_state", 0),
            inplace=True,
            compute_scib=compute_scib,
            cell_type_key=cell_type_key,
        )
    elif method == "totalvi":
        print(f"[integration] totalVI  batch_key={batch_key}  max_epochs={params.get('max_epochs', 400)}")
        mdata, metrics = run_totalvi(
            mdata,
            batch_key=batch_key,
            max_epochs=params.get("max_epochs", 400),
            random_state=params.get("random_state", 0),
            inplace=True,
            compute_scib=compute_scib,
            cell_type_key=cell_type_key,
        )
    elif method == "both":
        print(f"[integration] MOFA+ + totalVI  batch_key={batch_key}")
        mdata, metrics = run_both(
            mdata,
            batch_key=batch_key,
            n_factors=params.get("n_factors", 15),
            max_epochs=params.get("max_epochs", 400),
            random_state=params.get("random_state", 0),
            inplace=True,
            compute_scib=compute_scib,
            cell_type_key=cell_type_key,
        )
    else:
        raise ValueError(
            f"[integration] Unknown method '{method}'. "
            f"Use 'mofa', 'totalvi', or 'both'."
        )

    mdata.uns["dataset_name"] = cfg["dataset"].get("name", cfg["dataset"]["id"])

    # mudata >=0.4 emits FutureWarnings from internal .update() during write
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        mdata.write(out)
    print(f"[integration] {method.upper()}  {mdata.n_obs:,} cells → {out}")

    # Per-step HTML report
    try:
        from reports.templates.cite.cite_integration_report import run_cite_integration_report
        run_cite_integration_report(
            mdata=mdata,
            metrics=metrics,
            report_path=str(reports_dir / "cite_06_integration_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
            color_keys=params.get("umap_color_keys") or None,
        )
    except Exception as e:
        print(f"[integration] WARNING: report failed: {e}")

    return out


# ── Step 7: deg_cite ──────────────────────────────────────────────────────────

def run_deg_cite(input_path: Path, out: Path, reports_dir: Path,
                 params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[deg_cite] cached → {out}")
        return out

    print("[deg_cite] running …")
    import mudata as mu
    from pipeline.modules.cite.cite_deg import cite_deg

    # cite_07 reads the MuData from cite_06 (RNA + ADT + labels)
    mdata = mu.read(str(input_path))

    groupby          = params.get("groupby", "adt_celltype_manual")
    groupby_fallback = params.get("groupby_fallback", "adt_celltype_score")

    mdata_deg, deg_dict = cite_deg(
        mdata,
        groupby=groupby,
        groupby_fallback=groupby_fallback,
        method=params.get("method", "wilcoxon"),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        n_genes=params.get("n_genes", 200),
        exclude_protein_prefixes=params.get("exclude_protein_prefixes") or [],
        exclude_gene_prefixes=params.get("exclude_gene_prefixes") or [],
        use_raw_rna=params.get("use_raw_rna", False),
        inplace=False,
    )

    prov = deg_dict["provenance"]
    print(
        f"[deg_cite] DPE: {prov['n_dpe_significant']} significant proteins"
        f"  RNA cross-modal: {prov['n_rna_crossmodal_significant']} genes"
        f"  ({prov['n_cell_types']} cell types) → {out}"
    )

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        mdata_deg.write(out)

    try:
        from reports.templates.cite.cite_deg_report import run_cite_deg_report
        run_cite_deg_report(
            deg_dict=deg_dict,
            report_path=str(reports_dir / "cite_07_deg_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[deg_cite] WARNING: report failed: {e}")
    return out


# ── Step 8: gsea_cite ─────────────────────────────────────────────────────────

def run_gsea_cite(input_path: Path, out: Path, reports_dir: Path,
                  params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[gsea_cite] cached → {out}")
        return out

    print("[gsea_cite] running …")
    import mudata as mu
    from pipeline.modules.cite.cite_gsea import cite_gsea

    # cite_08 reads the MuData from cite_07 (DEG results in uns)
    mdata = mu.read(str(input_path))

    # Reconstruct cite_deg_dict from uns — the runner only persists the MuData,
    # so we rebuild the dict from the stored provenance and uns keys.
    deg_dict = _reconstruct_deg_dict(mdata)

    if not deg_dict["rna_crossmodal"]:
        print(
            "[gsea_cite] WARNING: no cross-modal RNA DEG results found in "
            "cite_07 MuData. GSEA requires MuData from cite_06. Skipping."
        )
        # Write input through as output so pipeline can continue
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            mdata.write(out)
        return out

    mdata_gsea, gsea_dict = cite_gsea(
        mdata,
        cite_deg_dict=deg_dict,
        gene_sets=params.get("gene_sets") or None,
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        top_n_genes=params.get("top_n_genes") or None,
        min_genes=params.get("min_genes", 5),
        organism=params.get("organism", "human"),
        direction=params.get("direction", "up"),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes") or [],
        inplace=False,
    )

    prov = gsea_dict["provenance"]
    n_tested  = prov.get("n_groups_tested", 0)
    n_skipped = prov.get("n_groups_skipped", 0)
    print(
        f"[gsea_cite] {n_tested} groups tested  {n_skipped} skipped → {out}"
    )

    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        mdata_gsea.write(out)

    try:
        from reports.templates.cite.cite_gsea_report import run_cite_gsea_report
        run_cite_gsea_report(
            gsea_dict=gsea_dict,
            report_path=str(reports_dir / "cite_08_gsea_report.html"),
            dataset_name=cfg["dataset"].get("name", cfg["dataset"]["id"]),
        )
    except Exception as e:
        print(f"[gsea_cite] WARNING: report failed: {e}")
    return out


def _reconstruct_deg_dict(mdata) -> dict:
    """
    Reconstruct a minimal cite_deg_dict from a persisted MuData (cite_07).

    cite_deg writes DEG results into uns["rank_genes_groups_rna_cm"] on
    mdata["rna"] and stores provenance in mdata.uns["omicsage_cite_deg"].
    This function rebuilds the dict that cite_gsea() expects.
    """
    import numpy as np

    provenance = mdata.uns.get("omicsage_cite_deg", {})
    rna = mdata["rna"] if "rna" in mdata.mod else None

    rna_crossmodal: dict = {}

    if rna is not None and "rank_genes_groups_rna_cm" in rna.uns:
        rgg = rna.uns["rank_genes_groups_rna_cm"]
        groups = list(rgg["names"].dtype.names)
        n_genes = rgg["names"].shape[0]

        for group in groups:
            def _safe(field, g=group):
                if field not in rgg:
                    return np.zeros(n_genes)
                arr = rgg[field]
                if arr.dtype.names and g in arr.dtype.names:
                    return arr[g]
                return np.zeros(n_genes)

            df_data = {
                "gene":     rgg["names"][group],
                "score":    _safe("scores"),
                "pval":     _safe("pvals"),
                "logfc":    _safe("logfoldchanges"),
                "pval_adj": _safe("pvals_adj"),
            }
            import pandas as pd
            df = pd.DataFrame(df_data).dropna(subset=["gene"])
            df = df[df["gene"] != ""]
            rna_crossmodal[group] = df.reset_index(drop=True)

    return {
        "rna_crossmodal": rna_crossmodal,
        "provenance": provenance,
        "input_type": provenance.get("input_type", "mudata"),
    }


# ══════════════════════════════════════════════════════════════════════════════
# STEP → RUNNER MAP
# ══════════════════════════════════════════════════════════════════════════════

STEP_RUNNERS = {
    "normalize_adt": run_normalize_adt,
    "doublets":      run_doublets,
    "reduce_adt":    run_reduce_adt,
    "harmony_adt":   run_harmony_adt,
    "annotate_adt":  run_annotate_adt,
    "deg_cite":      run_deg_cite,
    "gsea_cite":     run_gsea_cite,
    # integration handled separately (needs extra args)
}


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="OmicSage CITE-seq pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to YAML config file (e.g. config/runs/GSE194122_cite.yaml)",
    )
    parser.add_argument(
        "--step",
        choices=STEP_ORDER + ["all"],
        default=None,
        help="Run exactly one step (or 'all' to run everything enabled)",
    )
    parser.add_argument(
        "--from-step",
        choices=STEP_ORDER,
        dest="from_step",
        default=None,
        help="Start from this step (inclusive)",
    )
    parser.add_argument(
        "--to-step",
        choices=STEP_ORDER,
        dest="to_step",
        default=None,
        help="Stop at this step (inclusive)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Re-run steps even if cached output exists",
    )
    return parser.parse_args()


def resolve_step_window(from_step, to_step, step):
    if step:
        if from_step or to_step:
            print("Error: --step cannot be combined with --from-step / --to-step")
            sys.exit(1)
        from_step = to_step = step

    from_idx = 0 if from_step is None else STEP_ORDER.index(from_step)
    to_idx   = len(STEP_ORDER) - 1 if to_step is None else STEP_ORDER.index(to_step)

    if from_idx > to_idx:
        print(f"Error: --from-step '{from_step}' comes after --to-step '{to_step}'")
        sys.exit(1)

    return from_idx, to_idx


def main():
    args = parse_args()

    if args.step == "all":
        args.step = None

    cfg          = load_config(args.config)
    dataset_id   = cfg["dataset"]["id"]
    dataset_name = cfg["dataset"].get("name", dataset_id)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    reports_dir   = Path(cfg["paths"]["reports_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    # Build active step list
    from_idx, to_idx = resolve_step_window(args.from_step, args.to_step, args.step)
    window       = STEP_ORDER[from_idx : to_idx + 1]
    active_steps = [s for s in window if get_step_cfg(cfg, s)["enabled"]]

    if not active_steps:
        print(f"No enabled steps in the requested range: {window}")
        sys.exit(0)

    # Validate inputs before touching anything
    validate_plan(active_steps, cfg, processed_dir)

    # Banner
    start_time = datetime.now()
    print("=" * 60)
    print(f"OmicSage CITE-seq — {dataset_name}")
    print(f"Steps    : {' → '.join(active_steps)}")
    print(f"Started  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    for step in active_steps:
        step_cfg   = get_step_cfg(cfg, step)
        params     = step_cfg["params"]
        input_path = resolve_input(step, cfg, processed_dir)
        out        = processed_dir / STEP_OUTPUT[step]

        if step == "integration":
            run_integration(
                input_path, out, reports_dir, params, cfg,
                processed_dir=processed_dir,
                force=args.force,
            )
        else:
            STEP_RUNNERS[step](
                input_path, out, reports_dir, params, cfg,
                force=args.force,
            )

    # Combined report — assembles all per-step reports into one tabbed file
    try:
        from reports.templates.cite.cite_combined_report import generate_cite_combined_report
        dataset_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
        combined_path = generate_cite_combined_report(
            reports_dir=reports_dir,
            dataset_name=dataset_name,
            output_path=reports_dir / "cite_00_combined_report.html",
        )
        if combined_path:
            print(f"[combined]   → {combined_path}")
    except Exception as e:
        print(f"[combined] WARNING: combined report failed: {e}")

    # Footer
    end_time = datetime.now()
    elapsed  = end_time - start_time
    hours, remainder = divmod(int(elapsed.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    print()
    print("=" * 60)
    print("CITE-seq pipeline complete.")
    print(f"Reports   : {reports_dir}")
    print(f"Elapsed   : {hours:02d}h {minutes:02d}m {seconds:02d}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
