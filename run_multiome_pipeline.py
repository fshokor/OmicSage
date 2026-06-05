"""
OmicSage — Multiome Pipeline Runner
=====================================
Runs the Phase 5 multiome pipeline from the terminal, with the same
step-control, checkpointing, and config-driven behaviour as run_cite_pipeline.py.

Usage
-----
  # Run all enabled steps
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml

  # Run a single step
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml --step atac_reduce

  # Stop before integration
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml --to-step atac_annotate

  # Resume from a specific step
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml --from-step atac_annotate

  # Run a range of steps
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml --from-step atac_reduce --to-step atac_annotate

  # Force re-run even if cached output exists
  python run_multiome_pipeline.py --config config/runs/GSE194122_atac_multiome.yaml --step atac_reduce --force

Step order
----------
  atac_qc → atac_reduce → atac_annotate → multiome_integration → multiome_deg

Checkpointing
-------------
  Each step writes its output to processed_dir/multiome_NN_<step>.<ext>.
  AnnData steps write .h5ad; MuData steps write .h5mu.
  If the output already exists the step is skipped (cached).
  Use --force to override.

Input files (set in config under paths)
----------------------------------------
  atac_input  : raw ATAC AnnData, e.g. data/processed/GSE194122/01_qc_atac.h5ad
                (the ATAC modality extracted from the original multiome h5 file)
  rna_input   : annotated RNA AnnData, e.g. data/processed/GSE194122/05_annotated.h5ad
                (needed for atac_annotate — label transfer — and multiome_integration)
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
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")

logging.getLogger("streamlit").setLevel(logging.ERROR)

# ── Repo root ──────────────────────────────────────────────────────────────────
root = Path(__file__).resolve().parent
while not (root / "pipeline").exists() and root != root.parent:
    root = root.parent
os.chdir(root)
sys.path.insert(0, str(root))


# ══════════════════════════════════════════════════════════════════════════════
# STEP REGISTRY
# ══════════════════════════════════════════════════════════════════════════════

STEP_ORDER = [
    "atac_qc",
    "atac_reduce",
    "atac_annotate",
    "multiome_integration",
    "multiome_deg",
    "multiome_grn",
]

# Output filename per step
STEP_OUTPUT = {
    "atac_qc":              "multiome_01_qc_atac.h5ad",
    "atac_reduce":          "multiome_02_reduce_atac.h5ad",
    "atac_annotate":        "multiome_03_annotate_atac.h5ad",
    "multiome_integration": "multiome_04_integration.h5mu",
    "multiome_deg":         "multiome_05_deg.h5mu",
    "multiome_grn":         "multiome_06_grn.h5mu",
}

# Predecessor for each step (None = reads from paths.atac_input)
STEP_PREDECESSOR = {
    "atac_qc":              None,
    "atac_reduce":          "atac_qc",
    "atac_annotate":        "atac_reduce",          # also needs rna_input
    "multiome_integration": "atac_annotate",         # also needs rna_input
    "multiome_deg":         "multiome_integration",
    "multiome_grn":         "multiome_deg",
}

# Steps that require rna_input in addition to their predecessor
_STEPS_NEEDING_RNA = {"atac_annotate", "multiome_integration"}


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
    Return the primary input path for a step.
      - atac_qc reads from paths.atac_input
      - all other steps read the predecessor's output
    Raises FileNotFoundError / ValueError with a clear message if not found.
    """
    pred = STEP_PREDECESSOR[step]

    if pred is None:
        raw = cfg.get("paths", {}).get("atac_input")
        if not raw:
            raise ValueError(
                f"[{step}] paths.atac_input is not set in the config."
            )
        p = Path(raw)
        if not p.exists():
            raise FileNotFoundError(
                f"[{step}] atac_input not found: {p}\n"
                f"  Extract the ATAC modality from the multiome h5 file first."
            )
        return p

    pred_out = processed_dir / STEP_OUTPUT[pred]
    if pred_out.exists():
        return pred_out
    raise FileNotFoundError(
        f"[{step}] requires output of '{pred}' at {pred_out}\n"
        f"  Run '{pred}' first, or use --from-step to include it."
    )


def resolve_rna_input(cfg: dict, step: str, processed_dir: Path) -> Path:
    """Resolve the annotated RNA AnnData (needed for atac_annotate + integration)."""
    rna_path = cfg.get("paths", {}).get("rna_input")
    if not rna_path:
        raise ValueError(
            f"[{step}] paths.rna_input is not set in the config.\n"
            "  Set it to the annotated RNA file, e.g.:\n"
            "    paths:\n"
            "      rna_input: data/processed/GSE194122/05_annotated.h5ad"
        )
    p = Path(rna_path)
    if not p.exists():
        raise FileNotFoundError(
            f"[{step}] rna_input not found: {p}\n"
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

    for step in active_steps:
        if step in _STEPS_NEEDING_RNA:
            try:
                resolve_rna_input(cfg, step, processed_dir)
            except (FileNotFoundError, ValueError) as e:
                errors.append(str(e))
            break  # only need to check once — same rna_input for all

    if errors:
        print("\n[OmicSage Multiome] Validation failed — fix these before running:\n")
        for e in errors:
            print(f"  ✗  {e}\n")
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

def _load_atac(path: Path):
    import anndata as ad
    return ad.read_h5ad(path)


def _load_rna(path: Path):
    import anndata as ad
    return ad.read_h5ad(path)


def _load_mdata(path: Path):
    import mudata as mu
    return mu.read_h5mu(path)


def _save_mdata(mdata, path: Path) -> None:
    mdata.write_h5mu(path)


def _dataset_name(cfg: dict) -> str:
    return cfg["dataset"].get("name", cfg["dataset"]["id"])


# ══════════════════════════════════════════════════════════════════════════════
# STEP RUNNERS
# ══════════════════════════════════════════════════════════════════════════════

# ── Step 1: atac_qc ───────────────────────────────────────────────────────────

def run_atac_qc(input_path: Path, out: Path, reports_dir: Path,
                params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[atac_qc] cached → {out}")
        return out

    print("[atac_qc] running …")
    from pipeline.modules.multiome.atac_qc import atac_qc

    atac = _load_atac(input_path)
    atac_qcd, metrics = atac_qc(
        atac,
        min_peaks=params.get("min_peaks", 750),
        max_peaks=params.get("max_peaks", 500_000),
        min_peak_counts=params.get("min_peak_counts", 1_500),
        max_peak_counts=params.get("max_peak_counts", 100_000),
        max_nucleosome_signal=params.get("max_nucleosome_signal", 2.0),
        min_cells=params.get("min_cells", 15),
        run_scrublet=params.get("run_scrublet", True),
        filter_cells=params.get("filter_cells", False),
        inplace=False,
    )

    atac_qcd.write_h5ad(out)
    print(
        f"[atac_qc] {metrics['n_cells_after']:,} cells"
        f"  {metrics['n_peaks_after']:,} peaks → {out}"
    )
    try:
        from reports.templates.multiome.atac_qc_report import run_atac_qc_report
        run_atac_qc_report(
            atac=atac_qcd,
            metrics=metrics,
            report_path=str(reports_dir / "multiome_01_qc_report.html"),
            dataset_name=_dataset_name(cfg),
        )
    except Exception as e:
        print(f"[atac_qc] WARNING: report failed: {e}")
    return out


# ── Step 2: atac_reduce ───────────────────────────────────────────────────────

def run_atac_reduce(input_path: Path, out: Path, reports_dir: Path,
                    params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[atac_reduce] cached → {out}")
        return out

    print("[atac_reduce] running …")
    from pipeline.modules.multiome.atac_reduce import atac_reduce

    atac = _load_atac(input_path)
    atac_reduced, metrics = atac_reduce(
        atac,
        n_components=params.get("n_components", 50),
        n_neighbors=params.get("n_neighbors", 15),
        leiden_resolution=params.get("leiden_resolution", 0.5),
        use_raw_counts=params.get("use_raw_counts", True),
        random_state=params.get("random_state", 0),
        inplace=False,
    )

    atac_reduced.write_h5ad(out)
    print(
        f"[atac_reduce] {metrics['n_lsi_components_used']} LSI components"
        f"  {metrics['n_leiden_clusters']} clusters → {out}"
    )
    try:
        from reports.templates.multiome.atac_reduce_report import run_atac_reduce_report
        run_atac_reduce_report(
            atac=atac_reduced,
            metrics=metrics,
            report_path=str(reports_dir / "multiome_02_reduce_report.html"),
            dataset_name=_dataset_name(cfg),
        )
    except Exception as e:
        print(f"[atac_reduce] WARNING: report failed: {e}")
    return out


# ── Step 3: atac_annotate ─────────────────────────────────────────────────────

def run_atac_annotate(input_path: Path, out: Path, reports_dir: Path,
                      params: dict, cfg: dict, processed_dir: Path,
                      force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[atac_annotate] cached → {out}")
        return out

    print("[atac_annotate] running …")
    from pipeline.modules.multiome.atac_annotate import annotate_atac

    atac    = _load_atac(input_path)
    rna_path = resolve_rna_input(cfg, "atac_annotate", processed_dir)
    rna     = _load_rna(rna_path)
    print(f"[atac_annotate] RNA loaded: {rna.shape}  ATAC: {atac.shape}")

    atac_annotated, metrics = annotate_atac(
        atac,
        rna=rna,
        promoter_upstream_bp=params.get("promoter_upstream_bp", 2000),
        min_peaks_per_gene=params.get("min_peaks_per_gene", 1),
        leiden_key=params.get("leiden_key", "atac_leiden"),
        rna_label_key=params.get("rna_label_key", "cell_type_vote"),
        inplace=True,
    )

    try:
        from reports.templates.multiome.atac_annotate_report import run_atac_annotate_report
        run_atac_annotate_report(
            adata=atac_annotated,
            metrics=metrics,
            report_path=str(reports_dir / "multiome_03_annotate_report.html"),
            dataset_name=_dataset_name(cfg),
        )
    except Exception as e:
        print(f"[atac_annotate] WARNING: report failed: {e}")
    
    # Strip heavy matrices — report already done, not needed downstream
    if "tf_idf" in atac_annotated.layers:
        del atac_annotated.layers["tf_idf"]
    if "gene_activity" in atac_annotated.obsm:
        del atac_annotated.obsm["gene_activity"]
    if "gene_activity_var_names" in atac_annotated.uns:
        del atac_annotated.uns["gene_activity_var_names"]
    if "counts" in atac_annotated.layers:
        atac_annotated.X = atac_annotated.layers["counts"]

    # Write lean checkpoint
    atac_annotated.write_h5ad(out)
    print(
        f"[atac_annotate] {metrics['n_genes_activity']:,} gene activity genes"
        f"  {metrics['n_rna_barcodes_matched']:,} barcodes matched → {out}"
    )
    
    return out


# ── Step 4: multiome_integration ──────────────────────────────────────────────

def run_multiome_integration(input_path: Path, out: Path, reports_dir: Path,
                             params: dict, cfg: dict, processed_dir: Path,
                             force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[multiome_integration] cached → {out}")
        return out

    print("[multiome_integration] running …")
    import scanpy as sc
    import mudata as mu
    from pipeline.modules.multiome.multiome_integration import run_mofa, run_multivi

    atac    = _load_atac(input_path)
    rna_path = resolve_rna_input(cfg, "multiome_integration", processed_dir)
    rna     = _load_rna(rna_path)
    print(f"[multiome_integration] RNA: {rna.shape}  ATAC: {atac.shape}")

    # Align barcodes — keep intersection
    common = atac.obs_names.intersection(rna.obs_names)
    if len(common) == 0:
        raise ValueError(
            "[multiome_integration] RNA and ATAC share no cell barcodes."
        )
    atac = atac[common].copy()
    rna  = rna[common].copy()

    mdata     = mu.MuData({"rna": rna, "atac": atac})
    method    = params.get("method", "mofa").lower()
    batch_key = params.get("batch_key",
                cfg.get("multiome", {}).get("batch_key", "batch"))

    # ── Subset to variable features before building MuData ─────────────
    # RNA — use HVGs already computed by the RNA pipeline
    if "highly_variable" in rna.var.columns:
        rna = rna[:, rna.var["highly_variable"]].copy()
        print(f"[multiome_integration] RNA after HVG subset: {rna.shape}")

    # ATAC — compute and subset to top variable peaks
    n_top_peaks = params.get("n_top_peaks", 20_000)
    sc.pp.highly_variable_genes(
        atac,
        n_top_genes=n_top_peaks,
        flavor="cell_ranger",
        inplace=True,
    )
    atac = atac[:, atac.var["highly_variable"]].copy()
    print(f"[multiome_integration] ATAC after HVP subset: {atac.shape}")
    # ───────────────────────────────────────────────────────────────────

    mdata = mu.MuData({"rna": rna, "atac": atac})
    # ───────────────────────────────────────────────────────────────────

    if method == "mofa":
        print(f"[multiome_integration] MOFA+  batch_key={batch_key}"
              f"  n_factors={params.get('n_factors', 15)}")
        mdata, metrics = run_mofa(
            mdata,
            batch_key=batch_key,
            n_factors=params.get("n_factors", 15),
            random_state=params.get("random_state", 0),
            inplace=True,
        )
    elif method == "multivi":
        print(f"[multiome_integration] MultiVI  batch_key={batch_key}"
              f"  max_epochs={params.get('max_epochs', 500)}")
        mdata, metrics = run_multivi(
            mdata,
            batch_key=batch_key,
            n_latent=params.get("n_latent", 20),
            max_epochs=params.get("max_epochs", 500),
            random_state=params.get("random_state", 0),
            inplace=True,
        )
    else:
        raise ValueError(
            f"[multiome_integration] Unknown method '{method}'. "
            "Use 'mofa' or 'multivi'."
        )

    _save_mdata(mdata, out)
    print(f"[multiome_integration] method={method}  {mdata.n_obs:,} cells → {out}")

    try:
        from reports.templates.multiome.multiome_integration_report import (
            run_multiome_integration_report,
        )
        run_multiome_integration_report(
            mdata=mdata,
            metrics=metrics,
            report_path=str(reports_dir / "multiome_04_integration_report.html"),
            dataset_name=_dataset_name(cfg),
            batch_key=batch_key,
        )
    except Exception as e:
        print(f"[multiome_integration] WARNING: report failed: {e}")
    return out


# ── Step 5: multiome_deg ──────────────────────────────────────────────────────

def run_multiome_deg(input_path: Path, out: Path, reports_dir: Path,
                     params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[multiome_deg] cached → {out}")
        return out

    print("[multiome_deg] running …")
    from pipeline.modules.multiome.multiome_deg import multiome_deg

    mdata = _load_mdata(input_path)
    mdata, deg_dict = multiome_deg(
        mdata,
        groupby=params.get("groupby", "atac_celltype"),
        leiden_fallback=params.get("leiden_fallback", "atac_leiden"),
        method=params.get("method", "wilcoxon"),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        n_genes=params.get("n_genes", 200),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes") or None,
        exclude_peak_prefixes=params.get("exclude_peak_prefixes") or None,
        inplace=True,
    )

    _save_mdata(mdata, out)
    prov = deg_dict["provenance"]
    print(
        f"[multiome_deg] {prov['n_cell_types']} cell types"
        f"  {prov['n_rna_significant']} RNA DEG"
        f"  {prov['n_dca_significant']} DCA hits → {out}"
    )
    try:
        from reports.templates.multiome.multiome_deg_report import run_multiome_deg_report
        run_multiome_deg_report(
            deg_dict=deg_dict,
            report_path=str(reports_dir / "multiome_05_deg_report.html"),
            dataset_name=_dataset_name(cfg),
        )
    except Exception as e:
        print(f"[multiome_deg] WARNING: report failed: {e}")
    return out


# ── Step 6: multiome_grn ──────────────────────────────────────────────────────
 
def run_multiome_grn(input_path: Path, out: Path, reports_dir: Path,
                     params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[multiome_grn] cached → {out}")
        return out
 
    print("[multiome_grn] running …")
    from pipeline.modules.multiome.multiome_grn import multiome_grn
    from pipeline.modules.multiome.multiome_deg  import multiome_deg as _deg   # noqa: F401
 
    mdata = _load_mdata(input_path)
 
    # deg_dict is not cached separately; reconstruct a minimal version
    # from the checkpoint (provenance is in mdata.uns["omicsage_deg"])
    deg_dict = {"provenance": mdata.uns.get("omicsage_deg", {})}
 
    mdata, grn_dict = multiome_grn(
        mdata,
        deg_dict=deg_dict,
        motif_db=params.get("motif_db", "jaspar"),
        groupby=params.get("groupby", "atac_celltype"),
        n_top_peaks=params.get("n_top_peaks", 500),
        min_cells=params.get("min_cells", 10),
        random_state=params.get("random_state", 0),
        inplace=True,
    )
 
    _save_mdata(mdata, out)
    prov = grn_dict["provenance"]["outputs"]
    print(
        f"[multiome_grn] {prov['n_tfs_rna']} RNA TFs"
        f"  {prov['n_tfs_atac']} ATAC TFs"
        f"  {prov['n_grn_edges']} edges → {out}"
    )
    try:
        from reports.templates.multiome.multiome_grn_report import run_multiome_grn_report
        run_multiome_grn_report(
            result=grn_dict,
            report_path=str(reports_dir / "multiome_06_grn_report.html"),
            dataset_name=_dataset_name(cfg),
        )
    except Exception as e:
        print(f"[multiome_grn] WARNING: report failed: {e}")
    return out

# ══════════════════════════════════════════════════════════════════════════════
# STEP → RUNNER MAP
# ══════════════════════════════════════════════════════════════════════════════

STEP_RUNNERS = {
    "atac_qc":     run_atac_qc,
    "atac_reduce": run_atac_reduce,
    "multiome_deg": run_multiome_deg,
    "multiome_grn": run_multiome_grn,
}

# Steps that need special handling (extra inputs beyond predecessor)
_SPECIAL_STEPS = {"atac_annotate", "multiome_integration"}


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="OmicSage Multiome pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to YAML config file (e.g. config/runs/GSE194122_multiome.yaml)",
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

    cfg           = load_config(args.config)
    dataset_id    = cfg["dataset"]["id"]
    dataset_name  = cfg["dataset"].get("name", dataset_id)
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
    print(f"OmicSage Multiome — {dataset_name}")
    print(f"Steps    : {' → '.join(active_steps)}")
    print(f"Started  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    for step in active_steps:
        step_cfg   = get_step_cfg(cfg, step)
        params     = step_cfg["params"]
        input_path = resolve_input(step, cfg, processed_dir)
        out        = processed_dir / STEP_OUTPUT[step]

        if step == "atac_annotate":
            run_atac_annotate(
                input_path, out, reports_dir, params, cfg,
                processed_dir=processed_dir,
                force=args.force,
            )
        elif step == "multiome_integration":
            run_multiome_integration(
                input_path, out, reports_dir, params, cfg,
                processed_dir=processed_dir,
                force=args.force,
            )
        else:
            STEP_RUNNERS[step](
                input_path, out, reports_dir, params, cfg,
                force=args.force,
            )

    # Combined report
    try:
        from reports.templates.multiome.multiome_combined_report import (
            generate_multiome_combined_report,
        )
        combined_path = generate_multiome_combined_report(
            reports_dir=reports_dir,
            dataset_name=dataset_name,
            output_path=reports_dir / "multiome_00_combined_report.html",
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
    print("Multiome pipeline complete.")
    print(f"Reports   : {reports_dir}")
    print(f"Elapsed   : {hours:02d}h {minutes:02d}m {seconds:02d}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
