"""
OmicSage — Generic Pipeline Runner
====================================
Usage
-----
  # Run all enabled steps (two equivalent forms)
  python run_pipeline.py --config config/runs/GSE194122.yaml
  python run_pipeline.py --config config/runs/GSE194122.yaml --step all

  # Stop at a checkpoint (inclusive)
  python run_pipeline.py --config config/runs/GSE194122.yaml --to-step cluster

  # Resume from a checkpoint (inclusive)
  python run_pipeline.py --config config/runs/GSE194122.yaml --from-step annotate

  # Run a specific range
  python run_pipeline.py --config config/runs/GSE194122.yaml --from-step normalize --to-step reduce

  # Run exactly one step
  python run_pipeline.py --config config/runs/GSE194122.yaml --step normalize

  # Run all steps with AI features enabled
  python run_pipeline.py --config config/runs/GSE166635.yaml --step all --ai

  # Run only the AI layer (pipeline steps already cached)
  python run_pipeline.py --config config/runs/GSE166635.yaml --ai-only

  # Run only one specific AI module (all pipeline steps already cached)
  python run_pipeline.py --config config/runs/GSE166635.yaml --ai-only --ai-module coherence_reviewer

  # Run AI layer starting from a specific module
  python run_pipeline.py --config config/runs/GSE166635.yaml --ai-only --ai-from-module narrative_generator

  # Run AI layer up to a specific module (inclusive)
  python run_pipeline.py --config config/runs/GSE166635.yaml --ai-only --ai-to-module deg_validator

  # Run a range of AI modules
  python run_pipeline.py --config config/runs/GSE166635.yaml --ai-only --ai-from-module cluster_annotator --ai-to-module coherence_reviewer

  # Disable AI entirely even when --ai is passed (useful for testing)
  python run_pipeline.py --config config/runs/GSE166635.yaml --step all --ai --ai-module off

AI module order
---------------
  clustering_advisor → cluster_annotator → pipeline_advisor → deg_validator
  → coherence_reviewer → downstream_suggester → narrative_generator

Checkpointing
-------------
  Every step writes its output to processed_dir/NN_<step>.h5ad.
  If the file already exists the step is skipped (cached).
  Use --from-step to force re-execution from a given step onward.

Resolution override
-------------------
  After inspecting clusters, add to the config:
      cluster:
        params:
          best_resolution_override: 0.8
  Then re-run:
      python run_pipeline.py --config config/runs/GSE194122.yaml --step cluster
      python run_pipeline.py --config config/runs/GSE194122.yaml --from-step annotate
"""

import sys
import io
import logging

# Force line-buffered stdout so print() appears immediately even when piped
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace", line_buffering=True)

import argparse
import os
import re
import warnings
from datetime import datetime
from pathlib import Path
from typing import Optional

import scanpy as sc
import yaml

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Suppress streamlit "missing ScriptRunContext" noise from biochatter imports
logging.getLogger("streamlit").setLevel(logging.ERROR)
logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context").setLevel(logging.ERROR)

# ── repo root ──────────────────────────────────────────────────────────────────
root = Path(__file__).resolve().parent
while not (root / "pipeline").exists() and root != root.parent:
    root = root.parent
os.chdir(root)
sys.path.insert(0, str(root))

# ── step order (canonical) ─────────────────────────────────────────────────────
STEP_ORDER = [
    "qc",
    "normalize",
    "reduce",
    "cluster",
    "annotate",
    "deg",
    "gsea",
    "harmony",
    "cluster_harmony",
    "pseudobulk",
]

# ── AI module execution order (canonical) ──────────────────────────────────────
AI_MODULE_ORDER = [
    "clustering_advisor",
    "cluster_annotator",
    "pipeline_advisor",
    "deg_validator",
    "coherence_reviewer",
    "downstream_suggester",
    "narrative_generator",
]

# Output filename for each step
STEP_OUTPUT = {
    "qc":              "01_qc.h5ad",
    "normalize":       "02_normalized.h5ad",
    "reduce":          "03_reduced.h5ad",
    "cluster":         "04_clustered.h5ad",
    "annotate":        "05_annotated.h5ad",
    "deg":             "06_deg.h5ad",
    "gsea":            "07_gsea.h5ad",
    "harmony":         "08_harmony.h5ad",
    "cluster_harmony": "09_harmony_clustered.h5ad",
    "pseudobulk":      "10_pseudobulk_deg.h5ad",
}

# Each step's required predecessor (None = needs raw input)
STEP_PREDECESSOR = {
    "qc":              None,
    "normalize":       "qc",
    "reduce":          "normalize",
    "cluster":         "reduce",
    "annotate":        "cluster",
    "deg":             "annotate",
    "gsea":            "deg",
    "harmony":         "gsea",
    "cluster_harmony": "harmony",
    "pseudobulk":      "annotate",   # runs off annotated, not harmony
}


# ══════════════════════════════════════════════════════════════════════════════
# CONFIG
# ══════════════════════════════════════════════════════════════════════════════

def load_config(config_path: str) -> dict:
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    return cfg


def get_step_cfg(cfg: dict, step: str) -> dict:
    """Return the step sub-config, defaulting enabled=True and empty params."""
    step_cfg = cfg.get("steps", {}).get(step, {})
    if isinstance(step_cfg, bool):
        # shorthand: step: false
        return {"enabled": step_cfg, "params": {}}
    return {
        "enabled": step_cfg.get("enabled", True),
        "params":  step_cfg.get("params", {}),
        "input_override": step_cfg.get("input_override"),
    }


# ══════════════════════════════════════════════════════════════════════════════
# INPUT RESOLUTION
# ══════════════════════════════════════════════════════════════════════════════

def resolve_input(step: str, cfg: dict, processed_dir: Path) -> Path:
    """
    Return the input path for a step.
    Priority:
      1. input_override in step config
      2. output of predecessor step (cached file)
      3. raw_input from paths (for qc only)
    Raises if nothing is found.
    """
    step_cfg = get_step_cfg(cfg, step)

    # 1. explicit override
    if step_cfg.get("input_override"):
        p = Path(step_cfg["input_override"])
        if not p.exists():
            raise FileNotFoundError(
                f"[{step}] input_override path not found: {p}"
            )
        return p

    # 2. predecessor output
    pred = STEP_PREDECESSOR[step]
    if pred is not None:
        pred_out = processed_dir / STEP_OUTPUT[pred]
        if pred_out.exists():
            return pred_out
        raise FileNotFoundError(
            f"[{step}] requires output of '{pred}' at {pred_out}\n"
            f"  Options:\n"
            f"    • Run '{pred}' first\n"
            f"    • Add 'input_override' under steps.{step} in your config"
        )

    # 3. raw input (qc only)
    raw = cfg.get("paths", {}).get("raw_input")
    if raw:
        p = Path(raw)
        if not p.exists():
            raise FileNotFoundError(
                f"[{step}] raw_input not found: {p}"
            )
        return p

    raise ValueError(
        f"[{step}] No input source found. "
        f"Set paths.raw_input or steps.{step}.input_override in your config."
    )


# ══════════════════════════════════════════════════════════════════════════════
# VALIDATION
# ══════════════════════════════════════════════════════════════════════════════

def validate_plan(active_steps: list[str], cfg: dict, processed_dir: Path) -> None:
    """
    Check that every active step will have an input before we start running.
    Skips the check for steps whose predecessor is also active in this run
    (their input will be produced during the run itself).
    Fails loudly with a clear message if anything is missing.
    """
    errors = []
    for step in active_steps:
        pred = STEP_PREDECESSOR[step]
        # If the predecessor is also running in this batch, input will be
        # produced during the run — no need to validate it exists now.
        if pred is not None and pred in active_steps:
            continue
        try:
            resolve_input(step, cfg, processed_dir)
        except (FileNotFoundError, ValueError) as e:
            errors.append(str(e))

    if errors:
        print("\n[OmicSage] Validation failed — fix these before running:\n")
        for e in errors:
            print(f"  ✗  {e}\n")
        sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# STEP RUNNERS
# ══════════════════════════════════════════════════════════════════════════════

def run_qc(input_path: Path, out: Path, reports_dir: Path,
           params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[qc] cached → {out}")
        return out

    print("[qc] running …")
    from pipeline.modules.qc.ingest import load_dataset
    from pipeline.modules.qc.qc import run_qc
    from pipeline.modules.qc.data_report import run_data_report

    sample_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
    adata = load_dataset(input_path, sample_name=sample_name)

    # Data intake report — generated from raw data, before any QC filtering
    _dataset_id    = cfg["dataset"]["id"]
    _geo_accession = _dataset_id if re.match(r"^GSE\d+$", _dataset_id, re.IGNORECASE) else None
    _data_report_path = reports_dir / "00_data_report.html"
    try:
        run_data_report(
            adata,
            input_path=input_path,
            output_path=_data_report_path,
            geo_accession=_geo_accession,
        )
        print(f"[data_report] → {_data_report_path}")
    except Exception as _exc:
        print(f"[data_report] WARNING: could not generate data intake report: {_exc}")

    mdata, metrics = run_qc(
        adata,
        min_genes=params.get("min_genes", 200),
        max_genes=params.get("max_genes", 2500),
        max_mt_pct=params.get("max_mt_pct", 5.0),
        remove_doublets=params.get("remove_doublets", True),
        generate_report=True,
        report_path=str(reports_dir / "01_qc_report.html"),
        sample_name=sample_name,
    )

    rna = mdata["rna"]
    rna.write_h5ad(out)
    if "adt" in mdata.mod:
        mdata["adt"].write_h5ad(out.parent / "01_qc_adt.h5ad")

    pass_rate = 100 * metrics["n_cells_output"] / metrics["n_cells_input"]
    print(f"[qc] {metrics['n_cells_output']:,} cells kept ({pass_rate:.1f}%) → {out}")
    return out


def run_normalize(input_path: Path, out: Path, reports_dir: Path,
                  params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[normalize] cached → {out}")
        return out

    print("[normalize] running …")
    from pipeline.modules.qc.normalize import normalize
    from pipeline.modules.qc.normalization_report import run_normalization_report

    adata = sc.read_h5ad(input_path)
    adata_norm, metrics = normalize(
        adata,
        batch_key=params.get("batch_key", "batch"),
        target_sum=params.get("target_sum", 1e4),
        n_top_genes=params.get("n_top_genes", 2000),
        hvg_flavor=params.get("hvg_flavor", "seurat"),
        inplace=False,
    )

    dataset_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
    run_normalization_report(
        adata_norm=adata_norm,
        metrics=metrics,
        report_path=str(reports_dir / "02_normalization_report.html"),
        dataset_name=dataset_name,
    )

    adata_norm.write_h5ad(out)
    print(f"[normalize] HVGs={metrics['n_hvg_selected']} → {out}")
    return out


def run_reduce(input_path: Path, out: Path, reports_dir: Path,
               params: dict, cfg: dict, force: bool = False, batch_key: Optional[str] = None) -> Path:
    if out.exists() and not force:
        print(f"[reduce] cached → {out}")
        return out

    print("[reduce] running …")
    from pipeline.modules.qc.reduce import reduce
    from pipeline.modules.qc.reduce_report import run_reduce_report

    adata = sc.read_h5ad(input_path)
    adata_reduced, metrics = reduce(
        adata,
        n_comps=params.get("n_comps", 50),
        n_pcs=params.get("n_pcs"),
        n_pcs_method=params.get("n_pcs_method", "elbow"),
        n_neighbors=params.get("n_neighbors", 15),
        inplace=False,
    )

    dataset_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
    run_reduce_report(
        adata_reduced=adata_reduced,
        metrics=metrics,
        report_path=str(reports_dir / "03_reduce_report.html"),
        dataset_name=dataset_name,
        batch_key=batch_key
    )

    adata_reduced.write_h5ad(out)
    prov = adata_reduced.uns["omicsage_reduce"]
    print(
        f"[reduce] {prov['n_pcs_used']} PCs"
        f" ({prov['cumulative_variance_explained_by_selected_pcs']*100:.1f}% var)"
        f" → {out}"
    )
    return out


def run_cluster(input_path: Path, out: Path, reports_dir: Path,
                params: dict, cfg: dict, force: bool = False) -> Path:
    # Always re-run cluster if best_resolution_override is set
    # (user inspected and made a decision — honour it)
    best_resolution_override = params.get("best_resolution_override")
    if out.exists() and best_resolution_override is None and not force:
        print(f"[cluster] cached → {out}")
        return out
    if out.exists() and best_resolution_override is not None:
        print(f"[cluster] best_resolution_override={best_resolution_override}, re-running …")

    print("[cluster] running …")
    from pipeline.modules.clustering.cluster import cluster
    from pipeline.modules.clustering.cluster_report import run_cluster_report

    adata = sc.read_h5ad(input_path)
    adata_clustered, metrics = cluster(
        adata,
        resolution_range=params.get("resolution_range", [0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2]),
        best_resolution_override=best_resolution_override,
        inplace=False,
    )

    dataset_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
    run_cluster_report(
        adata_clustered=adata_clustered,
        metrics=metrics,
        report_path=str(reports_dir / "04_cluster_report.html"),
        dataset_name=dataset_name,
    )

    adata_clustered.write_h5ad(out)
    print(
        f"[cluster] res={metrics['best_resolution']}"
        f"  {metrics['best_n_clusters']} clusters"
        f"  ({metrics['selection_reason']}) → {out}"
    )
    return out


def run_annotate(input_path: Path, out: Path, reports_dir: Path,
                 params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[annotate] cached → {out}")
        return out

    print("[annotate] running …")
    from pipeline.modules.annotation.annotate import annotate
    from pipeline.modules.annotation.annotate_report import run_annotate_report

    adata = sc.read_h5ad(input_path)
    adata_annotated, annotation_dict = annotate(
        adata,
        methods=params.get("methods", ["celltypist", "markers", "sctype", "singler", "vote"]),
        leiden_col=params.get("leiden_col", "leiden"),
        celltypist_models=params.get("celltypist_models", ["Immune_All_High.pkl"]),
        tissue=params.get("tissue", "Immune system"),
        sctype_db_path=params.get("sctype_db_path"),
        scanvi_model=params.get("scanvi_model"),
        singler_ref=params.get("singler_ref"),
        singler_ref_label_col=params.get("singler_ref_label_col", "cell_type"),
        inplace=False,
    )

    dataset_name = cfg["dataset"].get("name", cfg["dataset"]["id"])
    run_annotate_report(
        adata_annotated=adata_annotated,
        annotation_dict=annotation_dict,
        report_path=str(reports_dir / "05_annotate_report.html"),
        dataset_name=dataset_name,
    )

    adata_annotated.write_h5ad(out)
    n_types = adata_annotated.obs["cell_type_vote"].nunique()
    print(f"[annotate] {n_types} consensus types → {out}")
    return out


def run_deg(input_path: Path, out: Path, reports_dir: Path,
            params: dict, cfg: dict, force: bool = False) -> tuple[Path, dict]:
    print("[deg] running …")
    from pipeline.modules.downstream.deg import deg
    from pipeline.modules.downstream.deg_report import generate_deg_report

    adata = sc.read_h5ad(input_path)
    adata_deg, deg_dict = deg(
        adata,
        groupby=params.get("groupby", "cell_type_vote"),
        method=params.get("method", "wilcoxon"),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        n_genes=params.get("n_genes", 500),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes", []),
        inplace=False,
    )

    generate_deg_report(
        adata=adata_deg,
        deg_dict=deg_dict,
        output_path=str(reports_dir / "06_deg_report.html"),
        top_n_volcano=params.get("top_n_volcano", 10),
        top_n_dotplot=params.get("top_n_dotplot", 5),
        max_volcano_groups=params.get("max_volcano_groups", 9),
    )

    adata_deg.write_h5ad(out)
    total_sig = sum(len(df) for df in deg_dict["results"].values())
    print(f"[deg] {len(deg_dict['results'])} groups  {total_sig:,} DEGs → {out}")
    return out, deg_dict


def _reload_deg_dict(processed_dir: Path, params: dict) -> tuple[Path, dict]:
    """Reconstruct deg_dict from cached file (needed when gsea runs after skipping deg)."""
    from pipeline.modules.downstream.deg import deg

    deg_path = processed_dir / STEP_OUTPUT["deg"]
    adata = sc.read_h5ad(deg_path)
    _, deg_dict = deg(
        adata,
        groupby=params.get("groupby", "cell_type_vote"),
        method=params.get("method", "wilcoxon"),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        n_genes=params.get("n_genes", 500),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes", []),
        inplace=False,
    )
    return deg_path, deg_dict


def run_gsea(input_path: Path, out: Path, reports_dir: Path,
             params: dict, cfg: dict, deg_dict: dict, force: bool = False) -> Path:
    print("[gsea] running …")
    from pipeline.modules.downstream.gsea import gsea
    from pipeline.modules.downstream.gsea_report import generate_gsea_report

    adata = sc.read_h5ad(input_path)
    adata_gsea, gsea_dict = gsea(
        adata,
        deg_dict=deg_dict,
        gene_sets=params.get("gene_sets", [
            "GO_Biological_Process_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
        ]),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        top_n_genes=params.get("top_n_genes"),
        min_genes=params.get("min_genes", 5),
        organism=params.get("organism", cfg["dataset"].get("organism", "human")),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes", []),
        inplace=False,
    )

    generate_gsea_report(
        gsea_dict=gsea_dict,
        output_path=str(reports_dir / "07_gsea_report.html"),
        top_n_table=params.get("top_n_table", 5),
        top_n_bar=params.get("top_n_bar", 10),
    )

    adata_gsea.write_h5ad(out)
    prov = gsea_dict["provenance"]
    print(f"[gsea] {prov['n_groups_tested']} groups tested  {prov['n_groups_skipped']} skipped → {out}")
    return out


def run_harmony(input_path: Path, out: Path, reports_dir: Path,
                params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[harmony] cached → {out}")
        return out

    print("[harmony] running …")
    from pipeline.modules.integration.harmony_correct import harmony_correct
    from pipeline.modules.integration.harmony_report import generate_harmony_report

    adata = sc.read_h5ad(input_path)
    adata = harmony_correct(
        adata,
        batch_key=params.get("batch_key", "batch"),
        n_pcs=params.get("n_pcs", 50),
        n_neighbors=params.get("n_neighbors", 15),
        umap_min_dist=params.get("umap_min_dist", 0.3),
        random_state=params.get("random_state", 0),
        copy=False,
    )

    generate_harmony_report(
        adata=adata,
        output_path=str(reports_dir / "08_harmony_report.html"),
    )

    adata.write_h5ad(out)
    prov = adata.uns["omicsage_harmony"]
    print(f"[harmony] {prov['n_batches']} batches corrected  {prov['elapsed_seconds']:.1f}s → {out}")
    return out


def run_cluster_harmony(input_path: Path, out: Path, reports_dir: Path,
                        params: dict, cfg: dict, force: bool = False) -> Path:
    if out.exists() and not force:
        print(f"[cluster_harmony] cached → {out}")
        return out

    print("[cluster_harmony] running …")
    from pipeline.modules.clustering.cluster import cluster, compute_ari

    adata = sc.read_h5ad(input_path)

    adata, metrics_pre = cluster(
        adata,
        resolution_range=params.get("resolution_range", [0.2, 0.4, 0.6, 0.8, 1.0]),
        cluster_key="leiden",
        inplace=True,
    )
    adata, metrics_post = cluster(
        adata,
        resolution_range=params.get("resolution_range", [0.2, 0.4, 0.6, 0.8, 1.0]),
        neighbors_key="neighbors_harmony",
        cluster_key="leiden_harmony",
        inplace=True,
    )

    ari = compute_ari(adata, "leiden", "leiden_harmony")
    print(
        f"[cluster_harmony] pre={metrics_pre['best_n_clusters']} clusters"
        f"  post={metrics_post['best_n_clusters']} clusters"
        f"  ARI={ari:.4f} → {out}"
    )

    adata.write_h5ad(out)
    return out


def run_pseudobulk(input_path: Path, out: Path, reports_dir: Path,
                   params: dict, cfg: dict, force: bool = False) -> Path:
    print("[pseudobulk] running …")
    from pipeline.modules.downstream.pseudobulk_deg import pseudobulk_deg
    from pipeline.modules.downstream.pseudobulk_deg_report import generate_pseudobulk_deg_report

    adata = sc.read_h5ad(input_path)
    adata_pb, pb_dict = pseudobulk_deg(
        adata,
        groupby=params.get("groupby", "cell_type_vote"),
        donor_key=params.get("donor_key", "batch"),
        counts_layer=params.get("counts_layer", "counts"),
        min_cells=params.get("min_cells", 10),
        min_samples=params.get("min_samples", 3),
        min_logfc=params.get("min_logfc", 0.25),
        max_pval_adj=params.get("max_pval_adj", 0.05),
        exclude_gene_prefixes=params.get("exclude_gene_prefixes", []),
        inplace=False,
    )

    generate_pseudobulk_deg_report(
        adata=adata_pb,
        pb_dict=pb_dict,
        output_path=str(reports_dir / "10_pseudobulk_deg_report.html"),
        top_n_volcano=params.get("top_n_volcano", 10),
        top_n_dotplot=params.get("top_n_dotplot", 5),
        max_volcano_groups=params.get("max_volcano_groups", 11),
    )

    adata_pb.write_h5ad(out)
    prov = pb_dict["provenance"]
    total_sig = sum(len(df) for df in pb_dict["results"].values())
    print(
        f"[pseudobulk] {prov['n_groups']} groups tested"
        f"  {prov['n_skipped']} skipped"
        f"  {total_sig:,} DEGs → {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# AI LAYER
# ══════════════════════════════════════════════════════════════════════════════

def _run_ai_layer(
    cfg: dict,
    study_context: dict,
    processed_dir: Path,
    reports_dir: Path,
    log_dir: str,
    active_ai_modules: "set[str] | None" = None,
) -> None:
    """
    Run Phase 3 AI modules in sequence after the analysis pipeline.

    Called only when --ai is passed on the CLI.  Every module call is wrapped
    in try/except so a single failure never aborts the rest of the layer.
    Results are threaded forward so later modules benefit from earlier ones.

    Parameters
    ----------
    active_ai_modules
        If None, all modules run (default).
        If a set, only modules whose names are in the set run.
        Controlled at runtime via --ai-module / --ai-from-module / --ai-to-module.

    Module order
    ------------
    A2  clustering_advisor   — resolution recommendation (reads cluster adata)
    B1  cluster_annotator    — LLM cluster labels    (reads annotated adata)
    A1  pipeline_advisor     — overall QC advice     (reads final adata)
    B2  deg_validator        — DEG literature links  (reads deg/gsea adata)
    B3  coherence_reviewer   — analysis coherence + analysis_summary.json
    A3  downstream_suggester — NEXT_STEPS.md
    C1  narrative_generator  — ai_narrative.md + AI tab in combined HTML
    """

    def _should_run(module_name: str) -> bool:
        """Return True if this module should execute given the active_ai_modules filter."""
        if active_ai_modules is None:
            return True
        return module_name in active_ai_modules

    # Build a display string of which modules will run
    if active_ai_modules is None:
        modules_label = "all"
    else:
        modules_label = ", ".join(
            m for m in AI_MODULE_ORDER if m in active_ai_modules
        ) or "none"

    print()
    print("=" * 60)
    print(f"AI Layer — starting  (modules: {modules_label})")
    print("=" * 60)

    def _safe(label: str, fn):
        """
        Call fn(); return result or None on any error.
        Catches BaseException (including KeyboardInterrupt) so a single
        module failure or user interrupt never kills the rest of the layer.
        """
        print(f"[ai/{label}] starting …", flush=True)
        try:
            return fn()
        except BaseException as exc:
            print(f"[ai/{label}] ERROR: {type(exc).__name__}: {exc}", flush=True)
            return None

    # ------------------------------------------------------------------
    # Locate adata files
    # ------------------------------------------------------------------
    cluster_path = processed_dir / STEP_OUTPUT["cluster"]
    ann_path     = processed_dir / STEP_OUTPUT["annotate"]
    deg_path     = processed_dir / STEP_OUTPUT["deg"]
    gsea_path    = processed_dir / STEP_OUTPUT["gsea"]

    # Richest available adata for modules that just need analysis context
    final_adata = None
    for p in [gsea_path, deg_path, ann_path, cluster_path]:
        if p.exists():
            print(f"[ai] loading {p.name} as primary adata …", flush=True)
            final_adata = sc.read_h5ad(p)
            break

    if final_adata is None:
        print("[ai] no processed adata found — skipping AI layer")
        return

    # Adata whose rank_genes_groups is grouped by cell type
    # (used by deg_validator — needs cell-type-level DEGs)
    rgg_adata = None
    for p in [gsea_path, deg_path, ann_path]:
        if p.exists():
            tmp = sc.read_h5ad(p)
            if "rank_genes_groups" in tmp.uns:
                rgg_adata = tmp
                break

    # Adata whose rank_genes_groups is grouped by leiden cluster IDs
    # (used by cluster_annotator — needs per-cluster marker genes)
    # The annotate step stores leiden-based rank_genes_groups in 05_annotated.h5ad,
    # but the deg step overwrites it with cell_type_vote grouping.
    # We reload 05_annotated.h5ad and re-run rank_genes_groups if needed.
    def _leiden_rgg_adata():
        """Return adata with rank_genes_groups grouped by leiden IDs."""
        leiden_col = None
        # Determine cluster column name from annotated adata
        for col in ("leiden", "louvain"):
            if ann_path.exists():
                tmp = sc.read_h5ad(ann_path)
                if col in tmp.obs.columns:
                    leiden_col = col
                    ann_adata = tmp
                    break

        if leiden_col is None:
            return None

        rgg = ann_adata.uns.get("rank_genes_groups", {})
        groupby = (rgg.get("params") or {}).get("groupby", "")
        if groupby == leiden_col:
            # Already grouped by leiden — use as-is
            return ann_adata

        # Re-run rank_genes_groups with leiden groupby (in-memory only)
        print(f"[ai/cluster_annotator] re-running rank_genes_groups "
              f"with groupby='{leiden_col}' for marker extraction …", flush=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")   # suppress PerformanceWarning flood
            sc.tl.rank_genes_groups(
                ann_adata,
                groupby=leiden_col,
                method="wilcoxon",
                n_genes=20,
                key_added="rank_genes_groups",
            )
        return ann_adata

    # ------------------------------------------------------------------
    # A2 — Clustering advisor
    # ------------------------------------------------------------------
    cluster_advice = None
    if cluster_path.exists() and _should_run("clustering_advisor"):
        def _clustering():
            from ai import clustering_advisor
            adata_cl = sc.read_h5ad(cluster_path)
            sweep_raw = adata_cl.uns.get("omicsage_cluster", {}).get("silhouette_scores", {})
            sweep = {float(k): float(v) for k, v in sweep_raw.items()}
            result = clustering_advisor.run(
                adata_cl, cfg, study_context,
                resolution_sweep_results=sweep,
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print(f"[ai/clustering_advisor] done — "
                      f"suggested resolution={result.suggested_resolution}", flush=True)
            return result
        cluster_advice = _safe("clustering_advisor", _clustering)
    elif not _should_run("clustering_advisor"):
        print("[ai/clustering_advisor] skipped", flush=True)

    # ------------------------------------------------------------------
    # B1 — Cluster annotator
    # ------------------------------------------------------------------
    annotations = None
    if ann_path.exists() and _should_run("cluster_annotator"):
        def _annotator():
            from ai import cluster_annotator
            leiden_adata = _leiden_rgg_adata()
            if leiden_adata is None:
                print("[ai/cluster_annotator] no leiden column found — skipping", flush=True)
                return None
            result = cluster_annotator.run(
                leiden_adata, cfg, study_context,
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print(f"[ai/cluster_annotator] done — annotated {len(result)} clusters", flush=True)
                # Propagate ai_cell_type into final_adata so downstream modules
                # (coherence_reviewer, narrative_generator) can read it.
                # leiden_adata and final_adata share the same cell barcodes.
                if "ai_cell_type" in leiden_adata.obs.columns:
                    final_adata.obs["ai_cell_type"] = (
                        leiden_adata.obs["ai_cell_type"]
                        .reindex(final_adata.obs.index)
                    )
                    print("[ai/cluster_annotator] ai_cell_type propagated to final_adata", flush=True)
            return result
        annotations = _safe("cluster_annotator", _annotator)
    elif not _should_run("cluster_annotator"):
        print("[ai/cluster_annotator] skipped", flush=True)

    # ------------------------------------------------------------------
    # A1 — Pipeline advisor
    # ------------------------------------------------------------------
    pipeline_advice = None
    if _should_run("pipeline_advisor"):
        def _pipeline():
            from ai import pipeline_advisor
            result = pipeline_advisor.run(
                final_adata, cfg, study_context,
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print(f"[ai/pipeline_advisor] done — "
                      f"{len(result.recommendations)} recommendations", flush=True)
            return result
        pipeline_advice = _safe("pipeline_advisor", _pipeline)
    else:
        print("[ai/pipeline_advisor] skipped", flush=True)

    # ------------------------------------------------------------------
    # B2 — DEG validator
    # ------------------------------------------------------------------
    deg_validation = None
    if rgg_adata is not None and _should_run("deg_validator"):
        def _deg_val():
            from ai import deg_validator
            result = deg_validator.run(
                rgg_adata, cfg, study_context,
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print(f"[ai/deg_validator] done — validated {len(result)} groups", flush=True)
            return result
        deg_validation = _safe("deg_validator", _deg_val)
    elif not _should_run("deg_validator"):
        print("[ai/deg_validator] skipped", flush=True)

    # ------------------------------------------------------------------
    # B3 — Coherence reviewer + analysis_summary.json
    # ------------------------------------------------------------------
    coherence_review = None
    if _should_run("coherence_reviewer"):
        def _coherence():
            from ai import coherence_reviewer
            summary_path = str(reports_dir / "analysis_summary.json")
            result = coherence_reviewer.run(
                final_adata, cfg, study_context,
                summary_path=summary_path,
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print("[ai/coherence_reviewer] done — analysis_summary.json written", flush=True)
            return result
        coherence_review = _safe("coherence_reviewer", _coherence)
    else:
        print("[ai/coherence_reviewer] skipped", flush=True)

    # ------------------------------------------------------------------
    # A3 — Downstream suggester + NEXT_STEPS.md
    # ------------------------------------------------------------------
    if _should_run("downstream_suggester"):
        def _downstream():
            from ai import downstream_suggester
            result = downstream_suggester.run(
                final_adata, cfg, study_context,
                coherence_review=coherence_review,
                output_path=str(reports_dir / "NEXT_STEPS.md"),
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print("[ai/downstream_suggester] done — NEXT_STEPS.md written", flush=True)
            return result
        _safe("downstream_suggester", _downstream)
    else:
        print("[ai/downstream_suggester] skipped", flush=True)

    # ------------------------------------------------------------------
    # C1 — Narrative generator + ai_narrative.md
    # ------------------------------------------------------------------
    narrative_result = None
    if _should_run("narrative_generator"):
        def _narrative():
            from ai import narrative_generator
            result = narrative_generator.run(
                final_adata, cfg, study_context,
                pipeline_advice=pipeline_advice,
                cluster_annotations=annotations,
                deg_validation=deg_validation,
                coherence_review=coherence_review,
                output_path=str(reports_dir / "ai_narrative.md"),
                log_dir=log_dir, runtime_ai=True,
            )
            if result:
                print("[ai/narrative_generator] done — ai_narrative.md written", flush=True)
            return result
        narrative_result = _safe("narrative_generator", _narrative)
    else:
        print("[ai/narrative_generator] skipped", flush=True)

    print()
    print("=" * 60)
    print("AI Layer — complete")
    print("=" * 60)


# Map step name → runner function
STEP_RUNNERS = {
    "qc":              run_qc,
    "normalize":       run_normalize,
    "reduce":          run_reduce,
    "cluster":         run_cluster,
    "annotate":        run_annotate,
    "deg":             run_deg,
    "gsea":            run_gsea,
    "harmony":         run_harmony,
    "cluster_harmony": run_cluster_harmony,
    "pseudobulk":      run_pseudobulk,
}


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    parser = argparse.ArgumentParser(
        description="OmicSage — Generic Pipeline Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True,
                        help="Path to dataset YAML config (e.g. config/runs/GSE194122.yaml)")
    parser.add_argument("--from-step", metavar="STEP", default=None,
                        help="Start execution from this step (inclusive)")
    parser.add_argument("--to-step", metavar="STEP", default=None,
                        help="Stop execution at this step (inclusive)")
    parser.add_argument("--step", metavar="STEP", default=None,
                        help="Run exactly one step (shorthand for --from-step X --to-step X)")
    parser.add_argument("--ai", action="store_true", default=False,
                        help="Enable AI features (requires ai section in config and LLM provider)")
    parser.add_argument("--ai-only", action="store_true", default=False,
                        help="Skip all pipeline steps and run only the AI layer "
                             "(implies --ai; requires cached processed data)")
    parser.add_argument("--ai-module", metavar="MODULE", default=None,
                        help="Run exactly one AI module (e.g. coherence_reviewer). "
                             "Use 'all' to run all modules (default), 'off' to disable AI entirely. "
                             f"Valid: all, off, {', '.join(AI_MODULE_ORDER)}")
    parser.add_argument("--force", action="store_true", default=False,
                        help="Re-run steps even if their output file already exists (bypass cache)")
    parser.add_argument("--ai-from-module", metavar="MODULE", default=None,
                        help="Start the AI layer from this module (inclusive), skip earlier ones. "
                             f"Valid: {', '.join(AI_MODULE_ORDER)}")
    parser.add_argument("--ai-to-module", metavar="MODULE", default=None,
                        help="Stop the AI layer at this module (inclusive). "
                             f"Valid: {', '.join(AI_MODULE_ORDER)}")
    return parser.parse_args()


def resolve_ai_module_window(
    ai_module,       # str | None
    ai_from_module,  # str | None
    ai_to_module,    # str | None
):
    """
    Return the set of AI module names to run.

    Returns
    -------
    None
        No filter — all modules run (default behaviour).
    set[str]
        Only these module names will execute.
        An empty set means AI is disabled (--ai-module off).
    """
    # Special values
    if ai_module == "off":
        return set()
    if ai_module == "all":
        return None   # same as not specifying → run all

    # --ai-module and --ai-from/to-module are mutually exclusive
    if ai_module and (ai_from_module or ai_to_module):
        print("Error: --ai-module cannot be combined with --ai-from-module / --ai-to-module")
        sys.exit(1)

    # Single module
    if ai_module:
        if ai_module not in AI_MODULE_ORDER:
            print(f"Error: unknown AI module '{ai_module}' for --ai-module")
            print(f"Valid: all, off, {', '.join(AI_MODULE_ORDER)}")
            sys.exit(1)
        return {ai_module}

    # Range (--ai-from-module and/or --ai-to-module)
    from_idx = 0
    to_idx   = len(AI_MODULE_ORDER) - 1

    if ai_from_module:
        if ai_from_module not in AI_MODULE_ORDER:
            print(f"Error: unknown AI module '{ai_from_module}' for --ai-from-module")
            print(f"Valid: {', '.join(AI_MODULE_ORDER)}")
            sys.exit(1)
        from_idx = AI_MODULE_ORDER.index(ai_from_module)

    if ai_to_module:
        if ai_to_module not in AI_MODULE_ORDER:
            print(f"Error: unknown AI module '{ai_to_module}' for --ai-to-module")
            print(f"Valid: {', '.join(AI_MODULE_ORDER)}")
            sys.exit(1)
        to_idx = AI_MODULE_ORDER.index(ai_to_module)

    if from_idx > to_idx:
        print(
            f"Error: --ai-from-module '{ai_from_module}' comes after "
            f"--ai-to-module '{ai_to_module}' in the module order"
        )
        sys.exit(1)

    # If neither flag was given, no filtering needed
    if ai_from_module is None and ai_to_module is None:
        return None

    return set(AI_MODULE_ORDER[from_idx : to_idx + 1])


def resolve_step_window(from_step, to_step, step) -> tuple[int, int]:
    """Return (from_idx, to_idx) into STEP_ORDER."""
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

    # ── "all" shorthand — treat --step all as no step filter (run everything) ───
    if args.step == "all":
        args.step = None

    # ── validate step names ────────────────────────────────────────────────────
    for arg_name, val in [("--from-step", args.from_step),
                           ("--to-step",   args.to_step),
                           ("--step",      args.step)]:
        if val and val not in STEP_ORDER:
            print(f"Error: unknown step '{val}' for {arg_name}")
            print(f"Valid steps: all, {', '.join(STEP_ORDER)}")
            sys.exit(1)

    cfg = load_config(args.config)
    dataset_id   = cfg["dataset"]["id"]
    dataset_name = cfg["dataset"].get("name", dataset_id)

    processed_dir = Path(cfg["paths"]["processed_dir"])
    reports_dir   = Path(cfg["paths"]["reports_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    # ── AI setup ───────────────────────────────────────────────────────────────
    runtime_ai = args.ai or args.ai_only   # --ai-only implies --ai
    log_dir = "logs/llm"
    study_context: dict = {}

    # Resolve which AI modules to run (may set runtime_ai=False if 'off')
    active_ai_modules = resolve_ai_module_window(
        getattr(args, "ai_module", None),
        getattr(args, "ai_from_module", None),
        getattr(args, "ai_to_module", None),
    )
    # --ai-module off disables AI entirely even if --ai was passed
    if active_ai_modules is not None and len(active_ai_modules) == 0:
        runtime_ai = False

    if runtime_ai:
        from ai.pipeline_advisor import load_study_context
        study_context_path = Path("config/runs") / dataset_id / "study_context.yaml"
        study_context = load_study_context(study_context_path)
        Path(log_dir).mkdir(parents=True, exist_ok=True)
        provider = cfg.get("ai", {}).get("provider", "ollama")
        if active_ai_modules is None:
            modules_info = "all modules"
        else:
            modules_info = f"modules: {', '.join(m for m in AI_MODULE_ORDER if m in active_ai_modules)}"
        print(f"[ai] AI features ENABLED  (provider={provider}, {modules_info})")
    else:
        if getattr(args, "ai_module", None) == "off":
            print("[ai] AI features OFF  (disabled via --ai-module off)")
        else:
            print("[ai] AI features OFF  (pass --ai to enable)")

    # ── --ai-only: skip pipeline steps entirely ───────────────────────────────
    if args.ai_only:
        start_time = datetime.now()
        print("=" * 60)
        print(f"OmicSage — {dataset_name}")
        print("Mode     : AI layer only (pipeline steps skipped)")
        print(f"Started  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 60)
        _run_ai_layer(cfg, study_context, processed_dir, reports_dir, log_dir,
                      active_ai_modules=active_ai_modules)
        end_time = datetime.now()
        elapsed  = end_time - start_time
        hours, remainder = divmod(int(elapsed.total_seconds()), 3600)
        minutes, seconds = divmod(remainder, 60)
        print()
        print("=" * 60)
        print("AI layer complete.")
        print(f"Reports   : {reports_dir}")
        print(f"Elapsed   : {hours:02d}h {minutes:02d}m {seconds:02d}s")
        print("=" * 60)
        return

    # ── build active step list ─────────────────────────────────────────────────
    from_idx, to_idx = resolve_step_window(args.from_step, args.to_step, args.step)
    window = STEP_ORDER[from_idx : to_idx + 1]
    active_steps = [s for s in window if get_step_cfg(cfg, s)["enabled"]]

    if not active_steps:
        print(f"No enabled steps in the requested range: {window}")
        sys.exit(0)

    # ── validate inputs before touching anything ───────────────────────────────
    validate_plan(active_steps, cfg, processed_dir)

    # ── banner ─────────────────────────────────────────────────────────────────
    start_time = datetime.now()
    print("=" * 60)
    print(f"OmicSage — {dataset_name}")
    print(f"Steps    : {' → '.join(active_steps)}")
    print(f"Started  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    # ── run ────────────────────────────────────────────────────────────────────
    force     = args.force
    deg_dict  = None   # carried from deg → gsea

    sc.settings.verbosity = 1

    for step in active_steps:
        step_cfg  = get_step_cfg(cfg, step)
        params    = step_cfg["params"]
        input_path = resolve_input(step, cfg, processed_dir)
        out        = processed_dir / STEP_OUTPUT[step]

        if step == "deg":
            out, deg_dict = run_deg(input_path, out, reports_dir, params, cfg, force=force)

        elif step == "gsea":
            # deg_dict may not have been computed this run if deg was cached/skipped
            if deg_dict is None:
                deg_params = get_step_cfg(cfg, "deg")["params"]
                _, deg_dict = _reload_deg_dict(processed_dir, deg_params)
            run_gsea(input_path, out, reports_dir, params, cfg, deg_dict, force=force)

        elif step == "reduce":
            batch_key = params.get("batch_key")
            STEP_RUNNERS[step](input_path, out, reports_dir, params, cfg,
                               force=force, batch_key=batch_key)

        else:
            STEP_RUNNERS[step](input_path, out, reports_dir, params, cfg, force=force)

    # ── combined report ────────────────────────────────────────────────────
    from reports.combined_report import generate_combined_report
    generate_combined_report(
        reports_dir=reports_dir,
        dataset_name=dataset_name,
        output_path=reports_dir / "00_combined_report.html",
    )

    # ── AI layer (runs after combined report so reviewer can inspect it) ───────
    if runtime_ai:
        _run_ai_layer(cfg, study_context, processed_dir, reports_dir, log_dir,
                      active_ai_modules=active_ai_modules)

    # ── footer ─────────────────────────────────────────────────────────────────
    end_time = datetime.now()
    elapsed  = end_time - start_time
    hours, remainder = divmod(int(elapsed.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    print()
    print("=" * 60)
    print("Pipeline complete.")
    print(f"Reports   : {reports_dir}")
    print(f"Elapsed   : {hours:02d}h {minutes:02d}m {seconds:02d}s")
    print("=" * 60)


if __name__ == "__main__":
    main()