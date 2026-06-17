"""
OmicSage — Spatial Transcriptomics Pipeline Runner
====================================================
Usage
-----
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --step all
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --to-step cluster
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step reduce
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step qc --to-step cluster
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --step cluster
  python run_spatial_pipeline.py --config config/runs/kuppe_heart.yaml --from-step reduce --force

Step order: ingest -> qc -> reduce -> cluster -> deconvolve -> downstream -> impute

Checkpointing: every step writes output_dir/NN_<step>.h5ad.
If the file exists the step is skipped. Use --force to override.
"""

from __future__ import annotations

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

import yaml

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
logging.getLogger("streamlit").setLevel(logging.ERROR)

root = Path(__file__).resolve().parent
while not (root / "pipeline").exists() and root != root.parent:
    root = root.parent
os.chdir(root)
sys.path.insert(0, str(root))


STEP_ORDER = ["ingest", "qc", "reduce", "cluster", "deconvolve", "downstream", "impute"]

STEP_OUTPUT = {
    "ingest":     "01_ingested.h5ad",
    "qc":         "02_qc.h5ad",
    "reduce":     "03_reduced.h5ad",
    "cluster":    "04_clustered.h5ad",
    "deconvolve": "05_deconvolved.h5ad",
    "downstream": "06_downstream.h5ad",
    "impute":     "07_imputed.h5ad",
}

STEP_PREDECESSOR = {
    "ingest":     None,
    "qc":         "ingest",
    "reduce":     "qc",
    "cluster":    "reduce",
    "deconvolve": "cluster",
    # downstream prefers the deconvolve checkpoint (has cell-type abundance
    # columns); falls back to cluster checkpoint when deconvolve was skipped.
    # resolve_input() implements the fallback logic.
    "downstream": "deconvolve",
    # impute reads from cluster checkpoint (does not require deconvolution).
    "impute":     "cluster",
}

STEP_REPORT = {
    "qc":         "spatial_qc_report.html",
    "reduce":     "spatial_reduce_report.html",
    "cluster":    "spatial_cluster_report.html",
    "deconvolve": "spatial_deconvolve_report.html",
    "downstream": "spatial_downstream_report.html",
    "impute":     "spatial_impute_report.html",
}


def load_config(config_path):
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def get_step_cfg(cfg, step):
    step_cfg = cfg.get("spatial", {}).get(step, {})
    if isinstance(step_cfg, bool):
        return {"enabled": step_cfg, "params": {}}
    return {
        "enabled": step_cfg.get("enabled", True) if isinstance(step_cfg, dict) else True,
        "params":  step_cfg if isinstance(step_cfg, dict) else {},
    }


# Steps that have a preferred predecessor but can fall back to an earlier one.
STEP_PREDECESSOR_FALLBACK = {
    "downstream": "cluster",  # deconvolve preferred; cluster is the fallback
}


def resolve_input(step, cfg, output_dir):
    pred = STEP_PREDECESSOR[step]
    if pred is None:
        return None
    pred_out = output_dir / STEP_OUTPUT[pred]
    if pred_out.exists():
        return pred_out

    # Try fallback predecessor (e.g. downstream: deconvolve -> cluster)
    fallback = STEP_PREDECESSOR_FALLBACK.get(step)
    if fallback is not None:
        fallback_out = output_dir / STEP_OUTPUT[fallback]
        if fallback_out.exists():
            print(
                f"  [{step}] deconvolve checkpoint not found — "
                f"falling back to '{fallback}' checkpoint ({fallback_out.name})"
            )
            return fallback_out

    raise FileNotFoundError(
        f"[{step}] requires output of '{pred}' at {pred_out}\n"
        f"  Run '{pred}' first or use --from-step {pred}"
    )


def validate_plan(active_steps, cfg, output_dir):
    for step in active_steps:
        if step == "ingest":
            source = cfg.get("spatial", {}).get("source")
            if not source:
                raise ValueError("[ingest] spatial.source not set in config.")
        else:
            pred = STEP_PREDECESSOR[step]
            # Only require the predecessor checkpoint when the predecessor
            # is NOT itself in the active window (won't be run this session)
            if pred not in active_steps:
                resolve_input(step, cfg, output_dir)


def run_ingest(cfg, output_dir, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_ingest import spatial_ingest
    out_path = output_dir / STEP_OUTPUT["ingest"]
    if out_path.exists() and not force:
        print(f"  [ingest] cached -> {out_path}")
        return out_path
    spatial_cfg = cfg.get("spatial", {})
    source = spatial_cfg.get("source", "benchmark")
    print(f"  [ingest] loading from: {source!r}")
    adata, _ = spatial_ingest(
        source=source,
        spatial_type=spatial_cfg.get("spatial_type", "auto"),
        counts_file=spatial_cfg.get("counts_file", "filtered_feature_bc_matrix.h5"),
        library_id=spatial_cfg.get("library_id", None),
        library_key=spatial_cfg.get("library_key", None),
        load_images=spatial_cfg.get("load_images", True),
    )
    print(f"  [ingest] {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    adata.write_h5ad(out_path)
    print(f"  [ingest] -> {out_path}")
    return out_path


def run_qc(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_qc import spatial_qc
    from reports.templates.spatial.spatial_qc_report import generate_spatial_qc_report
    out_path    = output_dir / STEP_OUTPUT["qc"]
    report_path = reports_dir / STEP_REPORT["qc"]
    dataset_id  = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [qc] cached -> {out_path}")
        return out_path
    adata  = sc.read_h5ad(input_path)
    qc_cfg = cfg.get("spatial", {}).get("qc", {})
    adata, params = spatial_qc(
        adata,
        min_counts=qc_cfg.get("min_counts",   500),
        max_counts=qc_cfg.get("max_counts",   100_000),
        min_genes=qc_cfg.get("min_genes",     200),
        max_genes=qc_cfg.get("max_genes",     10_000),
        max_mt_pct=qc_cfg.get("max_mt_pct",  20.0),
        mt_prefix=qc_cfg.get("mt_prefix",    "MT-"),
        filter_spots=qc_cfg.get("filter_spots", True),
        inplace=True,
    )
    out = params["outputs"]
    print(f"  [qc] {out['n_spots_before']:,} -> {out['n_spots_after']:,} spots "
          f"({out['n_spots_removed']:,} removed)")
    generate_spatial_qc_report(adata, str(report_path), dataset_id=dataset_id)
    print(f"  [qc] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [qc] -> {out_path}")
    return out_path


def run_reduce(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_reduce import spatial_reduce
    from reports.templates.spatial.spatial_reduce_report import generate_spatial_reduce_report
    out_path    = output_dir / STEP_OUTPUT["reduce"]
    report_path = reports_dir / STEP_REPORT["reduce"]
    dataset_id  = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [reduce] cached -> {out_path}")
        return out_path
    adata      = sc.read_h5ad(input_path)
    reduce_cfg = cfg.get("spatial", {}).get("reduce", {})
    adata, params = spatial_reduce(
        adata,
        n_top_genes=reduce_cfg.get("n_top_genes",  3000),
        n_comps=reduce_cfg.get("n_comps",           50),
        n_neighbors=reduce_cfg.get("n_neighbors",   6),
        coord_type=reduce_cfg.get("coord_type",     None),
        normalize_total=reduce_cfg.get("normalize_total", True),
        target_sum=reduce_cfg.get("target_sum",     1e4),
        log1p=reduce_cfg.get("log1p",               True),
        flavor=reduce_cfg.get("flavor",             "seurat"),
        inplace=True,
    )
    out = params["outputs"]
    print(f"  [reduce] {out['n_hvgs']:,} HVGs, {out['n_comps_computed']} PCs, "
          f"{out['spatial_graph_n_edges']:,} spatial edges")
    generate_spatial_reduce_report(adata, str(report_path), dataset_id=dataset_id)
    print(f"  [reduce] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [reduce] -> {out_path}")
    return out_path


def run_cluster(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_cluster import spatial_cluster
    from reports.templates.spatial.spatial_cluster_report import generate_spatial_cluster_report
    out_path     = output_dir / STEP_OUTPUT["cluster"]
    report_path  = reports_dir / STEP_REPORT["cluster"]
    dataset_id   = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [cluster] cached -> {out_path}")
        return out_path
    adata        = sc.read_h5ad(input_path)
    cluster_cfg  = cfg.get("spatial", {}).get("cluster", {})
    adata, params = spatial_cluster(
        adata,
        resolution=cluster_cfg.get("resolution",   0.5),
        n_neighbors=cluster_cfg.get("n_neighbors", 15),
        n_pcs=cluster_cfg.get("n_pcs",             30),
        random_state=cluster_cfg.get("random_state", 0),
        run_svg=cluster_cfg.get("run_svg",          True),
        svg_n_genes=cluster_cfg.get("svg_n_genes",  None),
        annotation_map=cluster_cfg.get("annotation_map", None),
        inplace=True,
    )
    out = params["outputs"]
    print(f"  [cluster] {out['n_clusters']} clusters found")
    if "n_significant_fdr05" in out:
        print(f"  [cluster] {out['n_significant_fdr05']} SVGs (FDR<0.05) "
              f"from {out['n_genes_tested']} tested")
    generate_spatial_cluster_report(adata, str(report_path), dataset_id=dataset_id)
    print(f"  [cluster] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [cluster] -> {out_path}")
    return out_path


def run_deconvolve(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_deconvolve import spatial_deconvolve
    from reports.templates.spatial.spatial_deconvolve_report import generate_spatial_deconvolve_report
    out_path    = output_dir / STEP_OUTPUT["deconvolve"]
    report_path = reports_dir / STEP_REPORT["deconvolve"]
    dataset_id  = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [deconvolve] cached -> {out_path}")
        return out_path
    adata      = sc.read_h5ad(input_path)
    deconv_cfg = cfg.get("spatial", {}).get("deconvolve", {})
    method     = deconv_cfg.get("method", "nnls")
    ref_path   = deconv_cfg.get("ref_path", None)
    ref_adata  = None
    if ref_path:
        print(f"  [deconvolve] loading reference: {ref_path!r}")
        ref_adata = sc.read_h5ad(ref_path)
        # NOTE: do NOT pre-assign ref_adata.X here — spatial_deconvolve
        # handles this internally for both methods.
        print(f"  [deconvolve] reference: {ref_adata.n_obs:,} cells x {ref_adata.n_vars:,} genes")
    else:
        print("  [deconvolve] no ref_path — will be skipped")
    print(f"  [deconvolve] method: {method}")
    adata, params = spatial_deconvolve(
        adata,
        ref_adata=ref_adata,
        # method selection
        method=method,
        per_sample=deconv_cfg.get("per_sample",             False),
        library_key=deconv_cfg.get("library_key",           None),
        # shared
        cell_type_key=deconv_cfg.get("cell_type_key",       "cell_type_original"),
        layer_ref=deconv_cfg.get("layer_ref",               "counts"),
        # NNLS-specific
        n_jobs=deconv_cfg.get("n_jobs",                     4),
        target_sum=deconv_cfg.get("target_sum",             1e4),
        # cell2location-specific (ignored when method=nnls)
        batch_key_ref=deconv_cfg.get("batch_key_ref",       "donor_id"),
        batch_key_st=deconv_cfg.get("batch_key_st",         "patient"),
        covariate_keys=deconv_cfg.get("covariate_keys",     None),
        N_cells_per_location=deconv_cfg.get("N_cells_per_location", 8),
        detection_alpha=deconv_cfg.get("detection_alpha",   20),
        max_epochs_ref=deconv_cfg.get("max_epochs_ref",     250),
        max_epochs_st=deconv_cfg.get("max_epochs_st",       30000),
        batch_size_ref=deconv_cfg.get("batch_size_ref",     2500),
        batch_size_st=deconv_cfg.get("batch_size_st",       None),
        num_samples_posterior=deconv_cfg.get("num_samples_posterior", 1000),
        cell_count_cutoff=deconv_cfg.get("cell_count_cutoff",         5),
        cell_percentage_cutoff2=deconv_cfg.get("cell_percentage_cutoff2", 0.03),
        nonz_mean_cutoff=deconv_cfg.get("nonz_mean_cutoff", 1.12),
        inplace=True,
    )
    if params["skipped"]:
        print(f"  [deconvolve] skipped: {params['skip_reason']}")
    else:
        out = params["outputs"]
        print(f"  [deconvolve] {out['n_cell_types']} cell types, "              f"{out['n_spots']:,} spots ({out['method']})")
        if out.get("per_sample"):
            print(f"  [deconvolve] per-sample loop over '{out['library_key']}'")
    generate_spatial_deconvolve_report(
        adata, str(report_path), dataset_id=dataset_id,
        img_key=cfg.get("spatial", {}).get("report", {}).get("img_key", "hires"),
    )
    print(f"  [deconvolve] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [deconvolve] -> {out_path}")
    return out_path


def run_downstream(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_downstream import spatial_downstream
    from reports.templates.spatial.spatial_downstream_report import generate_spatial_downstream_report
    out_path    = output_dir / STEP_OUTPUT["downstream"]
    report_path = reports_dir / STEP_REPORT["downstream"]
    dataset_id  = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [downstream] cached -> {out_path}")
        return out_path
    print(f"  [downstream] loading: {input_path.name}", flush=True)
    adata          = sc.read_h5ad(input_path)
    print(
        f"  [downstream] loaded: {adata.n_obs:,} spots x {adata.n_vars:,} genes  "
        f"| deconvolved: {'omicsage_spatial_deconvolve' in adata.uns}  "
        f"| cell types: {len(list(adata.uns.get('omicsage_spatial_deconvolve', {}).get('outputs', {}).get('cell_type_names', [])))}",
        flush=True,
    )
    ds_cfg         = cfg.get("spatial", {}).get("downstream", {})
    dominant_key   = ds_cfg.get("dominant_celltype_key", "dominant_cell_type")
    adata, params = spatial_downstream(
        adata,
        run_region_clustering=ds_cfg.get("run_region_clustering", True),
        region_resolution=ds_cfg.get("region_resolution", 0.5),
        region_n_neighbors=ds_cfg.get("region_n_neighbors", 15),
        run_celltype_expression=ds_cfg.get("run_celltype_expression", True),
        n_marker_genes=ds_cfg.get("n_marker_genes", 20),
        run_celltype_svg=ds_cfg.get("run_celltype_svg", True),
        svg_n_genes=ds_cfg.get("svg_n_genes", None),
        run_co_occurrence=ds_cfg.get("run_co_occurrence", True),
        run_nhood_enrichment=ds_cfg.get("run_nhood_enrichment", True),
        n_perms_nhood=ds_cfg.get("n_perms_nhood", 1000),
        run_ligrec=ds_cfg.get("run_ligrec", True),
        ligrec_n_perms=ds_cfg.get("ligrec_n_perms", 1000),
        ligrec_organism=ds_cfg.get("ligrec_organism", "human"),
        run_svg_gsea=ds_cfg.get("run_svg_gsea", True),
        svg_gsea_gene_sets=ds_cfg.get("svg_gsea_gene_sets", "GO_Biological_Process_2023"),
        svg_gsea_organism=ds_cfg.get("svg_gsea_organism", "Human"),
        dominant_celltype_key=dominant_key,
        n_jobs=ds_cfg.get("n_jobs", 1),
        inplace=True,
    )
    analyses = params.get("analyses", {})
    ran = [k for k, v in analyses.items() if not v.get("skipped")]
    skipped = [k for k, v in analyses.items() if v.get("skipped")]
    print(f"  [downstream] ran: {', '.join(ran) or 'none'}")
    if skipped:
        print(f"  [downstream] skipped: {', '.join(skipped)}")
    generate_spatial_downstream_report(
        adata, str(report_path), dataset_id=dataset_id, dominant_celltype_key=dominant_key
    )
    print(f"  [downstream] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [downstream] -> {out_path}")
    return out_path


def run_impute(input_path, output_dir, reports_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.scripts.spatial.spatial_impute import spatial_impute
    from reports.templates.spatial.spatial_impute_report import generate_spatial_impute_report
    out_path    = output_dir / STEP_OUTPUT["impute"]
    report_path = reports_dir / STEP_REPORT["impute"]
    dataset_id  = cfg.get("dataset_id", "spatial")

    if out_path.exists() and not force:
        print(f"  [impute] cached -> {out_path}")
        return out_path

    impute_cfg   = cfg.get("spatial", {}).get("impute", {})
    sc_ref_path  = impute_cfg.get("sc_reference_path", None)
    enabled      = impute_cfg.get("enabled", False)

    # Skip cleanly when disabled or no reference path configured
    if not enabled or not sc_ref_path:
        reason = "disabled in config" if not enabled else "sc_reference_path not set"
        print(f"  [impute] skipped: {reason}")
        # Passthrough: copy cluster checkpoint as impute checkpoint so
        # downstream steps can still resolve their predecessor.
        adata = sc.read_h5ad(input_path)
        adata.uns["omicsage_spatial_impute"] = {
            "module":      "spatial_impute",
            "timestamp":   "",
            "method":      impute_cfg.get("method", "tangram"),
            "skipped":     True,
            "skip_reason": reason,
            "outputs":     {},
        }
        generate_spatial_impute_report(
            adata, str(report_path), dataset_id=dataset_id, sc_ref_label=""
        )
        adata.write_h5ad(out_path)
        return out_path

    print(f"  [impute] loading spatial: {input_path.name}")
    adata = sc.read_h5ad(input_path)

    print(f"  [impute] loading sc reference: {sc_ref_path!r}")
    sc_adata = sc.read_h5ad(sc_ref_path)
    print(f"  [impute] sc reference: {sc_adata.n_obs:,} cells x {sc_adata.n_vars:,} genes")

    method            = impute_cfg.get("method", "tangram")
    n_top_genes       = impute_cfg.get("n_top_genes", 2000)
    cell_type_key     = impute_cfg.get("cell_type_key", "cell_type")
    device            = impute_cfg.get("device", "cpu")
    tangram_mode      = impute_cfg.get("tangram_mode", "clusters")
    max_cells_per_type = impute_cfg.get("max_cells_per_type", 500)

    print(f"  [impute] method: {method}, mode: {tangram_mode}, "
          f"n_top_genes: {n_top_genes}, device: {device}")
    adata, params = spatial_impute(
        adata,
        adata_sc=sc_adata,
        method=method,
        cell_type_key=cell_type_key,
        n_top_genes=n_top_genes,
        device=device,
        tangram_mode=tangram_mode,
        max_cells_per_type=max_cells_per_type,
        inplace=True,
    )

    out = params.get("outputs", {})
    if params.get("skipped"):
        print(f"  [impute] skipped: {params.get('skip_reason')}")
    else:
        print(
            f"  [impute] {out.get('n_genes_imputed', 0):,} genes imputed, "
            f"{out.get('n_spots', 0):,} spots, "
            f"mean mapping score: {out.get('mean_mapping_score', float('nan')):.3f}"
        )

    sc_ref_label = Path(sc_ref_path).name
    generate_spatial_impute_report(
        adata, str(report_path), dataset_id=dataset_id, sc_ref_label=sc_ref_label
    )
    print(f"  [impute] report -> {report_path}")

    # obsm["imputed_expression"] is already a float32 numpy array —
    # no extra serialization needed before h5ad checkpoint.
    adata.write_h5ad(out_path)
    print(f"  [impute] -> {out_path}")
    return out_path


STEP_RUNNERS = {
    "ingest":     lambda inp, out, cfg, force: run_ingest(cfg, out, force),
    "qc":         run_qc,
    "reduce":     run_reduce,
    "cluster":    run_cluster,
    "deconvolve": run_deconvolve,
    "downstream": run_downstream,
    "impute":     run_impute,
}


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


def parse_args():
    parser = argparse.ArgumentParser(
        description="OmicSage — Spatial Transcriptomics Pipeline"
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--step", default=None,
                        help=f"Run one step. One of: all, {', '.join(STEP_ORDER)}")
    parser.add_argument("--from-step", dest="from_step", default=None,
                        help="Start from this step (inclusive)")
    parser.add_argument("--to-step", dest="to_step", default=None,
                        help="Stop at this step (inclusive)")
    parser.add_argument("--force", action="store_true",
                        help="Re-run steps even if checkpoint exists")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.step == "all":
        args.step = None

    for arg_name, val in [("--from-step", args.from_step),
                           ("--to-step",   args.to_step),
                           ("--step",      args.step)]:
        if val and val not in STEP_ORDER:
            print(f"Error: unknown step '{val}' for {arg_name}")
            print(f"Valid steps: all, {', '.join(STEP_ORDER)}")
            sys.exit(1)

    cfg          = load_config(args.config)
    dataset_id   = cfg.get("dataset_id", "spatial")
    dataset_name = cfg.get("dataset_name", dataset_id)

    output_dir = Path(cfg["paths"]["output_dir"])
    reports_dir   = Path(cfg["paths"]["reports_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)
    print(output_dir)
    print(reports_dir)

    from_idx, to_idx = resolve_step_window(args.from_step, args.to_step, args.step)
    window = STEP_ORDER[from_idx : to_idx + 1]
    active_steps = [s for s in window
                    if get_step_cfg(cfg, s).get("enabled", True)]

    if not active_steps:
        print(f"No enabled steps in the requested range: {window}")
        sys.exit(0)

    validate_plan(active_steps, cfg, output_dir)

    start_time = datetime.now()
    print("=" * 60)
    print(f"OmicSage Spatial -- {dataset_name}")
    print(f"Steps    : {' -> '.join(active_steps)}")
    print(f"Started  : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    for step in active_steps:
        print(f"\n[{step}]")
        if step == "ingest":
            STEP_RUNNERS[step](None, output_dir, cfg, args.force)
        else:
            input_path = resolve_input(step, cfg, output_dir)
            STEP_RUNNERS[step](input_path, output_dir, reports_dir, cfg, args.force)

    from reports.templates.spatial.spatial_combined_report import generate_spatial_combined_report
    combined_path = reports_dir / "00_spatial_combined_report.html"
    generate_spatial_combined_report(
        reports_dir=reports_dir,
        dataset_name=dataset_name,
        output_path=combined_path,
    )

    end_time = datetime.now()
    elapsed  = end_time - start_time
    hours, remainder = divmod(int(elapsed.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)
    print()
    print("=" * 60)
    print("Spatial pipeline complete.")
    print(f"Output dir      : {output_dir}")
    print(f"Combined report : {combined_path}")
    print(f"Elapsed         : {hours:02d}h {minutes:02d}m {seconds:02d}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
