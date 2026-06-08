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

Step order: ingest -> qc -> reduce -> cluster -> deconvolve

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


STEP_ORDER = ["ingest", "qc", "reduce", "cluster", "deconvolve"]

STEP_OUTPUT = {
    "ingest":     "01_ingested.h5ad",
    "qc":         "02_qc.h5ad",
    "reduce":     "03_reduced.h5ad",
    "cluster":    "04_clustered.h5ad",
    "deconvolve": "05_deconvolved.h5ad",
}

STEP_PREDECESSOR = {
    "ingest":     None,
    "qc":         "ingest",
    "reduce":     "qc",
    "cluster":    "reduce",
    "deconvolve": "cluster",
}

STEP_REPORT = {
    "qc":         "spatial_qc_report.html",
    "reduce":     "spatial_reduce_report.html",
    "cluster":    "spatial_cluster_report.html",
    "deconvolve": "spatial_deconvolve_report.html",
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


def resolve_input(step, cfg, output_dir):
    pred = STEP_PREDECESSOR[step]
    if pred is None:
        return None
    pred_out = output_dir / STEP_OUTPUT[pred]
    if pred_out.exists():
        return pred_out
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
    from pipeline.modules.spatial.spatial_ingest import spatial_ingest
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
        load_images=spatial_cfg.get("load_images", True),
    )
    print(f"  [ingest] {adata.n_obs:,} spots x {adata.n_vars:,} genes")
    adata.write_h5ad(out_path)
    print(f"  [ingest] -> {out_path}")
    return out_path


def run_qc(input_path, output_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.spatial.spatial_qc import spatial_qc
    from reports.templates.spatial.spatial_qc_report import generate_spatial_qc_report
    out_path    = output_dir / STEP_OUTPUT["qc"]
    report_path = output_dir / STEP_REPORT["qc"]
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


def run_reduce(input_path, output_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.spatial.spatial_reduce import spatial_reduce
    from reports.templates.spatial.spatial_reduce_report import generate_spatial_reduce_report
    out_path    = output_dir / STEP_OUTPUT["reduce"]
    report_path = output_dir / STEP_REPORT["reduce"]
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


def run_cluster(input_path, output_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.spatial.spatial_cluster import spatial_cluster
    from reports.templates.spatial.spatial_cluster_report import generate_spatial_cluster_report
    out_path     = output_dir / STEP_OUTPUT["cluster"]
    report_path  = output_dir / STEP_REPORT["cluster"]
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


def run_deconvolve(input_path, output_dir, cfg, force=False):
    import scanpy as sc
    from pipeline.modules.spatial.spatial_deconvolve import spatial_deconvolve
    from reports.templates.spatial.spatial_deconvolve_report import generate_spatial_deconvolve_report
    out_path    = output_dir / STEP_OUTPUT["deconvolve"]
    report_path = output_dir / STEP_REPORT["deconvolve"]
    dataset_id  = cfg.get("dataset_id", "spatial")
    if out_path.exists() and not force:
        print(f"  [deconvolve] cached -> {out_path}")
        return out_path
    adata       = sc.read_h5ad(input_path)
    deconv_cfg  = cfg.get("spatial", {}).get("deconvolve", {})
    ref_path    = deconv_cfg.get("ref_path", None)
    ref_adata   = None
    if ref_path:
        print(f"  [deconvolve] loading reference: {ref_path!r}")
        ref_adata = sc.read_h5ad(ref_path)
        layer_ref = deconv_cfg.get("layer_ref", "counts")
        ref_adata.X = ref_adata.layers[layer_ref].copy()
        print(f"  [deconvolve] reference: {ref_adata.n_obs:,} cells x {ref_adata.n_vars:,} genes")
    else:
        print("  [deconvolve] no ref_path — will be skipped")
    adata, params = spatial_deconvolve(
        adata,
        ref_adata=ref_adata,
        cell_type_key=deconv_cfg.get("cell_type_key",        "cell_type_original"),
        batch_key_ref=deconv_cfg.get("batch_key_ref",        "donor_id"),
        batch_key_st=deconv_cfg.get("batch_key_st",          "patient"),
        covariate_keys=deconv_cfg.get("covariate_keys",      None),
        layer_ref=deconv_cfg.get("layer_ref",                "counts"),
        N_cells_per_location=deconv_cfg.get("N_cells_per_location", 8),
        detection_alpha=deconv_cfg.get("detection_alpha",    20),
        max_epochs_ref=deconv_cfg.get("max_epochs_ref",      250),
        max_epochs_st=deconv_cfg.get("max_epochs_st",        30000),
        batch_size_ref=deconv_cfg.get("batch_size_ref",      2500),
        inplace=True,
    )
    if params["skipped"]:
        print(f"  [deconvolve] skipped: {params['skip_reason']}")
    else:
        out = params["outputs"]
        print(f"  [deconvolve] {out['n_cell_types']} cell types, {out['n_spots']:,} spots")
    generate_spatial_deconvolve_report(adata, str(report_path), dataset_id=dataset_id)
    print(f"  [deconvolve] report -> {report_path}")
    adata.write_h5ad(out_path)
    print(f"  [deconvolve] -> {out_path}")
    return out_path


STEP_RUNNERS = {
    "ingest":     lambda inp, out, cfg, force: run_ingest(cfg, out, force),
    "qc":         run_qc,
    "reduce":     run_reduce,
    "cluster":    run_cluster,
    "deconvolve": run_deconvolve,
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
    output_dir   = Path(cfg.get("output_dir", f"outputs/{dataset_id}"))
    output_dir.mkdir(parents=True, exist_ok=True)

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
            STEP_RUNNERS[step](input_path, output_dir, cfg, args.force)

    from reports.templates.spatial.spatial_combined_report import generate_spatial_combined_report
    combined_path = output_dir / "00_spatial_combined_report.html"
    generate_spatial_combined_report(
        reports_dir=output_dir,
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
