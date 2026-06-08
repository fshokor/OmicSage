"""
run_spatial_pipeline.py — OmicSage Phase 7 spatial pipeline runner

Usage:
    python run_spatial_pipeline.py --config config/runs/visium_lymphnode.yaml
    python run_spatial_pipeline.py --source benchmark --dataset_id mouse_brain

Steps added per session:
    Session 1 (done): ingest + QC
    Session 2 (done): reduce (HVG, PCA, spatial neighbours)
    Session 3:        cluster + spatially variable genes
    Session 4:        deconvolution + combined report
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


def load_config(config_path: str) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def run_spatial_pipeline(config: dict) -> None:
    from pipeline.modules.spatial.spatial_ingest import spatial_ingest
    from pipeline.modules.spatial.spatial_qc import spatial_qc
    from pipeline.modules.spatial.spatial_reduce import spatial_reduce
    from reports.templates.spatial.spatial_qc_report import (
        generate_spatial_qc_report,
    )
    from reports.templates.spatial.spatial_reduce_report import (
        generate_spatial_reduce_report,
    )

    spatial_cfg = config.get("spatial", {})
    dataset_id  = config.get("dataset_id", "spatial")
    source      = spatial_cfg.get("source", "benchmark")
    output_dir  = config.get("output_dir", f"outputs/{dataset_id}")

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Step 1 — Ingest
    # ------------------------------------------------------------------ #
    print(f"[1/4] Ingesting spatial data from: {source!r}")
    adata, ingest_params = spatial_ingest(
        source=source,
        spatial_type=spatial_cfg.get("spatial_type", "auto"),
        counts_file=spatial_cfg.get("counts_file", "filtered_feature_bc_matrix.h5"),
        library_id=spatial_cfg.get("library_id", None),
        load_images=spatial_cfg.get("load_images", True),
    )
    print(f"      {adata.n_obs:,} spots × {adata.n_vars:,} genes loaded")

    # ------------------------------------------------------------------ #
    # Step 2 — QC
    # ------------------------------------------------------------------ #
    qc_cfg = spatial_cfg.get("qc", {})
    print("[2/4] Running spatial QC ...")
    adata, qc_params = spatial_qc(
        adata,
        min_counts=qc_cfg.get("min_counts", 500),
        max_counts=qc_cfg.get("max_counts", 100_000),
        min_genes=qc_cfg.get("min_genes", 200),
        max_genes=qc_cfg.get("max_genes", 10_000),
        max_mt_pct=qc_cfg.get("max_mt_pct", 20.0),
        mt_prefix=qc_cfg.get("mt_prefix", "MT-"),
        filter_spots=qc_cfg.get("filter_spots", True),
        inplace=True,
    )
    outputs = qc_params["outputs"]
    print(
        f"      {outputs['n_spots_before']:,} → {outputs['n_spots_after']:,} spots "
        f"({outputs['n_spots_removed']:,} removed)"
    )

    qc_report_path = str(Path(output_dir) / f"{dataset_id}_spatial_qc_report.html")
    print(f"      Writing QC report → {qc_report_path}")
    generate_spatial_qc_report(adata, qc_report_path, dataset_id=dataset_id)

    # ------------------------------------------------------------------ #
    # Step 3 — Reduce (normalize, HVG, PCA, spatial neighbours)
    # ------------------------------------------------------------------ #
    reduce_cfg = spatial_cfg.get("reduce", {})
    print("[3/4] Running spatial reduction (HVG, PCA, spatial neighbours) ...")
    adata, reduce_params = spatial_reduce(
        adata,
        n_top_genes=reduce_cfg.get("n_top_genes", 3000),
        n_comps=reduce_cfg.get("n_comps", 50),
        n_neighbors=reduce_cfg.get("n_neighbors", 6),
        coord_type=reduce_cfg.get("coord_type", None),
        normalize_total=reduce_cfg.get("normalize_total", True),
        target_sum=reduce_cfg.get("target_sum", 1e4),
        log1p=reduce_cfg.get("log1p", True),
        flavor=reduce_cfg.get("flavor", "seurat"),
        inplace=True,
    )
    reduce_outputs = reduce_params["outputs"]
    print(
        f"      {reduce_outputs['n_hvgs']:,} HVGs selected, "
        f"{reduce_outputs['n_comps_computed']} PCs computed, "
        f"{reduce_outputs['spatial_graph_n_edges']:,} spatial graph edges "
        f"(mean {reduce_outputs['spatial_graph_mean_neighbors']} neighbours/spot)"
    )
    if reduce_params["params"].get("skipped_normalization"):
        print("      ⚠ Normalization skipped — benchmark dataset detected.")

    reduce_report_path = str(
        Path(output_dir) / f"{dataset_id}_spatial_reduce_report.html"
    )
    print(f"      Writing reduce report → {reduce_report_path}")
    generate_spatial_reduce_report(adata, reduce_report_path, dataset_id=dataset_id)

    print("\n✓ Spatial QC + Reduction complete.")
    print(f"  QC report:    {qc_report_path}")
    print(f"  Reduce report: {reduce_report_path}")

    # TODO Session 3: spatial_cluster + spatially_variable_genes
    # TODO Session 4: spatial_deconvolve + combined report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="OmicSage — Spatial Transcriptomics Pipeline"
    )
    parser.add_argument("--config", type=str, default=None,
                        help="Path to YAML config file")
    parser.add_argument("--source", type=str, default=None,
                        help="Override config: data source path or 'benchmark'")
    parser.add_argument("--spatial_type", type=str, default="auto",
                        help="Override config: 'auto', 'visium', 'h5ad', 'benchmark'")
    parser.add_argument("--dataset_id", type=str, default="spatial",
                        help="Dataset identifier used in output paths")
    parser.add_argument("--output_dir", type=str, default=None,
                        help="Override config: output directory")
    args = parser.parse_args()

    if args.config:
        config = load_config(args.config)
    else:
        config = {}

    # CLI overrides
    if args.source:
        config.setdefault("spatial", {})["source"] = args.source
        config.setdefault("spatial", {})["spatial_type"] = args.spatial_type
    if args.dataset_id:
        config["dataset_id"] = args.dataset_id
    if args.output_dir:
        config["output_dir"] = args.output_dir

    if "spatial" not in config or "source" not in config.get("spatial", {}):
        print("Error: no data source specified. "
              "Pass --source or provide a config with spatial.source",
              file=sys.stderr)
        sys.exit(1)

    run_spatial_pipeline(config)


if __name__ == "__main__":
    main()
