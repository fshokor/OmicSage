"""
OmicSage — Phase 1 Pipeline
Dataset : GSE194122 BMMC CITE-seq (NeurIPS 2021)
Run     : python run_GSE194122.py

Steps
-----
1.  QC              — cell filtering, doublet removal
2.  Normalize       — CP10K, log1p, HVG selection
3.  Reduce          — PCA, kNN graph, UMAP
4.  Cluster         — Leiden sweep, best-resolution auto-selection
5.  Annotate        — CellTypist + marker scoring + majority vote
6.  DEG             — Wilcoxon one-vs-rest
7.  GSEA            — ORA via Enrichr (GO BP, KEGG, Reactome)
8.  Harmony         — batch correction, recompute UMAP
9.  Cluster harmony — Leiden on Harmony graph
10. Pseudobulk DEG  — DESeq2 Wald test via pydeseq2

Processed files → data/processed/GSE194122/
Reports         → reports/GSE194122/
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")
from datetime import datetime
import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ── repo root ──────────────────────────────────────────────────────────────────
root = Path(__file__).resolve().parent
while not (root / "pipeline").exists() and root != root.parent:
    root = root.parent
os.chdir(root)
sys.path.insert(0, str(root))

# ── output dirs ───────────────────────────────────────────────────────────────
PROCESSED = Path("data/processed/GSE194122")
REPORTS   = Path("reports/GSE194122")
PROCESSED.mkdir(parents=True, exist_ok=True)
REPORTS.mkdir(parents=True, exist_ok=True)

# ── input ──────────────────────────────────────────────────────────────────────
RAW = Path("data/benchmark/GSE194122_cite_raw_only.h5ad")

sc.settings.verbosity = 1


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — QC
# ══════════════════════════════════════════════════════════════════════════════
def step1_qc() -> Path:
    out = PROCESSED / "01_qc.h5ad"
    if out.exists():
        print(f"[step 1] QC — cached {out}")
        return out

    print("[step 1] QC — running …")
    from pipeline.modules.qc.ingest import load_dataset
    from pipeline.modules.qc.qc import run_qc

    adata = load_dataset(RAW, sample_name="GSE194122_BMMC_CITE")
    mdata, metrics = run_qc(
        adata,
        min_genes=200,
        max_genes=2500,
        max_mt_pct=5.0,
        remove_doublets=True,
        generate_report=True,
        report_path=str(REPORTS / "01_qc_report.html"),
        sample_name="GSE194122_BMMC_CITE",
    )

    rna = mdata["rna"]
    rna.write_h5ad(out)
    if "adt" in mdata.mod:
        mdata["adt"].write_h5ad(PROCESSED / "01_qc_adt.h5ad")

    pass_rate = 100 * metrics["n_cells_output"] / metrics["n_cells_input"]
    print(
        f"[step 1] QC — {metrics['n_cells_output']:,} cells kept"
        f" ({pass_rate:.1f}%)  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — NORMALIZE
# ══════════════════════════════════════════════════════════════════════════════
def step2_normalize(qc_path: Path) -> Path:
    out = PROCESSED / "02_normalized.h5ad"
    if out.exists():
        print(f"[step 2] Normalize — cached {out}")
        return out

    print("[step 2] Normalize — running …")
    from pipeline.modules.qc.normalize import normalize
    from pipeline.modules.qc.normalization_report import run_normalization_report

    adata = sc.read_h5ad(qc_path)
    adata_norm, metrics = normalize(
        adata,
        batch_key="batch",
        target_sum=1e4,
        n_top_genes=2000,
        hvg_flavor="seurat",
        inplace=False,
    )

    run_normalization_report(
        adata_norm=adata_norm,
        metrics=metrics,
        report_path=str(REPORTS / "02_normalization_report.html"),
        dataset_name="GSE194122_CITE",
    )

    adata_norm.write_h5ad(out)
    print(f"[step 2] Normalize — HVGs={metrics['n_hvg_selected']}  ->  {out}")
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — DIMENSIONALITY REDUCTION
# ══════════════════════════════════════════════════════════════════════════════
def step3_reduce(norm_path: Path) -> Path:
    out = PROCESSED / "03_reduced.h5ad"
    if out.exists():
        print(f"[step 3] Reduce — cached {out}")
        return out

    print("[step 3] Reduce — running …")
    from pipeline.modules.qc.reduce import reduce
    from pipeline.modules.qc.reduce_report import run_reduce_report

    adata = sc.read_h5ad(norm_path)
    adata_reduced, metrics = reduce(
        adata,
        n_comps=50,
        n_pcs=None,
        n_pcs_method="elbow",
        n_neighbors=15,
        inplace=False,
    )

    run_reduce_report(
        adata_reduced=adata_reduced,
        metrics=metrics,
        report_path=str(REPORTS / "03_reduce_report.html"),
        dataset_name="GSE194122_CITE",
    )

    adata_reduced.write_h5ad(out)
    prov = adata_reduced.uns["omicsage_reduce"]
    print(
        f"[step 3] Reduce — {prov['n_pcs_used']} PCs"
        f" ({prov['cumulative_variance_explained_by_selected_pcs']*100:.1f}% var)"
        f"  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — CLUSTERING
# ══════════════════════════════════════════════════════════════════════════════
def step4_cluster(reduced_path: Path) -> Path:
    out = PROCESSED / "04_clustered.h5ad"
    if out.exists():
        print(f"[step 4] Cluster — cached {out}")
        return out

    print("[step 4] Cluster — running …")
    from pipeline.modules.clustering.cluster import cluster
    from pipeline.modules.clustering.cluster_report import run_cluster_report

    adata = sc.read_h5ad(reduced_path)
    adata_clustered, metrics = cluster(
        adata,
        resolution_range=[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5],
        inplace=False,
    )

    run_cluster_report(
        adata_clustered=adata_clustered,
        metrics=metrics,
        report_path=str(REPORTS / "04_cluster_report.html"),
        dataset_name="GSE194122_CITE",
    )

    adata_clustered.write_h5ad(out)
    print(
        f"[step 4] Cluster — res={metrics['best_resolution']}"
        f"  {metrics['best_n_clusters']} clusters  ({metrics['selection_reason']})"
        f"  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — ANNOTATION
# ══════════════════════════════════════════════════════════════════════════════
def step5_annotate(clustered_path: Path) -> Path:
    out = PROCESSED / "05_annotated.h5ad"
    if out.exists():
        print(f"[step 5] Annotate — cached {out}")
        return out

    print("[step 5] Annotate — running …")
    from pipeline.modules.annotation.annotate import annotate
    from pipeline.modules.annotation.annotate_report import run_annotate_report

    adata = sc.read_h5ad(clustered_path)
    adata_annotated, annotation_dict = annotate(
        adata,
        methods=["celltypist", "markers", "vote"],
        leiden_col="leiden",
        celltypist_models=["Immune_All_High.pkl", "Immune_All_Low.pkl"],
        inplace=False,
    )

    run_annotate_report(
        adata_annotated=adata_annotated,
        annotation_dict=annotation_dict,
        report_path=str(REPORTS / "05_annotate_report.html"),
        dataset_name="GSE194122_CITE",
    )

    adata_annotated.write_h5ad(out)
    prov = annotation_dict["provenance"]
    n_types = adata_annotated.obs["cell_type_vote"].nunique()
    print(
        f"[step 5] Annotate — {n_types} consensus types"
        f"  methods={prov['methods_run']}"
        f"  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — DEG (Wilcoxon)
# ══════════════════════════════════════════════════════════════════════════════
def step6_deg(annotated_path: Path) -> tuple[Path, dict]:
    out = PROCESSED / "06_deg.h5ad"

    print("[step 6] DEG — running …")
    from pipeline.modules.downstream.deg import deg
    from pipeline.modules.downstream.deg_report import generate_deg_report

    adata = sc.read_h5ad(annotated_path)
    adata_deg, deg_dict = deg(
        adata,
        groupby="cell_type_vote",
        method="wilcoxon",
        min_logfc=0.25,
        max_pval_adj=0.05,
        n_genes=500,
        inplace=False,
    )

    generate_deg_report(
        adata=adata_deg,
        deg_dict=deg_dict,
        output_path=str(REPORTS / "06_deg_report.html"),
        top_n_volcano=10,
        top_n_dotplot=5,
        max_volcano_groups=9,
    )

    adata_deg.write_h5ad(out)
    total_sig = sum(len(df) for df in deg_dict["results"].values())
    print(
        f"[step 6] DEG — {len(deg_dict['results'])} groups"
        f"  {total_sig:,} significant DEGs  ->  {out}"
    )
    return out, deg_dict


# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — GSEA
# ══════════════════════════════════════════════════════════════════════════════
def step7_gsea(deg_path: Path, deg_dict: dict) -> Path:
    out = PROCESSED / "07_gsea.h5ad"

    print("[step 7] GSEA — running …")
    from pipeline.modules.downstream.gsea import gsea
    from pipeline.modules.downstream.gsea_report import generate_gsea_report

    adata = sc.read_h5ad(deg_path)
    adata_gsea, gsea_dict = gsea(
        adata,
        deg_dict=deg_dict,
        gene_sets=[
            "GO_Biological_Process_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
        ],
        min_logfc=0.25,
        max_pval_adj=0.05,
        top_n_genes=None,
        min_genes=5,
        organism="human",
        inplace=False,
    )

    generate_gsea_report(
        gsea_dict=gsea_dict,
        output_path=str(REPORTS / "07_gsea_report.html"),
        top_n_table=5,
        top_n_bar=10,
    )

    adata_gsea.write_h5ad(out)
    prov = gsea_dict["provenance"]
    print(
        f"[step 7] GSEA — {prov['n_groups_tested']} groups tested"
        f"  {prov['n_groups_skipped']} skipped  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — HARMONY BATCH CORRECTION
# ══════════════════════════════════════════════════════════════════════════════
def step8_harmony(gsea_path: Path) -> Path:
    out = PROCESSED / "08_harmony.h5ad"
    if out.exists():
        print(f"[step 8] Harmony — cached {out}")
        return out

    print("[step 8] Harmony — running …")
    from pipeline.modules.integration.harmony_correct import harmony_correct
    from pipeline.modules.integration.harmony_report import generate_harmony_report

    adata = sc.read_h5ad(gsea_path)
    adata = harmony_correct(
        adata,
        batch_key="batch",
        n_pcs=50,
        n_neighbors=15,
        umap_min_dist=0.3,
        random_state=0,
        copy=False,
    )

    generate_harmony_report(
        adata=adata,
        output_path=str(REPORTS / "08_harmony_report.html"),
    )

    adata.write_h5ad(out)
    prov = adata.uns["omicsage_harmony"]
    print(
        f"[step 8] Harmony — {prov['n_batches']} batches corrected"
        f"  {prov['elapsed_seconds']:.1f}s  ->  {out}"
    )
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — CLUSTERING ON HARMONY GRAPH
# ══════════════════════════════════════════════════════════════════════════════
def step9_cluster_harmony(harmony_path: Path) -> Path:
    out = PROCESSED / "09_harmony_clustered.h5ad"
    if out.exists():
        print(f"[step 9] Cluster (Harmony) — cached {out}")
        return out

    print("[step 9] Cluster (Harmony) — running …")
    from pipeline.modules.clustering.cluster import cluster, compute_ari

    adata = sc.read_h5ad(harmony_path)

    # pre-correction clustering (original kNN)
    adata, metrics_pre = cluster(
        adata,
        resolution_range=[0.2, 0.4, 0.6, 0.8, 1.0],
        cluster_key="leiden",
        inplace=True,
    )

    # post-correction clustering (Harmony kNN)
    adata, metrics_post = cluster(
        adata,
        resolution_range=[0.2, 0.4, 0.6, 0.8, 1.0],
        neighbors_key="neighbors_harmony",
        cluster_key="leiden_harmony",
        inplace=True,
    )

    ari = compute_ari(adata, "leiden", "leiden_harmony")
    print(
        f"[step 9] Cluster (Harmony)"
        f"  pre={metrics_pre['best_n_clusters']} clusters"
        f"  post={metrics_post['best_n_clusters']} clusters"
        f"  ARI={ari:.4f}  ->  {out}"
    )

    adata.write_h5ad(out)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 — PSEUDOBULK DEG
# ══════════════════════════════════════════════════════════════════════════════
def step10_pseudobulk(annotated_path: Path) -> Path:
    out = PROCESSED / "10_pseudobulk_deg.h5ad"

    print("[step 10] Pseudobulk DEG — running …")
    from pipeline.modules.downstream.pseudobulk_deg import pseudobulk_deg
    from pipeline.modules.downstream.pseudobulk_deg_report import (
        generate_pseudobulk_deg_report,
    )

    adata = sc.read_h5ad(annotated_path)
    adata_pb, pb_dict = pseudobulk_deg(
        adata,
        groupby="cell_type_vote",
        donor_key="batch",
        counts_layer="counts",
        min_cells=10,
        min_samples=3,
        min_logfc=0.25,
        max_pval_adj=0.05,
        inplace=False,
    )

    generate_pseudobulk_deg_report(
        adata=adata_pb,
        pb_dict=pb_dict,
        output_path=str(REPORTS / "10_pseudobulk_deg_report.html"),
        top_n_volcano=10,
        top_n_dotplot=5,
        max_volcano_groups=11,
    )

    adata_pb.write_h5ad(out)
    prov = pb_dict["provenance"]
    total_sig = sum(len(df) for df in pb_dict["results"].values())
    print(
        f"[step 10] Pseudobulk DEG — {prov['n_groups']} groups tested"
        f"  {prov['n_skipped']} skipped"
        f"  {total_sig:,} significant DEGs  ->  {out}"
    )
    return out

def _load_deg(annotated_path: Path) -> tuple[Path, dict]:
    """Re-run deg() to recover deg_dict when skipping step 6."""
    from pipeline.modules.downstream.deg import deg
    adata = sc.read_h5ad(annotated_path)
    _, deg_dict = deg(
        adata,
        groupby="cell_type_vote",
        method="wilcoxon",
        min_logfc=0.25,
        max_pval_adj=0.05,
        n_genes=500,
        inplace=False,
    )
    return PROCESSED / "06_deg.h5ad", deg_dict


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    
    start_time = datetime.now()
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--from-step", type=int, default=1, metavar="N",
        help="Resume from step N (1-10). Steps before N must already be cached."
    )
    args = parser.parse_args()
    start = args.from_step
    print("=" * 60)
    print("OmicSage — GSE194122 CITE-seq BMMC pipeline")
    print("=" * 60)

    qc_path        = step1_qc()        if start <= 1  else PROCESSED / "01_qc.h5ad"
    norm_path      = step2_normalize(qc_path)  if start <= 2  else PROCESSED / "02_normalized.h5ad"
    reduced_path   = step3_reduce(norm_path)   if start <= 3  else PROCESSED / "03_reduced.h5ad"
    clustered_path = step4_cluster(reduced_path) if start <= 4 else PROCESSED / "04_clustered.h5ad"
    annotated_path = step5_annotate(clustered_path) if start <= 5 else PROCESSED / "05_annotated.h5ad"
    deg_path, deg_dict = step6_deg(annotated_path) if start <= 6 else _load_deg(annotated_path)
    gsea_path      = step7_gsea(deg_path, deg_dict) if start <= 7 else PROCESSED / "07_gsea.h5ad"
    harmony_path   = step8_harmony(gsea_path)  if start <= 8  else PROCESSED / "08_harmony.h5ad"
    _              = step9_cluster_harmony(harmony_path) if start <= 9 else None
    _              = step10_pseudobulk(annotated_path)   if start <= 10 else None

    end_time = datetime.now()
    elapsed  = end_time - start_time
    hours, remainder = divmod(int(elapsed.total_seconds()), 3600)
    minutes, seconds = divmod(remainder, 60)

    print()
    print("=" * 60)
    print("Pipeline complete.")
    print(f"Processed files : {PROCESSED}")
    print(f"Reports         : {REPORTS}")
    print(f"Started         : {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Finished        : {end_time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Elapsed         : {hours:02d}h {minutes:02d}m {seconds:02d}s")
    print("=" * 60)


if __name__ == "__main__":
    main()
