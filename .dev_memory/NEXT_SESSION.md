## Session Context
Date: 2026-05-07
Phase: 1 — Core scRNA Pipeline
Last thing completed: Annotation module Session B — REPLANNED
                      Decided against rpy2/R integration permanently.
                      ScType-py + SingleR-py deferred to dedicated future session.
                      Full decision spec saved in docs/ANNOTATION_PLAN.md
File last worked on: pipeline/modules/qc/annotate.py (no changes made this session)

## Today's Goal
DEG module — implement deg.py (Wilcoxon rank-sum differential expression,
per cluster vs rest and pairwise, with volcano plot output).
ONE goal — do not start GSEA until DEG is fully tested.

## Step 1 — Verify annotation module still works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_annotate.py -v
# Expected: 18 passed, 1 skipped
```

## Step 2 — Implement deg.py
File to create: pipeline/modules/qc/deg.py

DEG requirements:
- Input: annotated AnnData with obs['cell_type_vote'] and obs[leiden_col]
- Method: Wilcoxon rank-sum via sc.tl.rank_genes_groups()
- Mode 1: each cluster vs all other cells (one-vs-rest)
- Mode 2: pairwise between specified groups (optional)
- Multiple testing correction: Benjamini-Hochberg (FDR)
- Minimum thresholds: log2FC and adjusted p-value (configurable)
- Output: deg_dict with results DataFrame per group
- Write results to adata.uns['omicsage_deg']
- Graceful handling: warn if fewer than 10 cells in any group

## Step 3 — Write test_deg.py
File to create: tests/test_deg.py

Tests to write:
- test_deg_returns_anndata_and_dict()
- test_deg_uns_provenance_keys()
- test_deg_output_columns()         — names, scores, pvals, logfoldchanges, pvals_adj
- test_deg_pval_range()             — all pvals between 0 and 1
- test_deg_every_cluster_has_results()
- test_deg_logfc_threshold_filters()
- test_deg_pval_threshold_filters()
- test_deg_inplace_false()          — caller object unchanged
- test_deg_small_group_warning()    — <10 cells triggers UserWarning

## Step 4 — Implement deg_report.py
File to create: reports/deg_report.py

Report sections:
- Summary table: cluster → top 5 DEGs (gene, log2FC, adj_pval)
- Volcano plot per cluster (log2FC vs -log10 adj_pval)
- Dot plot: top N genes across all clusters
- HTML output matching annotate_report.py style

## Step 5 — Update notebook
File to edit: notebooks/phase1_qc.ipynb

Add Step 6 section:
- Run deg() on GSE194122_cite_annotated.h5ad
- Sanity check: known marker genes appear as top DEGs
  (e.g. CD3D in T cells, CD14 in Monocytes, MS4A1 in B cells)
- Save GSE194122_cite_deg.h5ad

## Known Issues Carried Forward
- Docker images still not built locally (intentional)
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless
- Notebook must be opened from OmicSage root — os.chdir('/home/shoko/OmicSage') in cell 1
- seurat_v3 HVG flavor needs ≥500 cells
- numpy bool comparison: use bool() cast before `is True` checks in tests
- float dict keys must be stringified before storing in uns
- obs['cell_type'] in GSE194122 is publication ground-truth — preserved as obs['cell_type_groundtruth']
- Consensus vote column is obs['cell_type_vote'] (not obs['cell_type'])
- rpy2 NOT installed — do not attempt R integration

## Deferred Work (do not start this session)
- ScType-py + SingleR-py annotation methods → see docs/ANNOTATION_PLAN.md
- GSEA module → blocked on DEG completion
- ADT and ATAC QC modules
- Batch correction
- Clustering module improvements

## Files Modified This Session
- (none — planning session only)
- docs/ANNOTATION_PLAN.md   ← CREATED: full spec for ScType-py + SingleR-py

## Verify This Session's Work
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_annotate.py -v
# Expected: 18 passed, 1 skipped
```

## Relevant Context — deg.py API (proposed)

```python
from pipeline.modules.qc.deg import deg

adata_deg, deg_dict = deg(
    adata_annotated,
    groupby="cell_type_vote",      # obs column to group by
    leiden_col="leiden",           # fallback if groupby missing
    method="wilcoxon",             # wilcoxon | t-test | logreg
    min_logfc=0.25,                # log2 fold-change threshold
    max_pval_adj=0.05,             # adjusted p-value threshold
    n_genes=200,                   # top N genes to return per group
    use_raw=False,                 # use adata.layers['logcounts']
    inplace=False,
)
```

obs columns: unchanged (DEG writes to uns only)
uns keys after deg():
  adata.uns['omicsage_deg']
    groupby, method, min_logfc, max_pval_adj,
    n_genes, n_groups, results_keys,
    scanpy_version, omicsage_module, timestamp

deg_dict keys:
  'results'       — dict of {group: DataFrame(gene, score, pval, logfc, pval_adj)}
  'summary_df'    — wide DataFrame: top 5 genes per group
  'provenance'    — same as uns

## Relevant Context — GSE194122 Data Structure
Processed files:
  - data/processed/GSE194122_cite_rna_qc.h5ad
  - data/processed/GSE194122_cite_normalized.h5ad
  - data/processed/GSE194122_cite_reduced.h5ad
  - data/processed/GSE194122_cite_clustered.h5ad
  - data/processed/GSE194122_cite_annotated.h5ad   ← input for DEG
  - data/processed/GSE194122_cite_deg.h5ad          ← output of next session
  - data/processed/GSE194122_cite_adt_qc.h5ad
  - data/processed/GSE194122_multiome_rna_qc.h5ad
  - data/processed/GSE194122_multiome_atac_qc.h5ad
  - data/processed/GSE166635_HCC1_qc.h5ad

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   rpy2 3.6.7 (installed but NOT used — R methods deferred)
