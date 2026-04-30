## Session Context
Date: 2026-05-01 (or whenever next session starts)
Phase: 1 — Core scRNA Pipeline
Last thing completed: Data ingestion module — ingest.py — 30/30 tests passing
File last worked on: tests/test_ingest.py

## Today's Goal
Build the QC module — pipeline/modules/qc/qc.py
ONE goal only — do not start normalization until QC is tested and working.

## Step 1 — Verify last session still works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_phase0_structure.py tests/test_ingest.py -v
# Expected: 52 + 30 = 82 passed
```

## Step 2 — Implement qc.py
File to create: pipeline/modules/qc/qc.py

Key requirements:
- Input: AnnData with raw counts in adata.X (output of ingest.py)
- Compute per-cell metrics:
    - n_genes_by_counts    (genes detected per cell)
    - total_counts         (UMI sum per cell)
    - pct_counts_mt        (mitochondrial %)
- Detect MT genes automatically (prefix MT- or mt-)
- Run Scrublet for doublet detection → add obs['doublet_score'] + obs['predicted_doublet']
- Apply filters with configurable thresholds:
    - min_genes (default 200)
    - max_genes (default 6000)
    - max_mt_pct (default 20)
    - remove_doublets (default True)
- Return filtered AnnData + QC metrics dict
- Validate against GSE194122 ground truth: compare our pct_counts_mt to obs['GEX_pct_counts_mt']

## Step 3 — Write tests
File to create: tests/test_qc.py

Tests to write:
- test_mt_genes_detected()
- test_metrics_computed()
- test_filter_removes_low_quality_cells()
- test_doublet_scores_added()
- test_filtered_adata_has_correct_shape()
- test_qc_metrics_match_ground_truth()  <- compare to GEX_pct_counts_mt

## Step 4 — Run against real data
```bash
python3 -c "
from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc
adata = load_dataset('data/benchmark/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad')
adata_filtered, metrics = run_qc(adata)
print('Cells before QC:', adata.n_obs)
print('Cells after QC:', adata_filtered.n_obs)
print('Cells removed:', adata.n_obs - adata_filtered.n_obs)
"
```

## Known Issues From Last Session
- Docker images still not built locally (intentional)
- Always use `python -m pytest` not bare `pytest` (system Python issue)
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless, file written with old anndata version
- var_names not unique in CITE h5ad (RNA + ADT mixed) — handle in QC with var_names_make_unique()

## Files Modified This Session
- pipeline/modules/qc/ingest.py          <- CREATED, 30/30 tests passing
- tests/test_ingest.py                   <- CREATED
- scripts/download_test_data.py          <- CREATED (downloads GSE166635 for MTX testing)
- scripts/__init__.py                    <- CREATED
- pipeline/modules/qc/data_report.py    <- CREATED (earlier this session)
- scripts/download_benchmark.py         <- CREATED (earlier this session)
- .gitignore                            <- data/ fully excluded
- .dev_memory/NEXT_SESSION.md           <- UPDATED
- .dev_memory/PROGRESS.md              <- UPDATED

## Verify This Session's Work
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_ingest.py -v
# Expected: 30 passed
```

## Relevant Context — GSE194122 Data Structure
CITE-seq file (our primary benchmark):
  - adata.X                = NORMALIZED values
  - adata.layers['counts'] = RAW integer counts <- ingest.py moves this to X automatically
  - adata.obs['cell_type'] = 8+ annotated cell types (ground truth for validation)
  - adata.obs['GEX_pct_counts_mt'] = precomputed MT% (compare against our QC output)
  - adata.obs['batch'], obs['Site'], obs['DonorID'] = batch info for correction
  - adata.obsm['GEX_X_umap'], obsm['GEX_X_pca'] = ground truth embeddings
  - adata.var has mixed RNA + ADT features — call var_names_make_unique() in QC

Validation strategy:
  Our QC MT%   -> compare to obs['GEX_pct_counts_mt']  (should correlate > 0.99)
  Our clusters -> compare to obs['cell_type']            (Phase 1 milestone)
  Our UMAP     -> compare to obsm['GEX_X_umap']         (Phase 1 milestone)

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3
