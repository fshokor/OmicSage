## Session Context
Date: 2026-05-04
Phase: 1 — Core scRNA Pipeline
Last thing completed: QC module — qc.py + qc_report.py + 33/33 tests passing + notebook running
File last worked on: notebooks/phase1_qc.ipynb

## Today's Goal
Build the normalization module — pipeline/modules/qc/normalize.py
ONE goal only — do not start PCA/UMAP until normalization is tested and working.

## Step 1 — Verify last session still works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_phase0_structure.py tests/test_ingest.py tests/test_qc.py -v
# Expected: 82 + 33 = 115 passed, 2 skipped
```

## Step 2 — Implement normalize.py
File to create: pipeline/modules/qc/normalize.py

Key requirements:
- Input: AnnData with raw counts in adata.X (output of qc.py)
- Save raw counts to adata.layers['counts'] before modifying X
- Normalize with scran (scanpy's normalize_total as fallback)
- Log1p transform
- HVG selection: top 2000 highly variable genes (flavor='seurat_v3')
- Return normalized AnnData + normalization metrics dict
- Store normalization params in adata.uns['omicsage_normalization']

## Step 3 — Write tests
File to create: tests/test_normalize.py

Tests to write:
- test_raw_counts_preserved_in_layer()
- test_x_is_normalized_after_run()
- test_log1p_applied()
- test_hvg_selected()
- test_hvg_count_correct()
- test_normalization_params_in_uns()
- test_original_adata_not_mutated()

## Step 4 — Add normalization to notebook
Add a new section to notebooks/phase1_qc.ipynb:
- Load from data/processed/GSE194122_cite_qc.h5ad
- Run normalize.py
- Show HVG plot
- Save to data/processed/GSE194122_cite_normalized.h5ad

## Known Issues From Last Session
- Docker images still not built locally (intentional)
- Always use `python -m pytest` not bare `pytest` (system Python issue)
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless, file written with old anndata version
- CITE-seq h5ad has mixed RNA + ADT features (feature_types = 'GEX' / 'ADT')
  → always subset to GEX features before RNA-specific analysis
- Notebook must be opened from OmicSage root — add os.chdir('/home/shoko/OmicSage') to cell 1
- QC tests subsample CITE-seq to 5000 cells to stay within 7.6 GB RAM

## Files Modified This Session
- pipeline/modules/qc/qc.py              <- CREATED, 33/33 tests passing
- pipeline/modules/qc/qc_report.py       <- CREATED (HTML report generator)
- tests/test_qc.py                       <- CREATED, 33 passed / 2 skipped
- notebooks/phase1_qc.ipynb              <- CREATED (runs full QC pipeline)
- docs/MODULE_DOCS.md                    <- CREATED (documents all QC scripts)
- pipeline/modules/qc/data_report.py    <- UPDATED (backed='r' for large files)
- .dev_memory/NEXT_SESSION.md           <- UPDATED
- .dev_memory/PROGRESS.md              <- UPDATED

## Verify This Session's Work
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_qc.py -v
# Expected: 33 passed, 2 skipped
```

## Relevant Context — GSE194122 Data Structure
CITE-seq file (our primary benchmark):
  - adata.X                = NORMALIZED values in raw file
  - adata.layers['counts'] = RAW integer counts <- ingest.py moves this to X
  - adata.var['feature_types'] = 'GEX' (RNA) or 'ADT' (protein) — NOT 'Gene Expression'
  - adata.obs['cell_type'] = 8+ annotated cell types (ground truth for validation)
  - adata.obs['GEX_pct_counts_mt'] = precomputed MT% (validated r > 0.99 ✓)
  - adata.obs['batch'], obs['Site'], obs['DonorID'] = batch info for correction
  - adata.obsm['GEX_X_umap'], obsm['GEX_X_pca'] = ground truth embeddings
  - var_names NOT unique — RNA + ADT mixed, call var_names_make_unique() in QC

Processed files (output of QC, input for normalization):
  - data/processed/GSE194122_cite_qc.h5ad    <- filtered AnnData, raw counts in X
  - data/processed/GSE166635_HCC1_qc.h5ad   <- filtered AnnData, raw counts in X

Validation strategy:
  Our QC MT%   -> compare to obs['GEX_pct_counts_mt']  (done ✓, r > 0.99)
  Our clusters -> compare to obs['cell_type']            (Phase 1 milestone)
  Our UMAP     -> compare to obsm['GEX_X_umap']         (Phase 1 milestone)

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, ipykernel, jupyter
