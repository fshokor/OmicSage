## Session Context
Date: 2026-05-05
Phase: 1 — Core scRNA Pipeline
Last thing completed: Multi-modal QC fix — modality-aware run_qc() returning MuData, 42/42 tests passing, notebook updated
File last worked on: notebooks/phase1_qc.ipynb

## Today's Goal
Build the normalization module — pipeline/modules/qc/normalize.py
ONE goal only — do not start PCA/UMAP until normalization is tested and working.

## Step 1 — Verify last session still works
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_phase0_structure.py tests/test_ingest.py tests/test_qc.py -v
# Expected: 42 passed (test_qc.py), 2 skipped
```

## Step 2 — Implement normalize.py
File to create: pipeline/modules/qc/normalize.py

Key requirements:
- Input: AnnData with raw counts in adata.X (from mdata["rna"], output of qc.py)
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
- Load from data/processed/GSE194122_cite_rna_qc.h5ad  ← note: new filename from MuData fix
- Run normalize.py
- Show HVG plot
- Save to data/processed/GSE194122_cite_normalized.h5ad

## Known Issues From Last Session
- Docker images still not built locally (intentional)
- Always use `python -m pytest` not bare `pytest` (system Python issue)
- Always `conda activate omicsage` before running anything
- OldFormatWarning from GSE194122 — harmless, file written with old anndata version
- Notebook must be opened from OmicSage root — add os.chdir('/home/shoko/OmicSage') to cell 1
- QC tests subsample CITE-seq to 5000 cells to stay within 7.6 GB RAM

## Files Modified This Session
- pipeline/modules/qc/qc.py              <- UPDATED: modality-aware, MuData return
- tests/test_qc.py                       <- UPDATED: MuData access pattern, 42 tests
- notebooks/phase1_qc.ipynb             <- UPDATED: MuData API, no manual GEX subsetting
- .dev_memory/NEXT_SESSION.md           <- UPDATED
- .dev_memory/PROGRESS.md              <- UPDATED
- docs/MODULE_DOCS.md                   <- UPDATED

## Verify This Session's Work
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/test_qc.py -v
# Expected: 42 passed, 2 skipped
```

## Relevant Context — QC API (v2, updated this session)

run_qc() now returns (MuData, dict) — not (AnnData, dict).

    mdata, metrics = run_qc(adata)
    mdata["rna"]   # filtered RNA AnnData — QC metrics in .obs
    mdata["adt"]   # CITE-seq only — ADT features, same filtered cells
    mdata["atac"]  # Multiome only — ATAC peaks, same filtered cells

Modality is auto-detected from adata.var['feature_types'].
Pass the full mixed AnnData — no manual GEX subsetting needed.

For normalization, always pass mdata["rna"]:
    adata_rna = mdata["rna"]
    adata_norm = normalize(adata_rna, ...)

## Relevant Context — GSE194122 Data Structure
CITE-seq file (our primary benchmark):
  - adata.X                = NORMALIZED values in raw file
  - adata.layers['counts'] = RAW integer counts <- ingest.py moves this to X
  - adata.var['feature_types'] = 'GEX' (RNA) or 'ADT' (protein)
  - adata.obs['cell_type'] = 8+ annotated cell types (ground truth for validation)
  - adata.obs['GEX_pct_counts_mt'] = precomputed MT% (validated r > 0.99 ✓)
  - adata.obs['batch'], obs['Site'], obs['DonorID'] = batch info for correction
  - adata.obsm['GEX_X_umap'], obsm['GEX_X_pca'] = ground truth embeddings

Processed files (output of QC, input for normalization):
  - data/processed/GSE194122_cite_rna_qc.h5ad    <- RNA only, filtered, raw counts in X
  - data/processed/GSE194122_cite_adt_qc.h5ad    <- ADT only, filtered, ready for CLR
  - data/processed/GSE194122_multiome_rna_qc.h5ad   <- RNA only, filtered
  - data/processed/GSE194122_multiome_atac_qc.h5ad  <- ATAC peaks, filtered, Phase 4
  - data/processed/GSE166635_HCC1_qc.h5ad        <- filtered, raw counts in X

Validation strategy:
  Our QC MT%   -> compare to obs['GEX_pct_counts_mt']  (done ✓, r > 0.99)
  Our clusters -> compare to obs['cell_type']            (Phase 1 milestone)
  Our UMAP     -> compare to obsm['GEX_X_umap']         (Phase 1 milestone)

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata, ipykernel, jupyter
