# OmicSage — Current Status
> Last updated: 2026-05-05

## Phase
Phase 1 — Core scRNA Pipeline

## What Is Built and Tested Right Now

### ✅ Phase 0 — Foundation
- Repo structure scaffolded
- Docker base images (Python + R) — defined, not yet built locally
- CI/CD via GitHub Actions — configured
- YAML config schema — defined
- .dev_memory/ system — initialized

### ✅ Ingestion (pipeline/modules/qc/ingest.py)
- Auto-detects 10x MEX, H5, AnnData formats
- Moves normalized values out of .X, puts raw counts in .X
- Handles GSE194122 CITE-seq and multiome, GSE166635 HCC

### ✅ QC (pipeline/modules/qc/qc.py)
- Modality-aware: auto-detects GEX / ADT / ATAC from var['feature_types']
- Returns MuData with mdata["rna"], mdata["adt"], mdata["atac"] slots
- MT%, genes/cell, Scrublet doublets, SoupX ambient RNA
- 42 tests passing in tests/test_qc.py

### ✅ Normalization (pipeline/modules/qc/normalize.py)  ← NEW THIS SESSION
- Input: AnnData with raw counts in .X (from mdata["rna"])
- Saves raw counts to layers['counts']
- Normalizes to CP10K (normalize_total)
- log1p transform
- Saves log1p values to layers['logcounts']  (Seurat convention)
- HVG selection — top 2000 genes, flavor='seurat_v3' (pre-log1p on counts)
- batch_key support for per-batch HVG selection
- Input validation — rejects already-normalized data
- Provenance stored in uns['omicsage_normalization']
- 12 tests passing in tests/test_normalize.py

### ✅ Normalization Report (reports/normalization_report.py)  ← NEW THIS SESSION
- Self-contained HTML report (no Quarto dependency yet)
- Figures: HVG scatter, library size violin, top 20 HVGs, gene detection rate
- Summary table + provenance table from uns
- Callable from notebook or CLI
- CLI: python reports/normalization_report.py --input ... --output ... --report ...

## Total Tests Passing
138 (pre-session) + 12 (normalize) = 150 expected after this session

## What Is NOT Built Yet
- reduce.py     — PCA + UMAP + neighbors (next session)
- cluster.py    — Leiden clustering
- annotate.py   — SingleR + LLM cell type annotation
- DEG module
- GSEA module
- Trajectory
- scATAC module (Phase 4)
- Spatial module (Phase 5)
- Multiome module (Phase 6)
- Streamlit UI (Phase 7)
- CLI (Phase 7)
- Quarto reports (Phase 2)

## Processed Data Files
- data/processed/GSE194122_cite_rna_qc.h5ad        ← QC-filtered RNA, raw counts in .X
- data/processed/GSE194122_cite_adt_qc.h5ad        ← QC-filtered ADT
- data/processed/GSE194122_multiome_rna_qc.h5ad    ← QC-filtered multiome RNA
- data/processed/GSE194122_multiome_atac_qc.h5ad   ← QC-filtered ATAC (Phase 4)
- data/processed/GSE166635_HCC1_qc.h5ad            ← QC-filtered HCC RNA
- data/processed/GSE194122_cite_normalized.h5ad    ← TO CREATE in notebook
