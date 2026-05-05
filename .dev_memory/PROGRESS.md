# OmicSage — Progress Tracker
> Updated: 2026-05-05 (end of session 4)
> Update this file at the end of every session.

---

## Phase 0 — Foundation COMPLETE

- [x] GitHub repo created (https://github.com/fshokor/OmicSage)
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [x] Full repo structure scaffolded
- [x] Docker base images defined (Dockerfile.python + Dockerfile.r)
- [x] CI/CD via GitHub Actions
- [x] YAML config schema defined
- [x] .dev_memory/ system initialized
- [x] README with project description and badges
- [x] .gitignore — entire data/ directory excluded

---

## Phase 1 — Core scRNA Pipeline IN PROGRESS

### Environment
- [x] Conda environment created (omicsage, Python 3.11.15)
- [x] All Phase 1 packages installed and verified (scanpy 1.11.5)
- [x] GEOparse installed
- [x] ipykernel + jupyter installed
- [x] mudata installed
- [x] Always use: conda activate omicsage && python -m pytest

### Benchmark Dataset
- [x] Dataset selected: GSE194122 NeurIPS 2021 BMMC (switched from GSE166635)
- [x] CITE-seq file downloaded and unzipped (~2.5 GB)
- [ ] Multiome file — large, deferred
- [x] Data structure fully understood (see NEXT_SESSION.md)

### Utility Scripts
- [x] scripts/download_benchmark.py — downloads GSE194122 (CITE + multiome)
- [x] scripts/download_test_data.py — downloads GSE166635 into data/test/
- [x] scripts/__init__.py

### Data Intake Report
- [x] pipeline/modules/qc/data_report.py
      - Works for any .h5ad (public or personal)
      - Uses backed='r' mode — works on files of any size (fixed session 3)
      - Fetches GEO metadata via NCBI eutils API if --geo passed
      - Static matplotlib QC plots embedded in HTML
      - Inventories modalities, obs/var columns, embeddings, layers
- [x] Report generated and verified for GSE194122 CITE file

### Data Ingestion
- [x] pipeline/modules/qc/ingest.py — COMPLETE, ALL TESTS PASSING
      - Auto-detects .h5ad / .h5 / MTX directory formats
      - Extracts raw counts from layers['counts'] for processed h5ad
      - Falls back to adata.raw, then layer scan
      - Populates adata.uns['omicsage_source'] provenance metadata
      - Ensures sparse CSR output
- [x] tests/test_ingest.py — 30/30 tests passing

### QC Module — COMPLETE (v2)
- [x] pipeline/modules/qc/qc.py — COMPLETE, ALL TESTS PASSING
      - Modality auto-detection: rna / cite / multiome via feature_types
      - GEX subsetting handled internally — callers pass the full mixed AnnData
      - MT gene detection, per-cell metrics, Scrublet doublet detection
      - Configurable filters: min_genes, max_genes, max_mt_pct, remove_doublets
      - Returns MuData: mdata["rna"] always; mdata["adt"] for CITE; mdata["atac"] for Multiome
      - QC obs columns (n_genes_by_counts, total_counts, pct_counts_mt, doublet_score)
        live on mdata["rna"] only — other modalities start clean
      - generate_report=True triggers HTML report generation
      - Validated: MT% correlation vs GEX_pct_counts_mt r > 0.99 ✓
- [x] pipeline/modules/qc/qc_report.py — COMPLETE
      - Self-contained HTML report, no external dependencies
      - Violin plots before/after, UMI vs genes scatter, doublet histogram
      - Summary cards, filter tables, ground-truth MT% correlation plot
- [x] tests/test_qc.py — 42 passed, 2 skipped
      - Real data only — no synthetic fixtures (except synthetic Multiome AnnData
        in test_detect_modality_multiome_synthetic)
      - Session-scoped fixtures, CITE-seq subsampled to 5000 cells
      - Ground-truth MT% validation (r > 0.99, MAE < 0.5%)
      - Multi-modal tests: modality detection, MuData structure, ADT preservation,
        total_counts RNA-only validation, barcode alignment across modalities

### Documentation
- [x] docs/MODULE_DOCS.md — documents all 4 QC module scripts
      (ingest.py, qc.py, qc_report.py, data_report.py)
      Updated session 4: MuData API, modality parameter, new usage examples

### Notebook
- [x] notebooks/phase1_qc.ipynb — UPDATED (session 4)
      - Passes full mixed AnnData to run_qc() — no manual GEX subsetting
      - Accesses results via mdata["rna"], mdata["adt"], mdata["atac"]
      - Saves separate .h5ad per modality
      - MT% validation uses permissive run_qc() directly — no manual RNA subset

### Processed Data (output of QC, input for normalization)
- [x] data/processed/GSE194122_cite_rna_qc.h5ad     ← RNA only, ready for normalization
- [x] data/processed/GSE194122_cite_adt_qc.h5ad     ← ADT only, ready for CLR (future)
- [x] data/processed/GSE194122_multiome_rna_qc.h5ad ← RNA only, ready for normalization
- [x] data/processed/GSE194122_multiome_atac_qc.h5ad← ATAC peaks, ready for Phase 4
- [x] data/processed/GSE166635_HCC1_qc.h5ad         ← RNA only, ready for normalization

### Processing — NEXT SESSION
- [ ] pipeline/modules/qc/normalize.py  ← NEXT
      - Input: mdata["rna"] from qc.py (raw counts in X)
      - scran normalization (normalize_total as fallback)
      - log1p transform
      - HVG selection (top 2000, seurat_v3)
      - Store params in adata.uns['omicsage_normalization']
- [ ] tests/test_normalize.py
- [ ] PCA + UMAP + t-SNE
- [ ] Batch correction (Harmony + scVI)

### Clustering — NOT STARTED
- [ ] Leiden clustering

### Annotation — NOT STARTED
- [ ] SingleR annotation

### DEG + Pathway — NOT STARTED

### MILESTONE — Phase 1 Complete
- [ ] Our clusters match GSE194122 cell_type labels (>80% agreement)
- [ ] Our UMAP matches GEX_X_umap structure visually

---

## Phase 2 — Report Engine NOT STARTED
## Phase 3 — AI Layer NOT STARTED
## Phase 4 — scATAC Module NOT STARTED
## Phase 5 — Spatial Module NOT STARTED
## Phase 6 — Multiome Integration NOT STARTED
## Phase 7 — User Interfaces NOT STARTED
## Phase 8 — Benchmarking + Paper NOT STARTED

---

## Dataset Registry

| File | Modality | Status | Use |
|------|----------|--------|-----|
| GSE194122 CITE BMMC processed.h5ad | RNA + ADT | Ready | Phase 1 dev + validation |
| GSE194122 multiome BMMC processed.h5ad | RNA + ATAC | Deferred (large) | Phase 1 + Phase 4 |
| GSE166635 HCC1+HCC2 MTX | RNA | Ready (data/test/) | MTX ingestion + QC testing |
| PBMC 10k (10x Genomics) | RNA | Not downloaded | Validation |
| 10x Visium human brain | Spatial | Not downloaded | Phase 5 |

---

## Key Decisions Log

| Date | Decision | Reason |
|------|----------|--------|
| 2026-04-30 | Switched benchmark from GSE166635 to GSE194122 | Already analyzed GSE166635; GSE194122 is multi-modal and a standard benchmark |
| 2026-04-30 | Exclude entire data/ from git | Simpler and safer than per-extension rules |
| 2026-04-30 | data_report.py outputs static HTML | No Quarto dependency for intake report |
| 2026-04-30 | GSE166635 kept as test fixture in data/test/ | Good real MTX data for testing ingestion MTX path |
| 2026-04-30 | ingest.py checks layers then adata.raw for raw counts | GSE194122 uses layers['counts']; other datasets may use adata.raw |
| 2026-05-04 | data_report.py uses backed='r' | Large h5ad files (2.9 GB) were killing the process |
| 2026-05-04 | QC tests use real data only, no synthetic fixtures | Synthetic fixtures produced unrealistic data and failed; real data is more meaningful |
| 2026-05-04 | CITE-seq tests subsample to 5000 cells | 90k cells × 14k genes exceeds 7.6 GB RAM when running Scrublet |
| 2026-05-04 | QC report embedded in qc.py via generate_report=True | Every step should produce a readable output — core OmicSage design principle |
| 2026-05-04 | Notebook-first workflow for running pipeline | Better visibility per step, good portfolio piece |
| 2026-05-05 | run_qc() returns MuData instead of AnnData | Consistent API for all modalities; RNA/ADT/ATAC always accessible by key; no manual splitting in callers |
| 2026-05-05 | QC metrics computed on RNA only, other modalities start clean | Biologically correct — ADT counts must not inflate total_counts; ATAC peaks must not inflate n_genes |
| 2026-05-05 | ADT and ATAC QC deferred to separate functions | Each modality has its own QC criteria (CLR/isotype for ADT; TSS/FRiP for ATAC); adding them to run_qc() today would have been out of scope |
| 2026-05-05 | Processed files renamed to *_rna_qc.h5ad and *_adt_qc.h5ad | Modality-explicit naming prevents confusion now that outputs are split per modality |
