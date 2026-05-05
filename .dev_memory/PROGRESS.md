# OmicSage — Progress Tracker
> Updated: 2026-05-04 (end of session 3)
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
      - Uses backed='r' mode — works on files of any size (fixed this session)
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

### QC Module — COMPLETE
- [x] pipeline/modules/qc/qc.py — COMPLETE, ALL TESTS PASSING
      - MT gene auto-detection (MT- human / mt- mouse, gene_ids fallback)
      - Per-cell metrics: n_genes_by_counts, total_counts, pct_counts_mt
      - Scrublet doublet detection (graceful failure handling)
      - Configurable filters: min_genes, max_genes, max_mt_pct, remove_doublets
      - generate_report=True triggers HTML report generation
      - Validated: MT% correlation vs GEX_pct_counts_mt r > 0.99 ✓
- [x] pipeline/modules/qc/qc_report.py — COMPLETE
      - Self-contained HTML report, no external dependencies
      - Violin plots before/after, UMI vs genes scatter, doublet histogram
      - Summary cards, filter tables, ground-truth MT% correlation plot
- [x] tests/test_qc.py — 33 passed, 2 skipped
      - Real data only — no synthetic fixtures
      - Session-scoped fixtures, CITE-seq subsampled to 5000 cells
      - Ground-truth MT% validation (r > 0.99, MAE < 0.5%)

### Documentation
- [x] docs/MODULE_DOCS.md — documents all 4 QC module scripts
      (ingest.py, qc.py, qc_report.py, data_report.py)

### Notebook
- [x] notebooks/phase1_qc.ipynb — CREATED
      - Runs full QC pipeline on CITE-seq + HCC MTX
      - Generates HTML reports to reports/
      - Saves filtered AnnData to data/processed/

### Processed Data (output of QC)
- [x] data/processed/GSE194122_cite_qc.h5ad
- [x] data/processed/GSE166635_HCC1_qc.h5ad

### Processing — NEXT SESSION
- [ ] pipeline/modules/qc/normalize.py  ← NEXT
      - scran normalization
      - log1p transform
      - HVG selection (top 2000, seurat_v3)
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
| 2026-05-04 | Ground-truth MT% tests subset to feature_types='GEX' | CITE file mixes RNA+ADT; var_names_make_unique renames MT genes and breaks detection |
| 2026-05-04 | QC report embedded in qc.py via generate_report=True | Every step should produce a readable output — core OmicSage design principle |
| 2026-05-04 | Notebook-first workflow for running pipeline | Better visibility per step, good portfolio piece |
