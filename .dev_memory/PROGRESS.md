# OmicSage — Progress Tracker
> Updated: 2026-04-30 (end of session 2)
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
- [x] Always use: conda activate omicsage && python -m pytest

### Benchmark Dataset
- [x] Dataset selected: GSE194122 NeurIPS 2021 BMMC (switched from GSE166635)
- [x] CITE-seq file downloaded and unzipped (~2.5 GB)
- [ ] Multiome file — downloading/unzipping (in progress)
- [x] Data structure fully understood (see NEXT_SESSION.md)

### Utility Scripts
- [x] scripts/download_benchmark.py — downloads GSE194122 (CITE + multiome)
- [x] scripts/download_test_data.py — downloads GSE166635 into data/test/
- [x] scripts/__init__.py

### Data Intake Report
- [x] pipeline/modules/qc/data_report.py
      - Works for any .h5ad (public or personal)
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
- [x] tests/test_ingest.py — 30/30 tests passing (98s)
      - Unit tests: format detection, integer matrix check, layer detection
      - Contract tests: shape, dtype, uns, sample column, normalized layer
      - Real data: GSE166635 HCC1 MTX + GSE194122 CITE h5ad

### QC Module — NEXT SESSION
- [ ] pipeline/modules/qc/qc.py
- [ ] tests/test_qc.py

### Processing — NOT STARTED
- [ ] Normalization (scran)
- [ ] HVG selection
- [ ] PCA + UMAP
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
| GSE194122 multiome BMMC processed.h5ad | RNA + ATAC | Downloading | Phase 1 + Phase 4 |
| GSE166635 HCC1+HCC2 MTX | RNA | Ready (data/test/) | Test fixture for MTX ingestion |
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
