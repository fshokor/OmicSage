## Session Context
Date: 2026-04-30
Phase: 1 — Core scRNA Pipeline
Last thing completed: Phase 0 fully done — repo, Docker, CI/CD all green on GitHub
File last worked on: .github/workflows/ci.yml

## Today's Goal
Download GSE166635 benchmark dataset and implement the data ingestion module.
ONE goal only — do not start QC until ingestion is tested and working.

## Step 1 — Download benchmark data
```bash
cd ~/OmicSage
pip install GEOparse --break-system-packages
python3 -c "
import GEOparse
gse = GEOparse.get_GEO('GSE166635', destdir='data/benchmark/', silent=False)
print('Done:', list(gse.gsms.keys())[:5])
"
```
Expected output: 10x MEX files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)
in data/benchmark/GSE166635/

## Step 2 — Implement data ingestion
File to create: pipeline/modules/qc/ingest.py
Should auto-detect format and return AnnData:
- 10x MEX directory → sc.read_10x_mtx()
- .h5 file → sc.read_10x_h5()
- .h5ad file → sc.read_h5ad()

## Known Issues From Last Session
- Docker images still not built locally (intentional — too heavy, do manually when needed)
- `python` must be called as `python3` in this WSL2 setup
- CI Docker build only runs on push to main branch (not dev)
- Dockerfile timezone fix (DEBIAN_FRONTEND=noninteractive) applied to both images

## Files Modified Last Session
- docker/Dockerfile.python   ← added ENV DEBIAN_FRONTEND + tzdata fix
- docker/Dockerfile.r        ← added ENV DEBIAN_FRONTEND + tzdata fix
- .github/workflows/ci.yml   ← Docker build main-only, removed pip cache
- .gitignore                 ← added !**/.gitkeep
- reports/templates/.gitkeep ← added
- reports/slides/.gitkeep    ← added
- data/benchmark/.gitkeep    ← added (force-added, overrides .gitignore)
- data/references/.gitkeep   ← added (force-added, overrides .gitignore)
- docs/.gitkeep              ← added

## Verify Last Session Works
```bash
cd ~/OmicSage
python3 -m pytest tests/test_phase0_structure.py -v
# Expected: 52 passed
```

## Relevant Context
- Benchmark paper: Wang et al. 2025 HCC, DOI: 10.1038/s41698-025-00952-3
- GEO accession: GSE166635 (scRNA-seq, hepatocellular carcinoma, multiple samples)
- Phase 1 Python stack: scanpy, scrublet, harmonypy, scvi-tools, rpy2 (for scran/SoupX)
- Start from count matrices — upstream alignment is out of scope
- Do NOT start Phase 2 until Phase 1 milestone confirmed: reproduce Wang et al. findings
