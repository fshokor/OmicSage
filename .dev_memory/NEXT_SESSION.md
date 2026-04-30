## Session Context
Date: 2026-04-30
Phase: 0 — Foundation (COMPLETE) → Phase 1 — Core scRNA Pipeline
Last thing we completed: Full Phase 0 scaffold — repo structure, Docker images, CI/CD, .dev_memory system
File we were working on last: tests/test_phase0_structure.py

## Today's Goal
[CHOOSE ONE]:
A) Download benchmark dataset (GSE166635 HCC) and write the data ingestion module
   → Start: `python cli/omicsage.py create-project hcc_benchmark --modality scrna`
B) Implement the scRNA QC module (pipeline/modules/qc/scrna_qc.nf + Python script)
   → Prerequisite: benchmark data must exist at data/benchmark/

Recommended: Start with A if GSE166635 is not yet downloaded. Otherwise B.

## Known Issues From Last Session
- Docker images are defined but NOT YET BUILT locally — run these first:
  ```
  docker build -f docker/Dockerfile.python -t omicsage/python:latest .
  docker build -f docker/Dockerfile.r -t omicsage/r:latest .
  ```
  Note: R image takes 20-40 min. Start it and work on Python tasks while it builds.
- GitHub Actions CI will fail on nextflow --help step until Nextflow is confirmed installed in CI runner — check ci.yml if first push fails
- Nextflow not yet verified working locally — run: `nextflow run pipeline/main.nf --help`

## Files Modified Last Session
- config/schema.yaml           ← master config, all parameters documented
- pipeline/main.nf             ← entry point with modality router
- pipeline/workflows/scrna.nf  ← stub, imports pending Phase 1 modules
- pipeline/workflows/scatac.nf ← stub
- pipeline/workflows/spatial.nf← stub
- pipeline/workflows/integration.nf ← stub
- pipeline/modules/qc/scrna_qc.nf   ← stub process, ready for Phase 1
- nextflow.config              ← profiles: local, docker, singularity, slurm, ci, test
- docker/Dockerfile.python     ← mambaforge + scanpy + scvi + BioChatter + Quarto
- docker/Dockerfile.r          ← Bioconductor + Seurat v5 + Signac + ArchR + SingleR
- docker-compose.yml           ← python, streamlit, r, nextflow services
- environment.yml              ← conda env for local dev without Docker
- .github/workflows/ci.yml     ← lint + pytest + nextflow smoke + docker build
- .github/workflows/docker-publish.yml ← push to GHCR on version tags
- ai/biochatter_client.py      ← typed stub, Phase 3 implementation pending
- cli/omicsage.py              ← Click CLI: create-project, run, list, compare
- ui/app.py                    ← Streamlit landing page stub
- tests/test_phase0_structure.py ← smoke tests: dir structure, config, memory files
- pyproject.toml               ← pytest config
- .dev_memory/ (all 5 files)   ← initialized today

## Verify Last Session Works
```bash
# From repo root in WSL2:
pip install pytest pyyaml
pytest tests/test_phase0_structure.py -v
```
All tests should pass (structure tests check file existence).

## Relevant Context
- GitHub: https://github.com/fshokor/OmicSage
- Benchmark paper: Wang et al. 2025 HCC (DOI: 10.1038/s41698-025-00952-3)
- Benchmark dataset: GSE166635 on GEO (scRNA-seq, hepatocellular carcinoma)
- Expected Phase 1 milestone: reproduce key findings of Wang et al. from raw counts
- Phase 1 Python stack: Scanpy, scran (via rpy2), Scrublet, SoupX, Harmony, scvi-tools
- DO NOT start Phase 2 (reports) until Phase 1 milestone is confirmed working
