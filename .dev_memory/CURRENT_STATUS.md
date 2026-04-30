# OmicSage — Current Status
Last updated: 2026-04-30

## Phase
Phase 0 — Foundation: COMPLETE
Next: Phase 1 — Core scRNA Pipeline

## Built and Verified
- Repo structure: all dirs and files per spec
- config/schema.yaml: complete, all parameters documented
- pipeline/main.nf: modality router
- pipeline/workflows/*.nf: stubs for all 4 modalities
- nextflow.config: profiles local/docker/singularity/slurm/ci/test
- docker/Dockerfile.python: mambaforge + scanpy + scvi + BioChatter + Quarto
- docker/Dockerfile.r: Bioconductor + Seurat v5 + Signac + ArchR + SingleR
- docker-compose.yml: python, streamlit, r, nextflow services
- environment.yml: conda env
- .github/workflows/ci.yml: lint + pytest + docker build
- .github/workflows/docker-publish.yml: GHCR on version tags
- ai/biochatter_client.py: typed stub (Phase 3 pending)
- cli/omicsage.py: create-project, run, list, compare
- ui/app.py: Streamlit landing page stub
- tests/test_phase0_structure.py: 54 tests, all passing
- .dev_memory/: all 5 files initialized
- README.md, LICENSE, .gitignore

## Not Yet Done
- Docker images not built locally
- GSE166635 not downloaded
- First git push pending
- DVC not initialized
