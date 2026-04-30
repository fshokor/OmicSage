# OmicSage — Blockers

## OPEN

### B001 — Docker images not built
Run:
  docker build -f docker/Dockerfile.python -t omicsage/python:latest .
  docker build -f docker/Dockerfile.r -t omicsage/r:latest .
R image: 20-40 min.

### B002 — GSE166635 not downloaded
pip install GEOparse
python -c "import GEOparse; GEOparse.get_GEO('GSE166635', destdir='data/benchmark/')"

### B003 — First git push pending
git add . && git commit -m "feat: Phase 0 complete" && git push origin main

### B004 — DVC not initialized
cd ~/OmicSage && dvc init

## RESOLVED
(none yet)

## Open Questions
Q001 — Default LLM provider: claude vs ollama? Decide before Phase 3.
Q002 — Quarto HTML vs PDF as primary? Plan: HTML primary, PDF secondary.
Q003 — pip-installable CLI? Add entry_points before Phase 2.
