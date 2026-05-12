## Session Context
Date: 2026-05-12
Phase: 2 — Report Engine
Last thing completed: Phase 1 MILESTONE — full pipeline run on GSE166635 HCC
                      Wang et al. 2025 key findings reproduced.
                      Phase 1 is complete.
File last worked on: config/runs/GSE166635.yaml

## Today's Goal
Phase 2 kickoff — decide the report engine approach and build the first
Quarto report template (QC report).

The question to answer at session start:
  Do we replace the existing hand-coded HTML reports with Quarto,
  or do we keep HTML for interactive use and add Quarto as a separate
  PDF/pptx layer on top?

Recommended: keep HTML reports (they work, they're fast, no dependency),
add Quarto as an additional rendering layer that reads the same data and
produces PDF + PowerPoint. One command triggers both.

ONE goal — decide the approach and build ONE template (QC or normalization).
Do not build all templates in one session.

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~231 passed, 1-2 skipped
```

## Step 2 — Quarto setup check
```bash
quarto check
# If not installed: https://quarto.org/docs/get-started/
# Quarto needs to be installed separately (not via conda)
```

## Step 3 — Decide report engine approach before writing any code
Options:
  A. Quarto .qmd templates that read .h5ad files directly via Python chunks
     → most flexible, best PDF output, requires Quarto + pandoc
  B. python-pptx only — generate slides from the existing HTML figures
     → simpler, no new dependencies, slides only (no PDF)
  C. Both: Quarto for PDF, python-pptx for slides, HTML stays as-is
     → most complete, most work

Recommended: C — but build Quarto QC template first this session,
python-pptx in the next session.

## Phase 2 Plan (across multiple sessions)
Session 1 (this): Quarto QC report template
Session 2: Quarto analysis report template (clustering + annotation + DEG)
Session 3: python-pptx slide deck generator
Session 4: auto-methods text from provenance metadata
Milestone: biologist receives PDF + slides from one command

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- Use `2>&1 | tee logs/<dataset>.log` to capture all output
- run_pipeline.py must be at repo root
- obs['cell_type_vote'] is the consensus column (not obs['cell_type'])
- rpy2 NOT used — do not attempt R integration
- gseapy.enrichr returns NO Overlap column in v1.1.3
- harmonypy ho.Z_corr is (n_cells, n_pcs) — do NOT transpose
- obs[batch_key] must be cast .astype(str) before pandas operations
- pseudobulk requires layers['counts'] (raw integers) — NOT layers['logcounts']
- GSE166635 pseudobulk is DISABLED — only 2 samples

## Files Modified Last Session
- run_pipeline.py                          ← generic runner (validation bug fixed)
- pipeline/modules/qc/ingest.py            ← multi-sample MTX support
- config/schema.yaml                       ← updated
- config/runs/GSE194122.yaml               ← new
- config/runs/GSE166635.yaml               ← new
- config/runs/GSE194122_multiome.yaml      ← new
- .dev_memory/MODULE_DOCS.md               ← section 21 replaced
- .dev_memory/CURRENT_STATUS.md            ← updated
- .dev_memory/NEXT_SESSION.md              ← updated (this file)
- .dev_memory/PROGRESS.md                  ← updated

## Verify Last Session Works
```bash
cd ~/OmicSage
conda activate omicsage
python run_pipeline.py --config config/runs/GSE194122.yaml --step qc
# Should use cached 01_qc.h5ad and print: [qc] cached → ...
```

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2
New dependency needed: quarto (install separately from quarto.org)
