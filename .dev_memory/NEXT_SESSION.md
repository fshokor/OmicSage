## Session Context
Date: 2026-05-17
Phase: 1 — Core scRNA Pipeline (annotation step)
Last thing we completed: Phase 3 AI layer — all 10 sessions done (466+ tests passing).
                          Decision made: stop AI pipeline development, keep manual version.
                          The AI layer code is complete and tested but will NOT be extended further
                          for now. Manual pipeline is the primary path going forward.
File we were working on last: ai/report_reviewer.py + tests/test_groundedness.py

## Decision Logged
After completing the full AI layer (Sessions 0-10 of Phase 3), we evaluated the added
complexity vs. value at this stage of the project and decided:
  → Keep the manual pipeline as the single primary path
  → Pause all AI pipeline development
  → AI layer code stays in the repo, intact and tested, but no further extension
  → The architecture still supports ai_features: true — just not being developed further now
  → Next focus: complete the scRNA pipeline's annotation module

See .dev_memory/DECISIONS.md for full rationale.

## Today's Goal
Build the **Annotation Module** — manual cell type annotation using SingleR + marker review.

This is the final major piece of the Phase 1 scRNA pipeline before the milestone run.

## Annotation Plan — What to Build

The annotation module lives in `pipeline/modules/annotation/`.

### Step 1 — SingleR automatic annotation (R)
File: `pipeline/modules/annotation/singler_annotate.R`
- Input: clustered AnnData saved as RDS or read via zellkonverter
- Reference: `celldex::HumanPrimaryCellAtlasData()` (download on first run)
- Output:
  - Per-cell predictions (cell_type, singler_score, delta_next)
  - Per-cluster consensus label (majority vote)
  - Flag low-confidence cells: delta_next < 0.1
  - Save: `results/{project}/annotation/singler_predictions.csv`

### Step 2 — Marker review summary (Python)
File: `pipeline/modules/annotation/marker_review.py`
- Input: AnnData with clustering + rank_genes_groups results
- For each cluster: extract top 10 markers (logFC + p-val)
- Output: `results/{project}/annotation/marker_summary.csv`
  Columns: cluster_id, top_markers (comma-sep), known_cell_type (manual fill), confidence

### Step 3 — Nextflow process
File: `pipeline/modules/annotation/scrna_annotation.nf`
- Wraps singler_annotate.R as a Nextflow process
- Connects to scrna.nf workflow after CLUSTERING step
- Outputs annotated h5ad + annotation CSVs into results dir

### Step 4 — Annotated AnnData output
- Write cell_type column to adata.obs from SingleR consensus
- Write singler_score and delta_next columns
- Save as `results/{project}/annotated.h5ad`

### Step 5 — Annotation section in HTML report
- Add annotation tab to the combined HTML report (00_combined_report.html)
- Include: UMAP colored by cell_type, annotation table, low-confidence cell count

## Expected File Outputs
```
pipeline/modules/annotation/
  ├── scrna_annotation.nf
  ├── singler_annotate.R
  └── marker_review.py

results/{project}/annotation/
  ├── singler_predictions.csv
  ├── marker_summary.csv
  └── annotation_umap.png
```

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- rpy2 NOT used — R scripts are called via subprocess or Nextflow, not rpy2
- Docker images may still need building (see BLOCKERS.md B001)
- obs['cell_type_vote'] is the consensus column naming convention

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 466+ passed, 1 skipped
```

## Relevant Context
- Primary benchmark dataset: GSE166635 HCC (Wang et al. 2025)
- Annotation reference: celldex HumanPrimaryCellAtlasData (download on first run in R)
- SingleR is available in the R environment (Seurat v5 + Bioconductor stack)
- Phase milestone: reproduce Wang et al. 2025 cell type assignments from raw counts
- The combined HTML report system (Phase 2) is complete — annotation output plugs into it
