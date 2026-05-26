## Session Context
Date: 2026-05-26
Phase: Phase 4 — CITE-seq Module (Session 2 of 4)
Last thing we completed: adt_normalize.py — CLR normalization module
                          37 tests written and passing (392 total, 1 skipped)

Files added this session:
  - pipeline/modules/cite/adt_normalize.py   — CLR normalization for ADT
  - pipeline/modules/cite/__init__.py        — package init (create if missing)
  - tests/test_adt_normalize.py              — 37 tests, all passing

## Modality Roadmap (corrected — carry forward every session)
  Phase 4 — CITE-seq Module      ← IN PROGRESS
  Phase 5 — Spatial Module       ← Visium, MERFISH, Xenium
  Phase 6 — Multiome Integration ← RNA + ATAC jointly

scATAC-seq as a standalone modality is deferred — covered in Phase 6.

## Today's Goal
Build the **ADT doublet detection module** — Phase 4, Session 2.

One clearly scoped deliverable:
  `pipeline/modules/cite/adt_doublets.py`

ADT doublet detection catches doublets that RNA-only Scrublet misses,
using the protein signal. It must run on CLR-normalized ADT values
(output of adt_normalize.py), not raw counts.

This module takes `mdata["adt"].layers["adt_clr"]` and performs:
  1. Doublet detection on CLR-normalized ADT using scrublet or
     an ADT-appropriate method (check reference for recommendation)
  2. Adds obs columns to mdata["adt"]:
       obs["adt_doublet_score"]     — continuous score (0–1)
       obs["adt_predicted_doublet"] — boolean flag
  3. Optionally filters doublets from both mdata["adt"] AND mdata["rna"]
     (cross-modal filtering — a cell doublet in protein space should be
     removed from RNA too)
  4. Returns updated AnnData + metrics dict

Do NOT start dimensionality reduction this session — doublets only.

## Correct Phase 4 session order
  Session 1 — adt_normalize.py      ✓ done (392 tests, 1 skipped)
  Session 2 — adt_doublets.py       ← today
  Session 3 — adt_reduce.py         ← PCA + neighbor graph + UMAP
  Session 4 — wnn.py                ← RNA + protein joint embedding

## Reference for this session
ADT doublet detection:
  https://www.sc-best-practices.org/surface_protein/doublet_detection.html

ADT dimensionality reduction (Session 3):
  https://www.sc-best-practices.org/surface_protein/dimensionality_reduction.html

ADT batch correction (relevant for Session 4 WNN):
  https://www.sc-best-practices.org/surface_protein/batch_correction.html

Full surface protein chapter index (bookmark):
  https://www.sc-best-practices.org/surface_protein/normalization.html

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- rpy2 NOT used — R scripts called via subprocess or Nextflow, not rpy2
- `seurat_v3` HVG flavor numerically unstable on small fixtures — use
  `flavor='seurat'` in small-fixture tests
- obs['cell_type_vote'] is the consensus column naming convention
- `cell_type_confidence` = weighted fraction of methods agreeing (0.0–1.0)
  NOT AI self-reported confidence
- muon CLR FutureWarning on __version__ is internal to muon — not our code,
  suppress with warnings.catch_warnings() as already done in adt_normalize.py

## adt_normalize.py API (carry forward — adt_reduce.py depends on this)
  normalize_adt(adata, clr_axis=0, dsb_empty_adata=None,
                isotype_controls=None, inplace=False)
    → (AnnData, dict)

  Input:  mdata["adt"].X  — raw integer ADT counts
  Output layers:
    layers["counts"]   — raw counts (preserved)
    layers["adt_clr"]  — CLR-normalized values
    .X                 — CLR-normalized values (same as adt_clr)
  Provenance: adata.uns["omicsage_adt_normalize"]
  clr_axis=0 default  (per-protein across cells, muon default)

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 392 passed, 1 skipped
# If count is wrong, stop and investigate before writing any new code
```

## Relevant Context
- Benchmark dataset: GSE194122 BMMC CITE-seq (NeurIPS 2021)
  → BioLegend TotalSeq B Universal Human Panel (~140 antibodies)
  → mdata["adt"] post-normalize has layers["adt_clr"] with CLR values
- RNA doublets already caught by Scrublet in qc.py (obs["predicted_doublet"])
  ADT doublet detection is complementary — catches doublets RNA missed
- Cross-modal filtering rule: a doublet flagged in ADT space must also be
  removed from mdata["rna"] — they share barcodes, so both must stay in sync
- Phase 4 milestone: WNN UMAP showing RNA + protein integrated embedding
  (Session 4 goal — do not start until Sessions 2 and 3 are complete)
- UMAP key naming convention (establish now, enforce across all sessions):
    obsm["X_umap"]      → RNA UMAP  (already exists from reduce.py)
    obsm["X_umap_adt"]  → ADT UMAP  (Session 3)
    obsm["X_umap_wnn"]  → WNN UMAP  (Session 4)
