## Session Context
Date: 2026-05-26
Phase: Phase 4 — CITE-seq Module (Session 3 of 4)
Last thing we completed: adt_doublets.py — ADT doublet detection module
                          65 tests written and passing (457 total, 1 skipped)

Files added this session:
  - pipeline/modules/cite/adt_doublets.py  — marker-pair doublet detection
  - tests/test_adt_doublets.py             — 65 tests, all passing

## Modality Roadmap (corrected — carry forward every session)
  Phase 4 — CITE-seq Module      ← IN PROGRESS
  Phase 5 — Spatial Module       ← Visium, MERFISH, Xenium
  Phase 6 — Multiome Integration ← RNA + ATAC jointly

scATAC-seq as a standalone modality is deferred — covered in Phase 6.

## Today's Goal
Build the **ADT dimensionality reduction module** — Phase 4, Session 3.

One clearly scoped deliverable:
  `pipeline/modules/cite/adt_reduce.py`

This module takes `mdata["adt"].layers["adt_clr"]` (CLR-normalized ADT,
output of adt_normalize.py) and performs:
  1. PCA on CLR-normalized ADT values
  2. Neighbor graph on ADT PCA embedding
  3. UMAP (required — always computed)
  4. Writes results to:
       obsm["X_pca_adt"]   — ADT PCA embedding
       obsm["X_umap_adt"]  — ADT UMAP (UMAP key naming convention)
  5. Returns updated AnnData + metrics dict

Do NOT start WNN this session — dimensionality reduction only.

## Correct Phase 4 session order
  Session 1 — adt_normalize.py      ✓ done (392 tests, 1 skipped)
  Session 2 — adt_doublets.py       ✓ done (457 total, 1 skipped)
  Session 3 — adt_reduce.py         ← today
  Session 4 — wnn.py                ← RNA + protein joint embedding

## Reference for this session
ADT dimensionality reduction:
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
- _resolve_marker is a closure inside detect_adt_doublets — do not test it
  as a module-level symbol; test marker resolution via the public API instead

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

## adt_doublets.py API (carry forward — wnn.py depends on clean cells)
  detect_adt_doublets(mdata, marker_pairs=None, threshold=2.5,
                      filter_doublets=False, inplace=False)
    → (MuData, dict)

  Input:  mdata["adt"].layers["adt_clr"]
  Default marker pairs: [("CD3", "CD19"), ("CD3", "CD14")]
  Marker matching: case-insensitive prefix match on var_names
  Output obs columns on mdata["adt"]:
    obs["adt_doublet_score"]      — float, fraction of pairs co-expressed (0–1)
    obs["adt_predicted_doublet"]  — bool
  Cross-modal: filter_doublets=True removes cells from both adt AND rna
  Provenance: adata.uns["omicsage_adt_doublets"]
  Threshold: strict greater-than (value exactly at threshold is NOT flagged)

## UMAP key naming convention (enforce across all sessions)
  obsm["X_umap"]      → RNA UMAP  (exists — from reduce.py)
  obsm["X_umap_adt"]  → ADT UMAP  (Session 3 — today)
  obsm["X_umap_wnn"]  → WNN UMAP  (Session 4)

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 457 passed, 1 skipped
# If count is wrong, stop and investigate before writing any new code
```

## Relevant Context
- Benchmark dataset: GSE194122 BMMC CITE-seq (NeurIPS 2021)
  → BioLegend TotalSeq B Universal Human Panel (~140 antibodies)
  → mdata["adt"] post-normalize has layers["adt_clr"] with CLR values
  → ~140 antibodies = low-dimensional — PCA on ADT typically uses fewer
    components than RNA (15–30 is common, not 50)
- Check the reference page before implementing — it may recommend a specific
  number of PCA components or neighbor graph parameters for ADT data
- RNA PCA already exists in obsm["X_pca"] from reduce.py — ADT PCA goes
  into obsm["X_pca_adt"] to avoid collision
- Phase 4 milestone: WNN UMAP showing RNA + protein integrated embedding
  (Session 4 goal — do not start until Sessions 2 and 3 are complete)
