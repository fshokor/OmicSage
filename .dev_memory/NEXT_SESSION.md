## Session Context
Date: 2026-05-26
Phase: Phase 4 — CITE-seq Module (Session 4 of 6)
Last thing we completed: adt_harmony.py — Abatch correction
                          66 tests written and passing (736 total, 1 skipped)

Files added last session:
  - pipeline/modules/cite/adt_harmony.py  
  - tests/test_adt_harmony.py             — 66 tests, all passing

## Modality Roadmap (corrected — carry forward every session)
  Phase 4 — CITE-seq Module      ← IN PROGRESS
  Phase 5 — Spatial Module       ← Visium, MERFISH, Xenium
  Phase 6 — Multiome Integration ← RNA + ATAC jointly

scATAC-seq as a standalone modality is deferred — covered in Phase 6.

## Revised Phase 4 Session Order (6 sessions total)
  Session 1 — adt_normalize.py      ✓ done (392 tests, 1 skipped)
  Session 2 — adt_doublets.py       ✓ done (457 total, 1 skipped)
  Session 3 — adt_reduce.py         ✓ done (516 total, 1 skipped)
  Session 4 — adt_harmony.py        ✓ done (736 total, 1 skipped)
  Session 5 — adt_annotate.py       ← TODAY Leiden clustering + marker annotation
  Session 6 — wnn.py + cite_pipeline.py + reports  ← Phase 4 milestone

## Today's Goal
Build the **ADT Annotation module** — Phase 4, Session 5.

One clearly scoped deliverable:
  `pipeline/modules/cite/adt_annotate.py`


## Session 5 Preview — adt_annotate.py
  Reference: https://www.sc-best-practices.org/surface_protein/annotation.html
  Steps (from notebook source, confirmed):
    1. Leiden clustering at low resolution (resolution=0.1, flavor="igraph",
       n_iterations=2, directed=False) on harmony-corrected neighbor graph
    2. rank_genes_groups on leiden clusters (for marker dotplot)
    3. Manual cluster → cell type mapping stored in obs["adt_celltype"]
       (OmicSage convention: "adt_celltype", not "celltype", to avoid
       collision with RNA annotation obs["cell_type_vote"])
    4. Returns updated AnnData + metrics dict (n_clusters, cluster_sizes,
       annotation_map)
  Annotation is marker-based (manual map) — no LLM this session.

## Session 6 Preview — wnn.py + cite_pipeline.py + reports
  Three deliverables in one session (they are tightly coupled):

  A. wnn.py
     - WNN joint embedding: RNA (X_pca) + ADT (X_pca_harmony_adt)
     - Uses muon.pp.neighbors(mdata, ...) — multi-modal WNN
     - Writes mdata.obsm["X_umap_wnn"]
     - Writes mdata.obsp["wnn_connectivities"], mdata.obsp["wnn_distances"]

  B. cite_pipeline.py
     - End-to-end orchestrator script for the full CITE-seq pipeline
     - Calls: normalize_adt → detect_adt_doublets → reduce_adt →
              run_harmony_adt → annotate_adt → run_wnn
     - Accepts MuData, returns fully processed MuData
     - One function: run_cite_pipeline(mdata, config=None) → MuData

  C. reports/
     - HTML report generator covering all Phase 4 steps grouped together:
         QC summary, normalization, doublet detection, dimensionality
         reduction, batch correction, annotation, WNN embedding
     - One report per pipeline run, not per module
     - Template: reports/templates/cite_report.py (or .qmd if Quarto)

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- rpy2 NOT used — R scripts called via subprocess or Nextflow, not rpy2
- `seurat_v3` HVG flavor numerically unstable on small fixtures — use
  `flavor='seurat'` in small-fixture tests
- obs['cell_type_vote'] is the consensus column naming convention for RNA
- `cell_type_confidence` = weighted fraction of methods agreeing (0.0–1.0)
  NOT AI self-reported confidence
- muon CLR FutureWarning on __version__ is internal to muon — not our code,
  suppress with warnings.catch_warnings() as already done in adt_normalize.py
- _resolve_marker is a closure inside detect_adt_doublets — do not test it
  as a module-level symbol; test marker resolution via the public API instead

## API Reference (carry forward every session)

### adt_normalize.py
  normalize_adt(adata, clr_axis=0, dsb_empty_adata=None,
                isotype_controls=None, inplace=False)
    → (AnnData, dict)
  Input:  mdata["adt"].X  — raw integer ADT counts
  Output layers:
    layers["counts"]   — raw counts (preserved)
    layers["adt_clr"]  — CLR-normalized values
    .X                 — CLR-normalized values (same as adt_clr)
  Provenance: adata.uns["omicsage_adt_normalize"]

### adt_doublets.py
  detect_adt_doublets(mdata, marker_pairs=None, threshold=2.5,
                      filter_doublets=False, inplace=False)
    → (MuData, dict)
  Input:  mdata["adt"].layers["adt_clr"]
  Output obs columns on mdata["adt"]:
    obs["adt_doublet_score"]      — float (0–1)
    obs["adt_predicted_doublet"]  — bool
  Threshold: strict greater-than
  Provenance: adata.uns["omicsage_adt_doublets"]

### adt_reduce.py
  reduce_adt(mdata, n_comps=50, n_pcs=20, n_neighbors=15,
             svd_solver="arpack", random_state=0, inplace=False)
    → (AnnData, dict)
  Input:  mdata["adt"].layers["adt_clr"]
  Output on mdata["adt"]:
    obsm["X_pca_adt"]   — ADT PCA (n_cells × n_comps_actual)
    obsm["X_umap_adt"]  — ADT UMAP before batch correction (n_cells × 2)
  n_pcs=20 default (sc-best-practices ch.37)
  Provenance: adata.uns["omicsage_adt_reduce"]

## Embedding Key Naming Convention (enforce across all sessions)
  RNA space (on mdata["rna"]):
    obsm["X_pca"]              — RNA PCA           (reduce.py)
    obsm["X_pca_harmony"]      — RNA Harmony PCA   (future RNA integration)
    obsm["X_umap"]             — RNA UMAP          (reduce.py)

  ADT space (on mdata["adt"]):
    obsm["X_pca_adt"]          — ADT PCA           (adt_reduce.py)
    obsm["X_pca_harmony_adt"]  — ADT Harmony PCA   (adt_harmony.py — today)
    obsm["X_umap_adt"]         — ADT UMAP          (overwritten by adt_harmony.py)

  Joint space (on mdata):
    obsm["X_umap_wnn"]         — WNN UMAP          (wnn.py — Session 6)

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 516 passed, 1 skipped
# If count is wrong, stop and investigate before writing any new code
```

## Relevant Context
- Benchmark dataset: GSE194122 BMMC CITE-seq (NeurIPS 2021)
  → BioLegend TotalSeq B Universal Human Panel (~140 antibodies)
  → batch key for Harmony is "donor" (multiple donors in GSE194122)
- Harmony operates on obsm["X_pca_adt"] and writes obsm["X_pca_harmony_adt"]
  via sc.external.pp.harmony_integrate
- After Harmony: recompute neighbors with use_rep="X_pca_harmony_adt",
  then recompute UMAP → obsm["X_umap_adt"] (overwrite uncorrected version)
- The batch_key parameter must be validated to exist in adata.obs before
  calling Harmony — raise a clear ValueError if missing
- harmonypy is already in the omicsage conda environment
