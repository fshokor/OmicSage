## Session Context
Date: 2026-05-26
Phase: Phase 4 — CITE-seq Module (Session 7 of 7)
Last thing we completed: cite_integration.py — MOFA+ and totalVI integration
                          50 tests written and passing (855 total, 1 skipped)

Files added last session:
  - pipeline/modules/cite/cite_integration.py  — MOFA+ and totalVI
  - tests/test_cite_integration.py             — 50 tests, all passing

## Modality Roadmap (corrected — carry forward every session)
  Phase 4 — CITE-seq Module      ← IN PROGRESS (final session)
  Phase 5 — Spatial Module       ← Visium, MERFISH, Xenium
  Phase 6 — Multiome Integration ← RNA + ATAC jointly

scATAC-seq as a standalone modality is deferred — covered in Phase 6.

## Revised Phase 4 Session Order (7 sessions total)
  Session 1 — adt_normalize.py      ✓ done (392 tests, 1 skipped)
  Session 2 — adt_doublets.py       ✓ done (457 total, 1 skipped)
  Session 3 — adt_reduce.py         ✓ done (516 total, 1 skipped)
  Session 4 — adt_harmony.py        ✓ done (736 total, 1 skipped)
  Session 5 — adt_annotate.py       ✓ done (805 total, 1 skipped)
  Session 6 — cite_integration.py   ✓ done (855 total, 1 skipped)
  Session 7 — cite_pipeline.py + reports  ← TODAY

## Today's Goal

### Session 7 — cite_pipeline.py + HTML report

  A. cite_pipeline.py
     - End-to-end orchestrator for the full CITE-seq pipeline
     - One function: run_cite_pipeline(mdata, config=None) → MuData
     - Calls in order:
         normalize_adt → detect_adt_doublets → reduce_adt →
         run_harmony_adt → annotate_adt → run_mofa (default integration)
     - config dict controls optional steps and parameters:
         config["batch_key"]        — default "donor"
         config["integration"]      — "mofa" | "totalvi", default "mofa"
         config["n_factors"]        — MOFA+ factors, default 15
         config["max_epochs"]       — totalVI epochs, default 400
         config["annotation_map"]   — passed to annotate_adt, default None
         config["filter_doublets"]  — passed to detect_adt_doublets, default False
     - Returns fully processed MuData with all embeddings populated
     - Provenance: mdata.uns["omicsage_cite_pipeline"] — dict of per-step
       metrics and timestamps

  B. reports/templates/cite_report.py
     - HTML report generator covering all Phase 4 steps in one document
     - Input: fully processed MuData (output of run_cite_pipeline)
     - Output: reports/cite_report_{timestamp}.html
     - Sections:
         1. Dataset summary (n_cells, n_proteins, n_donors)
         2. ADT normalization (CLR distribution plot)
         3. Doublet detection (score histogram, n_doublets)
         4. Dimensionality reduction (ADT UMAP before/after Harmony)
         5. Annotation (Leiden cluster sizes, cell type map if provided)
         6. Integration (MOFA/totalVI UMAP coloured by batch + cluster)
     - Uses matplotlib/scanpy plotting only — no Quarto dependency
     - One function: generate_cite_report(mdata, output_dir) → str (path)

## WNN — Deferred
  WNN via muon.pp.neighbors is deferred to a future standalone session.
  Root cause: muon's Jaccard+Euclidean NNDescent metric hangs on small
  fixtures due to pynndescent behaviour in CI/sandbox environments.
  The API contract is already documented in cite_integration.py.
  When resuming: write a @pytest.mark.slow WNN test class with n_cells=2000
  fixture and run separately from the main suite.
  Expected output keys (when implemented):
    mdata.obsm["X_umap_wnn"]         — WNN UMAP
    mdata.obsp["wnn_connectivities"] — WNN graph connectivities
    mdata.obsp["wnn_distances"]      — WNN graph distances

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
- MOFA+: muon namespaces batch_key as "rna:donor"/"adt:donor" in mdata.obs —
  fixed in cite_integration.py by pushing batch_key to mdata.obs before
  calling mu.tl.mofa. cite_pipeline.py must do the same if calling run_mofa.

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

### adt_harmony.py
  run_harmony_adt(mdata, batch_key, n_pcs=20, n_neighbors=15,
                  random_state=0, inplace=False)
    → (AnnData, dict)
  Input:  mdata["adt"].obsm["X_pca_adt"], mdata["adt"].obs[batch_key]
  Output on mdata["adt"]:
    obsm["X_pca_harmony_adt"] — Harmony-corrected ADT PCA
    obsm["X_umap_adt"]        — UMAP from harmony embedding (overwrites)
  Provenance: adata.uns["omicsage_adt_harmony"]

### adt_annotate.py
  annotate_adt(mdata, annotation_map=None, resolution=0.1,
               n_iterations=2, random_state=0, inplace=False)
    → (AnnData, dict)
  Input:  mdata["adt"].obsp["connectivities"] (from adt_harmony.py)
  Output on mdata["adt"]:
    obs["leiden"]       — Leiden cluster IDs (always)
    obs["adt_celltype"] — cell type labels (only if annotation_map provided)
  Provenance: adata.uns["omicsage_adt_annotate"]

### cite_integration.py
  run_mofa(mdata, batch_key, n_factors=15, use_layer=None,
           random_state=0, inplace=False)
    → (MuData, dict)
  Input:  mdata["rna"].X (log1p), mdata["adt"].X (CLR), obs[batch_key]
  Output on mdata:
    obsm["X_mofa"]      — MOFA+ latent factors
    obsm["X_umap_mofa"] — UMAP from MOFA+ embedding
  NOTE: pushes batch_key to mdata.obs before calling mu.tl.mofa
  Provenance: mdata.uns["omicsage_mofa"]

  run_totalvi(mdata, batch_key, max_epochs=400, random_state=0,
              inplace=False)
    → (MuData, dict)
  Input:  mdata["rna"].layers["counts"], mdata["adt"].layers["counts"],
          mdata["rna"].obs[batch_key]
  Output on mdata:
    obsm["X_totalVI"]      — totalVI latent representation
    obsm["X_umap_totalVI"] — UMAP from totalVI embedding
  Provenance: mdata.uns["omicsage_totalVI"]

## Embedding Key Naming Convention (enforce across all sessions)
  RNA space (on mdata["rna"]):
    obsm["X_pca"]              — RNA PCA           (reduce.py)
    obsm["X_pca_harmony"]      — RNA Harmony PCA   (future RNA integration)
    obsm["X_umap"]             — RNA UMAP          (reduce.py)

  ADT space (on mdata["adt"]):
    obsm["X_pca_adt"]          — ADT PCA           (adt_reduce.py)
    obsm["X_pca_harmony_adt"]  — ADT Harmony PCA   (adt_harmony.py)
    obsm["X_umap_adt"]         — ADT UMAP          (adt_harmony.py — overwritten)

  Joint space (on mdata):
    obsm["X_mofa"]             — MOFA+ factors     (cite_integration.py)
    obsm["X_umap_mofa"]        — MOFA+ UMAP        (cite_integration.py)
    obsm["X_totalVI"]          — totalVI latent    (cite_integration.py)
    obsm["X_umap_totalVI"]     — totalVI UMAP      (cite_integration.py)
    obsm["X_umap_wnn"]         — WNN UMAP          (deferred)

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/ -v
# Expected: 855 passed, 1 skipped
# If count is wrong, stop and investigate before writing any new code
```

## Relevant Context
- Benchmark dataset: GSE194122 BMMC CITE-seq (NeurIPS 2021)
  → BioLegend TotalSeq B Universal Human Panel (~140 antibodies)
  → batch key for Harmony and integration is "donor"
- cite_pipeline.py should default to MOFA+ (faster, no GPU needed)
  and allow switching to totalVI via config["integration"] = "totalvi"
- The HTML report uses matplotlib/scanpy only — no Quarto, no Jinja2
  templates, no external dependencies beyond what is already in the
  omicsage conda environment
- mofapy2 and scvi-tools must be added to environment.yml / CI deps
  if not already present (check before the session starts)
