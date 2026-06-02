## Session Context
Date: 2026-06-01
Phase: 4 — CITE-seq pipeline
Last thing we completed: Full audit of cite pipeline against sc-best-practices tutorial.
  Identified systematic issues in normalization, doublet detection, and annotation
  that need to be fixed to align with the tutorial results.
Last files modified:
  - pipeline/modules/cite/adt_annotate.py
  - reports/templates/cite/cite_annotate_report.py

## Today's Goal
Align every step of the CITE-seq pipeline and its reports with the
sc-best-practices tutorial (chapters: QC, normalization, doublet detection,
dimensionality reduction, batch correction, annotation).

Reference URLs (fetch and read before implementing each section):
  https://www.sc-best-practices.org/surface_protein/quality_control.html
  https://www.sc-best-practices.org/surface_protein/normalization.html
  https://www.sc-best-practices.org/surface_protein/doublet_detection.html
  https://www.sc-best-practices.org/surface_protein/dimensionality_reduction.html
  https://www.sc-best-practices.org/surface_protein/batch_correction.html
  https://www.sc-best-practices.org/surface_protein/annotation.html

One deliverable per step — do not implement multiple steps in one session.
Agree on the step order with Fatima before writing any code.

## What Needs To Change — Per Step

---

### STEP 1 — Normalization (HIGHEST PRIORITY — unblocks everything else)
File: pipeline/modules/cite/adt_normalize.py
Report: reports/templates/cite/cite_normalize_report.py
Checkpoint: data/processed/GSE194122/cite_01_normalized.h5ad

PROBLEM:
  Currently only CLR normalization is applied (DSB applied: No in report).
  CLR does not remove ambient protein background. All proteins show a
  floor near 0 with long right tails — no bimodal distributions.
  The tutorial explicitly recommends DSB as the primary normalization
  because it uses empty droplets to model and subtract ambient signal,
  producing interpretable background-corrected values where
  negative = background, positive = truly expressed.
  This affects every downstream step: doublet detection, PCA, annotation scoring.

WHAT TO CHANGE:
  1. Implement DSB normalization in adt_normalize.py:
     - Use the `dsb` Python package (pip install dsb) or implement via
       the DSB algorithm: background estimation from empty droplets +
       protein-specific denoising
     - Empty droplet barcodes for GSE194122 are available — the NeurIPS
       2021 dataset was specifically designed with empty droplets included.
       They are typically stored alongside the raw count matrix.
     - Store DSB-normalized values in adata.layers["adt_dsb"]
     - Keep CLR in adata.layers["adt_clr"] as a fallback
     - Set adata.X to the DSB layer after normalization (pipeline default)
  2. Add `dsb` to requirements-ci.txt and environment.yml
  3. Update the report:
     - Change stat card "DSB applied: No" → "DSB applied: Yes"
     - Add before/after DSB violin plots showing bimodal distributions
       (the key quality check from the tutorial)
     - Add ambient background level per protein (from DSB model params)
  4. Update config schema to add dsb_empty_droplets_path parameter

TUTORIAL REFERENCE RESULT:
  After DSB, each protein shows a clear bimodal distribution:
  left peak = background-expressing cells, right peak = truly expressing cells.
  This is Figure 1 in the normalization chapter.

---

### STEP 2 — Doublet Detection
File: pipeline/modules/cite/adt_doublets.py
Report: reports/templates/cite/cite_doublets_report.py
Checkpoint: data/processed/GSE194122/cite_02_doublets.h5ad

PROBLEM:
  Currently only 1 doublet flagged across 21,778 cells (0.0%).
  The tutorial expects ~4-5% doublets. The histogram shows all cells
  piled at score ≈ 0 — the doublet scoring is not producing a real
  distribution. Root cause: doublet scoring is running on CLR-normalized
  data where the background floor makes lineage marker co-expression
  invisible. With DSB normalization the bimodal distributions make it
  clear when a cell co-expresses two exclusive lineage markers.

WHAT TO CHANGE:
  1. Re-run doublet detection after DSB normalization is in place
     (implement Step 1 first — this step depends on it)
  2. The tutorial uses co-expression of lineage-exclusive marker pairs
     to score doublets: CD3/CD19 (T+B), CD3/CD14 (T+Mono), CD19/CD14 (B+Mono)
     A cell expressing both markers of a pair above the positive peak
     (as defined by DSB bimodal threshold) is a doublet candidate.
  3. The threshold should be set on the DSB-normalized values, not CLR.
     DSB values have a natural threshold near 0 (negative = background).
     Cells with DSB > 0 for both markers in a pair are doublets.
  4. Expected output: histogram with a real right tail, ~4-5% flagged.
  5. Update the report to show:
     - Real doublet score distribution (not a spike at 0)
     - Per-donor doublet rate (some donors may have higher rates)
     - Scatter plots of exclusive marker pairs with doublets highlighted
       (CD3 vs CD19, CD3 vs CD14) — this is the key figure in the tutorial

TUTORIAL REFERENCE RESULT:
  Doublet score histogram shows a clear bimodal or right-skewed distribution.
  Scatter of CD3 vs CD19 shows a main diagonal cluster of doublets
  co-expressing both markers above background.

---

### STEP 3 — Dimensionality Reduction
File: pipeline/modules/cite/adt_reduce.py
Report: reports/templates/cite/cite_reduce_report.py
Checkpoint: data/processed/GSE194122/cite_03_reduced.h5ad

PROBLEM:
  Minor — currently 20 PCs used, tutorial uses 18. Not critical.
  Main issue: PCA is running on CLR values — should run on DSB after Step 1.
  No PCA loadings figure showing which proteins drive each PC.

WHAT TO CHANGE:
  1. After DSB is in place, confirm PCA runs on layers["adt_dsb"] not CLR.
  2. Add PCA loadings heatmap to the report:
     - Top 10 proteins contributing to PC1–PC4
     - Confirms PCs are biologically meaningful
     - (PC1 should separate T vs B, PC2 myeloid, etc.)
  3. Add elbow plot with a vertical line at the selected n_components
     so the threshold decision is visible.
  4. Tutorial uses 18 PCs — consider changing default from 20 to 18,
     or making this data-driven (e.g. 80% variance explained).

TUTORIAL REFERENCE RESULT:
  PC1 separates lymphoid from myeloid. Loadings show CD3/CD4/CD8 on one
  side, CD14/CD16/HLA-DR on the other. Elbow at ~18 PCs.

---

### STEP 4 — Batch Correction
File: pipeline/modules/cite/adt_harmony.py
Report: reports/templates/cite/cite_harmony_report.py
Checkpoint: data/processed/GSE194122/cite_04_harmony.h5ad

PROBLEM:
  Harmony is working correctly — pre/post UMAP shows good batch mixing
  with biological structure preserved. This step aligns with the tutorial.
  Minor issues only.

WHAT TO CHANGE:
  1. Add quantitative batch mixing metrics to the report:
     - iLISI (integration LISI): measures batch mixing — higher is better.
       Computable from the post-Harmony UMAP embedding.
     - cLISI (cell-type LISI): measures cell type separation — higher is better.
       Requires adt_celltype to be present (run after annotation).
     - Show as a summary table: iLISI before / iLISI after / improvement %
  2. Add per-batch ADT library size violin plot:
     - Shows whether batches differ in library size before correction
     - Justifies why Harmony was needed
  These are nice-to-have additions — implement only if time allows after
  higher priority steps are done.

TUTORIAL REFERENCE RESULT:
  Post-Harmony UMAP shows all 12 donors intermixed within each cell type
  cluster. Quantitative LISI scores confirm mixing.

---

### STEP 5 — Annotation (SECOND PRIORITY)
File: pipeline/modules/cite/adt_annotate.py
Report: reports/templates/cite/cite_annotate_report.py
Checkpoint: data/processed/GSE194122/cite_05_annotated.h5ad

PROBLEM (three separate issues):

  A. BMMC_MARKER_PANEL is missing erythroid markers.
     Clusters 0 and 1 (largest clusters, ~36% of cells combined) show
     CD71, CD36, CD88 as top markers in the dotplot — these are erythroid
     markers. The panel has no "Erythroid" entry so these clusters get
     mislabelled. This is the biggest annotation error.

  B. Scoring assigns wrong labels to T cell clusters.
     Cluster 4 (4841 cells, 22%) is labelled "Treg" but shows CD3, CD4,
     CD5, CD45RA, CD45RO — this is clearly CD4 T. Treg and CD4 T share
     most markers; the fold-change scoring cannot discriminate them because
     the key discriminating markers (CD25-high, CD127-low for Treg) are not
     weighted differently from shared markers.
     Cluster 3 (3120 cells, 14%) is labelled "Platelet" but the dotplot
     shows CD11c, CD172a, CD41 — this is more consistent with monocytes
     or a mixed population.

  C. CD4 T is completely absent from the annotation output despite being
     the expected dominant population in BMMC.

WHAT TO CHANGE:

  1. Add Erythroid to BMMC_MARKER_PANEL immediately:
     "Erythroid": ["CD71", "CD36", "CD235a", "CD88"]
     Also add "HSPC" (haematopoietic stem/progenitor):
     "HSPC": ["CD34", "CD38", "CD90", "CD117", "CD133", "CD135"]
     Note: CD235a (Glycophorin A) may appear as "CD235a" or "GYPA"
     depending on the panel — check var_names.

  2. Use manual annotation_map for this run (bypass scoring):
     Based on the dotplot (rank_genes_groups on leiden clusters):
       "0": "Erythroid"      # CD71, CD36, CD88 — largest cluster 29%
       "1": "Erythroid"      # CD71, CD36, CD33 — second erythroid cluster
       "2": "B"              # CD19, CD72, CD9
       "3": "CD14 Mono"      # CD11c, CD172a — review after DSB
       "4": "CD4 T"          # CD5, CD4, CD45RA
       "5": "CD8 T"          # CD2, CD8, CD3
       "6": "Plasma"         # CD38, CD63, CD54
       "7": "NK"             # CD16, CD45RA, CD56
       "8": "DC"             # CD123, CD162, CD304
     Add this as the default annotation_map in the config YAML.
     NOTE: Verify cluster 3 (CD14 Mono vs Platelet) by checking
     CD41 expression — if CD41-high it is Platelet, CD14-high = Mono.

  3. After DSB is in place, re-run scoring with updated panel.
     DSB values will make the bimodal thresholds clearer and the
     fold-change scoring more accurate.

  4. No changes needed to report structure — the dotplot, marker UMAP
     grid, and cell type UMAP are all correct. Just fix the data.

TUTORIAL REFERENCE RESULT:
  Tutorial annotation for GSE194122 BMMC shows:
  CD4 T, CD8 T, NK, B, Plasma, CD14 Mono, CD16 Mono, DC, pDC,
  Erythroid progenitor, HSPC — 11 broad populations.
  Erythroid progenitors are the largest cluster in BMMC (~30%).

---

### STEP 6 — Integration (MOFA+ / totalVI)
File: pipeline/modules/cite/cite_integration.py
Report: reports/templates/cite/cite_integration_report.py
Checkpoint: data/processed/GSE194122/cite_06_integrated.h5ad

PROBLEM:
  Integration appears to run but the report does not show biological
  validation figures. MOFA+ variance decomposition is not reported.

WHAT TO CHANGE:
  1. Add MOFA+ variance decomposition bar chart to the report:
     - Per-factor variance explained broken down by modality (RNA % vs ADT %)
     - Shows which factors are driven by which modality
     - Key biological interpretation figure from the tutorial
  2. Add MOFA+ top weights heatmap:
     - Top 5 RNA genes + top 5 ADT proteins per factor (top 4 factors)
     - Two-panel heatmap (RNA genes | ADT proteins)
  3. Add totalVI latent space UMAP (if totalVI was run):
     - Coloured by cell type and by batch
     - Separate from the Harmony ADT UMAP
  These are deferred until Steps 1–5 are correct, since the integration
  quality depends on clean normalization and annotation.

TUTORIAL REFERENCE RESULT:
  MOFA+ Factor 1 is driven by ADT (lymphoid vs myeloid surface markers).
  Factor 2 separates by RNA (transcriptional programs within T cells).
  Integration UMAP shows cleaner separation than ADT-only UMAP.

---

## Known Issues From Last Session

1. BMMC_MARKER_PANEL missing Erythroid — clusters 0+1 mislabelled
2. DSB not implemented — CLR only, affects doublet detection and annotation
3. Only 1 doublet flagged — doublet detection broken without DSB
4. run_cite_pipeline.py was not forwarding `preset` param to annotate_adt()
   — fixed by adding params.get("preset") to the annotate_adt() call
5. adt_annotate.py rank_genes_groups runs on "leiden" (correct — diagnostic)
   and annotation scoring uses fold-change mean_in - mean_out (correct logic,
   but still misannotates without DSB and without Erythroid in panel)

## Files Modified Last Session
- pipeline/modules/cite/adt_annotate.py
  (preset system, BMMC_MARKER_PANEL, fold-change scoring, step ordering)
- reports/templates/cite/cite_annotate_report.py
  (scanpy import fix, removed concordance plot, added cell type UMAP,
   dotplot always uses leiden groupby)

## Verify Last Session Works
```bash
conda activate omicsage
python -m pytest tests/test_adt_annotate.py -v --tb=short
```
Expected: all tests passing (check count matches last known baseline).

## Recommended Session Order
Session A (today): Step 1 — DSB normalization (highest impact, unblocks all)
Session B (next):  Step 5 — Fix annotation map + Erythroid in panel
Session C:         Step 2 — Re-run doublets on DSB-normalized data
Session D:         Steps 3+4 — Minor additions (loadings heatmap, LISI)
Session E:         Step 6 — MOFA+ variance decomposition in report

## Relevant Context
- Dataset: GSE194122 NeurIPS 2021 BMMC CITE-seq
  134 proteins, BioLegend TotalSeq-A Human Universal Cocktail V1.0
  12 batches (s1d1 through s4d9), 21,778 cells after RNA QC
- Empty droplet barcodes needed for DSB — check how the data was ingested.
  In the NeurIPS dataset the raw (unfiltered) matrix contains empty droplets.
  Path likely: data/raw/GSE194122/ — look for barcodes_unfiltered.tsv or
  the unfiltered_feature_bc_matrix folder from cellranger.
- ADT data path: data/processed/GSE194122/01_qc_adt.h5ad
- RNA annotated path: data/processed/GSE194122/05_annotated.h5ad
- Batch key: "batch" throughout
- Pipeline convention: inplace=False, return (AnnData, dict),
  provenance in uns["omicsage_<module>"]
- Always run python -m pytest (not pytest) to use conda env
- Always verify baseline test count before writing new code
