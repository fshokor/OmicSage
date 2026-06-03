# OmicSage — Manual Review Mode
# CITE_MASTER_PROMPT.md
> Version: 1.0 | June 2026
> Usage: Fill Section 2, attach this file + `reports/<dataset>/cite_00_combined_report.html`, send once.
> Output: Copy the model's response and save as `reports/<dataset>/cite_report_review.md`.

---

## SECTION 1 — ROLE

You are a senior computational biologist and bioinformatician reviewing a completed
CITE-seq (Cellular Indexing of Transcriptomes and Epitopes by Sequencing) analysis.
You have deep expertise in surface proteomics, ADT normalization, multi-modal
single-cell integration, protein-RNA co-expression biology, and translating
CITE-seq results into immunological and cell-biological meaning.

You are reviewing a full pipeline report produced by the OmicSage CITE-seq pipeline.
Your job is to produce a structured interpretation document (`cite_report_review.md`)
that mirrors — and in several areas exceeds — what the automated AI pipeline would
have produced.

**Instructions:**
- Run all 10 tasks in Section 4 in sequence without stopping.
- Do not ask clarifying questions. If a value is missing from the report, write "not reported".
- Return the complete `cite_report_review.md` using the output template in Section 6 at the very end.
- Every factual claim must be grounded as described in Section 5.

---

## SECTION 2 — STUDY CONTEXT BLOCK

**Instructions for the model:**
Before running any task in Section 4, read the **Normalize ADT** tab (the first tab)
of the attached HTML report. It contains the dataset name, key metrics, and provenance
table for the ADT arm of the pipeline.

Extract the following fields directly from that tab and fill the study context block
below. Only fill a field from the report if the analyst has left it as a placeholder
(i.e. still shows the example text in brackets). Never overwrite a field the analyst
has already filled in.

Extraction rules per field:
- **Tissue**: read from the dataset title or description in the Normalize ADT tab header
- **Disease**: infer from the same description; if the word "healthy" or "PBMC/BMMC"
  appears and no disease is mentioned, write "healthy (no disease)"
- **Species**: infer from tissue/dataset context; normalise to: human / mouse / [other as written]
- **Conditions**: infer from the dataset title (e.g. healthy donors, tumour vs normal,
  stimulated vs control); if unclear write "not reported"
- **N cells (post-QC)**: read the "Cells" stat card from the Normalize ADT tab
- **N proteins (ADT panel)**: read the "Proteins" stat card from the Normalize ADT tab
- **N donors**: search the description or provenance table for a phrase like "N donors"
  or "N batches"; if not found write "not reported"
- **Batch key**: check the provenance table; prefer the column named: batch, donor,
  sample_id in that order; if no batch column is present write "null (single batch)"
- **Integration method used**: read from the Integration tab — look for the Run Summary
  provenance table (e.g. MOFA+, totalVI, WNN, or combinations)
- **Biological question**: leave blank if the analyst has not stated one; do not
  infer it here — you will infer it during Task 1 if needed
- **Known biology**: leave blank unless the analyst has filled it in

After filling the block, print it in full so the analyst can verify it before
you proceed to Section 4.

```
Tissue:                   [e.g. bone marrow / PBMC / liver]
Disease:                  [e.g. AML / null if healthy]
Species:                  [human / mouse / other]
Conditions:               [e.g. healthy donors / tumour vs normal / stimulated vs control]
N cells (post-QC):        [from Normalize ADT tab]
N proteins (ADT panel):   [from Normalize ADT tab]
N donors:                 [from report, or "not reported"]
Batch key:                [e.g. donor / null if single batch]
Integration method used:  [e.g. MOFA+ + totalVI / MOFA+ only / WNN]
Biological question:      [your question — or leave blank for the model to infer]
Known biology:            [optional — e.g. "NeurIPS 2021 BMMC dataset, 10 donors, 13 cell types expected"]
```

---

## SECTION 3 — REPORT STRUCTURE GUIDE

The attached HTML report (`cite_00_combined_report.html`) was produced by the OmicSage
CITE-seq pipeline. It contains the following 10 tabs. Read each tab carefully before
running the corresponding task.

| Tab | Contents |
|-----|----------|
| **Normalize ADT** | CLR normalization provenance (CLR axis, raw vs CLR distributions, DSB warning if empty droplets absent), N cells, N proteins, raw count statistics |
| **Doublets** | Marker-pair doublet scores (CD3/CD19, CD3/CD14), doublet rate, whether filter was applied or flagged only |
| **Reduce ADT** | ADT PCA (N PCs computed and used), variance explained per PC, UMAP colored by batch and by top proteins |
| **Harmony ADT** | Batch correction on ADT PCA embedding; before/after UMAP; iLISI / graph connectivity scores |
| **Annotate ADT** | Leiden clustering on ADT (resolution used, N clusters), top ADT markers per cluster (ranked by CLR score), optional `adt_celltype` / `adt_celltype_manual` labels, UMAP |
| **Integration** | MOFA+ and/or totalVI results; integration quality metrics (iLISI, cLISI, ASW, graph connectivity); joint UMAP; factor loadings (MOFA+) or latent space summary (totalVI) |
| **DEG / DPE** | Differential Protein Expression (DPE): Wilcoxon one-vs-rest on CLR ADT values, grouped by `adt_celltype_manual`; RNA DEG if run; N significant proteins and genes |
| **GSEA** | Pathway analysis on the RNA arm only (GO BP, KEGG, Reactome); NES, adjusted p-values; note if not run |
| **Protein-RNA** | Per-protein Pearson/Spearman correlation between ADT and RNA expression; high, low, and negative correlators flagged |
| **Epitope** | Protein-level cell type marker summary; epitope coverage table; isotype control QC |

If a tab is missing or empty, note "tab not present" in the relevant section of your output.

---

## SECTION 4 — TASK LIST
> Run all 10 tasks in sequence. Do not stop between tasks. Do not ask questions.
> Produce the final cite_report_review.md using the template in Section 6.

---

### Task 1 — ADT Normalization QC

Read the Normalize ADT tab. Evaluate the following:

1. **CLR distribution quality**: Is the CLR normalization appropriate?
   - CLR axis = 0 (per-protein across cells) is the standard for surface protein data.
     Flag if axis = 1 was used without justification.
   - The CLR distribution per protein should be roughly bimodal (negative vs positive
     population). A flat or unimodal distribution may indicate a failed staining or
     ambient protein signal.
2. **DSB check**: Was DSB normalization applied?
   - DSB (Denoised and Scaled by Background) is the recommended method when empty
     droplets are available. It removes ambient protein background and produces more
     interpretable distributions.
   - If DSB was not applied (no empty droplets provided), flag as **warning** and note
     that downstream results rely on CLR alone.
3. **Raw count statistics**: Is the raw median total counts/cell within a plausible range?
   - Typical CITE-seq: 1,000–5,000 raw ADT counts/cell for panels of 50–200 proteins.
     Flag if < 500 (under-staining) or > 20,000 (over-staining or doublets).
4. **Isotype controls**: Are isotype controls present in the panel?
   - Check the Epitope tab for isotype controls. If absent, note that background
     estimation is limited.
5. List any normalization flags as info / warning / critical.

---

### Task 2 — Doublet Assessment

Read the Doublets tab. Evaluate the following:

1. **Doublet rate**: Is the observed doublet rate within the expected range?
   - Expected: ~0.8% per 1,000 cells loaded (10x Genomics rule of thumb).
   - For 10,000–25,000 cells: expect ~8–20% doublet flag rate from computational tools.
   - Flag if the doublet rate is 0% (likely under-detection, especially for large datasets)
     or > 25% (possible over-detection or high loading density).
2. **Marker pairs evaluated**: Were the correct lineage-exclusive pairs used?
   - Standard pairs: CD3/CD19 (T cell vs B cell), CD3/CD14 (T cell vs Monocyte).
   - Additional pairs appropriate for bone marrow: CD19/CD14, CD56/CD3.
   - Flag as info if additional relevant pairs were not evaluated.
3. **Filter vs flag**: Was the doublet filter applied or only flagged?
   - If flagged only (filter not applied), check whether `obs[adt_predicted_doublet]`
     was used in downstream steps. Note if doublets were carried into clustering.
4. **Cross-modal concordance**: Does the ADT doublet rate align with the RNA Scrublet
   doublet rate from the RNA pipeline (if known from the context block)?
   Report as concordant / discordant / RNA rate not available.

---

### Task 3 — ADT Dimensionality Reduction & Clustering Assessment

Read the Reduce ADT and Annotate ADT tabs. Evaluate the following:

1. **N PCs used**: Is the number of PCs appropriate for this ADT panel size?
   - For panels of 50–100 proteins: 15–25 PCs typical.
   - For panels of 100–200 proteins: 20–40 PCs typical.
   - Flag if N PCs used is < 10 (under-representation) or > 50 (over-fitting noise).
2. **Variance explained**: Does the cumulative variance plateau reasonably?
   - Report the cumulative variance explained by the chosen PCs.
3. **Clustering resolution**: Is the resolution appropriate for the N cells and expected
   biology?
   - Low (< 5,000 cells): resolution 0.1–0.3 typical for ADT
   - Medium (5,000–30,000): resolution 0.1–0.4 typical
   - High (> 30,000): resolution 0.2–0.6 typical
   **State explicitly whether the chosen resolution is appropriate / too coarse /
   too fine for this cell count. Do not leave this implicit.**
4. **N ADT clusters vs N RNA clusters**: Are the two arms producing comparable
   resolution? Divergence may indicate one modality is dominating or that the ADT
   panel has limited discriminatory power for rare cell types.
5. Identify any ADT clusters that are candidates for sub-clustering (large, mixed
   marker profiles) or merging (near-identical marker profiles).

---

### Task 4 — Batch Correction Assessment (ADT Harmony)

Read the Harmony ADT tab. Evaluate the following:

1. **iLISI score**: Higher values (→ N batches) indicate good batch mixing.
   - Report the value and interpret: well mixed / partially mixed / poor mixing.
2. **Graph connectivity**: Values > 0.5 indicate good global connectivity post-correction.
3. **cLISI score**: Values close to 1.0 indicate cell type labels are not scrambled
   by batch correction (biological structure preserved).
4. **ASW (Average Silhouette Width)**: Values > 0.2 indicate that cell types are still
   separable after correction.
5. **Over-correction check**: If cLISI ≈ 1.0 but iLISI is very low (batches not mixed),
   Harmony may have failed to converge. If both iLISI is high and ASW is low, Harmony
   may have over-corrected (biological signal lost). Flag either scenario as warning.
6. **Visual check instruction**: Review the before/after UMAP figures in the Harmony tab.
   Note whether batch-specific clusters visible before correction are dispersed after.

---

### Task 5 — ADT Cluster Annotation

Read the marker protein table in the Annotate ADT tab.

For each ADT cluster, provide:
- **Predicted cell type** — your interpretation based on the top ADT marker proteins
- **Reviewer confidence** — high / medium / low, based on specificity of the protein markers
- **Supporting markers** — specific proteins from the table that support your call
  (use canonical surface marker names, e.g. CD4, CD8a, CD19, CD14, CD56)
- **Alternatives** — other plausible cell types if confidence is medium or low,
  and what additional proteins would distinguish them
- **Agreement with `adt_celltype`** — does your annotation agree with the automated
  `adt_celltype` or `adt_celltype_manual` label in the report? If not, explain the discrepancy.

**CITE-seq-specific note on confidence:**
Surface protein markers are generally more specific than gene markers for major immune
lineages. High confidence is appropriate when ≥ 2 lineage-defining markers are concordant
(e.g. CD3+ and CD4+ for CD4 T cells). Low confidence is warranted when only one marker
is present or when isotype controls are absent.

Present results as a table (see output template, Section 6).

---

### Task 6 — Protein-RNA Coherence

This is a CITE-seq-specific cross-modal task. Read both the Annotate ADT tab (ADT
cluster labels) and the RNA annotation from context (RNA `cell_type_vote` if mentioned
in the study context block or known from the benchmark dataset).

For each ADT cluster:
1. **Does the ADT-based annotation agree with the expected RNA-based annotation?**
   - Concordant: both arms predict the same cell type.
   - Discordant — biologically expected: e.g. ADT detects activated vs resting
     state that RNA clusters do not separate cleanly.
   - Discordant — technically concerning: e.g. ADT says Monocyte but RNA markers
     are T cell genes. Flag as warning or critical.
2. **Protein-RNA correlation tab**: Read the Protein-RNA tab.
   - List the top 5 highest-correlating proteins and their r values.
   - List the top 5 lowest-correlating (or negatively correlating) proteins and their r values.
   - For each low correlator: is the decoupling biologically expected
     (post-translational regulation, protein stability, receptor shedding) or potentially
     a technical artefact (staining failure, antibody cross-reactivity)?
3. Flag any protein with r < 0.1 as an annotation risk: it cannot be used to validate
   RNA-based annotation.

---

### Task 7 — DPE Validation (Differential Protein Expression)

Read the DEG / DPE tab. For each cell type comparison listed:

1. Classify each reported significant protein as:
   - **Expected** — a well-known surface marker for this cell type (e.g. CD19 up in B cells)
   - **Unexpected** — a protein not typically differentially expressed in this context;
     explain whether this is biologically interesting or potentially artefactual
   - **Panel artefact** — a protein that appears significant due to panel imbalance
     (CLR normalization is sensitive to the protein composition of the panel; if one
     lineage expresses 70% of the panel, all other proteins are artificially suppressed)
2. **Required check — isotype controls**: Are any isotype controls appearing in the
   DPE results? If yes, flag as critical: background signal is contaminating results.
3. **Required check — CLR panel imbalance**: For any cell type where > 30% of the
   top 10 DPE hits are from a single protein family (e.g. all CD3-complex subunits,
   or all HLA variants), flag as warning: CLR normalization may be amplifying
   compositional contrast rather than true differential expression.
4. **Identify discovery highlights**: proteins not typically reported as differentially
   expressed for this cell type in this tissue context, that warrant follow-up.

If RNA DEG results are also present in the tab, perform the same classification for
the top 10 RNA DEGs per comparison and note whether the RNA and protein arms
are coherent (e.g. CD4 protein up in a cluster AND CD4 RNA up in the same cluster).

---

### Task 8 — Integration Quality Review

Read the Integration tab. This tab contains first-class biological results — read it fully.

**MOFA+ (if present):**
1. **Factor loadings**: Read the top factors. Do they correspond to known biology?
   - Factor 1 should typically capture the dominant source of variation (cell type, donor,
     or cell cycle). State which source of variation each reported factor most likely captures.
   - A factor loaded primarily on a single donor or batch is a batch effect, not biology. Flag.
   - A factor loaded on both RNA and ADT for the same cell type signature is the ideal
     signal (multi-modal concordant factor).
2. **iLISI / graph connectivity / cLISI / ASW** (MOFA+ metrics): interpret using the
   same thresholds as Task 4. Report values and verdict.

**totalVI (if present):**
1. **iLISI / graph connectivity / cLISI / ASW** (totalVI metrics): interpret and compare
   to MOFA+ if both were run. Which method produced better batch mixing (iLISI)
   while preserving biological structure (cLISI)?
2. **Joint latent space**: Does the totalVI UMAP separate major cell types cleanly?
   Note any major cell types that appear merged in the joint embedding.

**Head-to-head comparison (if both methods were run):**
Report which method is preferred for downstream analysis based on the four metrics
and explain why. This recommendation will feed into Task 10.

---

### Task 9 — Literature Linking

For each unexpected protein from Task 7 and each integration factor interpretation
from Task 8, link supporting literature using the following rule:

**If you have a PubMed or bioRxiv search tool available:**
Use it directly. Report: protein/finding → PMID + title + one-sentence biological context.
Do not include abstract text. PMIDs and titles only.

**If you do not have a search tool:**
Provide the exact PubMed search string for each protein or finding so the analyst
can run them manually:
```
[PROTEIN or GENE] AND [tissue/disease] AND [cell type] AND ("CITE-seq" OR "surface proteomics" OR "ADT")
```
Label the output clearly as "Manual search required — strings provided below."

Format all literature links as: Protein/Finding | PMID | Title | Context (one sentence)

---

### Task 10 — Downstream Suggestions + Write cite_report_review.md

**Part A — Downstream Suggestions:**
Based on the cell types found, the DPE results, the integration quality, and the
biological question stated in Section 2, suggest the most valuable next analysis steps.

Prioritize suggestions that are CITE-seq-specific and not available from RNA alone:

| Priority | Step | Rationale | Recommended tool | Expected output |
|----------|------|-----------|-----------------|----------------|
| 1 | WNN joint embedding | If WNN was not run, it should be: it weights RNA and ADT by their local informativeness per cell | Seurat v5 WNN / muon wnn | Joint UMAP with higher resolution than either modality alone |
| 2 | (add from your analysis) | | | |

Additional suggestions to consider where relevant to the data:
- **TCR/BCR integration**: if T/B cell enrichment is seen, paired VDJ data may be available
- **Receptor-ligand with surface proteins**: use ADT-confirmed receptor expression for
  CellChat/LIANA instead of RNA proxies (higher confidence)
- **DSB re-normalization**: if DSB was not applied and empty droplets are obtainable,
  reprocessing with DSB will improve protein distribution interpretability
- **Sub-clustering of large ADT clusters**: use both RNA and ADT jointly
- **Epitope QC re-run with isotype control subtraction**: if isotype controls are present
  in the panel but were not used for background correction

**Part B — Assemble cite_report_review.md:**
Compile all findings from Tasks 1–10 into the output template in Section 6.

**Grounding requirement (Section 5 applies here):**
- Every factual claim must cite a metric from the report, a specific protein/gene name, or a PMID.
- Never invent numbers. If a value is not in the report, write "not reported".
- Aim for ≥ 85% of factual sentences to be explicitly grounded.

Return the complete document — do not truncate or summarise.

---

## SECTION 5 — GROUNDEDNESS RULE

Every factual claim in `cite_report_review.md` must be grounded by one of:

| Ground type | Example |
|-------------|---------|
| Metric from report | "134 proteins in panel" / "CLR axis = 0" / "iLISI = 0.054 (MOFA+)" |
| Specific protein or gene name | "CD19 upregulated log2FC = 2.4 in cluster 3" |
| PMID from Task 9 | "PMID:12345678" |
| Widely accepted and no search available | "literature-supported" (use sparingly) |

**Never invent numbers or protein names.**
If a value is not present in the HTML report, write "not reported" — do not estimate.
Aim for ≥ 85% of factual sentences to be explicitly grounded.

---

## SECTION 6 — OUTPUT TEMPLATE

Copy everything below this line into your response and fill it in completely.
Do not truncate any section. Return as a single markdown document.

---

```markdown
# CITE-seq AI Interpretation Report
**Dataset:** [dataset name or GEO accession]
**Date:** [today's date]
**Reviewer:** Manual Review Mode (CITE D1)
**Model used:** [e.g. Claude Sonnet 4.6 / GPT-4o / Claude Opus 4]
**Tissue:** [from Section 2]
**Disease:** [from Section 2]

---

## 1. Study Overview

- **Tissue:** [value]
- **Disease:** [value]
- **Species:** [value]
- **Conditions:** [value]
- **N cells (post-QC):** [value]
- **N proteins (ADT panel):** [value]
- **N donors:** [value]
- **Batch key:** [value or "none"]
- **Integration method used:** [value]
- **Biological question:** [stated by analyst / inferred by reviewer]

---

## 2. ADT Normalization Assessment

- **Normalization method:** [CLR only / DSB + CLR]
- **CLR axis:** [0 (per-protein) / 1 (per-cell)]
- **DSB applied:** [yes / no — reason if no]
- **Raw median counts/cell:** [value]
- **N proteins in panel:** [value]
- **Isotype controls present:** [yes / no / not reported]
- **Overall normalization assessment:** [appropriate / review recommended / issues found]

| Severity | Flag | Suggestion |
|----------|------|------------|
| [info/warning/critical] | [description] | [action] |

---

## 3. Doublet Assessment

- **Doublet rate:** [%] — [within expected / below expected / above expected] for [N cells] cells
- **Marker pairs evaluated:** [list]
- **Filter applied:** [yes / no — flagged only]
- **Cross-modal concordance with RNA doublets:** [concordant / discordant / not available]
- **Flags:**

| Severity | Flag | Suggestion |
|----------|------|------------|
| [info/warning/critical] | [description] | [action] |

---

## 4. ADT Dimensionality Reduction & Clustering

- **N PCs used:** [value] — [appropriate / too few / too many] for [N proteins] proteins
- **Clustering resolution:** [value] — [appropriate / too coarse / too fine] for [N cells] cells
- **N ADT clusters:** [value]
- **N RNA clusters (if known):** [value or "not reported"]
- **Sub-clustering candidates:** [cluster IDs and rationale, or "none identified"]
- **Merge candidates:** [cluster IDs and rationale, or "none identified"]

---

## 5. Batch Correction Assessment (Harmony ADT)

| Metric | Value | Interpretation |
|--------|-------|----------------|
| iLISI | | well mixed / partially mixed / poor mixing |
| Graph connectivity | | [> 0.5 = good / < 0.5 = poor] |
| cLISI | | [≈ 1.0 = biology preserved / < 0.9 = biology partially scrambled] |
| ASW (label) | | [> 0.2 = separable / < 0.2 = poor separation] |

- **Over-correction risk:** [yes / no — rationale]
- **Visual assessment of before/after UMAP:** [describe what changed]
- **Overall batch correction verdict:** [good / adequate / review recommended / failed]

---

## 6. ADT Cell Type Annotation

### 6a. Cluster-level predictions

| Cluster | N cells | Predicted cell type | Reviewer confidence | Supporting markers | Alternatives | Agreement with adt_celltype |
|---------|---------|--------------------|--------------------|-------------------|--------------|------------------------------|
| 0 | | | | | | |

### 6b. Surface protein marker reference (CITE-seq specific)

| Cell type | Canonical markers | Tissue-context-specific markers | Distinguishing markers | Counter-markers | Confidence |
|-----------|--------------------|--------------------------------|----------------------|-----------------|------------|
| | | | | | |

---

## 7. Protein-RNA Coherence

### 7a. Cross-modal annotation agreement

| ADT cluster | ADT annotation | Expected RNA annotation | Agreement | Discordance type | Action |
|-------------|---------------|------------------------|-----------|-----------------|--------|
| | | | Concordant / Discordant | Biological / Technical / N/A | |

### 7b. Protein-RNA correlation summary

**Top 5 highest-correlating proteins:**

| Protein | r value | Interpretation |
|---------|---------|----------------|
| | | |

**Top 5 lowest / negatively correlating proteins:**

| Protein | r value | Decoupling explanation | Annotation risk |
|---------|---------|----------------------|----------------|
| | | Biological / Technical / Unclear | High / Low |

---

## 8. DPE Validation (Differential Protein Expression)

**Comparison: [cell type] vs rest**

| Protein | log2FC | adj p-value | Classification | Rationale |
|---------|--------|-------------|---------------|-----------|
| | | | Expected / Unexpected / Panel artefact | |

**Isotype control check:** [clean / contaminated — describe]
**CLR panel imbalance check:** [pass / warning — describe]

**Discovery highlights:**
- [Protein]: [why this is noteworthy — grounded with log2FC or PMID]

---

## 9. Integration Quality Review

### MOFA+ metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| iLISI | | |
| Graph connectivity | | |
| cLISI | | |
| ASW (label) | | |

**Factor interpretation:**

| Factor | Dominant loadings | Likely biological source | Multi-modal concordance |
|--------|------------------|------------------------|------------------------|
| Factor 1 | | | yes / no |

### totalVI metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| iLISI | | |
| Graph connectivity | | |
| cLISI | | |
| ASW (label) | | |

**Head-to-head verdict:** [MOFA+ preferred / totalVI preferred / comparable — rationale]

**Recommended embedding for downstream analysis:** [MOFA+ / totalVI / RNA alone / ADT alone]

---

## 10. Literature Links

| Protein / Finding | PMID | Title | Context (one sentence) |
|-------------------|------|-------|------------------------|
| | | | |

*PMIDs and titles only. No abstract text.*

[If no search tool was available, list manual PubMed search strings here instead.]

---

## 11. Downstream Suggestions

| Priority | Step | Rationale | Recommended tool | Expected output |
|----------|------|-----------|-----------------|----------------|
| 1 | | | | |
| 2 | | | | |
| 3 | | | | |

---

## 12. Summary

**Key findings:**
1. [finding — grounded with metric / protein / PMID]
2. [finding — grounded with metric / protein / PMID]
3. [finding — grounded with metric / protein / PMID]

**Open questions:**
- [question raised by the CITE-seq analysis]

**Suggested validation experiments:**
- [wet-lab or computational experiment to validate a key finding]

---

## Appendix: Report Figures Checklist

Check each item against the HTML report and mark ✓ present or ✗ missing.

| # | Figure / element | Status | Note |
|---|-----------------|--------|------|
| 1 | CLR distribution plot (before vs after, per-protein or summary) | | |
| 2 | DSB normalization applied or DSB warning present | | |
| 3 | Doublet score distribution plot or summary table | | |
| 4 | ADT PCA variance explained plot | | |
| 5 | ADT UMAP colored by batch (before Harmony) | | |
| 6 | ADT UMAP colored by batch (after Harmony) | | |
| 7 | ADT UMAP colored by adt_celltype or Leiden cluster | | |
| 8 | Integration metrics table (iLISI / cLISI / ASW) present for each method | | |
| 9 | MOFA+ factor loading plot OR totalVI latent UMAP | | |
| 10 | DPE results exclude isotype controls | | |
| 11 | Protein-RNA correlation heatmap or ranked table | | |
| 12 | Epitope coverage table present | | |
```

---

*OmicSage CITE D1 Manual Review Mode — CITE_MASTER_PROMPT.md v1.0*
*No Python. No infrastructure. Attach + send.*
