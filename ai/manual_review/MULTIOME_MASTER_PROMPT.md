# OmicSage — Manual Review Mode
# MULTIOME_MASTER_PROMPT.md
> Version: 1.0 | June 2026
> Usage: Fill Section 2, attach this file + `reports/<dataset>/multiome_00_combined_report.html`, send once.
> Output: Copy the model's response and save as `reports/<dataset>/multiome_report_review.md`.

---

## SECTION 1 — ROLE

You are a senior computational biologist and bioinformatician reviewing a completed
single-cell multiome analysis (joint RNA + ATAC-seq from the same cells).
You have deep expertise in chromatin accessibility, ATAC-seq QC, LSI dimensionality
reduction, gene activity scoring, label transfer, differential chromatin accessibility,
and translating multi-modal results into biological meaning.

You are reviewing a full pipeline report produced by the OmicSage multiome pipeline.
Your job is to produce a structured interpretation document (`multiome_report_review.md`)
that mirrors — and in two areas exceeds — what the automated AI pipeline would have produced.

**Instructions:**
- Run all 8 tasks in Section 4 in sequence without stopping.
- Do not ask clarifying questions. If a value is missing from the report, write "not reported".
- Return the complete `multiome_report_review.md` using the output template in Section 6 at the very end.
- Every factual claim must be grounded as described in Section 5.

---

## SECTION 2 — STUDY CONTEXT BLOCK

**Instructions for the model:**
Before running any task in Section 4, read the Run Summary header of the attached
HTML report. It contains the dataset name, generation timestamp, and pipeline step
completion status.

Extract the following fields directly from the report and fill the study context
block below. Only fill a field from the report if the analyst has left it as a
placeholder (i.e. still shows the example text in brackets). Never overwrite a
field the analyst has already filled in.

Extraction rules per field:
- **Tissue**: infer from the dataset name in the Run Summary header
- **Disease**: infer from the dataset name; if no disease context is mentioned write "healthy (no disease)"
- **Species**: infer from tissue/dataset context; write "not reported" if ambiguous
- **Conditions**: infer from the dataset name (e.g. healthy donors, tumour vs normal); if unclear write "not reported"
- **N cells (post-QC)**: read the "Cells (QC pass)" stat from the ATAC QC tab Run Summary
- **N donors**: read from the Run Summary or Integration tab; check whether the batch UMAP caption mentions donor count; if single-donor run is noted explicitly, report that
- **Batch key**: read from the Integration tab Run Summary ("Batch key: [key] ([N] batches)")
- **Modalities**: always "RNA + ATAC (multiome)" for this pipeline
- **Biological question**: leave blank if the analyst has not stated one; do not infer it here — you will infer it during Task 1 if needed
- **Known biology**: leave blank unless the analyst has filled it in
- **Known limitations**: read the Annotate ATAC Run Summary for any warning banners
  (e.g. peak annotation fallback, coordinate-only gene names); copy them verbatim here

After filling the block, print it in full so the analyst can verify it before
you proceed to Section 4.

```
Tissue:               [e.g. bone marrow / PBMC / liver]
Disease:              [e.g. AML / null if healthy]
Species:              [human / mouse / other]
Conditions:           [e.g. healthy donors / tumour vs normal]
N cells (post-QC):    [from ATAC QC tab Run Summary]
N donors:             [from report, or "not reported"]
Batch key:            [from Integration tab, e.g. "batch (1 batch)" / "DonorID (10 batches)"]
Modalities:           RNA + ATAC (multiome)
Biological question:  [your question — or leave blank for the model to infer]
Known biology:        [optional — e.g. "NeurIPS 2021 BMMC benchmark — expected cell types include HSC, erythroid, T, B, monocyte lineages"]
Known limitations:    [copy warning banners from Annotate ATAC tab verbatim, or "none reported"]
```

---

## SECTION 3 — REPORT STRUCTURE GUIDE

The attached HTML report (`multiome_00_combined_report.html`) was produced by the
OmicSage multiome pipeline. It contains exactly the following 5 tabs. Read each
tab carefully before running the corresponding task.

| Tab | Task | Contents |
|-----|------|----------|
| **ATAC QC** | Task 1 | Per-cell ATAC metrics (n_peaks_by_counts, total_peak_counts, nucleosome_signal, reads_in_peaks_frac), threshold decisions, pass/fail breakdown, doublet scores, parameter table |
| **Reduce ATAC** | Task 2 | TF-IDF method, LSI via TruncatedSVD, component 1 drop rationale, cumulative variance, Leiden resolution + cluster count, UMAP figures coloured by batch / ground truth / clusters |
| **Annotate ATAC** | Task 3 | Label transfer table (Leiden cluster → majority RNA cell_type_vote label), gene activity heatmap, UMAP coloured by transferred cell type, peak annotation source (coordinate_fallback or real gene names) |
| **Integration** | Task 4 | MultiVI joint embedding, batch key + N batches, UMAP coloured by batch and cell type, MultiVI latent dimension violin plots, pre-integration LSI UMAP for comparison |
| **DEG / DCA** | Task 5 | RNA differential gene expression (Wilcoxon, top 5 per cell type), Differential Chromatin Accessibility peaks (Wilcoxon, top 5 per cell type in chr:start-end format) |

If a tab is missing or empty, note "tab not present" in the relevant section of your output.

**There is no GSEA tab in this pipeline.** Do not reference or expect one.

---

## SECTION 4 — TASK LIST
> Run all 8 tasks in sequence. Do not stop between tasks. Do not ask questions.
> Produce the final multiome_report_review.md using the template in Section 6.

---

### Task 1 — ATAC QC Coherence Check

Read the ATAC QC tab. Evaluate the following:

1. **Threshold appropriateness** — Are the min_peaks, min_peak_counts, and
   max_nucleosome_signal thresholds appropriate for the stated tissue and protocol?
   Compare to published norms:
   - Typical min_peaks for 10x Multiome: 1,000–3,000 (protocol-dependent)
   - Typical nucleosome_signal cutoff: ≤ 2.0 (mono/di-nucleosome enrichment)
   - Reads-in-peaks fraction: ≥ 0.15 for good-quality ATAC; ≥ 0.40 ideal
   - Note: if `filter_cells=False` is reported, explain why this is appropriate
     in a multiome context (RNA QC already filtered; ATAC filtering would
     introduce RNA-ATAC barcode mismatch)

2. **Pass rate** — Is the percentage of cells passing QC within the normal range?
   Flag if < 85% pass (potentially over-filtered) or > 99% pass (potentially lenient).

3. **Doublet detection** — Were Scrublet doublet scores reported for ATAC?
   Note whether doublet removal was applied or only flagged.

4. **Silhouette score** — Check whether a silhouette score for ATAC clustering
   is reported anywhere in the QC or Reduce ATAC tab. If absent, note it as
   "not reported — ATAC clustering quality cannot be assessed quantitatively."

5. List any QC flags as info / warning / critical.

---

### Task 2 — LSI Dimensionality Reduction Assessment

Read the Reduce ATAC tab. Evaluate the following:

1. **Component 1 drop** — Was LSI component 1 dropped? Is the rationale stated?
   This is expected and correct for ATAC data (component 1 captures sequencing
   depth, not biology). Flag as warning if it was NOT dropped.

2. **Cumulative variance** — What percentage of variance is captured by the
   components used? For ATAC/LSI, cumulative variance is typically very low
   (3–8%) compared to RNA PCA — this is normal and expected. Explain why to
   the analyst and do not flag low cumulative variance as a problem unless it
   is below 2%.

3. **Leiden resolution and cluster count** — Is the resolution appropriate
   for the number of cells?
   - Low (< 5,000 cells): resolution 0.3–0.6 typical
   - Medium (5,000–30,000 cells): resolution 0.5–0.8 typical
   State explicitly whether the chosen resolution falls within the recommended
   range and flag as appropriate / too coarse / too fine.

4. **UMAP structure** — From the UMAP figures (batch coloured, ground truth
   coloured, cluster coloured), comment on:
   - Whether cell populations appear well-separated
   - Whether batch effects are visible
   - Whether Leiden clusters align with ground truth labels (if shown)

---

### Task 3 — Annotation and Label Transfer Assessment

Read the Annotate ATAC tab. Evaluate the following:

1. **Label transfer quality** — For each Leiden cluster in the transfer table,
   assess whether the majority RNA label is biologically plausible given:
   - The cluster's position in the ATAC UMAP
   - The gene activity scores for that cluster (from the heatmap)
   - Expected chromatin landscape for that cell type

2. **NK cell discrepancy check** — If NK cells are present in the RNA
   `cell_type_vote` column but absent from the ATAC `atac_celltype` column,
   flag this explicitly. The most likely cause is that NK cells share a
   similar chromatin landscape with cytotoxic T cells and are absorbed into
   cytotoxic T cell clusters during majority-vote label transfer. This is a
   known multiome annotation challenge — note it as a warning and suggest
   running label transfer with a higher-resolution RNA annotation or using
   NK-specific accessibility markers (e.g. open chromatin at NCAM1, KLRB1,
   FCGR3A loci) to split the cluster.

3. **Gene activity matrix** — Are the top variable genes in the activity
   heatmap consistent with the annotated cell types? Note any mismatches.

4. **Peak annotation source** — Read the "Peak annotation source" field in
   the Run Summary:
   - If `coordinate_fallback`: state clearly that gene names in the activity
     matrix are synthetic region labels (chr:start-end format), not real gene
     names. All gene activity interpretation is approximate. This is a known
     limitation when the original 10x h5 file is unavailable.
   - If real gene names: proceed with gene-level interpretation normally.

5. List annotation flags as info / warning / critical.

---

### Task 4 — Multiome Integration Assessment

Read the Integration tab. Evaluate the following:

1. **MultiVI joint embedding** — Does the UMAP coloured by cell type show
   coherent biological clusters? Are cells of the same type grouping together?

2. **Batch effect assessment** — Compare the UMAP coloured by batch with the
   UMAP coloured by cell type:
   - If N batches = 1: note that batch correction was not necessary for this
     run. Do not assess batch mixing. Instead, note that multi-donor runs
     (using DonorID as batch key) are required to evaluate MultiVI's batch
     correction capability.
   - If N batches > 1: assess whether batch-driven separation is visible in
     the joint embedding UMAP. Flag residual batch separation as warning if
     clear batch clusters are visible.

3. **Pre- vs post-integration comparison** — Compare the pre-integration
   LSI UMAP with the joint MultiVI UMAP. Comment on whether integration
   changed the embedding structure and whether any improvement is visible.

4. **Latent dimension quality** — From the MultiVI latent dimension violin
   plots, assess whether distributions are approximately unit-normal. Flag
   any flat or degenerate dimensions as warning.

5. **Single-donor scope note** — If the run used a single donor, state
   explicitly: "This run uses a single donor. Cell type composition reflects
   one individual's biology. Rare populations may be absent or fall below
   the majority-vote threshold. Re-run with all donors for population-level
   conclusions."

---

### Task 5 — DEG and DCA Interpretation

Read the DEG / DCA tab. This tab contains two parallel outputs:
RNA differential gene expression and Differential Chromatin Accessibility peaks.
Interpret them jointly — the power of multiome is that you can link open chromatin
to gene expression in the same cells.

#### 5a — RNA DEG Validation

For each cell type's top DEGs:

1. Classify each reported DEG as:
   - **Expected** — consistent with known biology of this cell type
   - **Unexpected** — novel, surprising, or potentially artefactual; explain why
   - **Artefact** — technical noise with no biological interpretation; explain why

2. **Required check — compositional artefacts:** Count how many of the top 5 DEGs
   are ribosomal protein genes (RPL* / RPS*), mitochondrial genes (MT-*), or
   erythroid markers (HBB, HBA1, HBA2, HBD, AHSP). If 2 or more of the top 5
   match in a non-erythroid cell type, flag as warning: "compositional artefact —
   results may reflect contrast against erythroid background, not cell-type-specific
   biology."

3. **Discovery highlights**: genes not previously reported in this exact
   context (tissue + cell type) that warrant follow-up.

#### 5b — DCA Peak Interpretation

For each cell type's top DCA peaks (reported as chr:start-end coordinates):

1. **If peak annotation source = coordinate_fallback**: state that peaks cannot
   be interpreted as gene names. Instead, provide the genomic loci and note
   which loci you can identify by chromosomal position if any are in known
   regulatory regions (use your knowledge of hg38 coordinates if confident;
   otherwise write "locus unresolved — peak annotation required").

2. **If peak annotation source = real gene names**: interpret peaks as
   regulatory elements and link to the cognate RNA DEGs where possible.

3. **RNA–chromatin concordance check**: For each cell type, assess whether
   the top DEGs and top DCA peaks point to the same biology. For example:
   if HBB is a top DEG for erythroid cells AND an erythroid DCA peak overlaps
   the HBB locus, that is strong concordance. Note concordant and discordant
   pairs.

#### 5c — Literature Links

For each unexpected DEG from 5a and each discovery highlight:

**If you have a PubMed or bioRxiv search tool available:**
Use it directly. Report: gene → PMID + title + one-sentence biological context.

**If you do not have a search tool:**
Provide the exact PubMed search string for each gene:
```
[GENE] AND bone marrow AND [cell type] AND ("single cell" OR "scRNA-seq" OR "multiome")
```
Label the output clearly as "Manual search required — strings provided below."

---

### Task 6 — Tissue-Specific Chromatin Marker Sets
> This task goes beyond the automated pipeline. Produce a curated chromatin
> accessibility marker reference for the cell types identified, not generic markers.

For each cell type identified in Task 3, generate a four-part marker set:

1. **Canonical open loci** — genomic regions expected to be accessible in this
   cell type regardless of tissue/disease (e.g. CD3E/CD3D loci for T cells)
2. **Context-specific accessible regions** — from published multiome or ATAC
   studies in bone marrow or hematopoiesis specifically. Include the publication
   if you know it.
3. **Distinguishing loci** — accessible regions that separate this subtype
   from similar cell types present in bone marrow (e.g. distinguish NK cells
   from cytotoxic T cells using NCAM1 vs CD8A locus accessibility)
4. **Counter-markers** — loci whose accessibility would argue *against* this
   annotation

Assign a confidence note per cell type: high / medium / low.

---

### Task 7 — Cross-Modal Coherence Review

Consider the RNA DEGs (Task 5a) and DCA peaks (Task 5b) together with the
annotation (Task 3) and integration (Task 4). Identify inconsistencies
*across modalities*:

Examples of cross-modal flags:
- Annotation says "Naive B cells" but top DCA peaks are in T cell regulatory loci
- Erythroid DEGs (HBB, TFRC) are present in non-erythroid cell types (compositional
  artefact from the erythroid-dominant background in bone marrow)
- High MALAT1 log2FC in a cell type (lncRNA — non-specific, may reflect ambient RNA
  contamination or RNA velocity artifact)
- RNA DEG and DCA peak sets point to completely different biology for the same
  cell type (poor RNA-ATAC concordance — may indicate label transfer error)
- NK cells present in RNA cell_type_vote but absent from ATAC-derived cell types
  (label transfer absorption into cytotoxic T clusters)

For each flag:
- Describe the inconsistency
- Cite the specific values, genes, or peaks involved
- Assign severity: info / warning / critical
- Suggest a corrective action

Conclude with an overall cross-modal coherence assessment:
**pass** / **review recommended** / **issues found**

---

### Task 8 — Downstream Analysis Suggestions

Based on the cell types found, the DEGs, the DCA peaks, the joint embedding,
and the biological question stated in Section 2, suggest the most valuable
next analysis steps specific to multiome data.

For each suggestion provide:
- **Step name**
- **Rationale** — why this is the logical next step given these multiome results
- **Recommended tool** — specific tool or package
- **Expected output** — what the analyst will learn
- **Multiome advantage** — why this step benefits from having both RNA and ATAC
  in the same cells (vs. doing it on RNA alone)

Rank suggestions by priority (1 = highest).

---

## SECTION 5 — GROUNDEDNESS RULE

Every factual claim in `multiome_report_review.md` must be grounded by one of:

| Ground type | Example |
|-------------|---------|
| Metric from report | "95.3% pass rate" / "resolution = 0.5" / "13 Leiden clusters" / "logFC = 3.66" |
| Specific gene or peak | "GNLY upregulated logFC = 4.23 in Tem/Temra" / "chr17-82126430-82127342 logFC = 3.80" |
| PMID | "PMID:12345678" |
| Widely accepted and no search available | "literature-supported" (use sparingly) |

**Never invent numbers, gene names, or peak coordinates.**
If a value is not present in the HTML report, write "not reported" — do not estimate.
Aim for ≥ 85% of factual sentences to be explicitly grounded.

**Specific rules for DCA peaks:**
- Never convert a coordinate-only peak (chr:start-end) to a gene name unless
  you are confident in the hg38 locus. Write "locus unresolved" if uncertain.
- Never state that a peak "regulates" a gene without peak annotation confirming
  the overlap.

---

## SECTION 6 — OUTPUT TEMPLATE

Copy everything below this line into your response and fill it in completely.
Do not truncate any section. Return as a single markdown document.

---

```markdown
# Multiome AI Interpretation Report
**Dataset:** [dataset name or GEO accession]
**Date:** [today's date]
**Reviewer:** Manual Review Mode (Multiome D1)
**Model used:** [e.g. Claude Sonnet 4.6 / GPT-4o / Claude Opus 4]
**Tissue:** [from Section 2]
**Disease:** [from Section 2]
**Modalities:** RNA + ATAC (multiome)

---

## 1. Study Overview

- **Tissue:** [value]
- **Disease:** [value]
- **Species:** [value]
- **Conditions:** [value]
- **N cells (post-QC):** [value]
- **N donors:** [value]
- **Batch key:** [value or "null (single batch)"]
- **Modalities:** RNA + ATAC (multiome)
- **Biological question:** [stated by analyst / inferred by reviewer]
- **Known limitations:** [verbatim from report warning banners, or "none reported"]

---

## 2. ATAC QC Assessment

- **filter_cells applied:** [True / False] — [rationale if False]
- **Pass rate:** [%] — [within normal / above typical / below typical] for [tissue]
- **Thresholds:** min_peaks = [x], min_peak_counts = [y], max_nucleosome_signal = [z]
- **Doublet detection:** [run / not reported]
- **Silhouette score:** [value / not reported]
- **Overall QC assessment:** [appropriate / review recommended / issues found]

| Severity | Flag | Suggestion |
|----------|------|------------|
| [info/warning/critical] | [description] | [action] |

---

## 3. LSI Reduction Assessment

- **Component 1 dropped:** [Yes / No] — [appropriate / flag if No]
- **Cumulative variance (components used):** [%] — [normal for ATAC / below expected]
- **LSI components used:** [N]
- **Resolution:** [x] — [appropriate / too coarse / too fine] for [N cells] cells
- **Recommended range:** [x–y] for this tissue and cell count
- **N Leiden clusters:** [n] — [reasonable / too many / too few] for [tissue]
- **UMAP structure:** [comment on separation, batch effects, ground truth alignment]

---

## 4. Annotation and Label Transfer Assessment

### 4a. Cluster-level label transfer

| Cluster | N cells | % total | Majority RNA label | Reviewer assessment | Chromatin plausibility | Flag |
|---------|---------|---------|-------------------|---------------------|----------------------|------|
| 0 | | | | | | |

### 4b. NK cell discrepancy
[Present / Not applicable]
- **RNA cell_type_vote NK count:** [n cells]
- **ATAC atac_celltype NK count:** [0 / n cells]
- **Most likely absorbing cluster:** [cluster ID and label]
- **Suggested resolution:** [action]

### 4c. Peak annotation source
- **Source:** [coordinate_fallback / real gene names]
- **Impact on interpretation:** [statement]

### 4d. Tissue-specific chromatin marker sets

| Cell type | Canonical open loci | Context-specific (hematopoiesis) | Distinguishing loci | Counter-markers | Confidence |
|-----------|--------------------|---------------------------------|---------------------|-----------------|------------|
| | | | | | |

---

## 5. DEG and DCA Interpretation

### 5a. RNA DEG Validation

**[Cell type name]**

| Gene | logFC | adj p-value | Classification | Rationale |
|------|-------|-------------|---------------|-----------|
| | | | Expected / Unexpected / Artefact | |

*(repeat for each cell type)*

**Compositional artefact check:**
[Pass / Warning per cell type]

**Discovery highlights:**
- [Gene]: [why noteworthy — grounded with logFC or PMID]

### 5b. DCA Peak Interpretation

**Peak annotation source:** [coordinate_fallback / real gene names]

| Cell type | Peak | logFC | Locus / Gene | RNA–chromatin concordance |
|-----------|------|-------|-------------|--------------------------|
| | | | [gene name / locus unresolved] | [concordant / discordant / unknown] |

### 5c. Literature Links

| Gene | PMID | Title | Context (one sentence) |
|------|------|-------|------------------------|
| | | | |

*PMIDs and titles only. No abstract text.*

[If no search tool was available, list manual PubMed search strings here instead.]

---

## 6. Cross-Modal Coherence Review

| Severity | Modalities involved | Flag | Suggestion |
|----------|---------------------|------|------------|
| [info/warning/critical] | [e.g. Annotation + DEG / RNA + ATAC] | [description] | [action] |

**Overall cross-modal coherence:** [pass / review recommended / issues found]

---

## 7. Downstream Suggestions

| Priority | Step | Rationale | Recommended tool | Expected output | Multiome advantage |
|----------|------|-----------|-----------------|----------------|--------------------|
| 1 | | | | | |
| 2 | | | | | |
| 3 | | | | | |

---

## 8. Summary

**Key findings:**
1. [finding — grounded with metric / gene / PMID]
2. [finding — grounded with metric / gene / PMID]
3. [finding — grounded with metric / gene / PMID]

**Open questions:**
- [question raised by the analysis]

**Suggested validation experiments:**
- [wet-lab or computational experiment to validate a key finding]

---

## Appendix: Multiome Report Figures Checklist

Check each item against the HTML report and mark ✓ present or ✗ missing.

| # | Figure / element | Status | Note |
|---|-----------------|--------|------|
| 1 | ATAC QC metric distributions (peaks per cell, total counts, nucleosome signal, reads-in-peaks fraction) | | |
| 2 | Scrublet doublet score distribution for ATAC | | |
| 3 | LSI variance explained bar chart with component 1 highlighted | | |
| 4 | UMAP figures coloured by: batch, ground truth, Leiden clusters | | |
| 5 | Label transfer table (cluster → majority RNA label, N cells, % total) | | |
| 6 | Gene activity heatmap (top variable genes × clusters, z-scored) | | |
| 7 | MultiVI joint embedding UMAP coloured by batch and cell type | | |
| 8 | Pre-integration LSI UMAP for comparison with joint embedding | | |
| 9 | MultiVI latent dimension violin plots | | |
| 10 | RNA DEG table (top 5 per cell type, logFC, adj p-value) | | |
| 11 | DCA peak table (top 5 per cell type, logFC, adj p-value, chr:start-end format) | | |
| 12 | Peak annotation source reported (coordinate_fallback or real gene names) | | |
| 13 | Silhouette score for ATAC clustering | | |
```

---

*OmicSage Multiome D1 Manual Review Mode — MULTIOME_MASTER_PROMPT.md v1.0*
*No Python. No infrastructure. Attach + send.*
