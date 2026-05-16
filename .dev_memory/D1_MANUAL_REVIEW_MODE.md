# OmicSage — D1 Manual Review Mode: Complete Plan
> Written: May 2026
> Status: Approved design — ready to implement in one session
> Position in pipeline: final step, after C1 produces the combined HTML report
> Hard rule: one session = write MASTER_PROMPT.md + test on GSE166635

---

## Context and Motivation

The automated AI pipeline (A2 → C1) requires a running LLM API (Claude, OpenAI, or
local Ollama). D1 is a complete fallback for the entire pipeline AI — not just a
reviewer of it. It reads the combined HTML report (which always exists regardless of
AI mode) and produces `report_review.md`, a structured document that mirrors what
the automated pipeline would have produced.

**Two usage modes:**

| Situation | What you do |
|-----------|-------------|
| Pipeline AI ran successfully | Run D1 optionally for enrichment and independent validation |
| Pipeline AI disabled or failed | Run D1 — it produces everything A2–C1 would have produced |

**Why D1 is also a validation tool:** Run the pipeline AI, then run D1 independently,
compare outputs. Significant disagreement on annotation or coherence flags a problem
in the automated pipeline.

---

## Revised Phase 3 Pipeline Shape

```
Automated pipeline AI (runs with --ai flag):
  A2  clustering_advisor     → clustering_advice.json
  B1  cluster_annotator      → adata.obs[ai_cell_type, ai_confidence]
  A1  pipeline_advisor       → pipeline_advice.json
  B2  deg_validator          → deg_validation.json  (PMIDs + titles only)
  B3  coherence_reviewer     → analysis_summary.json + coherence_flags.json
  A3  downstream_suggester   → NEXT_STEPS.md
  C1  narrative_generator    → ai_narrative.md + AI tab injected into combined HTML
                                ↑
                reads: all JSON outputs above + analysis_summary.json

Manual fallback — D1 (zero cost, zero infrastructure):
  Input:  combined HTML report  (always exists)
  Output: report_review.md      (mirrors C1 AI tab content)
  Covers: everything A2–C1 would have produced
```

**C2 (report_writer + pptx) is removed.** The combined HTML report is the deliverable.
The AI tab (written by C1) is the automated AI interpretation layer.
`report_review.md` (written by D1) is the manual AI interpretation layer.

---

## Files to Create

```
ai/
└── manual_review/
    └── MASTER_PROMPT.md    ← single self-contained file (fill study context, attach HTML, send)

reports/<dataset>/
└── report_review.md        ← saved manually after the D1 session
```

**That is the entire D1 implementation.** No Python, no tests, no CI.
`MASTER_PROMPT.md` contains everything: role, study context block, report structure
guide, all 10 tasks written in full, PubMed instructions, and the output template
embedded at the end. You attach it + the HTML to a chatbot, send once, get
`report_review.md` back.

---

## MASTER_PROMPT.md — Structure

Single file. Self-contained. Six sections written once and maintained as the
report format evolves. You fill in Section 2 each time (2 minutes), attach the
HTML report, and send. The model runs all 10 tasks in sequence and returns the
complete `report_review.md` in one response.

### Section 1 — Role

```
You are a senior computational biologist and bioinformatician reviewing a completed
single-cell RNA-seq analysis. You have deep expertise in QC thresholds, clustering,
cell type annotation, differential expression, pathway analysis, and translating
results into biological meaning. You are reviewing a full pipeline report and will
produce a structured interpretation document that mirrors what an automated AI
pipeline would have produced. Run all tasks below in sequence without stopping.
Return the complete report_review.md at the end.
```

### Section 2 — Study Context Block (fill in each time)

```
Tissue: [e.g. liver]
Disease: [e.g. hepatocellular carcinoma / null if healthy]
Species: [human / mouse / other]
Conditions: [e.g. tumour vs adjacent normal]
N cells (post-QC): [from report]
N donors: [from report]
Batch key: [from report, or null if single batch]
Biological question: [your question, or leave blank for AI to infer]
Known biology: [optional — e.g. "Wang et al. 2025 found immunosuppressive macrophages"]
```

### Section 3 — Report Structure Guide

```
The attached HTML report contains the following tabs:
- QC: per-cell metrics (n_genes, MT%, n_counts), threshold decisions, pass/fail
      flags, cells removed percentage
- Normalization: provenance table (target_sum, HVG flavor, batch_key, n_HVGs)
- Dimensionality reduction: PCs computed, PCs used, cumulative variance, UMAP figure
- Clustering: resolution used, n_clusters, silhouette score, marker gene table
              per cluster (top 10 genes, log2FC)
- Annotation: SingleR votes per cluster, confidence scores, UMAP colored by cell type
- DEG: top DEGs per comparison (gene, log2FC, adj p-value), volcano plots
- GSEA: significant pathways per group (GO, KEGG, Reactome), adj p-values
- AI tab: present if pipeline AI ran — ignore for this review, produce independently
```

### Section 4 — Task List (run all in sequence, no stopping)

```
Task 1  — QC coherence check
          Read the QC tab. Were thresholds appropriate for this tissue and disease?
          Is the % cells removed within normal range? Flag any concerns.

Task 2  — Clustering advice
          Read the clustering tab. Is the resolution appropriate? Is the number of
          clusters reasonable for this tissue? Identify sub-clustering candidates.
          Cite comparable published studies if you know them.

Task 3  — Cluster annotation
          Read the marker gene table. Predict the cell type for each cluster.
          Assign confidence (high / medium / low). List supporting markers and
          alternative interpretations for ambiguous clusters.

Task 4  — Tissue-specific marker sets
          For each cell type identified in Task 3, generate a marker validation set
          specific to [tissue] + [disease]. For each cell type provide:
            1. Canonical markers (universally accepted)
            2. Context-specific markers (from published [tissue]+[disease] studies)
            3. Distinguishing markers (separates this subtype from similar types here)
            4. Counter-markers (high expression would change the annotation)
          Add a confidence note per cell type: high / medium / low.

Task 5  — DEG validation
          Read the DEG tab. For each comparison, classify genes as:
            - Expected (consistent with known biology of this cell type + disease)
            - Unexpected (novel or potentially artefactual — worth investigating)
          Flag discovery highlights: genes not previously reported in this context.

Task 6  — Literature linking
          For each unexpected gene from Task 5 and each key expected gene, search
          for supporting literature.

          If you have a PubMed or bioRxiv search tool available, use it directly
          and report: gene → PMID + title + one-sentence context.

          If you do not have a search tool, provide the exact PubMed search string
          for each gene so I can run them manually and paste back results.

          PMIDs and titles only — never include abstract text.

Task 7  — GSEA coherence
          Read the GSEA tab alongside the cell type annotations from Task 3.
          Are the enriched pathways consistent with the annotated cell types?
          Flag unexpected enrichments with severity: info / warning / critical.

Task 8  — Cross-module coherence review
          Consider everything from Tasks 1–7 together. Flag inconsistencies across
          modules (e.g. annotation contradicts markers, GSEA contradicts cell types,
          QC removal unusually high). Rate each flag: info / warning / critical.
          Give an overall coherence assessment: pass / review recommended / issues found.

Task 9  — Downstream analysis suggestions
          Based on the cell types found, DEGs, and biological question, suggest the
          most valuable next analysis steps. For each: step name, rationale, recommended
          tool, expected output. Rank by priority.

Task 10 — Write report_review.md
          Assemble all findings from Tasks 1–9 into the output format specified in
          Section 6. Every factual claim must cite a metric value from the report,
          a specific gene name, or a PMID from Task 6.
          Never invent numbers. If a value is not in the report, say "not reported".
```

### Section 5 — Groundedness Rule

```
Every factual claim in the output must be grounded by one of:
  - A metric value from the report (e.g. "12.3% cells removed")
  - A specific gene name (e.g. "PDCD1 upregulated log2FC=2.1 in cluster 4")
  - A PMID from the literature search (e.g. "PMID:12345678")
  - The phrase "literature-supported" if widely accepted and no search was run

Never invent numbers or gene names. If something is not in the report, say
"not reported". Aim for ≥85% of factual sentences to be explicitly grounded.
```

### Section 6 — Output Template (embedded)

The model writes this directly. You copy the output and save it as
`reports/<dataset>/report_review.md`.

```markdown
# AI Interpretation Report
**Dataset:** [name]
**Date:** [date]
**Reviewer:** Manual Review Mode (D1)
**Model used:** [e.g. Claude Sonnet 4.6 / GPT-4o]

---

## 1. Study Overview
- Tissue, disease, conditions, n_cells, n_donors
- Biological question: [stated / inferred]

---

## 2. QC Assessment
- Thresholds: min_genes=[x], max_MT%=[y]
- Cells removed: [%] — [within normal / above typical / below typical] for [tissue]
- Assessment: [appropriate / review recommended]

| Severity | Flag | Suggestion |
|----------|------|------------|

---

## 3. Clustering Assessment
- Resolution: [x] — [appropriate / too coarse / too fine]
- Recommended range: [x–y]
- N clusters: [n] — [reasonable / too many / too few] for [tissue + disease]
- Sub-clustering candidates: [cluster IDs or "none"]
- Literature context: [comparable studies + resolutions used]

---

## 4. Cell Type Annotation

| Cluster | N cells | Predicted type | Confidence | Supporting markers | Alternatives |
|---------|---------|---------------|------------|-------------------|--------------|

**Tissue-specific marker sets:**

| Cell type | Canonical | [Tissue+Disease]-specific | Distinguishing | Counter-markers | Confidence |
|-----------|-----------|--------------------------|----------------|-----------------|------------|

---

## 5. DEG Validation

**Comparison: [name]**
- Expected genes: [list]
- Unexpected genes: [list + rationale]
- Discovery highlights: [list]

---

## 6. Literature Links

| Gene | PMID | Title | Context |
|------|------|-------|---------|

*PMIDs and titles only. No abstract text.*

---

## 7. GSEA Coherence
- Consistent pathways: [list]

| Unexpected pathway | Group | Adj p-value | Severity |
|--------------------|-------|-------------|----------|

---

## 8. Cross-Module Coherence

| Category | Severity | Flag | Suggestion |
|----------|----------|------|------------|

**Overall coherence:** [pass / review recommended / issues found]

---

## 9. Downstream Suggestions

| Priority | Step | Rationale | Tool | Expected output |
|----------|------|-----------|------|----------------|

---

## 10. Summary

**Key findings:**
1. [finding — grounded]
2. [finding — grounded]
3. [finding — grounded]

**Open questions:**
- [question]

**Suggested validation experiments:**
- [experiment]
```

---

## What D1 Covers That the Automated Pipeline Does Not

Two genuine additions beyond mirroring the pipeline:

| Capability | Pipeline AI | D1 |
|-----------|-------------|-----|
| QC coherence | ✅ A1 | ✅ Task 1 |
| Clustering advice | ✅ A2 | ✅ Task 2 |
| Cluster annotation | ✅ B1 | ✅ Task 3 |
| Tissue-specific marker sets | ❌ | ✅ Task 4 (new) |
| DEG validation | ✅ B2 | ✅ Task 5 |
| Literature linking | ✅ B2 (RAG) | ✅ Task 6 (auto if tool available, manual loop fallback) |
| GSEA coherence | ❌ | ✅ Task 7 (new) |
| Cross-module coherence | ✅ B3 | ✅ Task 8 |
| Downstream suggestions | ✅ A3 | ✅ Task 9 |
| Written interpretation | ✅ C1 | ✅ Task 10 |
| Image analysis (UMAP, plots) | ❌ (text only) | ✅ Claude.ai sees embedded images |

---

## Implementation Session Plan

**Single session. No code. No tests. No CI.**
Deliverable: `MASTER_PROMPT.md` written, tested on GSE166635, `report_review.md` saved.

```
Step 1 — Write MASTER_PROMPT.md                               (45 min)
           Fill in all 6 sections. Task instructions must be
           detailed enough that the model runs all 10 tasks
           without steering.

Step 2 — Test: attach MASTER_PROMPT.md + 00_combined_report.html
           to Claude.ai, fill study context, send once.        (45 min)
           Input:  reports/GSE166635/00_combined_report.html
           Output: reports/GSE166635/report_review.md

Step 3 — Review output quality                                 (15 min)
           - Are all 10 sections present?
           - Are factual claims grounded (metric / gene / PMID)?
           - Did PubMed search run automatically or fall back?
           - Did the marker sets go beyond what B1 would produce?

Step 4 — Iterate MASTER_PROMPT.md based on what needed steering (15 min)

Step 5 — Update PHASE3_PLAN.md and .dev_memory/               (10 min)
```

Total: ~2 hours. At the end you have a tested prompt and a real scientific
interpretation of the HCC analysis ready to share.

---

## Relationship to PHASE3_PLAN.md

### Modules removed
- **C2 (report_writer + pptx):** removed. Combined HTML report is the deliverable.

### Modules unchanged
- A1, A2, B1, B2, B3, A3, C1: as specified in PHASE3_PLAN.md

### C1 responsibility clarified
C1 (narrative_generator) is the assembly point for all upstream AI outputs.
It reads all JSON outputs (clustering_advice, deg_validation, coherence_flags, etc.)
plus analysis_summary.json from B3, and produces:
- `ai_narrative.md` — standalone markdown
- The AI tab injected into the combined HTML report

### D1 added
D1 is a complete fallback for the entire AI pipeline (A2–C1).
It lives in `ai/manual_review/` as a single prompt file: `MASTER_PROMPT.md`.
It is not part of the automated pipeline and has no Python code, no tests, no CI.
It is documented in the README under "AI Features — Manual Review Mode".

**How to use D1:**
1. Fill in the study context block in `MASTER_PROMPT.md` (2 min)
2. Open Claude.ai or ChatGPT
3. Attach `MASTER_PROMPT.md` and `reports/<dataset>/00_combined_report.html`
4. Send once
5. Copy output → save as `reports/<dataset>/report_review.md`

---

## Build Order (revised)

```
Session 0 — Infrastructure: _llm_client, _audit_log, _skill_loader, _config_gate, _base
Session 1 — A2: Clustering advisor
Session 2 — B1: Cluster annotator
Session 3 — A1: Pipeline advisor
Session 4 — B2: DEG validator + literature linker  ← completed
Session 5 — B3: Coherence reviewer (build analysis_summary.json here)
Session 6 — A3: Downstream analysis suggester
Session 7 — C1: Narrative generator + AI tab injection into combined HTML
Session 8 — D1: Manual Review Mode (prompt files + test on GSE166635)
Session 9 — Milestone validation: groundedness test + end-to-end GSE166635
```

Total: 9 sessions (down from 10 — C2 removed).

---

*Plan written: May 2026*
*Replaces C2 section of PHASE3_PLAN.md*
*D1 is new — not in original PHASE3_PLAN.md*
