# OmicSage — Architectural Decisions

## D001 — Nextflow DSL2
Industry standard, HPC/cloud ready, nf-core module reuse.

## D002 — BioChatter as AI middleware
Already handles RAG/KG/LLM routing. Nature Biotech paper = citable.
Pin to v0.8.5. Wrap in ai/biochatter_client.py to isolate from API changes.

## D003 — Quarto + python-pptx
Quarto: HTML+PDF from same source, parameterized rendering.
python-pptx: wet-lab groups live in PowerPoint, editable slides > static PDF.

## D004 — Streamlit for UI
Fastest path to biologist-facing interface. Pure Python. Replace with React for SaaS.

## D005 — Docker + Singularity
Docker for dev/CI. Singularity for HPC (most clusters block Docker).

## D006 — Separate Python and R images
Single image = 8-10GB. Separate = ~4GB each. Nextflow routes via labels.

## D007 — MIT License
Maximum adoption. Pharma/biotech can use commercially. Citable.

## D008 — AI-free mode mandatory
Full pipeline runs without any API key. Competitive differentiator.

## D009 — GSE166635 as benchmark
Known expected results. Multiple tumor/normal samples for batch correction testing.

## D010 — Start from count matrices
nf-core/scrnaseq handles upstream. OmicSage value is downstream + AI + reports.

## D011 — Click for CLI
Cleaner API than argparse. Natural sub-command hierarchy.

## D012 — GHCR for Docker registry
Free for public repos. Integrated via GITHUB_TOKEN. Co-located with code.

ClawBio (clawbio.ai): open-source bioinformatics skill library.
Potential future integration as optional skill backends (Phase 4-5).
Potential OpenClaw marketplace listing (Phase 7+).
Not a dependency for Phase 1-3.
Monitor: claw-spatial, claw-gwas, scRNA Orchestrator updates.

## D013 — Stop AI Pipeline Development; Keep Manual Version (2026-05-17)
After completing all 10 sessions of Phase 3 (AI layer) with 466+ tests passing, we evaluated
the trade-off between AI pipeline complexity and project stage priorities and decided:

**What we decided:**
- Stop all AI pipeline development immediately after Phase 3 milestone validation
- The manual pipeline (ai_features: false) is the sole primary path going forward
- All AI layer code stays in the repo — intact, tested, and importable — but is not extended
- No new AI modules will be added until Phase 8 (benchmarking) or later, if at all

**Why:**
- Phase 1 (annotation) is still incomplete — SingleR-based annotation not yet built
- The AI layer adds significant complexity with diminishing returns at this project stage
- The manual pipeline already reproduces Wang et al. 2025 HCC findings correctly
- For job interviews and freelance (Month 1-4 goals), a clean manual pipeline > a buggy AI one
- The architecture still allows ai_features: true — the option is preserved, not removed

**What this means in practice:**
- Default config: ai_features: false
- ai/ directory: frozen at current state, no new modules
- Next steps: complete Phase 1 annotation → Phase 4 scATAC → Phase 5 Spatial

---

Report Reviewer (D1) — text-only analysis, no figure inspection
Decision: report_reviewer.py strips HTML to plain text before the LLM call. Images referenced in the report are not extracted or passed to the model. The reviewer therefore evaluates narrative quality, methods completeness, and whether conclusions are supported by described results — it cannot visually inspect figures. This was a deliberate scope decision for Phase 3: the module is a scientific writing reviewer, not a figure quality checker. True figure analysis would require a vision-capable model call (Claude vision API or GPT-4o with images), per-figure base64 extraction from the HTML before tag stripping, and a separate prompt designed around visual interpretation. The complexity and token cost do not justify inclusion at this stage. If figure review becomes a priority in a later phase, it should be a standalone module (ai/figure_reviewer.py) rather than an extension of this one, so the two concerns remain independently testable and independently disableable.

---
 
## D014 — AI layer paused (2026-05-17)
**Decision:** Stop extending the AI pipeline. Manual pipeline is the primary path.
**Rationale:** After completing all 7 AI modules (466+ tests), the added complexity
vs. value at this stage was evaluated. The manual review mode (MASTER_PROMPT.md) delivers
better groundedness and transparency for now. AI layer stays in repo, intact and tested.
**Impact:** Phase 3 is complete but frozen. No new AI modules until further notice.
 
---
 
## D015 — Modality roadmap corrected: ATAC → CITE-seq (2026-05-21)
**Decision:** Phase 4 is CITE-seq (RNA + protein/ADT), not scATAC-seq standalone.
**Rationale:** The original roadmap listed "Phase 4 — scATAC-seq Module" but this
was an error. The benchmark dataset (GSE194122) is CITE-seq, not ATAC. scATAC-seq
as a standalone modality is not a priority. ATAC is best analysed jointly with RNA
in the Multiome context (Phase 6). CITE-seq is the logical next step because the
benchmark data is already processed through QC and mdata["adt"] is available.
 
**Corrected roadmap:**
  Phase 4: CITE-seq (RNA + protein/ADT)
  Phase 5: Multiome Integration (RNA + ATAC jointly)
  Phase 6: Spatial (Visium, MERFISH, Xenium)
  scATAC-seq standalone: deferred, covered within Phase 5
 
**Files to update:** OMICSAGE_PROJECT_CONTEXT.md Section 9 checkboxes.
 
---
 
## D016 — Report improvements: separate consensus score from reviewer confidence (2026-05-21)
**Decision:** The annotation report now has two separate confidence columns:
  - `Consensus score`: weighted fraction of automated methods agreeing (0.0–1.0)
  - `Confidence`: qualitative label (High / Medium / Low) from the consensus score
The MASTER_PROMPT.md reviewer task further splits these into:
  - `Consensus score`: read directly from the report (numeric)
  - `Reviewer confidence`: reviewer's own marker-based judgement (high/medium/low)
**Rationale:** The original single "Confidence" column conflated method agreement
(a statistical quantity) with AI self-assessment (poorly calibrated). A biologist
reading "Medium" could not tell which was meant. Review of GSE194122 confirmed this
caused ambiguity in the cluster 9 (pDC) and cluster 10 (CD4+ T) annotations.
---
 
## D017 — exclude_gene_prefixes: default empty, opt-in per dataset (2026-05-21)
**Decision:** `exclude_gene_prefixes` defaults to `[]` (no filtering) in deg.py,
gsea.py, and pseudobulk_deg.py. Filtering is opt-in via the YAML config.
**Rationale:** The default should preserve all data. Filtering is dataset-specific —
the RPL/RPS artefact is driven by the erythroid-heavy BMMC background in GSE194122
and would not apply to e.g. a pure T cell dataset. GSE194122.yaml sets
`["RPL", "RPS", "MT-"]` for deg and gsea, and additionally `["HBB", "HBA1", "HBA2",
"HBD", "AHSP"]` for pseudobulk where erythroid markers dominate downregulated lists.
 
---
 
## D018 — data_report.py: sample and batch as separate detected columns (2026-05-21)
**Decision:** `data_report.py` now independently detects a `sample` column and a
`batch` column, showing each as a separate metric card. When neither is found, a
single "Single batch" card is shown rather than rendering a broken card.
**Rationale:** In many datasets `sample` and `batch` are different things (e.g.
12 donors but 4 processing sites). Conflating them into one card lost information.
The no-batch case previously rendered `N/A` and `Batches (none)` which was unclear.