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