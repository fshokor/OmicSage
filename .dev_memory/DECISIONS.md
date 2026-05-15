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

Report Reviewer (D1) — text-only analysis, no figure inspection
Decision: report_reviewer.py strips HTML to plain text before the LLM call. Images referenced in the report are not extracted or passed to the model. The reviewer therefore evaluates narrative quality, methods completeness, and whether conclusions are supported by described results — it cannot visually inspect figures. This was a deliberate scope decision for Phase 3: the module is a scientific writing reviewer, not a figure quality checker. True figure analysis would require a vision-capable model call (Claude vision API or GPT-4o with images), per-figure base64 extraction from the HTML before tag stripping, and a separate prompt designed around visual interpretation. The complexity and token cost do not justify inclusion at this stage. If figure review becomes a priority in a later phase, it should be a standalone module (ai/figure_reviewer.py) rather than an extension of this one, so the two concerns remain independently testable and independently disableable.