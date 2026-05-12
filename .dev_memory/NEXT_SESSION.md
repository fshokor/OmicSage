## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Phase 2 MILESTONE — combined tabbed HTML report
                      (00_combined_report.html) generating automatically
                      at end of every pipeline run. Tested on GSE166635.
                      Phase 2 is complete.
File last worked on: reports/combined_report.py + run_pipeline.py

## Today's Goal
Phase 3 kickoff — BioChatter integration.

The question to answer at session start:
  What is the smallest useful AI feature to build first?

Recommended: QC threshold suggester — reads QC metrics from the h5ad,
calls an LLM, returns suggested min_genes / max_genes / max_mt_pct with
reasoning. Standalone function, no UI needed yet. One session, one function.

ONE goal — decide the Phase 3 build order and implement the first AI feature.
Do not build more than one AI feature per session.

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~231 passed, 1-2 skipped
```

## Step 2 — Verify Phase 2 milestone still works
```bash
python run_pipeline.py --config config/runs/GSE194122.yaml --step qc
# Should cache and print: [combined_report] N tabs → reports/GSE194122/00_combined_report.html
```

## Step 3 — Decide Phase 3 build order before writing any code
Suggested order:
  1. QC threshold suggester (ai/threshold_suggester.py)
     → reads QC metrics, returns suggested thresholds + reasoning
  2. Cluster interpreter (ai/cluster_interpreter.py)
     → marker genes → LLM → cell type label + evidence
  3. PubMed RAG tied to DEG results (via BioChatter)
  4. Narrative generator for combined report
  5. Multi-LLM support (Claude / GPT-4o / local Ollama)
  Milestone: AI narrative groundedness score > 0.85

## Phase 3 Notes
- Use BioChatter as the AI middleware — do NOT build LLM routing from scratch
- All LLM calls must be audit-logged to logs/llm/ in JSONL format
- ai_features: false in config must disable the AI layer completely
- Pipeline must still run and reports must still generate without any API key
- Start with Claude API (ANTHROPIC_API_KEY) — add Ollama support later

## Known Issues Carried Forward
- Always use `python -m pytest` not bare `pytest`
- Always `conda activate omicsage` before running anything
- Use `2>&1 | tee logs/<dataset>.log` to capture all output
- run_pipeline.py must be at repo root
- obs['cell_type_vote'] is the consensus column (not obs['cell_type'])
- rpy2 NOT used — do not attempt R integration
- gseapy.enrichr returns NO Overlap column in v1.1.3
- harmonypy ho.Z_corr is (n_cells, n_pcs) — do NOT transpose
- obs[batch_key] must be cast .astype(str) before pandas operations
- pseudobulk requires layers['counts'] (raw integers) — NOT layers['logcounts']
- GSE166635 pseudobulk is DISABLED — only 2 samples

## Files Modified Last Session
- reports/combined_report.py               ← new combined tabbed report generator
- run_pipeline.py                          ← wired in combined report at end of main()
- .dev_memory/CURRENT_STATUS.md            ← updated
- .dev_memory/NEXT_SESSION.md              ← updated (this file)
- .dev_memory/PROGRESS.md                  ← updated

## Verify Last Session Works
```bash
cd ~/OmicSage
conda activate omicsage
python run_pipeline.py --config config/runs/GSE166635.yaml --step qc
# Should cache and auto-generate reports/GSE166635/00_combined_report.html
```

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2
New dependency needed: biochatter (pip install biochatter)
