## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 0 — Shared AI infrastructure built and tested.
                      All 4 foundation modules passing. Encoding fix applied
                      to test_phase0_structure.py (cp1252 → utf-8).
File last worked on: tests/test_ai_infrastructure.py

## Today's Goal
Phase 3 — Session 1: A1 — Pipeline Advisor.

Build the first feature module. Nothing else.

Files to create:
  - ai/skills/pipeline_advisor.yaml   ← skill YAML (prompts live here, not in .py)
  - ai/pipeline_advisor.py            ← feature module
  - tests/test_pipeline_advisor.py    ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~251 passed, 1-2 skipped
```

## Step 2 — Read the pipeline advisor spec in PHASE3_PLAN.md
Section: "Session 1 — A1: Pipeline Advisor"
Read it fully before writing any code.

## Step 3 — Build ai/skills/pipeline_advisor.yaml
Follows the cluster_annotator.yaml pattern exactly.
Inputs the skill needs:
  - tissue, disease_context, experiment_design, biological_question
  - n_cells, n_genes, modalities (list), n_batches, n_donors, n_conditions
Output schema (JSON):
  - recommended_steps: list of {step_name, priority, rationale}
    priority values: "required" | "recommended" | "optional"
  - inferred_biological_question: str | null
  - warnings: list[str]
  - reasoning: str

## Step 4 — Build ai/pipeline_advisor.py

### Public API
```python
from dataclasses import dataclass, field
from ai._base import AiResult

@dataclass
class StepRecommendation:
    step_name: str
    priority: str   # required | recommended | optional
    rationale: str

@dataclass
class PipelineAdvice(AiResult):
    recommended_steps: list[StepRecommendation] = field(default_factory=list)
    inferred_biological_question: str | None = None
    warnings: list[str] = field(default_factory=list)

def run(
    adata_or_mdata,         # AnnData or MuData — reads properties, never modifies
    config: dict,
    study_context: dict,    # loaded from config/runs/<dataset>/study_context.yaml
    *,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> PipelineAdvice | None:
    ...
```

### Rule-based pre-checks (run BEFORE the LLM call — fast, no API cost)
These fire regardless of AI being enabled:
  - n_batches > 1 → add warning if batch_key not set in config
  - n_donors > 2 AND n_conditions > 1 → recommend pseudobulk over Wilcoxon
  - modalities includes "ADT" → recommend WNN (flag as Phase 6, not yet built)
  - n_cells < 500 → warn about unreliable clustering and doublet detection

The LLM adds rationale and literature context on top of these rule outputs.

### Usage pattern (identical to all feature modules)
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata_or_mdata, config, study_context, *, log_dir, runtime_ai=True):
    try:
        check_ai_enabled(config, module="pipeline_advisor", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
    # ... rest of the module
```

### study_context loading helper
```python
# ai/pipeline_advisor.py — at module level
import yaml
from pathlib import Path

def load_study_context(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}
```

## Step 5 — Tests (all without a real API key)
Mock pattern — same as test_ai_infrastructure.py:
  monkeypatch _build_conversation to return a MagicMock that returns
  a valid JSON string matching the pipeline_advisor output schema.

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns PipelineAdvice when mock LLM returns valid JSON
  - recommended_steps list is non-empty
  - StepRecommendation priority values are all valid ("required"|"recommended"|"optional")
  - Batch warning fires when n_batches > 1 and batch_key not in config
  - Pseudobulk recommendation fires when n_donors > 2 and n_conditions > 1
  - inferred_biological_question populated when blank in study_context
  - n_cells < 500 warning fires
  - Audit log written to log_dir after successful call
  - Graceful handling of malformed LLM JSON (returns None, logs warning)
  - AiResult base fields (timestamp, model, provider) populated correctly

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ← TODAY — A1: Pipeline advisor
Session 2  — A2: Clustering advisor (first PubMed RAG use)
Session 3  — B1: Cluster annotator
Session 4  — B2: DEG validator + literature linker
Session 5  — B3: Coherence reviewer (build analysis_summary.json here)
Session 6  — A3: Downstream analysis suggester
Session 7  — C1: Narrative generator
Session 8  — C2: Full report + PowerPoint
Session 9  — Milestone validation: groundedness test + end-to-end GSE166635
```
Full details: .dev_memory/PHASE3_PLAN.md

## AI Layer Design Principles (carry into every session)
- Skills are the prompts — no raw strings in Python, ever
- AI is optional at three levels: global / per-module / per-step runtime flag
- A disabled module always returns None silently — never errors, never logs
- Every LLM call is audit-logged — no exceptions
- Graceful degradation — a failed AI module returns None, never crashes
- Context scales quality — more study_context = better AI output
- analysis_summary.json is load-bearing — design carefully in Session 5
- Copyright enforced at report layer — PMIDs + titles only, no abstract text

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()
- ai/skills/cluster_annotator.yaml  ← reference pattern for all skill YAMLs

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
- test_phase0_structure.py: schema fixture uses encoding="utf-8" (cp1252 fix applied)
- biochatter==0.14.2 pinned — do NOT upgrade without re-running all AI tests

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2
