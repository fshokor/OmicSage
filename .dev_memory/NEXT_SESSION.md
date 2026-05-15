## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 2 — A2 Clustering Advisor built and tested.
                      22 tests passing. Fixed: call_llm returns str only
                      (not 4-tuple), write_audit_record uses token_usage=
                      dict (not separate kwargs), _advice_to_dict needed
                      for JSON-serializing tuple + nested dataclasses.
File last worked on: tests/test_clustering_advisor.py

## Today's Goal
Phase 3 — Session 3: B1 — Cluster Annotator.

First module that writes results back to adata.obs. Nothing else.

Files to create:
  - ai/skills/cluster_annotator.yaml   ← skill YAML
  - ai/cluster_annotator.py            ← feature module
  - tests/test_cluster_annotator.py    ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~307 passed, 1 skipped
```

## Step 2 — Read the cluster annotator spec in PHASE3_PLAN.md
Section: "Session 3 — B1: Cluster Annotator"
Read it fully before writing any code.

## Step 3 — Build ai/skills/cluster_annotator.yaml
Follow pipeline_advisor.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - cluster_id, n_cells
  - marker_genes: list of top N marker genes ranked by log2FC

Output schema (JSON) per cluster:
  - cell_type: str
  - confidence: "high" | "medium" | "low"
  - supporting_markers: list of str
  - alternative_types: list of str
  - recommended_db: str
  - manual_marker_set: list of str
  - reasoning: str

## Step 4 — Build ai/cluster_annotator.py

### Public API
```python
from dataclasses import dataclass, field
from ai._base import AiResult

@dataclass
class ClusterAnnotation(AiResult):
    cluster_id: str = ""
    cell_type: str = "unknown"
    confidence: str = "low"           # "high" | "medium" | "low"
    supporting_markers: list[str] = field(default_factory=list)
    alternative_types: list[str] = field(default_factory=list)
    recommended_db: str = ""
    manual_marker_set: list[str] = field(default_factory=list)

def run(
    adata,                  # AnnData after marker computation — WRITES to obs
    config: dict,
    study_context: dict,
    *,
    n_markers: int = 20,    # top N markers per cluster to pass to LLM
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> list[ClusterAnnotation] | None:
    ...
```

### What it writes to adata
- `adata.obs['ai_cell_type']`     — predicted label per cell
- `adata.obs['ai_confidence']`    — high | medium | low per cell
- `adata.obs['ai_alt_types']`     — comma-separated alternatives per cell

### Input it reads from adata
- `adata.uns['rank_genes_groups']` — top marker genes per cluster
  If missing → raise informative error (not silent None)
- cluster assignments from `adata.obs['leiden']` (or first available
  clustering column)

### One LLM call per cluster
Loop over clusters. Each call is independent. If one cluster fails to
parse, log a warning and skip that cluster — do not abort the full run.

### Usage pattern (identical to all feature modules)
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, study_context, *, n_markers, log_dir, runtime_ai):
    try:
        check_ai_enabled(config, module="cluster_annotator", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
```

### call_llm signature (confirmed in Session 2)
```python
# Returns str only — not a tuple
raw_response = call_llm(
    skill_name="cluster_annotator",
    inputs=skill_inputs,
    config=config,
    log_dir=log_dir,
    module="cluster_annotator",
)
```

### write_audit_record signature (confirmed in Session 1)
```python
write_audit_record(
    log_dir=log_dir,
    module="cluster_annotator",
    skill_version="1.0",
    model=model,          # from config.ai.model
    provider=provider,    # from config.ai.provider
    input_summary={...},
    token_usage=None,     # call_llm already logged tokens; pass None here
    raw_response=raw_response,
    parsed_output={...} | None,
    parse_success=True | False,
)
```

## Step 5 — Tests (all without a real API key)
Mock pattern: patch("ai.cluster_annotator.call_llm") — returns str only.

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns list[ClusterAnnotation] when mock LLM returns valid JSON
  - adata.obs columns written correctly (ai_cell_type, ai_confidence, ai_alt_types)
  - Missing rank_genes_groups raises informative error
  - Partial parse failure (1 of 3 clusters fails) — other 2 still written
  - AiResult base fields (timestamp, model, provider) populated correctly
  - Audit log written to log_dir after successful run
  - Graceful handling of malformed LLM JSON for one cluster (logs warning, skips)
  - confidence values are always one of: high | medium | low

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ✅ DONE — A2: Clustering advisor (22 tests passing)
Session 3  ← TODAY — B1: Cluster annotator (first module writing to adata.obs)
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

## Confirmed API Signatures (do not guess — use these exactly)

### call_llm (ai/_llm_client.py)
```python
call_llm(
    skill_name: str,       # name of skill YAML, e.g. "cluster_annotator"
    inputs: dict,          # fills user_prompt_template
    config: dict,          # full pipeline config (not just ai section)
    *,
    log_dir: str,
    module: str | None,
    runtime_ai: bool,
) -> str                   # raw response string ONLY — no tuple
```

### write_audit_record (ai/_audit_log.py)
```python
write_audit_record(
    *,
    log_dir: str | Path,
    module: str,
    skill_version: str,
    model: str,
    provider: str,
    input_summary: dict,
    token_usage: dict | None,   # dict with prompt_tokens/completion_tokens or None
    raw_response: str,
    parsed_output: dict | None,
    parse_success: bool,
) -> None
```

### check_ai_enabled (ai/_config_gate.py)
```python
check_ai_enabled(config: dict, module: str, runtime_ai: bool = True) -> None
# Raises AiDisabledError if: global off, module off, or runtime_ai=False
```

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-2 Modules (done, do not modify)
- ai/pipeline_advisor.py            ← 13 tests passing
- ai/skills/pipeline_advisor.yaml
- ai/clustering_advisor.py          ← 22 tests passing
- ai/skills/clustering_advisor.yaml

## Bugs Fixed This Session (carry as warnings)
- call_llm returns str only — never unpack as tuple
- write_audit_record uses token_usage=dict, not separate prompt_tokens=/completion_tokens=
- Nested dataclasses and tuples in result objects need explicit _to_dict()
  helper before passing to parsed_output= in write_audit_record
- Provider and model always read from config — never hardcoded anywhere
- Default provider: "ollama", default model: "llama3" (free, local, no key)

## ClawBio Note (added this session)
ClawBio (clawbio.ai) is an open-source bioinformatics skill library.
Not a competitor — no report engine, no multi-project, no multiome.
Potential future integration as optional skill backends (Phase 4-5).
Potential OpenClaw marketplace listing (Phase 7+).
Not a dependency for Phase 1-3.
Monitor: claw-spatial, scRNA Orchestrator updates.
Entry added to DECISIONS.md.

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
- YAML skill files: output_schema must always use block scalar (|)
  Never use bare type annotations: list[{...}], str | null, list[str]
  These are illegal YAML and will cause ScannerError at load time.

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2
