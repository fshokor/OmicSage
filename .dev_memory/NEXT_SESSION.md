## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 4 — B2 DEG Validator + Literature Linker built and tested.
                      25 new tests passing (355 total, 1 skipped).
                      Fixed: rank_genes_groups structured array rows must have length
                      n_groups (one value per group field), not n_genes_per_group.
                      Fixed: mock functions patching _query_pubmed must match full
                      signature (gene, tissue, disease_context) — not just (gene).
File last worked on: tests/test_deg_validator.py

## Today's Goal
Phase 3 — Session 5: B3 — Coherence Reviewer.

This session also builds analysis_summary.json — the load-bearing compressed
run summary consumed by B3, C1, and C2. Design it before writing any module code.
Nothing else.

Files to create:
  - ai/skills/coherence_reviewer.yaml     ← skill YAML
  - ai/coherence_reviewer.py              ← feature module
  - tests/test_coherence_reviewer.py      ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 355 passed, 1 skipped
```

## Step 2 — Read the Coherence Reviewer spec in PHASE3_PLAN.md
Section: "Session 5 — B3: Coherence Reviewer"
Read it fully before writing any code.
Pay special attention to the analysis_summary.json schema — it is load-bearing.

## Step 3 — Design analysis_summary.json first

This file is written by B3 and consumed by C1 (narrative generator) and C2
(report writer). It must compress the full run into ~2000 tokens.

Schema (write this to disk as part of the module, not hardcoded in tests):
```json
{
  "study_context": {
    "tissue": "str",
    "disease": "str or null",
    "n_cells": "int",
    "n_batches": "int"
  },
  "qc_decisions": {
    "min_genes": "int",
    "max_mt_pct": "float",
    "cells_removed_pct": "float"
  },
  "clustering": {
    "resolution": "float",
    "n_clusters": "int",
    "silhouette_score": "float or null"
  },
  "cell_types": [
    { "cluster": "str", "cell_type": "str", "confidence": "str", "n_cells": "int" }
  ],
  "top_degs": [
    { "comparison": "str", "cluster": "str", "gene": "str", "log2fc": "float" }
  ],
  "pathways": [
    { "cluster": "str", "pathway": "str", "padj": "float" }
  ]
}
```

Rules:
- top_degs: max 3 genes per comparison (keeps token budget)
- pathways: max 3 pathways per cluster (keeps token budget)
- All fields optional except study_context.tissue and clustering.n_clusters
- Missing fields written as null, never omitted entirely
- Written to: reports/<dataset>/analysis_summary.json

## Step 4 — Build ai/skills/coherence_reviewer.yaml
Follow cluster_annotator.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - analysis_summary: str (JSON-serialised analysis_summary.json contents)
  - tissue, disease_context (from study_context)

Output schema (JSON):
  - flags: list of dicts — each with category, severity, description, suggestion
  - sub_clustering_candidates: list of str (cluster IDs)
  - rare_cell_candidates: list of str (cluster IDs)
  - overall_assessment: str

Flag categories: qc | clustering | annotation | deg | pathway
Flag severities: info | warning | critical

## Step 5 — Build ai/coherence_reviewer.py

### Public API
```python
from dataclasses import dataclass, field
from ai._base import AiResult

@dataclass
class CoherenceFlag:
    category: str = ""       # qc | clustering | annotation | deg | pathway
    severity: str = ""       # info | warning | critical
    description: str = ""
    suggestion: str = ""

@dataclass
class CoherenceReview(AiResult):
    flags: list[CoherenceFlag] = field(default_factory=list)
    sub_clustering_candidates: list[str] = field(default_factory=list)
    rare_cell_candidates: list[str] = field(default_factory=list)
    overall_assessment: str = ""

def build_analysis_summary(adata, config: dict, study_context: dict) -> dict:
    """Build the analysis_summary dict from adata. Does NOT call LLM."""
    ...

def run(
    adata,
    config: dict,
    study_context: dict,
    *,
    summary_path: str | None = None,   # if set, write analysis_summary.json here
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> CoherenceReview | None:
    ...
```

### What it reads from adata
- adata.uns['rank_genes_groups'] — for top_degs (optional, null if missing)
- adata.uns['leiden_resolution'] or adata.uns['omicsage_cluster'] — resolution
- adata.obs['leiden'] — cluster labels and cell counts
- adata.obs['ai_cell_type'] + adata.obs['ai_confidence'] — from B1 (optional)
- adata.uns['omicsage_qc'] — QC decisions (optional)
- adata.uns['gsea_results'] or similar — pathway results (optional)

### Key design rule
build_analysis_summary() is a pure function — no LLM, no side effects.
It must be independently testable without mocking anything.
run() calls build_analysis_summary() then calls the LLM.

### Usage pattern (identical to all feature modules)
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, study_context, *, summary_path, log_dir, runtime_ai):
    try:
        check_ai_enabled(config, module="coherence_reviewer", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
```

## Step 6 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.coherence_reviewer.call_llm") — returns str only

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - build_analysis_summary returns correct structure from mock adata
  - build_analysis_summary handles missing rank_genes_groups gracefully (null)
  - build_analysis_summary handles missing ai_cell_type gracefully (null)
  - Returns CoherenceReview when mock LLM returns valid JSON
  - flags list populated correctly (category, severity, description, suggestion)
  - sub_clustering_candidates populated from mock response
  - rare_cell_candidates populated from mock response
  - analysis_summary.json written to summary_path when provided
  - analysis_summary.json is valid JSON and matches schema
  - AiResult base fields (timestamp, model, provider, skill_name) populated
  - Audit log written to log_dir after successful run
  - Degraded parse (invalid JSON from LLM) does not crash — returns CoherenceReview
    with empty flags

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ✅ DONE — A2: Clustering advisor (22 tests passing)
Session 3  ✅ DONE — B1: Cluster annotator (23 tests passing)
Session 4  ✅ DONE — B2: DEG validator + literature linker (25 tests passing)
Session 5  ← TODAY — B3: Coherence reviewer + analysis_summary.json
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
- analysis_summary.json is load-bearing — design carefully (done this session)
- Copyright enforced at report layer — PMIDs + titles only, no abstract text

## Confirmed API Signatures (do not guess — use these exactly)

### call_llm (ai/_llm_client.py)
```python
call_llm(
    skill_name: str,
    inputs: dict,
    config: dict,
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
    token_usage: dict | None,
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

### AiResult base fields (ai/_base.py)
```python
timestamp: str      # ISO-8601 UTC, set automatically
model: str
provider: str
skill_name: str     # ← NOTE: skill_name, NOT skill
skill_version: str
reasoning: str
```

### _query_pubmed (ai/deg_validator.py)
```python
_query_pubmed(gene: str, tissue: str, disease_context: str | None) -> list[GeneLitRef]
# Mock functions patching this MUST match all three arguments
```

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-4 Modules (done, do not modify)
- ai/pipeline_advisor.py            ← 13 tests passing
- ai/skills/pipeline_advisor.yaml
- ai/clustering_advisor.py          ← 22 tests passing
- ai/skills/clustering_advisor.yaml
- ai/cluster_annotator.py           ← 23 tests passing
- ai/skills/cluster_annotator.yaml
- ai/deg_validator.py               ← 25 tests passing
- ai/skills/deg_validator.yaml

## Bugs Fixed This Session (carry as warnings)
- rank_genes_groups structured array: rows must have length n_groups (one value
  per group field), NOT n_genes_per_group. The dtype field count is what matters.
  Correct: row rank = (group0_gene_rank, group1_gene_rank, ...)
- Mock functions patching _query_pubmed must match the full 3-argument signature
  (gene, tissue, disease_context). Partial signatures pass patch() silently but
  explode at call time with "takes 1 positional argument but 3 were given".

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
- test_phase0_structure.py: TestCLI asserts cfg.provider == "ollama" (not "claude")
- YAML skill files: output_schema must always use block scalar (|)
  Never use bare type annotations: list[{...}], str | null, list[str]
  These are illegal YAML and will cause ScannerError at load time.
- biochatter==0.14.2 pinned — do NOT upgrade without re-running all AI tests
- call_llm returns str only — never unpack as tuple
- write_audit_record uses token_usage=dict (or None), not separate kwargs
- Provider and model always read from config — never hardcoded anywhere
- Default provider: "ollama", default model: "llama3"

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2
