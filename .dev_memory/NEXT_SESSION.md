## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 1 — A1 Pipeline Advisor built and tested.
                      13 tests passing. YAML parse bug fixed (output_schema
                      must use block scalar | not bare type annotations like
                      list[{...}] or str | null — these are illegal in bare YAML).
File last worked on: tests/test_pipeline_advisor.py

## Today's Goal
Phase 3 — Session 2: A2 — Clustering Advisor.

First module that uses PubMed RAG. Nothing else.

Files to create:
  - ai/skills/clustering_advisor.yaml   ← skill YAML
  - ai/clustering_advisor.py            ← feature module
  - tests/test_clustering_advisor.py    ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: ~285 passed, 1 skipped
```

## Step 2 — Read the clustering advisor spec in PHASE3_PLAN.md
Section: "Session 2 — A2: Clustering Advisor"
Read it fully before writing any code.

## Step 3 — Build ai/skills/clustering_advisor.yaml
Follow the pipeline_advisor.yaml pattern exactly.
CRITICAL: output_schema must use a block scalar (|) — never bare type
annotations. This was the bug in Session 1.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - resolution_sweep_results: dict of resolution → silhouette_score
  - n_cells, n_highly_variable_genes

Output schema (JSON):
  - suggested_resolution: float
  - resolution_range: [float, float]   ← reasonable range to explore
  - expected_n_clusters: int
  - literature_context: list of {pmid, title, resolution_used, tissue, disease}
  - reasoning: str

## Step 4 — Build ai/clustering_advisor.py

### Public API
```python
from dataclasses import dataclass, field
from ai._base import AiResult

@dataclass
class PubMedRef:
    pmid: str
    title: str
    resolution_used: float | None
    tissue: str
    disease: str | None

@dataclass
class ClusteringAdvice(AiResult):
    suggested_resolution: float = 0.5
    resolution_range: tuple[float, float] = (0.3, 0.8)
    expected_n_clusters: int = 0
    literature_context: list[PubMedRef] = field(default_factory=list)

def run(
    adata,                  # AnnData after dim. reduction — reads only, never modifies
    config: dict,
    study_context: dict,
    resolution_sweep_results: dict[float, float],  # resolution → silhouette_score
    *,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ClusteringAdvice | None:
    ...
```

### PubMed RAG query pattern (from PHASE3_PLAN.md)
Query string: "Leiden clustering resolution {tissue} {disease} single-cell RNA-seq"
Use BioChatter's built-in retrieval — do NOT call PubMed API directly.
Return PMIDs + titles only. Never abstract text (copyright rule).

### Rule-based pre-checks (before LLM call)
  - Find resolution with highest silhouette score in sweep → anchor suggestion
  - resolution_sweep_results empty → warn and use default 0.5
  - n_cells < 1000 → bias toward lower resolution (fewer expected clusters)
  - n_cells > 50000 → bias toward higher resolution

### Usage pattern (identical to all feature modules)
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, study_context, resolution_sweep_results, *, log_dir, runtime_ai=True):
    try:
        check_ai_enabled(config, module="clustering_advisor", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
```

## Step 5 — Tests (all without a real API key)
Mock pattern: patch ai.clustering_advisor.call_llm (same as Session 1).
For PubMed RAG: patch the BioChatter retrieval call separately.

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns ClusteringAdvice when mock LLM returns valid JSON
  - suggested_resolution is within the resolution_sweep_results range
  - resolution_range is a 2-tuple with range[0] < range[1]
  - Empty resolution_sweep_results handled gracefully (no crash, default used)
  - PubMed query string constructed correctly from study_context fields
  - Empty PubMed result handled gracefully (literature_context = [])
  - AiResult base fields (timestamp, model, provider) populated correctly
  - Audit log written to log_dir after successful run
  - Graceful handling of malformed LLM JSON (returns None, logs warning)

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ← TODAY — A2: Clustering advisor (first PubMed RAG use)
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

## Session 1 Modules (done, do not modify)
- ai/pipeline_advisor.py            ← 13 tests passing
- ai/skills/pipeline_advisor.yaml   ← block scalar fix applied

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
