## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 5 — B3 Coherence Reviewer + analysis_summary.json
                      built and tested. 22 new tests passing (377 total, 1 skipped).
                      build_analysis_summary() is a pure function — no LLM, fully
                      tested without mocking. analysis_summary.json written to
                      reports/<dataset>/analysis_summary.json when summary_path set.
File last worked on: tests/test_coherence_reviewer.py

## Today's Goal
Phase 3 — Session 6: A3 — Downstream Analysis Suggester.
Nothing else.

Files to create:
  - ai/skills/downstream_suggester.yaml     ← skill YAML
  - ai/downstream_suggester.py              ← feature module
  - tests/test_downstream_suggester.py      ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 377 passed, 1 skipped
```

## Step 2 — Read the Downstream Suggester spec in PHASE3_PLAN.md
Section: "Session 6 — A3: Downstream Analysis Suggester"
Read it fully before writing any code.

## Step 3 — Build ai/skills/downstream_suggester.yaml
Follow coherence_reviewer.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - biological_question (from study_context)
  - analysis_summary: str (JSON-serialised analysis_summary.json)
  - coherence_flags: str (JSON-serialised flags from B3, optional)

Output schema (JSON):
  - suggestions: list of dicts — each with step_name, rationale,
    expected_output, relevant_tool
  - reasoning: str

## Step 4 — Build ai/downstream_suggester.py

### Public API
```python
@dataclass
class DownstreamSuggestion:
    step_name: str = ""
    rationale: str = ""
    expected_output: str = ""
    relevant_tool: str = ""

@dataclass
class DownstreamAdvice(AiResult):
    suggestions: list[DownstreamSuggestion] = field(default_factory=list)

def run(
    adata,
    config: dict,
    study_context: dict,
    coherence_review=None,       # CoherenceReview | None from B3
    *,
    output_path: str | None = None,   # if set, write NEXT_STEPS.md here
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> DownstreamAdvice | None:
    ...
```

### Rule-based logic (run BEFORE calling LLM — these always fire)
These fire based on what is present in adata.obs and adata.uns:
- Progenitor + mature cell types in obs['ai_cell_type'] → suggest trajectory
- Immune + non-immune cells both present → suggest cell-cell communication
- Clinical metadata column present in obs → suggest survival analysis
- sub_clustering_candidates non-empty in coherence_review → suggest sub-clustering
- n_conditions > 1 in study_context → suggest pseudobulk DEG if not already done

### Output written to NEXT_STEPS.md
Human-readable markdown. Each suggestion as a section:
  ## <step_name>
  **Rationale**: ...
  **Expected output**: ...
  **Relevant tool**: ...

### Usage pattern
```python
from ai._config_gate import check_ai_enabled, AiDisabledError

def run(adata, config, study_context, coherence_review=None, *, output_path, log_dir, runtime_ai):
    try:
        check_ai_enabled(config, module="downstream_suggester", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None
```

## Step 5 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.downstream_suggester.call_llm") — returns str only

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns DownstreamAdvice when mock LLM returns valid JSON
  - Trajectory suggestion fires when progenitor + mature cell types present
  - Communication suggestion fires when immune + non-immune cells present
  - Sub-clustering suggestion fires when coherence_review has sub_clustering_candidates
  - suggestions list populated correctly (step_name, rationale, expected_output, tool)
  - NEXT_STEPS.md written to output_path when provided
  - NEXT_STEPS.md is valid markdown with correct structure
  - Degraded parse (invalid JSON) does not crash — returns DownstreamAdvice with empty list
  - AiResult base fields populated (timestamp, model, provider, skill_name)
  - coherence_review=None handled gracefully (no crash)

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ✅ DONE — A2: Clustering advisor (22 tests passing)
Session 3  ✅ DONE — B1: Cluster annotator (23 tests passing)
Session 4  ✅ DONE — B2: DEG validator + literature linker (25 tests passing)
Session 5  ✅ DONE — B3: Coherence reviewer + analysis_summary.json (22 tests passing)
Session 6  ← TODAY — A3: Downstream analysis suggester
Session 7  — C1: Narrative generator
Session 8  — C2: Full report + PowerPoint
Session 9  — D1: Report reviewer (reads final HTML report, not adata)  ← NEW
Session 10 — Milestone validation: groundedness test + end-to-end GSE166635
```
Full details: .dev_memory/PHASE3_PLAN.md

## FUTURE DESIGN DECISION — Report Reviewer (Session 9)
A new module `ai/report_reviewer.py` to be built in Session 9, after C2.

What it does:
  Reads the final HTML report generated by C2 (not adata) and reviews it
  as a document — checking narrative consistency, figure-caption alignment,
  methods text accuracy, and whether the conclusions follow from the results.
  Complements B3 (which reviews adata) with a document-level review.

When it runs: after C2 has written the final HTML report to disk.

Input: path to the final HTML report file
Output:
  - report_flags: list of dicts (category, severity, description, suggestion)
    Categories: narrative | figures | methods | conclusions
  - overall_report_quality: str (one paragraph)
  - Written to: reports/<dataset>/report_review.md

Skill file: ai/skills/report_reviewer.yaml
Module file: ai/report_reviewer.py
Test file: tests/test_report_reviewer.py

Key design note: reads HTML as plain text — strips tags, extracts section
text, feeds to LLM. Does NOT require a browser or rendering engine.
BeautifulSoup4 (already in many envs) for tag stripping if available,
otherwise regex fallback. Check if bs4 is in omicsage env before building.

## AI Layer Design Principles (carry into every session)
- Skills are the prompts — no raw strings in Python, ever
- AI is optional at three levels: global / per-module / per-step runtime flag
- A disabled module always returns None silently — never errors, never logs
- Every LLM call is audit-logged — no exceptions
- Graceful degradation — a failed AI module returns None, never crashes
- Context scales quality — more study_context = better AI output
- analysis_summary.json is load-bearing — designed and built in Session 5
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

## Session 1-5 Modules (done, do not modify)
- ai/pipeline_advisor.py            ← 13 tests passing
- ai/skills/pipeline_advisor.yaml
- ai/clustering_advisor.py          ← 22 tests passing
- ai/skills/clustering_advisor.yaml
- ai/cluster_annotator.py           ← 23 tests passing
- ai/skills/cluster_annotator.yaml
- ai/deg_validator.py               ← 25 tests passing
- ai/skills/deg_validator.yaml
- ai/coherence_reviewer.py          ← 22 tests passing
- ai/skills/coherence_reviewer.yaml

## Bugs Fixed in Prior Sessions (carry as warnings)
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
