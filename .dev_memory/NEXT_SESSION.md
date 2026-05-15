## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 9 — D1 Report Reviewer built and tested.
                      18 new tests passing (466 total, 1 skipped).
                      Skill YAML with 4 review categories (narrative, figures,
                      methods, conclusions) and 3 severity levels.
                      HTML tag stripping via BeautifulSoup with regex fallback.
                      24000-char truncation before LLM call.
                      report_review.md written to report_dir.
                      Graceful degradation: missing HTML → None,
                      invalid JSON → result with empty flags,
                      missing report_dir → created automatically.
File last worked on: tests/test_report_reviewer.py

## Today's Goal
Phase 3 — Session 10: Milestone validation.
Nothing else.

Two parts:
  1. Groundedness test (tests/test_groundedness.py)
  2. End-to-end pipeline run on GSE166635 with ai_features: true and false

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 466 passed, 1 skipped
```

## Step 2 — Groundedness test
File to create: tests/test_groundedness.py

What it measures:
  Score = cited factual sentences / total factual sentences in ai_narrative.md
  Target: >= 0.85

A factual sentence is any sentence containing:
  - a numeric value (e.g. "14 clusters", "35% of cells")
  - a gene name (all-caps 2-6 char token, e.g. PDCD1, CD8A)
  - a PMID reference

A cited factual sentence is a factual sentence that also contains
a PMID, a gene name used as evidence, or an explicit metric value
that matches a value in analysis_summary.json.

Run against: reports/GSE166635/ai_narrative.md
             reports/GSE166635/analysis_summary.json

## Step 3 — End-to-end GSE166635 with AI on
```bash
python run_pipeline.py \
  --config config/runs/GSE166635.yaml \
  --step all --ai \
  2>&1 | tee logs/GSE166635_phase3_ai_on.log
```

Verify after run:
  - reports/GSE166635/ai_narrative.md exists and non-empty
  - reports/GSE166635/report_review.md exists and non-empty
  - reports/GSE166635/NEXT_STEPS.md exists and non-empty
  - reports/GSE166635/analysis_summary.json exists and valid JSON
  - reports/GSE166635/presentation.pptx exists (open and verify 8 slides)
  - logs/llm/*.jsonl all exist and non-empty
  - Groundedness score >= 0.85

## Step 4 — End-to-end GSE166635 with AI off
```bash
python run_pipeline.py \
  --config config/runs/GSE166635.yaml \
  --step all \
  2>&1 | tee logs/GSE166635_phase3_ai_off.log
```

Verify after run:
  - Pipeline completes without errors
  - No files written to logs/llm/
  - No ai_narrative.md, no report_review.md, no NEXT_STEPS.md
  - HTML report and figures still generated normally

## Step 5 — Ollama local provider test (manual)
```bash
# Requires ollama running locally with llama3 pulled
python run_pipeline.py \
  --config config/runs/GSE166635_ollama.yaml \
  --step annotate --ai
```
Verify: logs/llm/cluster_annotator.jsonl written, provider field = "ollama"

## Phase 3 Milestone Checklist

[ ] All prior tests still passing (466+)
[ ] All AI feature modules have passing tests
[ ] Full pipeline runs on GSE166635 with ai_features: true
[ ] Full pipeline runs on GSE166635 with ai_features: false
[ ] Groundedness score >= 0.85 on GSE166635 narrative
[ ] All LLM calls written to logs/llm/*.jsonl
[ ] Ollama provider works locally with llama3
[ ] PowerPoint generated with all 8 slides
[ ] NEXT_STEPS.md written to run output directory
[ ] analysis_summary.json written per run and valid JSON

## Step 5 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.report_reviewer.call_llm") — returns str only
  - Use tmp_path for all file output
  - Write a minimal fake HTML file to tmp_path for the happy-path tests

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns None when html_report_path does not exist
  - Returns ReportReviewerResult when mock LLM returns valid JSON
  - report_flags parsed correctly (3 flags from mock JSON)
  - overall_report_quality populated
  - report_review.md written to report_dir
  - report_review.md contains Overall Quality section
  - report_review.md contains at least one flag section
  - review_path set in result
  - AiResult base fields populated
  - One call_llm call made
  - HTML tags stripped before passing to LLM (assert "<" not in captured input)
  - Text truncated to max 24000 chars before LLM call
  - Invalid JSON from LLM handled gracefully (result not None, flags empty)
  - report_dir created if missing
  - bs4 absent falls back to regex (mock ImportError on bs4)

## Phase 3 Build Order (full reference)
```
Session 0  ✅ DONE — Infrastructure: _llm_client, _audit_log, _skill_loader,
                      _config_gate, _base
Session 1  ✅ DONE — A1: Pipeline advisor (13 tests passing)
Session 2  ✅ DONE — A2: Clustering advisor (22 tests passing)
Session 3  ✅ DONE — B1: Cluster annotator (23 tests passing)
Session 4  ✅ DONE — B2: DEG validator + literature linker (25 tests passing)
Session 5  ✅ DONE — B3: Coherence reviewer + analysis_summary.json (22 tests passing)
Session 6  ✅ DONE — A3: Downstream analysis suggester (20 tests passing)
Session 7  ✅ DONE — C1: Narrative generator (18 tests passing)
Session 8  ✅ DONE — C2: Full report + PowerPoint (20 tests passing)
Session 9  ✅ DONE — D1: Report reviewer (18 tests passing)
Session 10 ← TODAY — Milestone validatio
```
Full details: .dev_memory/PHASE3_PLAN.md

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

### ReportReviewerResult (ai/report_reviewer.py)
```python
@dataclass
class ReportFlag:
    category: str = ""      # narrative | figures | methods | conclusions
    severity: str = ""      # info | warning | critical
    description: str = ""
    suggestion: str = ""

@dataclass
class ReportReviewerResult(AiResult):
    report_flags: list[ReportFlag] = field(default_factory=list)
    overall_report_quality: str = ""
    review_path: str | None = None
```

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-9 Modules (done, do not modify)
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
- ai/downstream_suggester.py        ← 20 tests passing
- ai/skills/downstream_suggester.yaml
- ai/narrative_generator.py         ← 18 tests passing
- ai/skills/narrative_generator.yaml
- ai/report_writer.py               ← 20 tests passing
- ai/skills/report_writer.yaml
- ai/report_reviewer.py             ← 18 tests passing
- ai/skills/report_reviewer.yaml

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
- python-pptx is a soft dependency in report_writer.py — ImportError silently
  skips PowerPoint. Version confirmed: 1.0.2
- report_reviewer.py: text-only, no figure inspection — see DECISIONS.md

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2,
                   python-pptx 1.0.2, beautifulsoup4 4.14.3