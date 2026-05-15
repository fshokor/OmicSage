## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 8 — C2 Full Report + PowerPoint Generator
                      built and tested. 20 new tests passing (448 total, 1 skipped).
                      Two AI-generated report sections (conclusion.md,
                      perspectives.md) via separate LLM calls. Rule-based 8-slide
                      PowerPoint (presentation.pptx) via python-pptx — no LLM.
                      Speaker notes on all slides from C1 narrative blocks.
                      Graceful degradation: pptx_path=None if python-pptx missing,
                      placeholder text if figures missing, sections still written
                      if narrative_result=None.
                      python-pptx version confirmed: 1.0.2
File last worked on: tests/test_report_writer.py

## Today's Goal
Phase 3 — Session 9: D1 — Report Reviewer.
Nothing else.

Files to create:
  - ai/skills/report_reviewer.yaml     ← skill YAML
  - ai/report_reviewer.py              ← feature module
  - tests/test_report_reviewer.py      ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 448 passed, 1 skipped
```

## Step 2 — Check if bs4 is available
```bash
conda activate omicsage && python -c "import bs4; print(bs4.__version__)"
```
If present: use BeautifulSoup for HTML tag stripping.
If absent: use regex fallback. Do NOT install bs4 just for this module —
           the fallback must always work even when bs4 is present.

## Step 3 — Build ai/skills/report_reviewer.yaml
Follow narrative_generator.yaml / report_writer.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - report_text: str (plain text extracted from the HTML report, tags stripped)
  - tissue: str
  - disease_context: str | null
  - biological_question: str

Output schema (JSON):
  - report_flags: list of dicts, each with keys:
      category (narrative | figures | methods | conclusions),
      severity (info | warning | critical),
      description, suggestion
  - overall_report_quality: str (one paragraph)
  - reasoning: str

## Step 4 — Build ai/report_reviewer.py

### Public API
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
    review_path: str | None = None   # path to written report_review.md

def run(
    html_report_path: str,
    config: dict,
    study_context: dict,
    *,
    report_dir: str,         # required — output directory for this dataset run
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ReportReviewerResult | None:
    ...
```

### Key behaviour

#### HTML reading
- Open html_report_path as plain text (UTF-8, errors='replace')
- Strip HTML tags:
    - If bs4 available: BeautifulSoup(text, "html.parser").get_text(separator=" ")
    - Otherwise: re.sub(r"<[^>]+>", " ", text) then collapse whitespace
- Truncate to 6000 tokens (~24000 chars) — pass the truncated plain text to the LLM
- If html_report_path does not exist: log warning, return None

#### LLM call
- One call to report_reviewer skill
- Inputs: report_text (truncated plain text), tissue, disease_context,
          biological_question

#### Output written to disk
- reports/<dataset>/report_review.md
  Format:
    # Report Review

    ## Overall Quality
    <overall_report_quality paragraph>

    ## Flags
    ### <CATEGORY> — <SEVERITY>
    **Issue:** <description>
    **Suggestion:** <suggestion>
    (one section per flag)

- review_path in result = path to written file

### Graceful degradation
- html_report_path does not exist → log warning, return None
- LLM returns invalid JSON → log warning, return ReportReviewerResult with
  report_flags=[], overall_report_quality=raw_response[:500], review_path=None
- report_dir does not exist → create it

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
Session 9  ← TODAY — D1: Report reviewer
Session 10 — Milestone validation: groundedness test + end-to-end GSE166635
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

### NarrativeResult (ai/narrative_generator.py)
```python
@dataclass
class NarrativeBlock:
    block_name: str = ""
    narrative_text: str = ""
    cited_evidence: list[str] = field(default_factory=list)
    groundedness_score: float = 0.0

@dataclass
class NarrativeResult(AiResult):
    blocks: list[NarrativeBlock] = field(default_factory=list)
    overall_groundedness: float = 0.0
```

### ReportWriterResult (ai/report_writer.py)
```python
@dataclass
class ReportSection:
    section_name: str = ""
    section_text: str = ""
    key_findings: list[str] = field(default_factory=list)
    cited_evidence: list[str] = field(default_factory=list)

@dataclass
class ReportWriterResult(AiResult):
    sections: list[ReportSection] = field(default_factory=list)
    pptx_path: str | None = None
```

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-8 Modules (done, do not modify)
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

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2,
                   python-pptx 1.0.2
