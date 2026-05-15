## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 7 — C1 Biological Narrative Generator
                      built and tested. 18 new tests passing (418 total, 1 skipped).
                      Four narrative blocks (qc_rationale, cell_type_landscape,
                      differential_expression, interpretation) generated via
                      separate LLM calls. Graceful degradation per block when
                      upstream inputs missing. ai_narrative.md written to disk.
                      overall_groundedness = mean of block scores (target >= 0.85,
                      validated in Session 10).
File last worked on: tests/test_narrative_generator.py

## Today's Goal
Phase 3 — Session 8: C2 — Full Report + PowerPoint Generator.
Nothing else.

Files to create:
  - ai/skills/report_writer.yaml       ← skill YAML
  - ai/report_writer.py                ← feature module
  - tests/test_report_writer.py        ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 418 passed, 1 skipped
```

## Step 2 — Read the Report Writer spec in PHASE3_PLAN.md
Section: "Session 8 — C2: Full Report + PowerPoint Generator"
Read it fully before writing any code.

## Step 3 — Build ai/skills/report_writer.yaml
Follow narrative_generator.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - biological_question, objectives (from study_context)
  - analysis_summary: str (JSON-serialised analysis_summary.json from B3)
  - ai_narrative: str (contents of ai_narrative.md from C1, or "not available")
  - coherence_review: str (JSON-serialised CoherenceReview from B3, optional)
  - report_section: str (which section to generate — one of:
      "conclusion" | "perspectives")

Output schema (JSON):
  - section_text: str (the generated text for this report section)
  - key_findings: list of str (for conclusion: 3-5 bullet points summarising findings)
  - cited_evidence: list of str (metric values, gene names, or PMIDs cited)
  - reasoning: str

## Step 4 — Build ai/report_writer.py

### Public API
```python
@dataclass
class ReportSection:
    section_name: str = ""       # "conclusion" | "perspectives"
    section_text: str = ""
    key_findings: list[str] = field(default_factory=list)
    cited_evidence: list[str] = field(default_factory=list)

@dataclass
class ReportWriterResult(AiResult):
    sections: list[ReportSection] = field(default_factory=list)
    pptx_path: str | None = None    # path to written PowerPoint file

def run(
    adata,
    config: dict,
    study_context: dict,
    narrative_result=None,       # NarrativeResult | None from C1
    coherence_review=None,       # CoherenceReview | None from B3
    cluster_annotations=None,    # list[ClusterAnnotation] | None from B1
    deg_validation=None,         # DegValidation | None from B2
    *,
    report_dir: str,             # required — output directory for this dataset run
    figures_dir: str | None = None,  # directory containing saved figures
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> ReportWriterResult | None:
    ...
```

### Key behaviour — two deliverables

#### 1. HTML/PDF report sections (AI-generated)
Two LLM calls, one per section:
  - conclusion: 3-5 key findings as bullet points + one summary paragraph
  - perspectives: open questions + suggested follow-up experiments +
                  potential clinical relevance

Written to:
  - reports/<dataset>/conclusion.md
  - reports/<dataset>/perspectives.md

#### 2. PowerPoint (python-pptx, rule-based — NO LLM call)
Eight slides built programmatically from the data:
  1. Title slide — dataset name, tissue, disease, date
  2. Data overview — n_cells, modalities, QC metrics table
  3. UMAP — placeholder if figure not available
  4. Cell type proportions — bar chart from cluster_annotations
  5. Top DEGs — table from deg_validation
  6. Pathway enrichment — table from analysis_summary pathways
  7. Key findings — bullet points from conclusion section
  8. Conclusions + perspectives — from conclusion + perspectives sections

Speaker notes on every slide = corresponding C1 narrative paragraph
  (matched by slide topic to narrative block).

Written to: reports/<dataset>/presentation.pptx

### Graceful degradation
- If narrative_result is None, conclusion and perspectives are still
  generated from analysis_summary alone (C1 output is enrichment, not required).
- If a figure file is missing, that slide uses a text placeholder.
- If python-pptx is not installed, log a warning and skip PowerPoint only.
  HTML sections still written.
- pptx_path in result is None if PowerPoint was not written.

## Step 5 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.report_writer.call_llm") — returns str only
  - Use tmp_path for all file output

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns ReportWriterResult when mock LLM returns valid JSON
  - conclusion section generated from analysis_summary alone (no narrative_result)
  - perspectives section generated
  - Both sections written as .md files under report_dir
  - conclusion.md contains key findings as markdown bullets
  - PowerPoint written to report_dir/presentation.pptx
  - PowerPoint has correct number of slides (8)
  - Speaker notes present on all slides
  - pptx_path set in result when PowerPoint written
  - pptx_path is None when python-pptx unavailable (mock ImportError)
  - Missing figures_dir handled gracefully (placeholder text on slide)
  - AiResult base fields populated
  - Two separate call_llm calls made (conclusion + perspectives)
  - coherence_review=None and cluster_annotations=None handled gracefully

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
Session 8  ← TODAY — C2: Full report + PowerPoint
Session 9  — D1: Report reviewer (reads final HTML report, not adata)
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

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-7 Modules (done, do not modify)
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
- python-pptx must be a soft dependency in report_writer.py — wrap import
  in try/except ImportError, skip PowerPoint gracefully if not installed.
  Verify it is in the omicsage env before Session 8 starts:
  `conda activate omicsage && python -c "import pptx; print(pptx.__version__)"`

## Conda Environment
Name: omicsage
Activate: conda activate omicsage
Python: 3.11.15
Verified packages: scanpy 1.11.5, numpy 2.4.3, pytest 9.0.3, scrublet, mudata,
                   ipykernel, jupyter, scikit-misc, kneed>=0.8.5, celltypist,
                   gseapy 1.1.3, harmonypy, pydeseq2, biochatter==0.14.2
