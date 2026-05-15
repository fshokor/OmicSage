## Session Context
Date: next session
Phase: 3 — AI Layer
Last thing completed: Session 6 — A3 Downstream Analysis Suggester
                      built and tested. 20 new tests passing (400 total, 1 skipped).
                      Rule-based triggers (trajectory, communication, survival,
                      sub-clustering, pseudobulk) fire before LLM call. LLM
                      suggestions de-duplicated against rule-based by step_name.
                      NEXT_STEPS.md written as human-readable markdown.
File last worked on: tests/test_downstream_suggester.py

## Today's Goal
Phase 3 — Session 7: C1 — Biological Narrative Generator.
Nothing else.

Files to create:
  - ai/skills/narrative_generator.yaml     ← skill YAML
  - ai/narrative_generator.py              ← feature module
  - tests/test_narrative_generator.py      ← all tests without a real API key

## Step 1 — Verify all tests still pass
```bash
cd ~/OmicSage
conda activate omicsage
python -m pytest tests/ -v
# Expected: 400 passed, 1 skipped
```

## Step 2 — Read the Narrative Generator spec in PHASE3_PLAN.md
Section: "Session 7 — C1: Biological Narrative Generator"
Read it fully before writing any code.

## Step 3 — Build ai/skills/narrative_generator.yaml
Follow downstream_suggester.yaml pattern exactly.
CRITICAL: output_schema must use block scalar (|) — never bare type annotations.

Inputs the skill needs:
  - tissue, disease_context (from study_context)
  - biological_question (from study_context)
  - analysis_summary: str (JSON-serialised analysis_summary.json from B3)
  - pipeline_advice: str (JSON-serialised PipelineAdvice from A1, optional)
  - cluster_annotations: str (JSON-serialised per-cluster annotation from B1, optional)
  - deg_validation: str (JSON-serialised DegValidation from B2, optional)
  - coherence_review: str (JSON-serialised CoherenceReview from B3, optional)
  - narrative_block: str (which block to generate — one of:
      "qc_rationale" | "cell_type_landscape" | "differential_expression" | "interpretation")

Output schema (JSON):
  - narrative_text: str (the generated narrative paragraph(s) for this block)
  - cited_evidence: list of str (metric values, gene names, or PMIDs used)
  - groundedness_score: float (cited_evidence count / total factual sentences,
      self-reported by the model — cross-checked in test_groundedness.py)
  - reasoning: str (why the model wrote what it wrote)

## Step 4 — Build ai/narrative_generator.py

### Public API
```python
@dataclass
class NarrativeBlock:
    block_name: str = ""          # "qc_rationale" | "cell_type_landscape" | etc.
    narrative_text: str = ""
    cited_evidence: list[str] = field(default_factory=list)
    groundedness_score: float = 0.0

@dataclass
class NarrativeResult(AiResult):
    blocks: list[NarrativeBlock] = field(default_factory=list)
    overall_groundedness: float = 0.0   # mean across all blocks

def run(
    adata,
    config: dict,
    study_context: dict,
    pipeline_advice=None,        # PipelineAdvice | None from A1
    cluster_annotations=None,    # list[ClusterAnnotation] | None from B1
    deg_validation=None,         # DegValidation | None from B2
    coherence_review=None,       # CoherenceReview | None from B3
    *,
    output_path: str | None = None,   # if set, write ai_narrative.md here
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> NarrativeResult | None:
    ...
```

### Key behaviour
- Generates four narrative blocks in order:
    1. qc_rationale
    2. cell_type_landscape
    3. differential_expression
    4. interpretation
- Each block is a SEPARATE LLM call using the same skill YAML with
  narrative_block set to the block name. This keeps each call focused.
- If the upstream data for a block is missing (e.g. no deg_validation →
  skip differential_expression block), that block is skipped silently.
  The NarrativeResult still returns with the remaining blocks.
- overall_groundedness = mean of groundedness_score across generated blocks.
  Target ≥ 0.85 (validated in Session 10 test_groundedness.py).

### Output written to ai_narrative.md
```markdown
# AI Biological Narrative

*Generated: <timestamp>*
*Model: <model> (<provider>)*
*Overall groundedness score: <score>*

## QC Rationale
<narrative_text>

## Cell Type Landscape
<narrative_text>

## Differential Expression
<narrative_text>

## Interpretation and Perspectives
<narrative_text>
```

### Groundedness target
Every factual claim must cite a metric value, a gene name, or a PMID.
groundedness_score = cited factual sentences / total factual sentences.
Target ≥ 0.85. The model self-reports this — test_groundedness.py will
cross-check it on the real GSE166635 run in Session 10.

### Graceful degradation
If any upstream module output is None, that block is skipped.
If a block's LLM call fails to parse, that block is skipped.
The remaining blocks are returned normally.

## Step 5 — Tests (all without a real API key)
Mock pattern:
  - patch("ai.narrative_generator.call_llm") — returns str only

Required tests:
  - Returns None when ai_features=False
  - Returns None when runtime_ai=False
  - Returns NarrativeResult when mock LLM returns valid JSON
  - All four blocks generated when all upstream inputs provided
  - qc_rationale block skipped when analysis_summary missing qc data
  - differential_expression block skipped when deg_validation=None
  - NarrativeBlock fields populated (block_name, narrative_text, cited_evidence)
  - overall_groundedness calculated correctly (mean of block scores)
  - overall_groundedness = 0.0 when no blocks generated
  - ai_narrative.md written to output_path when provided
  - ai_narrative.md has correct markdown structure (headings per block)
  - Degraded parse (invalid JSON for one block) skips that block, others continue
  - AiResult base fields populated (timestamp, model, provider, skill_name)
  - pipeline_advice=None and coherence_review=None handled gracefully
  - Four separate call_llm calls made (one per block) when all inputs present

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
Session 7  ← TODAY — C1: Narrative generator
Session 8  — C2: Full report + PowerPoint
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

## Infrastructure Modules (Session 0 — all done, do not modify)
- ai/_base.py             ← AiResult base dataclass
- ai/_config_gate.py      ← check_ai_enabled() + AiDisabledError
- ai/_audit_log.py        ← write_audit_record()
- ai/_llm_client.py       ← call_llm() + _build_conversation()
- ai/_skill_loader.py     ← load_skill()

## Session 1-6 Modules (done, do not modify)
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
