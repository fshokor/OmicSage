"""
ai/narrative_generator.py
Phase 3 — Session 7 — C1: Biological Narrative Generator

Generates four narrative blocks for the analysis report, one LLM call each:
  1. qc_rationale
  2. cell_type_landscape
  3. differential_expression
  4. interpretation

Each block cites metric values, gene names, or PMIDs from upstream AI outputs.
Target groundedness_score >= 0.85 (validated in Session 10 test_groundedness.py).

Returns NarrativeResult | None.
Returns None when AI is disabled (global, per-module, or runtime flag).
A failed or missing block is skipped silently; the rest are returned normally.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Block definitions — order is fixed; presence depends on upstream inputs
# ---------------------------------------------------------------------------

_BLOCK_ORDER = [
    "qc_rationale",
    "cell_type_landscape",
    "differential_expression",
    "interpretation",
]

# Which upstream input must be non-None for each block to be generated.
# "interpretation" is only skipped if ALL upstream inputs are missing.
_BLOCK_REQUIREMENTS: dict[str, list[str]] = {
    "qc_rationale": ["analysis_summary"],
    "cell_type_landscape": ["cluster_annotations"],
    "differential_expression": ["deg_validation"],
    "interpretation": [],  # generated as long as analysis_summary is present
}


# ---------------------------------------------------------------------------
# Output dataclasses
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _serialise(obj) -> str:
    """Serialise an arbitrary object to a JSON string for prompt injection.
    Returns 'not available' for None."""
    if obj is None:
        return "not available"
    if isinstance(obj, str):
        return obj
    try:
        # dataclass or object with __dict__
        if hasattr(obj, "__dataclass_fields__"):
            return json.dumps(_dataclass_to_dict(obj), indent=2)
        return json.dumps(obj, indent=2, default=str)
    except Exception:
        return str(obj)


def _dataclass_to_dict(obj) -> dict:
    """Recursively convert a dataclass to a plain dict."""
    import dataclasses
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return {
            k: _dataclass_to_dict(v)
            for k, v in dataclasses.asdict(obj).items()
        }
    return obj


def _parse_block_response(raw: str, block_name: str) -> NarrativeBlock | None:
    """Parse a raw LLM response string into a NarrativeBlock.
    Returns None if parsing fails."""
    try:
        data = json.loads(raw)
        return NarrativeBlock(
            block_name=block_name,
            narrative_text=str(data.get("narrative_text", "")),
            cited_evidence=list(data.get("cited_evidence", [])),
            groundedness_score=float(data.get("groundedness_score", 0.0)),
        )
    except (json.JSONDecodeError, ValueError, TypeError) as exc:
        logger.warning(
            "narrative_generator: failed to parse block '%s': %s", block_name, exc
        )
        return None


def _block_inputs_present(
    block_name: str,
    analysis_summary: str | None,
    cluster_annotations,
    deg_validation,
) -> bool:
    """Return True if the required upstream inputs for this block are present."""
    reqs = _BLOCK_REQUIREMENTS[block_name]
    if not reqs:
        # interpretation: generate as long as analysis_summary exists
        return analysis_summary is not None
    for req in reqs:
        if req == "analysis_summary" and analysis_summary is None:
            return False
        if req == "cluster_annotations" and cluster_annotations is None:
            return False
        if req == "deg_validation" and deg_validation is None:
            return False
    return True


def _write_markdown(
    result: NarrativeResult,
    output_path: str,
) -> None:
    """Write the narrative result to ai_narrative.md."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    block_heading = {
        "qc_rationale": "QC Rationale",
        "cell_type_landscape": "Cell Type Landscape",
        "differential_expression": "Differential Expression",
        "interpretation": "Interpretation and Perspectives",
    }

    lines = [
        "# AI Biological Narrative",
        "",
        f"*Generated: {result.timestamp}*",
        f"*Model: {result.model} ({result.provider})*",
        f"*Overall groundedness score: {result.overall_groundedness:.2f}*",
        "",
    ]

    for block in result.blocks:
        heading = block_heading.get(block.block_name, block.block_name.replace("_", " ").title())
        lines.append(f"## {heading}")
        lines.append("")
        lines.append(block.narrative_text)
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")
    logger.info("narrative_generator: wrote %s", path)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run(
    adata,
    config: dict,
    study_context: dict,
    pipeline_advice=None,
    cluster_annotations=None,
    deg_validation=None,
    coherence_review=None,
    *,
    output_path: str | None = None,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> NarrativeResult | None:
    """Generate biological narrative blocks for the analysis report.

    Parameters
    ----------
    adata
        AnnData object (used only to extract analysis_summary from uns if present).
    config
        Pipeline config dict. Must contain ai section.
    study_context
        Dict loaded from config/study_context.yaml.
    pipeline_advice
        PipelineAdvice dataclass from A1, or None.
    cluster_annotations
        List of ClusterAnnotation dataclasses from B1, or None.
    deg_validation
        DegValidation dataclass from B2, or None.
    coherence_review
        CoherenceReview dataclass from B3, or None.
    output_path
        If set, write ai_narrative.md to this path.
    log_dir
        Directory for audit log JSONL files.
    runtime_ai
        Set False to disable AI at runtime regardless of config.

    Returns
    -------
    NarrativeResult | None
    """
    try:
        check_ai_enabled(config, module="narrative_generator", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- resolve model / provider from config ---
    ai_cfg = config.get("ai", {})
    model = ai_cfg.get("model", "llama3")
    provider = ai_cfg.get("provider", "ollama")

    # --- extract study context fields ---
    tissue = study_context.get("dataset", {}).get("tissue", "unknown tissue")
    disease_context = (
        study_context.get("disease", {}).get("context", None) or "not specified"
    )
    biological_question = study_context.get("biological_question", "not specified")

    # --- extract analysis_summary from adata.uns if available ---
    analysis_summary_str: str | None = None
    if adata is not None and hasattr(adata, "uns"):
        raw_summary = adata.uns.get("omicsage_analysis_summary")
        if raw_summary is not None:
            analysis_summary_str = (
                raw_summary
                if isinstance(raw_summary, str)
                else json.dumps(raw_summary, indent=2, default=str)
            )

    # Serialise all upstream inputs for prompt injection
    pipeline_advice_str = _serialise(pipeline_advice)
    cluster_annotations_str = _serialise(cluster_annotations)
    deg_validation_str = _serialise(deg_validation)
    coherence_review_str = _serialise(coherence_review)

    # --- generate blocks ---
    blocks: list[NarrativeBlock] = []

    for block_name in _BLOCK_ORDER:
        if not _block_inputs_present(
            block_name,
            analysis_summary_str,
            cluster_annotations,
            deg_validation,
        ):
            logger.info(
                "narrative_generator: skipping block '%s' — required inputs missing",
                block_name,
            )
            continue

        inputs = {
            "tissue": tissue,
            "disease_context": disease_context,
            "biological_question": biological_question,
            "analysis_summary": analysis_summary_str or "not available",
            "pipeline_advice": pipeline_advice_str,
            "cluster_annotations": cluster_annotations_str,
            "deg_validation": deg_validation_str,
            "coherence_review": coherence_review_str,
            "narrative_block": block_name,
        }

        try:
            raw = call_llm(
                skill_name="narrative_generator",
                inputs=inputs,
                config=config,
                log_dir=log_dir,
                module="narrative_generator",
                runtime_ai=runtime_ai,
            )
        except Exception as exc:
            logger.warning(
                "narrative_generator: LLM call failed for block '%s': %s",
                block_name,
                exc,
            )
            continue

        block = _parse_block_response(raw, block_name)
        if block is None:
            continue
        blocks.append(block)

    # --- assemble result ---
    overall_groundedness = (
        sum(b.groundedness_score for b in blocks) / len(blocks)
        if blocks
        else 0.0
    )

    # Extract skill version from skill YAML if loadable, else default
    skill_version = "1.0"
    try:
        from ai._skill_loader import load_skill
        _, _ = load_skill("narrative_generator", **{
            "tissue": tissue,
            "disease_context": disease_context,
            "biological_question": biological_question,
            "analysis_summary": "{}",
            "pipeline_advice": "not available",
            "cluster_annotations": "not available",
            "deg_validation": "not available",
            "coherence_review": "not available",
            "narrative_block": "qc_rationale",
        })
    except Exception:
        pass  # version stays "1.0"

    result = NarrativeResult(
        blocks=blocks,
        overall_groundedness=round(overall_groundedness, 4),
        # AiResult base fields
        timestamp=datetime.now(timezone.utc).isoformat(),
        model=model,
        provider=provider,
        skill_name="narrative_generator",
        skill_version=skill_version,
        reasoning=(
            f"Generated {len(blocks)} of {len(_BLOCK_ORDER)} narrative blocks."
        ),
    )

    if output_path is not None:
        try:
            _write_markdown(result, output_path)
        except Exception as exc:
            logger.warning("narrative_generator: failed to write markdown: %s", exc)

    return result
