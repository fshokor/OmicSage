"""
ai/pipeline_advisor.py
======================
Phase 3 — Session 1: A1 Pipeline Advisor

Runs once per project, before any analysis step.
Reads study_context + loaded data properties, then:
  1. Runs rule-based pre-checks (fast, no API cost)
  2. Calls the LLM for rationale + any remaining recommendations
  3. Returns a PipelineAdvice dataclass (or None if AI is disabled)

Public API
----------
run(adata_or_mdata, config, study_context, *, log_dir, runtime_ai) -> PipelineAdvice | None
load_study_context(path) -> dict
"""

from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from ai._audit_log import write_audit_record
from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm
from ai._skill_loader import load_skill

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class StepRecommendation:
    """A single recommended pipeline step."""
    step_name: str
    priority: str   # "required" | "recommended" | "optional"
    rationale: str


@dataclass
class PipelineAdvice(AiResult):
    """Full output of the Pipeline Advisor module."""
    recommended_steps: list[StepRecommendation] = field(default_factory=list)
    inferred_biological_question: str | None = None
    warnings: list[str] = field(default_factory=list)

    # inherited from AiResult: timestamp, model, provider, skill_version, reasoning


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_VALID_PRIORITIES = {"required", "recommended", "optional"}

_VALID_STEP_NAMES = {
    "qc", "normalization", "hvg_selection", "pca", "neighbor_graph",
    "umap", "clustering", "batch_correction", "cell_type_annotation",
    "wilcoxon_deg", "pseudobulk_deg", "gsea", "trajectory",
    "cell_communication", "wnn_integration",
}


def load_study_context(path: str | Path) -> dict:
    """Load and return a study_context.yaml as a plain dict.

    Returns an empty dict if the file is missing or empty.
    """
    p = Path(path)
    if not p.exists():
        logger.warning("study_context not found at %s — using empty context", path)
        return {}
    with open(p, encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def _extract_data_properties(adata_or_mdata: Any) -> dict:
    """Pull n_cells, n_genes, modalities, n_batches, n_donors, n_conditions
    from an AnnData or MuData object without modifying it.
    """
    try:
        import anndata
        import mudata
    except ImportError:
        pass

    props: dict[str, Any] = {}

    # MuData
    if hasattr(adata_or_mdata, "mod"):
        mdata = adata_or_mdata
        props["n_cells"] = mdata.n_obs
        # sum features across modalities
        props["n_genes"] = sum(m.n_vars for m in mdata.mod.values())
        props["modalities"] = [k.upper() for k in mdata.mod.keys()]
        # use first modality for batch/donor/condition sniffing
        primary = next(iter(mdata.mod.values()))
        obs = primary.obs
    else:
        adata = adata_or_mdata
        props["n_cells"] = adata.n_obs
        props["n_genes"] = adata.n_vars
        props["modalities"] = ["RNA"]
        obs = adata.obs

    # n_batches — count unique values of any plausible batch column
    batch_cols = ["batch", "sample", "sample_id", "donor", "donor_id", "library_id"]
    n_batches = 1
    for col in batch_cols:
        if col in obs.columns:
            n_batches = int(obs[col].nunique())
            break
    props["n_batches"] = n_batches

    # n_donors
    donor_cols = ["donor", "donor_id", "subject", "patient", "individual"]
    n_donors = 1
    for col in donor_cols:
        if col in obs.columns:
            n_donors = int(obs[col].nunique())
            break
    props["n_donors"] = n_donors

    # n_conditions
    cond_cols = ["condition", "group", "disease", "treatment", "label"]
    n_conditions = 1
    for col in cond_cols:
        if col in obs.columns:
            n_conditions = int(obs[col].nunique())
            break
    props["n_conditions"] = n_conditions

    return props


def _rule_based_checks(
    props: dict,
    config: dict,
    study_context: dict,
) -> tuple[list[str], list[str]]:
    """Fast rule-based pre-checks — no API cost.

    Returns
    -------
    warnings : list[str]
    flagged_steps : list[str]   step names the rules recommend
    """
    warnings_out: list[str] = []
    flagged_steps: list[str] = []

    n_cells = props.get("n_cells", 0)
    n_batches = props.get("n_batches", 1)
    n_donors = props.get("n_donors", 1)
    n_conditions = props.get("n_conditions", 1)
    modalities = [m.upper() for m in props.get("modalities", ["RNA"])]

    # Config-level batch_key
    batch_key = (
        config.get("preprocessing", {}).get("batch_key")
        or config.get("batch_key")
        or study_context.get("experiment", {}).get("batch_key")
    )

    # Rule 1: multiple batches detected but no batch_key configured
    if n_batches > 1:
        flagged_steps.append("batch_correction")
        if not batch_key:
            warnings_out.append(
                f"Detected {n_batches} batches but no batch_key is set in config "
                "or study_context — batch correction may not work correctly. "
                "Set experiment.batch_key in study_context.yaml."
            )

    # Rule 2: pseudobulk preferred when enough donors + conditions
    if n_donors > 2 and n_conditions > 1:
        flagged_steps.append("pseudobulk_deg")
        flagged_steps.append("gsea")
    else:
        if n_conditions > 1:
            flagged_steps.append("wilcoxon_deg")
            if n_donors <= 2:
                warnings_out.append(
                    f"Only {n_donors} donor(s) detected — pseudobulk DEG is unreliable "
                    "with fewer than 3 donors. Wilcoxon will be used instead."
                )

    # Rule 3: ADT modality → WNN (future phase)
    if "ADT" in modalities:
        flagged_steps.append("wnn_integration")
        warnings_out.append(
            "ADT modality detected. WNN integration (Seurat v5/Muon) is planned "
            "for Phase 6 and is not yet implemented. RNA-only steps will run now."
        )

    # Rule 4: very small dataset
    if n_cells < 500:
        warnings_out.append(
            f"Only {n_cells} cells detected — Leiden clustering and doublet detection "
            "are unreliable below 500 cells. Verify thresholds carefully."
        )

    # Always required steps
    for step in ["qc", "normalization", "hvg_selection", "pca",
                 "neighbor_graph", "umap", "clustering",
                 "cell_type_annotation"]:
        if step not in flagged_steps:
            flagged_steps.append(step)

    return warnings_out, flagged_steps


def _parse_llm_response(raw: str, module: str) -> dict | None:
    """Parse the raw LLM string as JSON.  Returns None on failure."""
    try:
        from ai._json_utils import extract_json_from_response
        return json.loads(extract_json_from_response(raw))
    except (json.JSONDecodeError, ValueError) as exc:
        logger.warning("[%s] Failed to parse LLM JSON: %s", module, exc)
        return None


def _build_advice_from_parsed(parsed: dict) -> tuple[list[StepRecommendation], str | None, list[str], str]:
    """Convert parsed JSON dict to typed fields.  Validates priorities."""
    steps = []
    for raw_step in parsed.get("recommended_steps", []):
        priority = raw_step.get("priority", "recommended")
        if priority not in _VALID_PRIORITIES:
            priority = "recommended"
        step_name = raw_step.get("step_name", "")
        if step_name not in _VALID_STEP_NAMES:
            logger.warning("LLM returned unknown step_name '%s' — skipping", step_name)
            continue
        steps.append(StepRecommendation(
            step_name=step_name,
            priority=priority,
            rationale=raw_step.get("rationale", ""),
        ))

    inferred_q = parsed.get("inferred_biological_question") or None
    extra_warnings = parsed.get("warnings") or []
    reasoning = parsed.get("reasoning", "")
    return steps, inferred_q, extra_warnings, reasoning


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run(
    adata_or_mdata: Any,
    config: dict,
    study_context: dict,
    *,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> PipelineAdvice | None:
    """Run the Pipeline Advisor.

    Parameters
    ----------
    adata_or_mdata : AnnData | MuData
        Loaded (but not yet QC-filtered) data object.  Read-only.
    config : dict
        Full pipeline config dict (loaded from config.yaml or similar).
    study_context : dict
        Loaded study_context.yaml content.
    log_dir : str
        Directory for audit log JSONL files.
    runtime_ai : bool
        False when the pipeline was invoked without the --ai flag.

    Returns
    -------
    PipelineAdvice | None
        None if AI is disabled at any level, or if the LLM response
        cannot be parsed. Never raises — real errors are logged.
    """
    # Gate — returns None at all three levels: global off, module off, runtime off
    try:
        check_ai_enabled(config, module="pipeline_advisor", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # Extract data properties
    props = _extract_data_properties(adata_or_mdata)

    # Rule-based pre-checks (no API cost)
    rule_warnings, rule_steps = _rule_based_checks(props, config, study_context)

    # Build skill inputs
    sc = study_context  # shorthand
    tissue = (sc.get("dataset") or {}).get("tissue") or sc.get("tissue") or "unknown"
    disease_context = (sc.get("disease") or {}).get("context") or sc.get("disease_context") or None
    experiment_design = (sc.get("experiment") or {}).get("design") or sc.get("experiment_design") or None
    biological_question = sc.get("biological_question") or None

    skill_inputs = {
        "tissue": tissue,
        "disease_context": disease_context or "not specified",
        "experiment_design": experiment_design or "not specified",
        "biological_question": biological_question or "",
        "n_cells": props["n_cells"],
        "n_genes": props["n_genes"],
        "modalities": props["modalities"],
        "n_batches": props["n_batches"],
        "n_donors": props["n_donors"],
        "n_conditions": props["n_conditions"],
        "rule_based_warnings": rule_warnings,
        "rule_based_steps": rule_steps,
    }

    # Load skill
    skill_dir = Path(__file__).parent / "skills"
    system_prompt, user_prompt = load_skill("pipeline_advisor", skill_dir=skill_dir, **skill_inputs)

    # Call LLM
    raw_response, token_usage, model_name, provider = call_llm(
        config=config,
        system_prompt=system_prompt,
        user_prompt=user_prompt,
    )

    # Parse response
    parsed = _parse_llm_response(raw_response, module="pipeline_advisor")
    parse_success = parsed is not None

    # Audit log — always, regardless of parse success
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    write_audit_record(
        log_dir=log_dir,
        module="pipeline_advisor",
        skill_version="1.0",
        model=model_name,
        provider=provider,
        input_summary={
            "tissue": tissue,
            "n_cells": props["n_cells"],
            "n_batches": props["n_batches"],
            "n_donors": props["n_donors"],
            "n_conditions": props["n_conditions"],
            "modalities": props["modalities"],
        },
        token_usage=token_usage,
        raw_response=raw_response,
        parsed_output=parsed,
        parse_success=parse_success,
    )

    if not parse_success:
        logger.warning("[pipeline_advisor] LLM response could not be parsed — returning None")
        return None

    # Build typed result
    steps, inferred_q, extra_warnings, reasoning = _build_advice_from_parsed(parsed)

    # If biological_question was blank and LLM didn't infer one, note it
    if not biological_question and not inferred_q:
        logger.info("[pipeline_advisor] biological_question was blank and LLM did not infer one")

    import datetime
    return PipelineAdvice(
        timestamp=datetime.datetime.utcnow().isoformat() + "Z",
        model=model_name,
        provider=provider,
        skill_version="1.0",
        reasoning=reasoning,
        recommended_steps=steps,
        inferred_biological_question=inferred_q,
        warnings=rule_warnings + extra_warnings,
    )
