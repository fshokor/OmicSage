"""
OmicSage — AI Cluster Annotator (Phase 3, Session 3 — B1)

Maps top marker genes per cluster to cell type labels using an LLM.
Writes results back to adata.obs:
  - adata.obs['ai_cell_type']    — predicted label
  - adata.obs['ai_confidence']   — high | medium | low
  - adata.obs['ai_alt_types']    — comma-separated alternative types

One LLM call per cluster. Partial failures are logged and skipped;
the remaining clusters are still written. A module-level failure
(missing rank_genes_groups, disabled AI) returns None cleanly.
"""

from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass, field
from datetime import datetime, timezone

from ai._base import AiResult
from ai._audit_log import write_audit_record
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm

logger = logging.getLogger(__name__)

_VALID_CONFIDENCE = {"high", "medium", "low"}
_REQUIRED_OUTPUT_KEYS = {
    "cell_type",
    "confidence",
    "supporting_markers",
    "alternative_types",
    "recommended_db",
    "manual_marker_set",
    "reasoning",
}


# ---------------------------------------------------------------------------
# Public result dataclass
# ---------------------------------------------------------------------------

@dataclass
class ClusterAnnotation(AiResult):
    """Annotation result for a single cluster."""
    cluster_id: str = ""
    cell_type: str = "unknown"
    confidence: str = "low"                         # "high" | "medium" | "low"
    supporting_markers: list[str] = field(default_factory=list)
    alternative_types: list[str] = field(default_factory=list)
    recommended_db: str = ""
    manual_marker_set: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_cluster_column(adata) -> str:
    """Return the name of the clustering obs column to use.

    Prefers 'leiden', then 'louvain', then the first column whose name
    contains 'cluster'. Raises ValueError if none found.
    """
    for preferred in ("leiden", "louvain"):
        if preferred in adata.obs.columns:
            return preferred
    for col in adata.obs.columns:
        if "cluster" in col.lower():
            return col
    raise ValueError(
        "No clustering column found in adata.obs. "
        "Expected 'leiden' or 'louvain'. Run clustering before cluster annotation."
    )


def _extract_markers(adata, cluster_id: str, n_markers: int) -> list[str]:
    """Extract top n_markers for a cluster from adata.uns['rank_genes_groups'].

    Returns a list of gene symbols ranked by score (log2FC proxy).
    """
    rgg = adata.uns["rank_genes_groups"]
    # rank_genes_groups stores structured arrays keyed by group name
    # Group names may be stored as str or int — normalise to str for lookup
    names = rgg["names"]
    # Find the column index matching this cluster_id
    group_names = list(names.dtype.names)
    # dtype.names gives the column labels (one per cluster/group)
    if cluster_id not in group_names:
        # Try int conversion (leiden stores groups as str "0", "1", …)
        logger.warning(
            "Cluster '%s' not found in rank_genes_groups. Available: %s",
            cluster_id,
            group_names,
        )
        return []
    markers = [str(g) for g in names[cluster_id][:n_markers]]
    return markers


def _parse_response(raw: str, cluster_id: str) -> dict | None:
    """Parse and validate LLM JSON response for one cluster.

    Returns the parsed dict on success, None on any parse/validation failure.
    Logs a warning (not an error) on failure so the caller can skip gracefully.
    """
    # Strip markdown fences if the model wrapped the JSON despite instructions
    text = raw.strip()
    if text.startswith("```"):
        lines = text.splitlines()
        # Drop opening fence (```json or ```) and closing fence
        text = "\n".join(
            line for line in lines
            if not line.strip().startswith("```")
        ).strip()

    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        logger.warning(
            "Cluster %s: JSON parse failed — %s. Raw response: %.200s",
            cluster_id, exc, raw,
        )
        return None

    if not isinstance(parsed, dict):
        logger.warning("Cluster %s: expected JSON object, got %s", cluster_id, type(parsed))
        return None

    missing = _REQUIRED_OUTPUT_KEYS - parsed.keys()
    if missing:
        logger.warning("Cluster %s: missing keys in LLM response: %s", cluster_id, missing)
        return None

    # Normalise confidence to one of the valid values
    confidence = str(parsed.get("confidence", "low")).strip().lower()
    if confidence not in _VALID_CONFIDENCE:
        logger.warning(
            "Cluster %s: invalid confidence '%s', defaulting to 'low'", cluster_id, confidence
        )
        confidence = "low"
    parsed["confidence"] = confidence

    return parsed


def _annotation_from_parsed(parsed: dict, cluster_id: str, config: dict) -> ClusterAnnotation:
    """Build a ClusterAnnotation dataclass from a validated parsed dict."""
    ai_cfg = config.get("ai", {})
    return ClusterAnnotation(
        # AiResult base fields
        timestamp=datetime.now(timezone.utc).isoformat(),
        model=ai_cfg.get("model", "llama3"),
        provider=ai_cfg.get("provider", "ollama"),
        skill_name="cluster_annotator",
        skill_version="1.0",
        reasoning=parsed.get("reasoning", ""),
        # ClusterAnnotation fields
        cluster_id=cluster_id,
        cell_type=str(parsed.get("cell_type", "unknown")),
        confidence=parsed["confidence"],
        supporting_markers=list(parsed.get("supporting_markers", [])),
        alternative_types=list(parsed.get("alternative_types", [])),
        recommended_db=str(parsed.get("recommended_db", "")),
        manual_marker_set=list(parsed.get("manual_marker_set", [])),
    )


def _write_obs_columns(adata, annotations: list[ClusterAnnotation], cluster_col: str) -> None:
    """Write ai_cell_type, ai_confidence, ai_alt_types to adata.obs.

    Cells in clusters that failed annotation receive "unknown", "low", "".
    """
    annotation_map = {ann.cluster_id: ann for ann in annotations}

    cell_types = []
    confidences = []
    alt_types = []

    for cluster_id in adata.obs[cluster_col].astype(str):
        ann = annotation_map.get(cluster_id)
        if ann is not None:
            cell_types.append(ann.cell_type)
            confidences.append(ann.confidence)
            alt_types.append(", ".join(ann.alternative_types))
        else:
            cell_types.append("unknown")
            confidences.append("low")
            alt_types.append("")

    adata.obs["ai_cell_type"] = cell_types
    adata.obs["ai_confidence"] = confidences
    adata.obs["ai_alt_types"] = alt_types


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run(
    adata,
    config: dict,
    study_context: dict,
    *,
    n_markers: int = 20,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> list[ClusterAnnotation] | None:
    """Annotate all clusters in adata using an LLM.

    Parameters
    ----------
    adata:
        AnnData object after marker gene computation. Must contain
        ``adata.uns['rank_genes_groups']`` and a clustering obs column.
        **This function writes to adata.obs in-place.**
    config:
        Full pipeline config dict (must include ``config['ai']`` section).
    study_context:
        Contents of study_context.yaml as a dict.
    n_markers:
        Number of top marker genes per cluster to pass to the LLM.
    log_dir:
        Directory for audit log JSONL files.
    runtime_ai:
        Set to False to disable AI at runtime regardless of config
        (equivalent to running without the ``--ai`` flag).

    Returns
    -------
    list[ClusterAnnotation] or None
        None if AI is disabled or rank_genes_groups is missing.
        Otherwise a list with one entry per successfully annotated cluster.
        Clusters that fail to parse are skipped and logged as warnings.
    """
    # --- Gate: return None silently if AI disabled at any level ---
    try:
        check_ai_enabled(config, module="cluster_annotator", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- Validate required input ---
    if "rank_genes_groups" not in adata.uns:
        raise ValueError(
            "adata.uns['rank_genes_groups'] not found. "
            "Run sc.tl.rank_genes_groups() before calling cluster_annotator.run()."
        )

    cluster_col = _get_cluster_column(adata)
    cluster_ids = [str(c) for c in adata.obs[cluster_col].unique()]

    tissue = study_context.get("dataset", {}).get("tissue", "unknown")
    disease_context = study_context.get("disease", {}).get("context", "healthy") or "healthy"
    ai_cfg = config.get("ai", {})
    model = ai_cfg.get("model", "llama3")
    provider = ai_cfg.get("provider", "ollama")

    annotations: list[ClusterAnnotation] = []

    for cluster_id in sorted(cluster_ids):
        markers = _extract_markers(adata, cluster_id, n_markers)
        if not markers:
            logger.warning("Cluster %s: no markers found, skipping.", cluster_id)
            continue

        n_cells = int((adata.obs[cluster_col].astype(str) == cluster_id).sum())

        skill_inputs = {
            "tissue": tissue,
            "disease_context": disease_context,
            "cluster_id": cluster_id,
            "n_cells": n_cells,
            "marker_genes": markers,
        }

        input_summary = {
            "tissue": tissue,
            "cluster_id": cluster_id,
            "n_cells": n_cells,
            "n_markers": len(markers),
        }

        # --- LLM call ---
        try:
            raw_response = call_llm(
                skill_name="cluster_annotator",
                inputs=skill_inputs,
                config=config,
                log_dir=log_dir,
                module="cluster_annotator",
                runtime_ai=runtime_ai,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning("Cluster %s: LLM call failed — %s", cluster_id, exc)
            write_audit_record(
                log_dir=log_dir,
                module="cluster_annotator",
                skill_version="1.0",
                model=model,
                provider=provider,
                input_summary=input_summary,
                token_usage=None,
                raw_response=str(exc),
                parsed_output=None,
                parse_success=False,
            )
            continue

        # --- Parse ---
        parsed = _parse_response(raw_response, cluster_id)
        parse_success = parsed is not None

        write_audit_record(
            log_dir=log_dir,
            module="cluster_annotator",
            skill_version="1.0",
            model=model,
            provider=provider,
            input_summary=input_summary,
            token_usage=None,
            raw_response=raw_response,
            parsed_output=parsed,
            parse_success=parse_success,
        )

        if not parse_success:
            logger.warning(
                "Cluster %s: parse failed, skipping this cluster. "
                "Other clusters are unaffected.",
                cluster_id,
            )
            continue

        ann = _annotation_from_parsed(parsed, cluster_id, config)
        annotations.append(ann)

    # --- Write obs columns for all successfully annotated clusters ---
    if annotations:
        _write_obs_columns(adata, annotations, cluster_col)
    else:
        warnings.warn(
            "cluster_annotator: no clusters were successfully annotated. "
            "adata.obs columns were not written.",
            UserWarning,
            stacklevel=2,
        )

    return annotations
