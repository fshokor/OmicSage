"""
ai/coherence_reviewer.py
B3 — Coherence Reviewer

Reads the completed analysis and acts as a senior bioinformatician reviewer.
Identifies inconsistencies, surprises, and sub-clustering candidates.

Also owns build_analysis_summary() — the pure function that compresses the
full run into ~2000 tokens for consumption by B3, C1, and C2.

Public API
----------
build_analysis_summary(adata, config, study_context) -> dict
run(adata, config, study_context, *, summary_path, log_dir, runtime_ai) -> CoherenceReview | None
"""

from __future__ import annotations

import json
import logging
import warnings
from dataclasses import dataclass, field
from pathlib import Path

from ai._base import AiResult
from ai._config_gate import AiDisabledError, check_ai_enabled
from ai._llm_client import call_llm

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class CoherenceFlag:
    category: str = ""    # qc | clustering | annotation | deg | pathway
    severity: str = ""    # info | warning | critical
    description: str = ""
    suggestion: str = ""


@dataclass
class CoherenceReview(AiResult):
    flags: list[CoherenceFlag] = field(default_factory=list)
    sub_clustering_candidates: list[str] = field(default_factory=list)
    rare_cell_candidates: list[str] = field(default_factory=list)
    overall_assessment: str = ""


# ---------------------------------------------------------------------------
# analysis_summary.json builder  (pure — no LLM, no side effects)
# ---------------------------------------------------------------------------

def build_analysis_summary(adata, config: dict, study_context: dict) -> dict:
    """
    Build the analysis_summary dict from adata.

    This is a pure function — no LLM calls, no file I/O, no side effects.
    It can be tested without mocking anything.

    All fields are optional except study_context.tissue and clustering.n_clusters.
    Missing fields are written as null, never omitted entirely.

    Token budget rules
    ------------------
    - top_degs : max 3 genes per comparison
    - pathways  : max 3 pathways per cluster
    """
    # ---- study_context ------------------------------------------------
    ds = study_context.get("dataset", {}) or {}
    disease_block = study_context.get("disease", {}) or {}
    exp_block = study_context.get("experiment", {}) or {}

    tissue = ds.get("tissue") or study_context.get("tissue") or "unknown"
    disease = disease_block.get("context") or study_context.get("disease_context") or None
    batch_key = exp_block.get("batch_key") or config.get("batch_key") or None

    n_batches = 1
    if batch_key and batch_key in adata.obs.columns:
        n_batches = int(adata.obs[batch_key].nunique())

    sc_block = {
        "tissue": tissue,
        "disease": disease,
        "n_cells": int(adata.n_obs),
        "n_batches": n_batches,
    }

    # ---- QC decisions -------------------------------------------------
    omicsage_qc = adata.uns.get("omicsage_qc") or {}
    qc_block = {
        "min_genes": omicsage_qc.get("min_genes", None),
        "max_mt_pct": omicsage_qc.get("max_mt_pct", None),
        "cells_removed_pct": omicsage_qc.get("cells_removed_pct", None),
    }

    # ---- clustering ---------------------------------------------------
    resolution = (
        adata.uns.get("leiden_resolution")
        or (adata.uns.get("omicsage_cluster") or {}).get("best_resolution")
        or (adata.uns.get("omicsage_cluster") or {}).get("resolution")  # fallback alias
        or None
    )
    leiden_labels = adata.obs.get("leiden")
    n_clusters = int(leiden_labels.nunique()) if leiden_labels is not None else None

    # silhouette score — stored by clustering module if computed
    silhouette = (adata.uns.get("omicsage_cluster") or {}).get("silhouette_score", None)

    cluster_block = {
        "resolution": float(resolution) if resolution is not None else None,
        "n_clusters": n_clusters,
        "silhouette_score": float(silhouette) if silhouette is not None else None,
    }

    # ---- cell types ---------------------------------------------------
    cell_types_list = []
    if leiden_labels is not None:
        ai_cell_type = adata.obs.get("ai_cell_type")
        ai_confidence = adata.obs.get("ai_confidence")
        for cluster_id in sorted(leiden_labels.unique(), key=lambda x: str(x)):
            mask = leiden_labels == cluster_id
            n_cells = int(mask.sum())
            ct = str(ai_cell_type[mask].mode()[0]) if (ai_cell_type is not None and mask.sum() > 0) else None
            conf = str(ai_confidence[mask].mode()[0]) if (ai_confidence is not None and mask.sum() > 0) else None
            cell_types_list.append({
                "cluster": str(cluster_id),
                "cell_type": ct,
                "confidence": conf,
                "n_cells": n_cells,
            })

    # ---- top DEGs (max 3 per comparison) ------------------------------
    top_degs = []
    rank = adata.uns.get("rank_genes_groups")
    if rank is not None:
        try:
            groups = list(rank["names"].dtype.names)
            for grp in groups:
                gene_names = [str(rank["names"][i][grp]) for i in range(len(rank["names"]))]
                logfc_vals = [float(rank["logfoldchanges"][i][grp]) for i in range(len(rank["logfoldchanges"]))]
                for gene, lfc in zip(gene_names[:3], logfc_vals[:3]):
                    top_degs.append({
                        "comparison": grp,
                        "cluster": grp,
                        "gene": gene,
                        "log2fc": round(lfc, 3),
                    })
        except Exception as exc:
            logger.warning("build_analysis_summary: could not parse rank_genes_groups: %s", exc)

    # ---- pathways (max 3 per cluster) ---------------------------------
    pathways = []
    gsea = adata.uns.get("gsea_results") or adata.uns.get("omicsage_gsea") or {}
    if isinstance(gsea, dict):
        for cluster_id, results in gsea.items():
            if isinstance(results, list):
                for entry in results[:3]:
                    pathways.append({
                        "cluster": str(cluster_id),
                        "pathway": entry.get("pathway") or entry.get("Term", ""),
                        "padj": entry.get("padj") or entry.get("Adjusted P-value", None),
                    })

    return {
        "study_context": sc_block,
        "qc_decisions": qc_block,
        "clustering": cluster_block,
        "cell_types": cell_types_list,
        "top_degs": top_degs,
        "pathways": pathways,
    }


# ---------------------------------------------------------------------------
# Main run function
# ---------------------------------------------------------------------------

def run(
    adata,
    config: dict,
    study_context: dict,
    *,
    summary_path: str | None = None,
    log_dir: str = "logs/llm",
    runtime_ai: bool = True,
) -> "CoherenceReview | None":
    """
    Run the coherence reviewer.

    Parameters
    ----------
    adata        : AnnData with completed analysis results
    config       : pipeline config dict (must contain ai section)
    study_context: loaded study_context.yaml as dict
    summary_path : if set, write analysis_summary.json to this path
    log_dir      : directory for audit JSONL logs
    runtime_ai   : False if --ai flag was not passed at runtime

    Returns
    -------
    CoherenceReview on success, None if AI is disabled or an unrecoverable
    error occurs.
    """
    # --- gate -----------------------------------------------------------
    try:
        check_ai_enabled(config, module="coherence_reviewer", runtime_ai=runtime_ai)
    except AiDisabledError:
        return None

    # --- build analysis summary (pure, no LLM) -------------------------
    summary = build_analysis_summary(adata, config, study_context)

    # --- optionally write analysis_summary.json ------------------------
    if summary_path is not None:
        try:
            out = Path(summary_path)
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(json.dumps(summary, indent=2))
            logger.info("analysis_summary.json written to %s", summary_path)
        except Exception as exc:
            logger.warning("Could not write analysis_summary.json: %s", exc)

    # --- extract context for prompts -----------------------------------
    ds = study_context.get("dataset", {}) or {}
    disease_block = study_context.get("disease", {}) or {}
    tissue = ds.get("tissue") or study_context.get("tissue") or "unknown"
    disease_context = (
        disease_block.get("context") or study_context.get("disease_context") or "none"
    )

    # --- call LLM -------------------------------------------------------
    raw_response = call_llm(
        skill_name="coherence_reviewer",
        inputs={
            "tissue": tissue,
            "disease_context": disease_context,
            "analysis_summary": json.dumps(summary),
        },
        config=config,
        log_dir=log_dir,
        module="coherence_reviewer",
        runtime_ai=runtime_ai,
    )

    # --- parse response -------------------------------------------------
    return _parse_response(raw_response, config)


# ---------------------------------------------------------------------------
# Response parser
# ---------------------------------------------------------------------------

def _parse_response(raw: str, config: dict) -> CoherenceReview:
    """
    Parse the LLM JSON response into a CoherenceReview.

    On any parse failure: log a warning and return a CoherenceReview with
    empty flags. Never raises — graceful degradation.
    """
    ai_cfg = config.get("ai", {})
    model = ai_cfg.get("model", "llama3")
    provider = ai_cfg.get("provider", "ollama")

    base_kwargs = dict(
        model=model,
        provider=provider,
        skill_name="coherence_reviewer",
        skill_version="1.0",
        reasoning="",
    )

    try:
        from ai._json_utils import extract_json_from_response
        data = json.loads(extract_json_from_response(raw))
    except (json.JSONDecodeError, ValueError) as exc:
        logger.warning("coherence_reviewer: JSON parse failed: %s", exc)
        return CoherenceReview(**base_kwargs)

    # ---- flags ---------------------------------------------------------
    flags = []
    for raw_flag in data.get("flags") or []:
        if not isinstance(raw_flag, dict):
            continue
        flags.append(CoherenceFlag(
            category=str(raw_flag.get("category", "")),
            severity=str(raw_flag.get("severity", "")),
            description=str(raw_flag.get("description", "")),
            suggestion=str(raw_flag.get("suggestion", "")),
        ))

    # ---- candidates ----------------------------------------------------
    sub = [str(x) for x in (data.get("sub_clustering_candidates") or [])]
    rare = [str(x) for x in (data.get("rare_cell_candidates") or [])]
    assessment = str(data.get("overall_assessment") or "")

    return CoherenceReview(
        **base_kwargs,
        flags=flags,
        sub_clustering_candidates=sub,
        rare_cell_candidates=rare,
        overall_assessment=assessment,
    )
