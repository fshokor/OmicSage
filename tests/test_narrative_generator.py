"""
tests/test_narrative_generator.py
Phase 3 — Session 7 — C1: Biological Narrative Generator

All tests run without a real API key.
call_llm is patched to return a controlled JSON string.
"""

from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pytest

from ai.narrative_generator import (
    NarrativeBlock,
    NarrativeResult,
    _parse_block_response,
    run,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def config_ai_on():
    return {
        "ai": {
            "features": True,
            "provider": "ollama",
            "model": "llama3",
            "modules": {},
        }
    }


@pytest.fixture
def config_ai_off():
    return {"ai": {"features": False, "provider": "ollama", "model": "llama3"}}


@pytest.fixture
def config_module_off():
    return {
        "ai": {
            "features": True,
            "provider": "ollama",
            "model": "llama3",
            "modules": {"narrative_generator": False},
        }
    }


@pytest.fixture
def study_context():
    return {
        "dataset": {"tissue": "liver", "species": "human"},
        "disease": {"context": "hepatocellular carcinoma"},
        "biological_question": "Characterise the TME of HCC.",
    }


@pytest.fixture
def adata_with_summary():
    adata = ad.AnnData(X=np.zeros((10, 5)))
    adata.uns["omicsage_analysis_summary"] = json.dumps({
        "study_context": {"tissue": "liver", "n_cells": 10},
        "qc_decisions": {"min_genes": 200, "max_mt_pct": 20, "cells_removed_pct": 8.5},
        "clustering": {"resolution": 0.5, "n_clusters": 4, "silhouette_score": 0.42},
        "cell_types": [{"cluster": 0, "cell_type": "CD8+ T cell", "n_cells": 5}],
        "top_degs": [{"comparison": "tumour_vs_normal", "cluster": 0, "gene": "PDCD1", "log2fc": 2.1}],
        "pathways": [{"cluster": 0, "pathway": "T cell exhaustion", "padj": 0.001}],
    })
    return adata


@pytest.fixture
def adata_no_summary():
    return ad.AnnData(X=np.zeros((10, 5)))


def _mock_block_response(block_name: str, score: float = 0.9) -> str:
    return json.dumps({
        "narrative_text": f"Narrative text for {block_name}. "
                          f"Cells removed: 8.5%. Gene PDCD1 log2FC=2.1. PMID:12345678.",
        "cited_evidence": ["cells_removed_pct=8.5", "PDCD1", "PMID:12345678"],
        "groundedness_score": score,
        "reasoning": f"Wrote {block_name} block citing available metrics.",
    })


def _make_mock_llm(score: float = 0.9):
    """Return a mock call_llm that infers block name from inputs."""
    def _mock(skill_name, inputs, config, *, log_dir, module, runtime_ai):
        return _mock_block_response(inputs.get("narrative_block", "unknown"), score)
    return _mock


# ---------------------------------------------------------------------------
# 1. Returns None when ai_features=False
# ---------------------------------------------------------------------------

def test_returns_none_when_ai_features_off(config_ai_off, study_context, adata_no_summary):
    result = run(adata_no_summary, config_ai_off, study_context)
    assert result is None


# ---------------------------------------------------------------------------
# 2. Returns None when runtime_ai=False
# ---------------------------------------------------------------------------

def test_returns_none_when_runtime_ai_false(config_ai_on, study_context, adata_no_summary):
    result = run(adata_no_summary, config_ai_on, study_context, runtime_ai=False)
    assert result is None


# ---------------------------------------------------------------------------
# 3. Returns NarrativeResult when mock LLM returns valid JSON
# ---------------------------------------------------------------------------

def test_returns_narrative_result_on_valid_mock(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"], "unexpected_genes": []}

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
        )

    assert isinstance(result, NarrativeResult)


# ---------------------------------------------------------------------------
# 4. All four blocks generated when all upstream inputs provided
# ---------------------------------------------------------------------------

def test_all_four_blocks_generated(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"], "unexpected_genes": []}
    mock_coherence_review = {"flags": [], "overall_assessment": "Looks good."}
    mock_pipeline_advice = {"recommended_steps": [], "warnings": []}

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            pipeline_advice=mock_pipeline_advice,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
            coherence_review=mock_coherence_review,
        )

    assert result is not None
    block_names = [b.block_name for b in result.blocks]
    assert "qc_rationale" in block_names
    assert "cell_type_landscape" in block_names
    assert "differential_expression" in block_names
    assert "interpretation" in block_names
    assert len(result.blocks) == 4


# ---------------------------------------------------------------------------
# 5. qc_rationale block skipped when analysis_summary missing from adata
# ---------------------------------------------------------------------------

def test_qc_rationale_skipped_when_no_analysis_summary(
    config_ai_on, study_context, adata_no_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_no_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
        )

    assert result is not None
    block_names = [b.block_name for b in result.blocks]
    assert "qc_rationale" not in block_names


# ---------------------------------------------------------------------------
# 6. differential_expression block skipped when deg_validation=None
# ---------------------------------------------------------------------------

def test_deg_block_skipped_when_deg_validation_none(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=None,
        )

    assert result is not None
    block_names = [b.block_name for b in result.blocks]
    assert "differential_expression" not in block_names


# ---------------------------------------------------------------------------
# 7. NarrativeBlock fields populated correctly
# ---------------------------------------------------------------------------

def test_narrative_block_fields_populated(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm(0.88)):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
        )

    assert result is not None
    ct_block = next(b for b in result.blocks if b.block_name == "cell_type_landscape")
    assert ct_block.narrative_text != ""
    assert isinstance(ct_block.cited_evidence, list)
    assert len(ct_block.cited_evidence) > 0
    assert 0.0 <= ct_block.groundedness_score <= 1.0


# ---------------------------------------------------------------------------
# 8. overall_groundedness calculated correctly (mean across blocks)
# ---------------------------------------------------------------------------

def test_overall_groundedness_is_mean_of_blocks(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"]}

    scores = [0.8, 0.9, 0.7, 0.85]
    call_count = [0]

    def _mock_with_varying_scores(skill_name, inputs, config, *, log_dir, module, runtime_ai):
        idx = call_count[0]
        call_count[0] += 1
        block = inputs.get("narrative_block", "unknown")
        return json.dumps({
            "narrative_text": f"Text for {block}.",
            "cited_evidence": ["gene_A"],
            "groundedness_score": scores[idx] if idx < len(scores) else 0.8,
            "reasoning": "test",
        })

    with patch("ai.narrative_generator.call_llm", side_effect=_mock_with_varying_scores):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
        )

    assert result is not None
    expected_mean = round(sum(b.groundedness_score for b in result.blocks) / len(result.blocks), 4)
    assert result.overall_groundedness == expected_mean


# ---------------------------------------------------------------------------
# 9. overall_groundedness = 0.0 when no blocks generated
# ---------------------------------------------------------------------------

def test_overall_groundedness_zero_when_no_blocks(
    config_ai_on, study_context, adata_no_summary
):
    # No analysis_summary, no cluster_annotations, no deg_validation
    # → no blocks can be generated
    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(adata_no_summary, config_ai_on, study_context)

    assert result is not None
    assert result.blocks == []
    assert result.overall_groundedness == 0.0


# ---------------------------------------------------------------------------
# 10. ai_narrative.md written to output_path when provided
# ---------------------------------------------------------------------------

def test_markdown_written_to_output_path(
    tmp_path, config_ai_on, study_context, adata_with_summary
):
    output_file = str(tmp_path / "reports" / "ai_narrative.md")
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            output_path=output_file,
        )

    assert Path(output_file).exists()


# ---------------------------------------------------------------------------
# 11. ai_narrative.md has correct markdown structure
# ---------------------------------------------------------------------------

def test_markdown_structure(
    tmp_path, config_ai_on, study_context, adata_with_summary
):
    output_file = str(tmp_path / "ai_narrative.md")
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"]}

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
            output_path=output_file,
        )

    content = Path(output_file).read_text(encoding="utf-8")
    assert "# AI Biological Narrative" in content
    assert "## QC Rationale" in content
    assert "## Cell Type Landscape" in content
    assert "## Differential Expression" in content
    assert "## Interpretation and Perspectives" in content
    assert "*Generated:" in content
    assert "*Model:" in content
    assert "*Overall groundedness score:" in content


# ---------------------------------------------------------------------------
# 12. Degraded parse (invalid JSON for one block) skips that block, others continue
# ---------------------------------------------------------------------------

def test_invalid_json_for_one_block_skips_it(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"]}

    call_count = [0]

    def _mock_partial_failure(skill_name, inputs, config, *, log_dir, module, runtime_ai):
        call_count[0] += 1
        block = inputs.get("narrative_block")
        if block == "cell_type_landscape":
            return "THIS IS NOT VALID JSON {{{"
        return _mock_block_response(block)

    with patch("ai.narrative_generator.call_llm", side_effect=_mock_partial_failure):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
        )

    assert result is not None
    block_names = [b.block_name for b in result.blocks]
    assert "cell_type_landscape" not in block_names
    # other blocks still present
    assert "qc_rationale" in block_names
    assert "differential_expression" in block_names


# ---------------------------------------------------------------------------
# 13. AiResult base fields populated
# ---------------------------------------------------------------------------

def test_airesult_base_fields_populated(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            cluster_annotations=mock_cluster_annotations,
        )

    assert result is not None
    assert result.timestamp != ""
    assert result.model == "llama3"
    assert result.provider == "ollama"
    assert result.skill_name == "narrative_generator"
    assert result.skill_version != ""


# ---------------------------------------------------------------------------
# 14. pipeline_advice=None and coherence_review=None handled gracefully
# ---------------------------------------------------------------------------

def test_none_optional_inputs_handled_gracefully(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]

    with patch("ai.narrative_generator.call_llm", side_effect=_make_mock_llm()):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            pipeline_advice=None,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=None,
            coherence_review=None,
        )

    assert result is not None
    # Should not crash; qc and cell_type blocks still generated


# ---------------------------------------------------------------------------
# 15. Four separate call_llm calls made (one per block) when all inputs present
# ---------------------------------------------------------------------------

def test_four_separate_llm_calls_when_all_inputs_present(
    config_ai_on, study_context, adata_with_summary
):
    mock_cluster_annotations = [{"cluster": 0, "cell_type": "CD8+ T cell"}]
    mock_deg_validation = {"expected_genes": ["PDCD1"]}
    mock_coherence_review = {"flags": [], "overall_assessment": "OK"}
    mock_pipeline_advice = {"recommended_steps": []}

    call_tracker = []

    def _tracking_mock(skill_name, inputs, config, *, log_dir, module, runtime_ai):
        call_tracker.append(inputs.get("narrative_block"))
        return _mock_block_response(inputs.get("narrative_block", "unknown"))

    with patch("ai.narrative_generator.call_llm", side_effect=_tracking_mock):
        result = run(
            adata_with_summary,
            config_ai_on,
            study_context,
            pipeline_advice=mock_pipeline_advice,
            cluster_annotations=mock_cluster_annotations,
            deg_validation=mock_deg_validation,
            coherence_review=mock_coherence_review,
        )

    assert result is not None
    assert len(call_tracker) == 4
    assert set(call_tracker) == {
        "qc_rationale",
        "cell_type_landscape",
        "differential_expression",
        "interpretation",
    }


# ---------------------------------------------------------------------------
# Unit test: _parse_block_response
# ---------------------------------------------------------------------------

def test_parse_block_response_valid():
    raw = json.dumps({
        "narrative_text": "CD8+ T cells were the dominant population.",
        "cited_evidence": ["CD8+ T cells", "n_cells=412"],
        "groundedness_score": 0.92,
        "reasoning": "Cited cluster stats.",
    })
    block = _parse_block_response(raw, "cell_type_landscape")
    assert block is not None
    assert block.block_name == "cell_type_landscape"
    assert block.groundedness_score == 0.92
    assert len(block.cited_evidence) == 2


def test_parse_block_response_invalid_json():
    block = _parse_block_response("NOT JSON", "qc_rationale")
    assert block is None


# ---------------------------------------------------------------------------
# Module disabled via modules key
# ---------------------------------------------------------------------------

def test_returns_none_when_module_disabled(
    config_module_off, study_context, adata_no_summary
):
    result = run(adata_no_summary, config_module_off, study_context)
    assert result is None
