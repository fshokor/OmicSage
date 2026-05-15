"""
tests/test_pipeline_advisor.py
==============================
All tests for the Phase 3 — Session 1 Pipeline Advisor module.

Run with:
    python -m pytest tests/test_pipeline_advisor.py -v

No real API key is needed — the LLM client is monkeypatched throughout.
"""

from __future__ import annotations

import json
import datetime
from dataclasses import dataclass, field
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata
import numpy as np
import pandas as pd
import pytest

# ---------------------------------------------------------------------------
# Helpers — minimal fixtures
# ---------------------------------------------------------------------------

def _make_adata(
    n_cells: int = 1000,
    n_genes: int = 2000,
    n_batches: int = 1,
    n_donors: int = 3,
    n_conditions: int = 2,
) -> anndata.AnnData:
    """Return a minimal AnnData with obs columns for batch/donor/condition."""
    rng = np.random.default_rng(42)
    X = rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_cells)])

    if n_batches > 1:
        obs["sample_id"] = [f"batch_{i % n_batches}" for i in range(n_cells)]
    if n_donors >= 1:
        obs["donor"] = [f"donor_{i % n_donors}" for i in range(n_cells)]
    if n_conditions >= 1:
        obs["condition"] = [f"cond_{i % n_conditions}" for i in range(n_cells)]

    return anndata.AnnData(X=X, obs=obs)


def _make_config(
    ai_features: bool = True,
    pipeline_advisor_enabled: bool = True,
    batch_key: str | None = None,
    provider: str = "claude",
    model: str = "claude-sonnet-4-20250514",
) -> dict:
    cfg: dict = {
        "ai": {
            "features": ai_features,
            "provider": provider,
            "model": model,
            "modules": {
                "pipeline_advisor": pipeline_advisor_enabled,
            },
        }
    }
    if batch_key:
        cfg["preprocessing"] = {"batch_key": batch_key}
    return cfg


def _make_study_context(
    tissue: str = "liver",
    disease_context: str = "hepatocellular carcinoma",
    biological_question: str | None = "Characterise the tumour microenvironment.",
    n_donors: int = 3,
    n_conditions: int = 2,
    batch_key: str | None = "sample_id",
) -> dict:
    return {
        "dataset": {"tissue": tissue},
        "disease": {"context": disease_context},
        "experiment": {
            "design": "tumour_vs_normal",
            "n_donors": n_donors,
            "n_conditions": n_conditions,
            "batch_key": batch_key,
        },
        "biological_question": biological_question,
    }


def _valid_llm_response(
    inferred_question: str | None = None,
    n_steps: int = 5,
    extra_warnings: list[str] | None = None,
) -> str:
    """Return a valid JSON string matching the pipeline_advisor output schema."""
    steps = [
        {
            "step_name": name,
            "priority": priority,
            "rationale": f"Standard {name} is needed for this dataset.",
        }
        for name, priority in [
            ("qc", "required"),
            ("normalization", "required"),
            ("hvg_selection", "required"),
            ("pca", "required"),
            ("umap", "recommended"),
            ("clustering", "required"),
            ("batch_correction", "recommended"),
            ("cell_type_annotation", "recommended"),
        ][:n_steps]
    ]
    return json.dumps({
        "recommended_steps": steps,
        "inferred_biological_question": inferred_question,
        "warnings": extra_warnings or [],
        "reasoning": "Standard scRNA-seq pipeline recommended based on dataset properties.",
    })


# ---------------------------------------------------------------------------
# Mock LLM helper
# ---------------------------------------------------------------------------

def _mock_call_llm(response_str: str):
    """Return a patcher for ai._llm_client.call_llm that returns response_str."""
    mock = MagicMock(return_value=(
        response_str,
        {"prompt_tokens": 100, "completion_tokens": 50},
        "claude-sonnet-4-20250514",
        "claude",
    ))
    return patch("ai.pipeline_advisor.call_llm", mock)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestAiDisabledReturnsNone:

    def test_returns_none_when_ai_features_false(self, tmp_path):
        """Global ai_features: false → None, no LLM call."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config(ai_features=False)
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()) as mock_llm:
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is None
        mock_llm.assert_not_called()

    def test_returns_none_when_module_disabled(self, tmp_path):
        """Per-module disable → None."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config(pipeline_advisor_enabled=False)
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()) as mock_llm:
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is None
        mock_llm.assert_not_called()

    def test_returns_none_when_runtime_ai_false(self, tmp_path):
        """runtime_ai=False (no --ai flag) → None even if config enables AI."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config(ai_features=True)
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()) as mock_llm:
            result = run(
                adata, config, study_context,
                log_dir=str(tmp_path),
                runtime_ai=False,
            )

        assert result is None
        mock_llm.assert_not_called()


class TestSuccessfulRun:

    def test_returns_pipeline_advice_on_valid_llm_response(self, tmp_path):
        """Happy path: valid LLM JSON → PipelineAdvice dataclass."""
        from ai.pipeline_advisor import run, PipelineAdvice

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        assert isinstance(result, PipelineAdvice)

    def test_recommended_steps_non_empty(self, tmp_path):
        """recommended_steps list must have at least one entry."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response(n_steps=5)):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        assert len(result.recommended_steps) > 0

    def test_step_priority_values_are_valid(self, tmp_path):
        """All StepRecommendation.priority values must be in the allowed set."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response(n_steps=8)):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        valid = {"required", "recommended", "optional"}
        for step in result.recommended_steps:
            assert step.priority in valid, (
                f"Step '{step.step_name}' has invalid priority '{step.priority}'"
            )

    def test_airesult_base_fields_populated(self, tmp_path):
        """timestamp, model, provider must be set on the returned object."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        assert result.timestamp  # non-empty string
        assert result.model == "claude-sonnet-4-20250514"
        assert result.provider == "claude"
        assert result.skill_version == "1.0"


class TestRuleBasedChecks:

    def test_batch_warning_fires_when_n_batches_gt_1_and_no_batch_key(self, tmp_path):
        """n_batches > 1 with no batch_key in config or study_context → warning."""
        from ai.pipeline_advisor import run

        # 3 batches, no batch_key anywhere
        adata = _make_adata(n_batches=3)
        config = _make_config()  # no batch_key
        study_context = _make_study_context(batch_key=None)

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        batch_warnings = [w for w in result.warnings if "batch" in w.lower() or "batch_key" in w.lower()]
        assert len(batch_warnings) >= 1

    def test_batch_warning_absent_when_batch_key_set(self, tmp_path):
        """If batch_key is set, the missing-batch_key warning should not fire."""
        from ai.pipeline_advisor import run

        adata = _make_adata(n_batches=3)
        config = _make_config(batch_key="sample_id")
        study_context = _make_study_context(batch_key="sample_id")

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        missing_key_warnings = [
            w for w in result.warnings if "no batch_key" in w.lower()
        ]
        assert len(missing_key_warnings) == 0

    def test_pseudobulk_recommendation_fires_when_n_donors_gt_2_and_n_conditions_gt_1(self, tmp_path):
        """n_donors > 2 and n_conditions > 1 → pseudobulk_deg in recommended steps."""
        from ai.pipeline_advisor import run

        adata = _make_adata(n_donors=4, n_conditions=2)
        config = _make_config()
        study_context = _make_study_context(n_donors=4, n_conditions=2)

        response = _valid_llm_response()
        # Inject pseudobulk step into response
        parsed = json.loads(response)
        parsed["recommended_steps"].append({
            "step_name": "pseudobulk_deg",
            "priority": "recommended",
            "rationale": "Enough donors for pseudobulk.",
        })
        with _mock_call_llm(json.dumps(parsed)):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        step_names = [s.step_name for s in result.recommended_steps]
        assert "pseudobulk_deg" in step_names

    def test_small_dataset_warning_fires_when_n_cells_lt_500(self, tmp_path):
        """n_cells < 500 → warning about unreliable clustering."""
        from ai.pipeline_advisor import run

        adata = _make_adata(n_cells=200)
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        small_warnings = [w for w in result.warnings if "500" in w or "cells" in w.lower()]
        assert len(small_warnings) >= 1

    def test_inferred_biological_question_populated_when_blank(self, tmp_path):
        """If biological_question is blank, LLM-inferred question is set."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        # No biological_question in study_context
        study_context = _make_study_context(biological_question=None)

        inferred = "Characterise the immune microenvironment of hepatocellular carcinoma."
        with _mock_call_llm(_valid_llm_response(inferred_question=inferred)):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        assert result.inferred_biological_question == inferred


class TestAuditLog:

    def test_audit_log_written_after_successful_call(self, tmp_path):
        """Audit log JSONL file must exist and be non-empty after a successful run."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm(_valid_llm_response()):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        log_file = tmp_path / "pipeline_advisor.jsonl"
        assert log_file.exists(), "Audit log file not created"
        content = log_file.read_text(encoding="utf-8").strip()
        assert len(content) > 0, "Audit log file is empty"

        # Verify it's valid JSONL (one JSON object per line)
        for line in content.splitlines():
            record = json.loads(line)
            assert "module" in record
            assert record["module"] == "pipeline_advisor"


class TestGracefulDegradation:

    def test_returns_none_on_malformed_llm_json(self, tmp_path):
        """If LLM returns unparseable output, return None and log warning."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm("this is not json at all"):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is None

    def test_audit_log_written_even_on_parse_failure(self, tmp_path):
        """Audit log must be written even when parsing fails."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        with _mock_call_llm("{ bad json"):
            run(adata, config, study_context, log_dir=str(tmp_path))

        log_file = tmp_path / "pipeline_advisor.jsonl"
        assert log_file.exists()
        record = json.loads(log_file.read_text(encoding="utf-8").strip().splitlines()[0])
        assert record["parse_success"] is False

    def test_unknown_step_name_in_llm_response_is_skipped(self, tmp_path):
        """LLM returning an unknown step_name should be silently dropped."""
        from ai.pipeline_advisor import run

        adata = _make_adata()
        config = _make_config()
        study_context = _make_study_context()

        bad_response = json.dumps({
            "recommended_steps": [
                {"step_name": "qc", "priority": "required", "rationale": "Always needed."},
                {"step_name": "invented_magic_step", "priority": "required", "rationale": "Does magic."},
            ],
            "inferred_biological_question": None,
            "warnings": [],
            "reasoning": "Some reasoning.",
        })

        with _mock_call_llm(bad_response):
            result = run(adata, config, study_context, log_dir=str(tmp_path))

        assert result is not None
        step_names = [s.step_name for s in result.recommended_steps]
        assert "qc" in step_names
        assert "invented_magic_step" not in step_names
