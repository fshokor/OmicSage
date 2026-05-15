"""
OmicSage — Tests for ai/cluster_annotator.py (Phase 3, Session 3 — B1)

All tests use a mock LLM — no real API key required.
Mock pattern: patch("ai.cluster_annotator.call_llm") returning str only.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_rank_genes_groups(cluster_ids: list[str], n_genes: int = 30):
    """Build a minimal adata.uns['rank_genes_groups'] structured array."""
    gene_names = [f"GENE{i:03d}" for i in range(n_genes)]
    # Each cluster gets a slightly different ordered list
    dtype = [(cid, "U20") for cid in cluster_ids]
    names_array = np.empty(n_genes, dtype=dtype)
    for cid in cluster_ids:
        names_array[cid] = gene_names
    return {
        "params": {"groupby": "leiden", "method": "wilcoxon"},
        "names": names_array,
    }


@pytest.fixture
def three_cluster_adata():
    """AnnData with 3 leiden clusters and rank_genes_groups populated."""
    n_cells = 90
    n_genes = 30
    X = sp.random(n_cells, n_genes, density=0.3, format="csr")
    obs = pd.DataFrame(
        {"leiden": (["0"] * 30) + (["1"] * 30) + (["2"] * 30)},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["rank_genes_groups"] = _make_rank_genes_groups(["0", "1", "2"], n_genes)
    return adata


@pytest.fixture
def single_cluster_adata():
    """AnnData with 1 leiden cluster."""
    n_cells = 30
    n_genes = 30
    X = sp.random(n_cells, n_genes, density=0.3, format="csr")
    obs = pd.DataFrame(
        {"leiden": ["0"] * n_cells},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)])
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["rank_genes_groups"] = _make_rank_genes_groups(["0"], n_genes)
    return adata


@pytest.fixture
def config_ai_enabled():
    return {
        "ai": {
            "features": True,
            "provider": "ollama",
            "model": "llama3",
            "modules": {"cluster_annotator": True},
        }
    }


@pytest.fixture
def config_ai_disabled():
    return {"ai": {"features": False}}


@pytest.fixture
def config_module_disabled():
    return {
        "ai": {
            "features": True,
            "provider": "ollama",
            "model": "llama3",
            "modules": {"cluster_annotator": False},
        }
    }


@pytest.fixture
def study_context():
    return {
        "dataset": {"tissue": "liver"},
        "disease": {"context": "hepatocellular carcinoma"},
    }


def _make_valid_response(cell_type: str = "CD8+ T cell") -> str:
    """Return a valid JSON string as the mock LLM would produce."""
    return json.dumps({
        "cell_type": cell_type,
        "confidence": "high",
        "supporting_markers": ["GENE000", "GENE001", "GENE002"],
        "alternative_types": [],
        "recommended_db": "CellTypist",
        "manual_marker_set": ["GENE000", "GENE001", "GENE002", "GENE003", "GENE004"],
        "reasoning": f"Cluster shows canonical {cell_type} markers.",
    })


# ---------------------------------------------------------------------------
# Gate tests — disabled AI always returns None
# ---------------------------------------------------------------------------

class TestAiGate:
    def test_returns_none_when_ai_features_false(
        self, three_cluster_adata, config_ai_disabled, study_context
    ):
        from ai.cluster_annotator import run
        result = run(three_cluster_adata, config_ai_disabled, study_context)
        assert result is None

    def test_returns_none_when_module_disabled(
        self, three_cluster_adata, config_module_disabled, study_context
    ):
        from ai.cluster_annotator import run
        result = run(three_cluster_adata, config_module_disabled, study_context)
        assert result is None

    def test_returns_none_when_runtime_ai_false(
        self, three_cluster_adata, config_ai_enabled, study_context
    ):
        from ai.cluster_annotator import run
        result = run(
            three_cluster_adata, config_ai_enabled, study_context, runtime_ai=False
        )
        assert result is None

    def test_no_obs_columns_written_when_disabled(
        self, three_cluster_adata, config_ai_disabled, study_context
    ):
        from ai.cluster_annotator import run
        run(three_cluster_adata, config_ai_disabled, study_context)
        for col in ("ai_cell_type", "ai_confidence", "ai_alt_types"):
            assert col not in three_cluster_adata.obs.columns


# ---------------------------------------------------------------------------
# Missing rank_genes_groups raises informative error
# ---------------------------------------------------------------------------

class TestMissingInput:
    def test_raises_on_missing_rank_genes_groups(
        self, three_cluster_adata, config_ai_enabled, study_context
    ):
        from ai.cluster_annotator import run
        del three_cluster_adata.uns["rank_genes_groups"]
        with pytest.raises(ValueError, match="rank_genes_groups"):
            with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
                run(three_cluster_adata, config_ai_enabled, study_context)

    def test_raises_on_missing_cluster_column(
        self, config_ai_enabled, study_context
    ):
        from ai.cluster_annotator import run
        n_cells = 30
        n_genes = 10
        X = sp.random(n_cells, n_genes, density=0.3, format="csr")
        obs = pd.DataFrame(
            {"sample": ["s1"] * n_cells},  # no leiden/louvain/cluster column
            index=[f"cell_{i}" for i in range(n_cells)],
        )
        var = pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["rank_genes_groups"] = _make_rank_genes_groups(["0"])
        with pytest.raises(ValueError, match="No clustering column"):
            with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
                run(adata, config_ai_enabled, study_context)


# ---------------------------------------------------------------------------
# Happy path — valid mock responses
# ---------------------------------------------------------------------------

class TestHappyPath:
    def test_returns_list_of_cluster_annotations(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run, ClusterAnnotation
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            result = run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert isinstance(result, list)
        assert len(result) == 3
        for ann in result:
            assert isinstance(ann, ClusterAnnotation)

    def test_obs_columns_written_correctly(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response("NK cell")):
            run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        adata = three_cluster_adata
        assert "ai_cell_type" in adata.obs.columns
        assert "ai_confidence" in adata.obs.columns
        assert "ai_alt_types" in adata.obs.columns
        assert (adata.obs["ai_cell_type"] == "NK cell").all()
        assert (adata.obs["ai_confidence"] == "high").all()

    def test_obs_columns_length_matches_n_cells(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert len(three_cluster_adata.obs["ai_cell_type"]) == 90

    def test_airesult_base_fields_populated(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            result = run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        ann = result[0]
        assert ann.timestamp != ""
        assert ann.model == "llama3"
        assert ann.provider == "ollama"
        assert ann.skill_name == "cluster_annotator"
        assert ann.skill_version == "1.0"

    def test_cluster_id_stored_on_annotation(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            result = run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert result[0].cluster_id == "0"

    def test_confidence_values_are_valid(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        valid = {"high", "medium", "low"}
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            result = run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        for ann in result:
            assert ann.confidence in valid
        assert all(
            v in valid
            for v in three_cluster_adata.obs["ai_confidence"].unique()
        )

    def test_alternative_types_empty_list_written_as_empty_string(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        # Empty alternative_types list → empty string in obs
        assert (single_cluster_adata.obs["ai_alt_types"] == "").all()

    def test_alternative_types_written_comma_separated(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        response = json.dumps({
            "cell_type": "T cell",
            "confidence": "medium",
            "supporting_markers": ["GENE000"],
            "alternative_types": ["NK cell", "ILC"],
            "recommended_db": "CellTypist",
            "manual_marker_set": ["GENE000", "GENE001", "GENE002", "GENE003", "GENE004"],
            "reasoning": "Ambiguous T/NK markers.",
        })
        with patch("ai.cluster_annotator.call_llm", return_value=response):
            run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert (single_cluster_adata.obs["ai_alt_types"] == "NK cell, ILC").all()


# ---------------------------------------------------------------------------
# Partial parse failure — other clusters still succeed
# ---------------------------------------------------------------------------

class TestPartialFailure:
    def test_one_cluster_parse_fail_others_succeed(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        responses = [
            _make_valid_response("CD8+ T cell"),   # cluster 0 — succeeds
            "NOT VALID JSON {{{{",                  # cluster 1 — fails
            _make_valid_response("Macrophage"),     # cluster 2 — succeeds
        ]
        call_count = {"n": 0}

        def mock_call_llm(**kwargs):
            resp = responses[call_count["n"]]
            call_count["n"] += 1
            return resp

        with patch("ai.cluster_annotator.call_llm", side_effect=mock_call_llm):
            result = run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )

        # Only 2 of 3 clusters annotated
        assert len(result) == 2
        cell_types = {ann.cell_type for ann in result}
        assert "CD8+ T cell" in cell_types
        assert "Macrophage" in cell_types

    def test_failed_cluster_cells_get_unknown_in_obs(
        self, three_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        responses = [
            _make_valid_response("CD8+ T cell"),
            "INVALID JSON",
            _make_valid_response("Macrophage"),
        ]
        call_count = {"n": 0}

        def mock_call_llm(**kwargs):
            resp = responses[call_count["n"]]
            call_count["n"] += 1
            return resp

        with patch("ai.cluster_annotator.call_llm", side_effect=mock_call_llm):
            run(
                three_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )

        # Cluster 1 cells (rows 30-59) should have "unknown" / "low"
        cluster1_mask = three_cluster_adata.obs["leiden"] == "1"
        assert (three_cluster_adata.obs.loc[cluster1_mask, "ai_cell_type"] == "unknown").all()
        assert (three_cluster_adata.obs.loc[cluster1_mask, "ai_confidence"] == "low").all()

    def test_malformed_json_logs_warning_not_error(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        """Parse failure must log a warning, not raise an exception."""
        import logging
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value="INVALID JSON"):
            with pytest.warns(UserWarning, match="no clusters were successfully annotated"):
                with patch.object(
                    logging.getLogger("ai.cluster_annotator"), "warning"
                ) as mock_warn:
                    result = run(
                        single_cluster_adata, config_ai_enabled, study_context,
                        log_dir=str(tmp_path),
                    )
        # Returns empty list (not None — AI was enabled, just nothing parsed)
        assert result == []


# ---------------------------------------------------------------------------
# Invalid confidence normalisation
# ---------------------------------------------------------------------------

class TestConfidenceNormalisation:
    def test_invalid_confidence_normalised_to_low(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        response = json.dumps({
            "cell_type": "B cell",
            "confidence": "VERY_HIGH",      # invalid value
            "supporting_markers": ["GENE000"],
            "alternative_types": [],
            "recommended_db": "PanglaoDB",
            "manual_marker_set": ["GENE000", "GENE001", "GENE002", "GENE003", "GENE004"],
            "reasoning": "Strong B cell markers.",
        })
        with patch("ai.cluster_annotator.call_llm", return_value=response):
            result = run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert result[0].confidence == "low"
        assert (single_cluster_adata.obs["ai_confidence"] == "low").all()

    def test_markdown_fenced_json_is_parsed(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        """Model sometimes wraps JSON in ```json fences despite instructions."""
        from ai.cluster_annotator import run
        inner = json.dumps({
            "cell_type": "Hepatocyte",
            "confidence": "high",
            "supporting_markers": ["GENE000"],
            "alternative_types": [],
            "recommended_db": "PanglaoDB",
            "manual_marker_set": ["GENE000", "GENE001", "GENE002", "GENE003", "GENE004"],
            "reasoning": "ALB and APOB expression.",
        })
        fenced = f"```json\n{inner}\n```"
        with patch("ai.cluster_annotator.call_llm", return_value=fenced):
            result = run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        assert len(result) == 1
        assert result[0].cell_type == "Hepatocyte"


# ---------------------------------------------------------------------------
# Audit log
# ---------------------------------------------------------------------------

class TestAuditLog:
    def test_audit_log_written_after_successful_run(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        log_file = tmp_path / "cluster_annotator.jsonl"
        assert log_file.exists()
        lines = log_file.read_text().strip().splitlines()
        assert len(lines) >= 1
        record = json.loads(lines[0])
        assert record["module"] == "cluster_annotator"
        assert record["parse_success"] is True

    def test_audit_log_written_on_parse_failure(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value="BAD JSON"):
            with pytest.warns(UserWarning):
                run(
                    single_cluster_adata, config_ai_enabled, study_context,
                    log_dir=str(tmp_path),
                )
        log_file = tmp_path / "cluster_annotator.jsonl"
        assert log_file.exists()
        record = json.loads(log_file.read_text().strip().splitlines()[0])
        assert record["parse_success"] is False

    def test_audit_log_contains_cluster_id_in_input_summary(
        self, single_cluster_adata, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            run(
                single_cluster_adata, config_ai_enabled, study_context,
                log_dir=str(tmp_path),
            )
        log_file = tmp_path / "cluster_annotator.jsonl"
        record = json.loads(log_file.read_text().strip().splitlines()[0])
        assert "cluster_id" in record["input_summary"]
        assert record["input_summary"]["cluster_id"] == "0"


# ---------------------------------------------------------------------------
# Louvain column fallback
# ---------------------------------------------------------------------------

class TestClusterColumnFallback:
    def test_louvain_column_used_when_leiden_absent(
        self, config_ai_enabled, study_context, tmp_path
    ):
        from ai.cluster_annotator import run
        n_cells = 30
        n_genes = 30
        X = sp.random(n_cells, n_genes, density=0.3, format="csr")
        obs = pd.DataFrame(
            {"louvain": ["0"] * n_cells},
            index=[f"cell_{i}" for i in range(n_cells)],
        )
        var = pd.DataFrame(index=[f"GENE{i:03d}" for i in range(n_genes)])
        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.uns["rank_genes_groups"] = _make_rank_genes_groups(["0"], n_genes)
        with patch("ai.cluster_annotator.call_llm", return_value=_make_valid_response()):
            result = run(adata, config_ai_enabled, study_context, log_dir=str(tmp_path))
        assert len(result) == 1
        assert "ai_cell_type" in adata.obs.columns
