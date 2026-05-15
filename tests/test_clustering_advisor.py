"""
tests/test_clustering_advisor.py
---------------------------------
Phase 3 — Session 2 — A2: Clustering Advisor tests.

All tests run without a real API key.
LLM call is mocked via: patch("ai.clustering_advisor.call_llm")
PubMed RAG is mocked via: patch("ai.clustering_advisor._fetch_pubmed_refs")
"""

from __future__ import annotations

import json
import os
import tempfile
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pytest

from ai.clustering_advisor import (
    ClusteringAdvice,
    PubMedRef,
    _best_from_sweep,
    _pubmed_query_string,
    _refs_to_prompt_str,
    _resolution_bias,
    _sweep_summary_str,
    run,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ENABLED_CONFIG = {
    "ai": {
        "features": True,
        "provider": "ollama",
        "model": "llama3",
        "modules": {
            "clustering_advisor": True,
        },
    }
}

DISABLED_GLOBAL_CONFIG = {
    "ai": {
        "features": False,
        "provider": "ollama",
        "model": "llama3",
    }
}

DISABLED_MODULE_CONFIG = {
    "ai": {
        "features": True,
        "provider": "ollama",
        "model": "llama3",
        "modules": {
            "clustering_advisor": False,
        },
    }
}

STUDY_CONTEXT = {
    "dataset": {"tissue": "liver"},
    "disease": {"context": "hepatocellular carcinoma"},
}

SWEEP = {0.3: 0.41, 0.5: 0.48, 0.8: 0.39, 1.0: 0.35}

VALID_LLM_RESPONSE = json.dumps({
    "suggested_resolution": 0.5,
    "resolution_range": [0.3, 0.8],
    "expected_n_clusters": 12,
    "literature_context": [
        {
            "pmid": "12345678",
            "title": "Single-cell RNA-seq reveals immune landscape of HCC",
            "resolution_used": 0.5,
            "tissue": "liver",
            "disease": "hepatocellular carcinoma",
        }
    ],
    "reasoning": "Resolution 0.5 yielded the highest silhouette score of 0.48 "
                 "and aligns with published HCC studies.",
})


def _make_adata(n_obs: int = 500, n_vars: int = 2000) -> ad.AnnData:
    rng = np.random.default_rng(42)
    X = rng.poisson(1, size=(n_obs, n_vars)).astype(np.float32)
    adata = ad.AnnData(X)
    # Mark first 500 vars as highly variable
    adata.var["highly_variable"] = False
    adata.var.iloc[:min(500, n_vars), adata.var.columns.get_loc("highly_variable")] = True
    return adata


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGating:
    """Tests for config gate — all must return None without touching LLM."""

    def test_returns_none_when_ai_features_false(self):
        adata = _make_adata()
        result = run(adata, DISABLED_GLOBAL_CONFIG, STUDY_CONTEXT, SWEEP, runtime_ai=True)
        assert result is None

    def test_returns_none_when_module_disabled(self):
        adata = _make_adata()
        result = run(adata, DISABLED_MODULE_CONFIG, STUDY_CONTEXT, SWEEP, runtime_ai=True)
        assert result is None

    def test_returns_none_when_runtime_ai_false(self):
        adata = _make_adata()
        result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, runtime_ai=False)
        assert result is None


class TestHappyPath:
    """Tests for successful ClusteringAdvice generation with mocked LLM."""

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_returns_clustering_advice(self, mock_llm, mock_pubmed):
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        assert isinstance(result, ClusteringAdvice)

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_suggested_resolution_within_sweep_range(self, mock_llm, mock_pubmed):
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        lo, hi = min(SWEEP), max(SWEEP)
        assert lo <= result.suggested_resolution <= hi

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_resolution_range_is_ordered(self, mock_llm, mock_pubmed):
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        assert result.resolution_range[0] < result.resolution_range[1]

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_base_fields_populated(self, mock_llm, mock_pubmed):
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        assert result.timestamp is not None
        assert result.model == "llama3"
        assert result.provider == "ollama"

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_audit_log_written(self, mock_llm, mock_pubmed):
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
            log_file = os.path.join(tmp, "clustering_advisor.jsonl")
            assert os.path.exists(log_file)
            with open(log_file) as f:
                record = json.loads(f.readline())
            assert record["module"] == "clustering_advisor"
            assert record["parse_success"] is True


class TestEdgeCases:
    """Tests for empty sweep, empty PubMed, malformed JSON."""

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_empty_sweep_handled_gracefully(self, mock_llm, mock_pubmed):
        """Empty sweep must not crash — default resolution used."""
        mock_llm.return_value = VALID_LLM_RESPONSE
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            # Should warn but not raise
            with pytest.warns(UserWarning, match="resolution_sweep_results is empty"):
                result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, {}, log_dir=tmp)
        assert result is not None  # LLM mock still returns valid response

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_empty_pubmed_returns_empty_literature_context(self, mock_llm, mock_pubmed):
        response = json.dumps({
            "suggested_resolution": 0.5,
            "resolution_range": [0.3, 0.8],
            "expected_n_clusters": 10,
            "literature_context": [],
            "reasoning": "No literature found.",
        })
        mock_llm.return_value = response
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        assert result.literature_context == []

    @patch("ai.clustering_advisor._fetch_pubmed_refs", return_value=[])
    @patch("ai.clustering_advisor.call_llm")
    def test_malformed_llm_json_returns_none(self, mock_llm, mock_pubmed):
        mock_llm.return_value = "this is not json at all }"
        adata = _make_adata()
        with tempfile.TemporaryDirectory() as tmp:
            result = run(adata, ENABLED_CONFIG, STUDY_CONTEXT, SWEEP, log_dir=tmp)
        assert result is None


class TestPubMedQuery:
    """Tests for PubMed query construction."""

    def test_query_string_includes_tissue_and_disease(self):
        q = _pubmed_query_string("liver", "hepatocellular carcinoma")
        assert "liver" in q
        assert "hepatocellular carcinoma" in q
        assert "Leiden" in q
        assert "single-cell RNA-seq" in q

    def test_query_string_without_disease(self):
        q = _pubmed_query_string("PBMC", None)
        assert "PBMC" in q
        assert "Leiden" in q


class TestHelpers:
    """Unit tests for rule-based helpers."""

    def test_best_from_sweep_returns_max_silhouette(self):
        best_res, best_sil = _best_from_sweep(SWEEP)
        assert best_res == 0.5
        assert abs(best_sil - 0.48) < 1e-6

    def test_best_from_sweep_empty_returns_default(self):
        best_res, best_sil = _best_from_sweep({})
        assert best_res == 0.5
        assert best_sil == 0.0

    def test_resolution_bias_low_cell_count(self):
        biased = _resolution_bias(500, 0.5)
        assert biased < 0.5

    def test_resolution_bias_high_cell_count(self):
        biased = _resolution_bias(100_000, 0.5)
        assert biased > 0.5

    def test_resolution_bias_normal_count_unchanged(self):
        biased = _resolution_bias(5_000, 0.5)
        assert biased == 0.5

    def test_sweep_summary_str_sorted(self):
        summary = _sweep_summary_str(SWEEP)
        assert "0.3" in summary
        assert "0.5" in summary

    def test_sweep_summary_str_empty(self):
        assert "No sweep" in _sweep_summary_str({})

    def test_refs_to_prompt_str_empty(self):
        assert _refs_to_prompt_str([]) == ""

    def test_refs_to_prompt_str_contains_pmid(self):
        refs = [{"pmid": "99999999", "title": "Test paper", "resolution_used": 0.5,
                 "tissue": "liver", "disease": "HCC"}]
        result = _refs_to_prompt_str(refs)
        assert "99999999" in result
        assert "Test paper" in result
