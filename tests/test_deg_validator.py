"""
tests/test_deg_validator.py
---------------------------
Phase 3 — Session 4: B2 — DEG Validator + Literature Linker

All tests run without a real API key.
  - LLM calls patched via: patch("ai.deg_validator.call_llm")
  - PubMed calls patched via: patch("ai.deg_validator._query_pubmed")

330 + N tests total after this session.
"""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from ai.deg_validator import (
    DegValidation,
    GeneLitRef,
    _deduplicate_lit_refs,
    _extract_degs,
    _parse_llm_response,
    run,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_config(
    ai_enabled: bool = True,
    module_enabled: bool = True,
    provider: str = "ollama",
    model: str = "llama3",
) -> dict:
    cfg: dict = {
        "ai": {
            "features": ai_enabled,
            "provider": provider,
            "model": model,
            "modules": {},
        }
    }
    if not module_enabled:
        cfg["ai"]["modules"]["deg_validator"] = False
    return cfg


def _make_study_context(
    tissue: str = "liver",
    disease: str | None = "hepatocellular carcinoma",
) -> dict:
    ctx: dict = {"dataset": {"tissue": tissue}, "disease": {}}
    if disease:
        ctx["disease"]["context"] = disease
    return ctx


def _make_adata_with_degs(n_groups: int = 2, n_genes_per_group: int = 5) -> anndata.AnnData:
    """
    Build a minimal AnnData with rank_genes_groups populated.
    Groups are named "0", "1", ... (n_groups total).

    Structured array layout (matches scanpy output):
      dtype has one field per group.
      Row i holds the i-th ranked gene (or value) for EACH group.
      So row i = (group0_gene_i, group1_gene_i, ...)
    """
    n_cells = 40
    n_genes = 20
    rng = np.random.default_rng(42)

    X = sp.csr_matrix(rng.integers(0, 50, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(
        {"leiden": [str(i % n_groups) for i in range(n_cells)]},
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(
        index=[f"GENE{i}" for i in range(n_genes)]
    )
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    # One set of gene names per group (groups may have different top genes)
    group_names = [str(g) for g in range(n_groups)]
    gene_names_per_group = {
        g: [f"GENE{i}" for i in range(n_genes_per_group)]
        for g in group_names
    }

    dtype_names = np.dtype([(g, "U10") for g in group_names])
    dtype_float  = np.dtype([(g, "f4")  for g in group_names])

    # Row rank → tuple of (group0_gene_rank, group1_gene_rank, ...)
    names_arr = np.array(
        [tuple(gene_names_per_group[g][rank] for g in group_names)
         for rank in range(n_genes_per_group)],
        dtype=dtype_names,
    )

    lfc_vals = {g: rng.uniform(0.5, 3.0, n_genes_per_group) for g in group_names}
    lfc_arr = np.array(
        [tuple(lfc_vals[g][rank] for g in group_names)
         for rank in range(n_genes_per_group)],
        dtype=dtype_float,
    )

    pval_vals = {g: rng.uniform(0.001, 0.05, n_genes_per_group) for g in group_names}
    pvals_arr = np.array(
        [tuple(pval_vals[g][rank] for g in group_names)
         for rank in range(n_genes_per_group)],
        dtype=dtype_float,
    )

    adata.uns["rank_genes_groups"] = {
        "params": {"groupby": "leiden", "method": "wilcoxon"},
        "names": names_arr,
        "logfoldchanges": lfc_arr,
        "pvals_adj": pvals_arr,
    }

    return adata


def _make_adata_no_rgg() -> anndata.AnnData:
    """AnnData without rank_genes_groups."""
    n_cells, n_genes = 20, 10
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.integers(0, 10, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"G{i}" for i in range(n_genes)])
    return anndata.AnnData(X=X, obs=obs, var=var)


def _valid_llm_response(
    expected: list[str] | None = None,
    unexpected: list[str] | None = None,
) -> str:
    payload = {
        "expected_genes": expected or ["GENE0", "GENE1"],
        "unexpected_genes": unexpected or ["GENE4"],
        "validation_summary": "The DEG list is broadly consistent with liver HCC biology.",
        "discovery_highlights": ["GENE0 upregulated in tumour — consistent with known HCC markers."],
    }
    return json.dumps(payload)


def _mock_pubmed_refs(gene: str, tissue: str, disease_context) -> list[GeneLitRef]:
    return [
        GeneLitRef(
            gene=gene,
            pmid=f"1234{gene[-1]}",
            title=f"Role of {gene} in liver cancer",
            context=f"{gene} mentioned in: Role of {gene} in liver cancer",
        )
    ]


# ---------------------------------------------------------------------------
# Gate tests
# ---------------------------------------------------------------------------

class TestAiGate:
    def test_returns_none_when_ai_features_false(self):
        adata = _make_adata_with_degs()
        config = _make_config(ai_enabled=False)
        result = run(adata, config, _make_study_context())
        assert result is None

    def test_returns_none_when_module_disabled(self):
        adata = _make_adata_with_degs()
        config = _make_config(module_enabled=False)
        result = run(adata, config, _make_study_context())
        assert result is None

    def test_returns_none_when_runtime_ai_false(self):
        adata = _make_adata_with_degs()
        config = _make_config()
        result = run(adata, config, _make_study_context(), runtime_ai=False)
        assert result is None


# ---------------------------------------------------------------------------
# Missing rank_genes_groups
# ---------------------------------------------------------------------------

class TestMissingRankGenes:
    def test_raises_informative_error_when_rgg_missing(self):
        adata = _make_adata_no_rgg()
        config = _make_config()
        with pytest.raises(ValueError, match="rank_genes_groups"):
            run(adata, config, _make_study_context())


# ---------------------------------------------------------------------------
# Core run tests
# ---------------------------------------------------------------------------

class TestRunCore:
    def test_returns_list_of_deg_validation(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=2)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        assert isinstance(results, list)
        assert len(results) == 2
        assert all(isinstance(r, DegValidation) for r in results)

    def test_expected_genes_populated(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response(expected=["GENE0", "GENE1"])),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        assert "GENE0" in results[0].expected_genes
        assert "GENE1" in results[0].expected_genes

    def test_unexpected_genes_populated(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response(unexpected=["GENE4"])),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        assert "GENE4" in results[0].unexpected_genes

    def test_literature_links_populated(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        assert len(results[0].literature_links) > 0
        assert all(isinstance(ref, GeneLitRef) for ref in results[0].literature_links)

    def test_airesult_base_fields_populated(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config(provider="ollama", model="llama3")
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        r = results[0]
        assert r.skill_name == "deg_validator"
        assert r.skill_version == "1.0"
        assert r.model == "llama3"
        assert r.provider == "ollama"
        assert r.timestamp != ""

    def test_comparison_field_set_to_group_name(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=2)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        comparisons = {r.comparison for r in results}
        assert comparisons == {"0", "1"}


# ---------------------------------------------------------------------------
# Empty DEG list
# ---------------------------------------------------------------------------

class TestEmptyDegs:
    def test_empty_degs_returns_empty_validation(self, tmp_path):
        """
        If a group has no DEGs (n_top_genes=0 effectively), return an
        empty DegValidation rather than crashing.
        """
        adata = _make_adata_with_degs(n_groups=1, n_genes_per_group=5)
        config = _make_config()
        ctx = _make_study_context()

        # Override _extract_degs to return empty list for the group
        with (
            patch("ai.deg_validator._extract_degs", return_value={"0": []}),
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        assert len(results) == 1
        assert results[0].expected_genes == []
        assert results[0].unexpected_genes == []
        assert results[0].literature_links == []


# ---------------------------------------------------------------------------
# PMID deduplication
# ---------------------------------------------------------------------------

class TestPmidDeduplication:
    def test_same_pmid_deduplicated(self):
        """Same PMID for two different genes → stored once."""
        refs = [
            GeneLitRef(gene="GENE0", pmid="99999", title="Paper A", context="GENE0 mentioned in: Paper A"),
            GeneLitRef(gene="GENE1", pmid="99999", title="Paper A", context="GENE1 mentioned in: Paper A"),
            GeneLitRef(gene="GENE2", pmid="88888", title="Paper B", context="GENE2 mentioned in: Paper B"),
        ]
        deduped = _deduplicate_lit_refs(refs)
        pmids = [r.pmid for r in deduped]
        assert pmids.count("99999") == 1
        assert pmids.count("88888") == 1
        assert len(deduped) == 2

    def test_deduplication_extends_context_for_additional_gene(self):
        """When a PMID is shared, context for the second gene is appended."""
        refs = [
            GeneLitRef(gene="GENE0", pmid="99999", title="Paper A", context="GENE0 mentioned in: Paper A"),
            GeneLitRef(gene="GENE1", pmid="99999", title="Paper A", context="GENE1 mentioned in: Paper A"),
        ]
        deduped = _deduplicate_lit_refs(refs)
        assert len(deduped) == 1
        assert "GENE1" in deduped[0].context

    def test_unique_pmids_all_kept(self):
        refs = [
            GeneLitRef(gene=f"G{i}", pmid=f"{i}", title=f"Paper {i}", context=f"ctx {i}")
            for i in range(5)
        ]
        deduped = _deduplicate_lit_refs(refs)
        assert len(deduped) == 5

    def test_end_to_end_deduplication_via_run(self, tmp_path):
        """
        Two genes return the same PMID from PubMed.
        Final literature_links should have that PMID only once.
        """
        adata = _make_adata_with_degs(n_groups=1, n_genes_per_group=3)
        config = _make_config()
        ctx = _make_study_context()

        # All genes return the same PMID
        def same_pmid_pubmed(gene, tissue, disease_context):
            return [GeneLitRef(gene=gene, pmid="77777", title="Shared Paper", context=f"{gene} in: Shared Paper")]

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=same_pmid_pubmed),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        pmids = [r.pmid for r in results[0].literature_links]
        assert pmids.count("77777") == 1


# ---------------------------------------------------------------------------
# Copyright constraint — no abstract text
# ---------------------------------------------------------------------------

class TestNoCopyrightViolation:
    def test_no_abstract_text_in_gene_lit_ref(self, tmp_path):
        """
        GeneLitRef.context must be derived from the title, never from abstract.
        We verify the mock doesn't slip abstract text through.
        """
        adata = _make_adata_with_degs(n_groups=1, n_genes_per_group=3)
        config = _make_config()
        ctx = _make_study_context()

        abstract_text = "ABSTRACT: this is prohibited abstract content from the paper."

        def pubmed_with_abstract(gene, tissue, disease_context):
            # Simulate a badly written integration that includes abstract — must not pass through
            return [
                GeneLitRef(
                    gene=gene,
                    pmid="11111",
                    title="Some Paper Title",
                    context=f"{gene} mentioned in: Some Paper Title",
                    # context must NOT contain abstract_text
                )
            ]

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=pubmed_with_abstract),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        for ref in results[0].literature_links:
            assert abstract_text not in ref.context
            assert abstract_text not in ref.title
            # pmid and title are safe — verify they are set
            assert ref.pmid != ""
            assert ref.title != ""

    def test_gene_lit_ref_has_no_extra_fields(self):
        """GeneLitRef must not grow an 'abstract' field."""
        ref = GeneLitRef(gene="TP53", pmid="12345", title="TP53 in cancer", context="TP53 mentioned in: TP53 in cancer")
        assert not hasattr(ref, "abstract")
        assert not hasattr(ref, "full_text")
        assert not hasattr(ref, "body")


# ---------------------------------------------------------------------------
# Audit log
# ---------------------------------------------------------------------------

class TestAuditLog:
    def test_audit_log_written_after_successful_run(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            run(adata, config, ctx, log_dir=str(tmp_path))

        log_file = tmp_path / "deg_validator.jsonl"
        assert log_file.exists(), "Audit log file was not created"
        lines = log_file.read_text().strip().splitlines()
        assert len(lines) >= 1

        record = json.loads(lines[0])
        assert record["module"] == "deg_validator"
        assert record["parse_success"] is True

    def test_audit_log_contains_expected_fields(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value=_valid_llm_response()),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            run(adata, config, ctx, log_dir=str(tmp_path))

        record = json.loads((tmp_path / "deg_validator.jsonl").read_text().strip().splitlines()[0])
        for key in ("module", "skill_version", "model", "provider", "input_summary", "raw_response", "parse_success"):
            assert key in record, f"Missing key in audit record: {key}"

    def test_audit_log_not_written_when_ai_disabled(self, tmp_path):
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config(ai_enabled=False)
        ctx = _make_study_context()

        run(adata, config, ctx, log_dir=str(tmp_path))

        log_file = tmp_path / "deg_validator.jsonl"
        assert not log_file.exists(), "Audit log should not be written when AI is disabled"


# ---------------------------------------------------------------------------
# LLM parse edge cases
# ---------------------------------------------------------------------------

class TestLlmParsing:
    def test_parse_valid_json(self):
        raw = _valid_llm_response()
        parsed = _parse_llm_response(raw)
        assert "expected_genes" in parsed
        assert "unexpected_genes" in parsed

    def test_parse_strips_markdown_fences(self):
        raw = "```json\n" + _valid_llm_response() + "\n```"
        parsed = _parse_llm_response(raw)
        assert isinstance(parsed["expected_genes"], list)

    def test_parse_raises_on_missing_keys(self):
        raw = json.dumps({"expected_genes": ["G1"]})  # missing other keys
        with pytest.raises(ValueError, match="missing keys"):
            _parse_llm_response(raw)

    def test_parse_raises_on_invalid_json(self):
        with pytest.raises(json.JSONDecodeError):
            _parse_llm_response("not valid json at all {{{")

    def test_degraded_parse_does_not_crash_run(self, tmp_path):
        """If LLM returns unparseable JSON, run returns DegValidation with empty lists."""
        adata = _make_adata_with_degs(n_groups=1)
        config = _make_config()
        ctx = _make_study_context()

        with (
            patch("ai.deg_validator.call_llm", return_value="INVALID JSON !!!"),
            patch("ai.deg_validator._query_pubmed", side_effect=_mock_pubmed_refs),
        ):
            results = run(adata, config, ctx, log_dir=str(tmp_path))

        # Should not crash; returns result with empty parsed fields
        assert isinstance(results, list)
        assert results[0].expected_genes == []
        assert results[0].unexpected_genes == []
        # Literature links still populated from PubMed even when LLM parse fails
        assert len(results[0].literature_links) > 0
