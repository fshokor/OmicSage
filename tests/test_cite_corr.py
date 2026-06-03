"""
tests/test_cite_corr.py — Tests for cite_corr (Phase 4, Step 9)

Coverage
--------
Input validation:
  - raises ValueError on AnnData input
  - raises ValueError on unsupported method

Matching:
  - exact match finds correct gene index
  - suffix-stripped match works ("CD3E-TotalSeqA" → "CD3E")
  - unmatched proteins recorded in corr_dict["unmatched"]
  - proteins below min_cells threshold recorded in corr_dict["skipped"]

Results DataFrame:
  - expected columns present (protein, gene, r, pval, pval_adj, n_cells, matched_by)
  - r values in [-1, 1]
  - sorted by r descending
  - matched_by values are "exact" or "stripped" only

Return contract:
  - returns (MuData, dict)
  - dict has keys: results, matched, unmatched, skipped, provenance
  - provenance written to mdata.uns["omicsage_cite_corr"]
  - inplace=False does not modify original
  - inplace=True modifies original
  - empty DataFrame returned when no pairs matched

Warnings:
  - warns when adt_clr layer absent
  - warns when logcounts layer absent
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

mudata = pytest.importorskip("mudata")
MuData = mudata.MuData


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Gene panel — includes names that match ADT proteins exactly and by stripping
_RNA_GENES = [
    "CD3E", "CD4", "CD8A", "CD19", "MS4A1",
    "CD14", "FCGR3A", "NKG7", "GNLY", "CD56",
    "IL7R", "CCR7", "GZMB", "PRF1", "FOXP3",
    "Gene15", "Gene16", "Gene17", "Gene18", "Gene19",
]

_ADT_PROTEINS = [
    "CD3E",              # exact match → CD3E
    "CD4",               # exact match → CD4
    "CD19-TotalSeqA",    # stripped → CD19
    "CD14",              # exact match → CD14
    "Mouse-IgG1",        # no match → unmatched
    "Rat-IgG2b",         # no match → unmatched
    "CD56",              # exact match → CD56
    "NoMatch-XYZ",       # no match → unmatched
]


def _make_rna(n_obs: int = 80, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(
        rng.poisson(2.0, size=(n_obs, len(_RNA_GENES))).astype(np.float32)
    )
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = _RNA_GENES
    adata.layers["logcounts"] = X.toarray().copy()
    return adata


def _make_adt(n_obs: int = 80, seed: int = 1) -> AnnData:
    rng = np.random.default_rng(seed)
    clr = rng.normal(0, 1, size=(n_obs, len(_ADT_PROTEINS))).astype(np.float32)
    # Make some proteins positively expressed (> 0)
    clr = np.abs(clr)
    adt = AnnData(X=clr.copy())
    adt.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adt.var_names = _ADT_PROTEINS
    adt.layers["adt_clr"] = clr.copy()
    return adt


@pytest.fixture
def mdata_fixture():
    rna = _make_rna()
    adt = _make_adt()
    return MuData({"rna": rna, "adt": adt})


@pytest.fixture
def adt_only_fixture():
    return _make_adt()


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------

from pipeline.modules.cite.cite_corr import cite_corr


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_raises_on_anndata_input(self, adt_only_fixture):
        with pytest.raises(ValueError, match="MuData"):
            cite_corr(adt_only_fixture)

    def test_raises_on_unsupported_method(self, mdata_fixture):
        with pytest.raises(ValueError, match="spearman"):
            cite_corr(mdata_fixture, method="pearson")


# ---------------------------------------------------------------------------
# Return contract
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self, mdata_fixture):
        result, corr_dict = cite_corr(mdata_fixture, inplace=False)
        assert isinstance(result, MuData)
        assert isinstance(corr_dict, dict)

    def test_dict_has_expected_keys(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, inplace=False)
        for key in ["results", "matched", "unmatched", "skipped", "provenance"]:
            assert key in corr_dict, f"Missing key: {key}"

    def test_provenance_in_uns(self, mdata_fixture):
        result, _ = cite_corr(mdata_fixture, inplace=False)
        assert "omicsage_cite_corr" in result.uns
        prov = result.uns["omicsage_cite_corr"]
        assert prov["module"] == "cite_corr"
        assert "timestamp" in prov
        assert "n_matched" in prov

    def test_inplace_false_does_not_modify_original(self, mdata_fixture):
        original_uns = set(mdata_fixture.uns.keys())
        cite_corr(mdata_fixture, inplace=False)
        assert set(mdata_fixture.uns.keys()) == original_uns

    def test_inplace_true_modifies_original(self, mdata_fixture):
        cite_corr(mdata_fixture, inplace=True)
        assert "omicsage_cite_corr" in mdata_fixture.uns


# ---------------------------------------------------------------------------
# Matching tests
# ---------------------------------------------------------------------------

class TestMatching:

    def test_exact_match_found(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, inplace=False)
        matched_proteins = [p for p, g in corr_dict["matched"]]
        assert "CD3E" in matched_proteins
        assert "CD4" in matched_proteins
        assert "CD14" in matched_proteins

    def test_suffix_stripped_match(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, inplace=False)
        # CD19-TotalSeqA should match CD19 after stripping
        matched_proteins = [p for p, g in corr_dict["matched"]]
        assert "CD19-TotalSeqA" in matched_proteins
        # Verify it matched to CD19 gene
        matched_genes = {p: g for p, g in corr_dict["matched"]}
        assert matched_genes.get("CD19-TotalSeqA") == "CD19"

    def test_unmatched_proteins_recorded(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, inplace=False)
        assert "Mouse-IgG1" in corr_dict["unmatched"]
        assert "Rat-IgG2b" in corr_dict["unmatched"]
        assert "NoMatch-XYZ" in corr_dict["unmatched"]

    def test_matched_by_values_valid(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        if not results_df.empty:
            valid_values = {"exact", "stripped"}
            assert set(results_df["matched_by"].unique()).issubset(valid_values)

    def test_provenance_counts_correct(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, inplace=False)
        prov = corr_dict["provenance"]
        assert prov["n_unmatched"] == 3  # Mouse-IgG1, Rat-IgG2b, NoMatch-XYZ
        assert prov["n_matched"] == 5    # CD3E, CD4, CD19-TotalSeqA, CD14, CD56


# ---------------------------------------------------------------------------
# Results DataFrame tests
# ---------------------------------------------------------------------------

class TestResultsDataFrame:

    def test_expected_columns(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        expected = {"protein", "gene", "r", "pval", "pval_adj", "n_cells", "matched_by"}
        assert expected.issubset(set(results_df.columns))

    def test_r_values_in_valid_range(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        if not results_df.empty:
            assert results_df["r"].between(-1, 1).all(), (
                "r values outside [-1, 1] range"
            )

    def test_sorted_descending_by_r(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        if len(results_df) > 1:
            assert results_df["r"].is_monotonic_decreasing, (
                "Results not sorted by r descending"
            )

    def test_n_cells_positive(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        if not results_df.empty:
            assert (results_df["n_cells"] > 0).all()

    def test_pval_in_valid_range(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        results_df = corr_dict["results"]
        if not results_df.empty:
            assert results_df["pval"].between(0, 1).all()

    def test_empty_dataframe_when_no_match(self):
        """All proteins unmatched → empty results DataFrame."""
        rng = np.random.default_rng(42)
        rna = AnnData(X=sp.csr_matrix(
            rng.poisson(1, size=(50, 10)).astype(np.float32)
        ))
        rna.obs_names = [f"cell_{i}" for i in range(50)]
        rna.var_names = [f"RNAGene{i}" for i in range(10)]

        adt = AnnData(X=sp.csr_matrix(
            np.abs(rng.normal(0, 1, size=(50, 3))).astype(np.float32)
        ))
        adt.obs_names = [f"cell_{i}" for i in range(50)]
        adt.var_names = ["Mouse-IgG1", "Rat-IgG2b", "FakeProtein"]
        adt.layers["adt_clr"] = adt.X.toarray().copy()

        mdata = MuData({"rna": rna, "adt": adt})
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, corr_dict = cite_corr(mdata, inplace=False)
        assert corr_dict["results"].empty


# ---------------------------------------------------------------------------
# min_cells threshold tests
# ---------------------------------------------------------------------------

class TestMinCells:

    def test_skipped_below_min_cells(self, mdata_fixture):
        """Setting a very high min_cells should skip all pairs."""
        _, corr_dict = cite_corr(mdata_fixture, min_cells=10_000, inplace=False)
        assert corr_dict["results"].empty
        assert len(corr_dict["skipped"]) > 0

    def test_skipped_dict_has_reason(self, mdata_fixture):
        _, corr_dict = cite_corr(mdata_fixture, min_cells=10_000, inplace=False)
        for item in corr_dict["skipped"]:
            assert "reason" in item
            assert "n_cells" in item

    def test_min_cells_one_includes_all_matched(self, mdata_fixture):
        """min_cells=1 should return results for all matched pairs."""
        _, corr_dict = cite_corr(mdata_fixture, min_cells=1, inplace=False)
        n_results  = len(corr_dict["results"])
        n_matched  = len(corr_dict["matched"])
        n_skipped  = len(corr_dict["skipped"])
        assert n_results + n_skipped == n_matched


# ---------------------------------------------------------------------------
# Warning tests
# ---------------------------------------------------------------------------

class TestWarnings:

    def test_warns_when_adt_clr_absent(self, mdata_fixture):
        del mdata_fixture["adt"].layers["adt_clr"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cite_corr(mdata_fixture, inplace=False)
        assert any("adt_clr" in str(warning.message) for warning in w)

    def test_warns_when_logcounts_absent(self, mdata_fixture):
        del mdata_fixture["rna"].layers["logcounts"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cite_corr(mdata_fixture, inplace=False)
        assert any("logcounts" in str(warning.message) for warning in w)

    def test_warns_when_no_pairs_matched(self):
        rng = np.random.default_rng(0)
        rna = AnnData(X=sp.csr_matrix(
            rng.poisson(1, size=(50, 5)).astype(np.float32)
        ))
        rna.obs_names = [f"cell_{i}" for i in range(50)]
        rna.var_names = [f"Gene{i}" for i in range(5)]

        adt = AnnData(X=sp.csr_matrix(
            np.abs(rng.normal(0, 1, size=(50, 3))).astype(np.float32)
        ))
        adt.obs_names = [f"cell_{i}" for i in range(50)]
        adt.var_names = ["IgG1", "IgG2", "IgG3"]
        adt.layers["adt_clr"] = adt.X.toarray().copy()

        mdata = MuData({"rna": rna, "adt": adt})
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cite_corr(mdata, inplace=False)
        assert any("no protein names matched" in str(warning.message).lower()
                   for warning in w)
