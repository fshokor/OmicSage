"""
tests/test_cite_deg.py — Tests for cite_deg (Phase 4, Step 7)

Coverage
--------
DPE (differential protein expression):
  - runs on MuData with adt_celltype_manual groupby
  - runs on AnnData fallback (DPE only, no cross-modal)
  - falls back to adt_celltype_score when manual labels absent
  - falls back to leiden when both celltype cols absent
  - raises ValueError when no valid groupby found
  - filters isotype control proteins via exclude_protein_prefixes
  - result DataFrames have expected columns (protein, score, pval, logfc, pval_adj)
  - provenance written to uns["omicsage_cite_deg"]
  - DPE result stored in adt.uns["rank_genes_groups_dpe"]

Cross-modal RNA DEG:
  - runs on MuData and produces rna_crossmodal dict
  - result DataFrames have gene column (not protein)
  - rank_genes_groups_rna_cm written to mdata["rna"].uns
  - skipped with UserWarning when AnnData fallback used
  - warns when logcounts layer absent

Return contract:
  - returns (MuData, dict) or (AnnData, dict)
  - dict has keys: dpe, rna_crossmodal, dpe_summary, rna_summary, provenance, input_type
  - input_type is "mudata" for MuData input, "anndata" for AnnData input
  - inplace=False returns a new object, does not modify original
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

mudata = pytest.importorskip("mudata")
MuData = mudata.MuData


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_rna(n_obs: int = 80, n_genes: int = 50, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32))
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"Gene{i}" for i in range(n_genes)]
    # logcounts layer (required for cross-modal DEG)
    adata.layers["logcounts"] = X.toarray().copy()
    return adata


def _make_adt(n_obs: int = 80, n_proteins: int = 20, seed: int = 1) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(3.0, size=(n_obs, n_proteins)).astype(np.float32))
    adt = AnnData(X=X)
    adt.obs_names = [f"cell_{i}" for i in range(n_obs)]
    protein_names = [f"CD{i}" for i in range(n_proteins - 2)] + ["Mouse-IgG1", "Rat-IgG2b"]
    adt.var_names = protein_names
    # CLR layer
    clr = rng.normal(0, 1, size=(n_obs, n_proteins)).astype(np.float32)
    adt.layers["adt_clr"] = clr
    adt.X = clr.copy()
    # Cell type labels — 4 groups of 20 cells each
    labels = np.repeat(["CD4 T", "CD8 T", "B", "Mono"], n_obs // 4)
    adt.obs["adt_celltype_manual"] = labels
    adt.obs["adt_celltype_score"]  = labels
    adt.obs["leiden"]              = [str(i % 4) for i in range(n_obs)]
    return adt


@pytest.fixture
def mdata_fixture():
    """MuData with RNA + ADT, shared barcodes."""
    rna = _make_rna()
    adt = _make_adt()
    return MuData({"rna": rna, "adt": adt})


@pytest.fixture
def adt_fixture():
    """Bare AnnData (ADT only) — cite_05 fallback path."""
    return _make_adt()


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------

from pipeline.modules.scripts.cite.cite_deg import cite_deg


# ---------------------------------------------------------------------------
# DPE tests
# ---------------------------------------------------------------------------

class TestDPE:

    def test_mudata_returns_tuple(self, mdata_fixture):
        result, deg_dict = cite_deg(mdata_fixture, inplace=False)
        assert isinstance(result, MuData)
        assert isinstance(deg_dict, dict)

    def test_anndata_fallback_returns_tuple(self, adt_fixture):
        result, deg_dict = cite_deg(adt_fixture, inplace=False)
        assert isinstance(result, AnnData)
        assert isinstance(deg_dict, dict)

    def test_dpe_dict_keys(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        for key in ["dpe", "rna_crossmodal", "dpe_summary", "rna_summary",
                    "provenance", "input_type"]:
            assert key in deg_dict, f"Missing key: {key}"

    def test_dpe_protein_column(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        for cell_type, df in deg_dict["dpe"].items():
            assert "protein" in df.columns, (
                f"DPE DataFrame for '{cell_type}' missing 'protein' column"
            )

    def test_dpe_expected_columns(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        expected = {"protein", "score", "pval", "logfc", "pval_adj"}
        for cell_type, df in deg_dict["dpe"].items():
            assert expected.issubset(df.columns), (
                f"DPE DataFrame for '{cell_type}' missing columns: "
                f"{expected - set(df.columns)}"
            )

    def test_dpe_groups_match_cell_types(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        expected_types = {"CD4 T", "CD8 T", "B", "Mono"}
        assert set(deg_dict["dpe"].keys()) == expected_types

    def test_dpe_written_to_uns(self, mdata_fixture):
        result, _ = cite_deg(mdata_fixture, inplace=False)
        assert "rank_genes_groups_dpe" in result["adt"].uns

    def test_provenance_in_uns(self, mdata_fixture):
        result, _ = cite_deg(mdata_fixture, inplace=False)
        assert "omicsage_cite_deg" in result.uns
        prov = result.uns["omicsage_cite_deg"]
        assert prov["module"] == "cite_deg"
        assert "timestamp" in prov
        assert "n_cell_types" in prov

    def test_input_type_mudata(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        assert deg_dict["input_type"] == "mudata"

    def test_input_type_anndata(self, adt_fixture):
        _, deg_dict = cite_deg(adt_fixture, inplace=False)
        assert deg_dict["input_type"] == "anndata"

    def test_exclude_protein_prefixes(self, mdata_fixture):
        _, deg_dict = cite_deg(
            mdata_fixture,
            exclude_protein_prefixes=["Mouse-IgG", "Rat-IgG"],
            inplace=False,
        )
        for cell_type, df in deg_dict["dpe"].items():
            if df.empty:
                continue
            for protein in df["protein"]:
                assert not protein.upper().startswith("MOUSE-IGG"), (
                    f"Isotype control not filtered: {protein}"
                )
                assert not protein.upper().startswith("RAT-IGG"), (
                    f"Isotype control not filtered: {protein}"
                )

    def test_inplace_false_does_not_modify_original(self, mdata_fixture):
        original_uns_keys = set(mdata_fixture.uns.keys())
        cite_deg(mdata_fixture, inplace=False)
        assert set(mdata_fixture.uns.keys()) == original_uns_keys

    def test_inplace_true_modifies_original(self, mdata_fixture):
        cite_deg(mdata_fixture, inplace=True)
        assert "omicsage_cite_deg" in mdata_fixture.uns

    def test_dpe_summary_dataframe(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        summary = deg_dict["dpe_summary"]
        assert isinstance(summary, __import__("pandas").DataFrame)
        if not summary.empty:
            assert "group" in summary.columns
            assert "protein" in summary.columns


# ---------------------------------------------------------------------------
# Groupby fallback tests
# ---------------------------------------------------------------------------

class TestGroupbyFallback:

    def test_falls_back_to_score_when_manual_absent(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual"], inplace=True
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        assert deg_dict["provenance"]["groupby"] == "adt_celltype_score"
        assert any("adt_celltype_score" in str(warning.message) for warning in w)

    def test_falls_back_to_leiden_when_both_absent(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual", "adt_celltype_score"], inplace=True
        )
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        assert deg_dict["provenance"]["groupby"] == "leiden"

    def test_raises_when_no_valid_groupby(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual", "adt_celltype_score", "leiden"],
            inplace=True,
        )
        with pytest.raises(ValueError, match="None of"):
            cite_deg(mdata_fixture, inplace=False)


# ---------------------------------------------------------------------------
# Cross-modal RNA DEG tests
# ---------------------------------------------------------------------------

class TestCrossModalRNADEG:

    def test_rna_crossmodal_populated_for_mudata(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        assert isinstance(deg_dict["rna_crossmodal"], dict)
        assert len(deg_dict["rna_crossmodal"]) > 0

    def test_rna_crossmodal_gene_column(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        for cell_type, df in deg_dict["rna_crossmodal"].items():
            assert "gene" in df.columns, (
                f"RNA cross-modal DataFrame for '{cell_type}' missing 'gene' column"
            )

    def test_rna_crossmodal_empty_for_anndata(self, adt_fixture):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, deg_dict = cite_deg(adt_fixture, inplace=False)
        assert deg_dict["rna_crossmodal"] == {}
        assert any("cross-modal" in str(warning.message).lower() for warning in w)

    def test_rank_genes_rna_cm_in_uns(self, mdata_fixture):
        result, _ = cite_deg(mdata_fixture, inplace=False)
        assert "rank_genes_groups_rna_cm" in result["rna"].uns

    def test_exclude_gene_prefixes(self, mdata_fixture):
        """Gene prefixes excluded from cross-modal results."""
        # Name some genes with RPL prefix so we can test exclusion
        mdata_fixture["rna"].var_names = (
            [f"RPL{i}" for i in range(10)]
            + [f"Gene{i}" for i in range(10, 50)]
        )
        _, deg_dict = cite_deg(
            mdata_fixture,
            exclude_gene_prefixes=["RPL"],
            inplace=False,
        )
        for cell_type, df in deg_dict["rna_crossmodal"].items():
            if df.empty:
                continue
            for gene in df["gene"]:
                assert not gene.upper().startswith("RPL"), (
                    f"RPL gene not filtered from cross-modal DEG: {gene}"
                )

    def test_rna_summary_dataframe(self, mdata_fixture):
        _, deg_dict = cite_deg(mdata_fixture, inplace=False)
        summary = deg_dict["rna_summary"]
        assert isinstance(summary, __import__("pandas").DataFrame)
        if not summary.empty:
            assert "gene" in summary.columns

    def test_warns_when_logcounts_absent(self, mdata_fixture):
        del mdata_fixture["rna"].layers["logcounts"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cite_deg(mdata_fixture, inplace=False)
        assert any("logcounts" in str(warning.message) for warning in w)
