"""
tests/test_cite_gsea.py — Tests for cite_gsea (Phase 4, Step 8)

Coverage
--------
- raises ValueError when rna_crossmodal is empty (AnnData fallback)
- returns (MuData, dict) with expected keys
- provenance written to uns["omicsage_cite_gsea"]
- provenance module field is "cite_gsea" not "gsea"
- inplace=False does not modify original
- inplace=True modifies original
- wraps existing gsea() — results dict has expected structure
- skipped list present in return dict
- warns when uns["omicsage_cite_deg"] absent
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from unittest.mock import patch, MagicMock

mudata = pytest.importorskip("mudata")
MuData = mudata.MuData
gseapy = pytest.importorskip("gseapy")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_rna(n_obs: int = 80, n_genes: int = 50, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32))
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    # Use real human gene symbols so Enrichr has a chance to match
    gene_names = [
        "CD3E", "CD4", "CD8A", "CD19", "MS4A1", "CD14", "FCGR3A",
        "NKG7", "GNLY", "CD56", "IL7R", "CCR7", "S100A4", "GZMB",
        "PRF1", "FOXP3", "IL2RA", "CD38", "IGHG1", "JCHAIN",
    ] + [f"Gene{i}" for i in range(30)]
    adata.var_names = gene_names[:n_genes]
    adata.layers["logcounts"] = X.toarray().copy()
    return adata


def _make_adt(n_obs: int = 80) -> AnnData:
    rng = np.random.default_rng(1)
    X = sp.csr_matrix(rng.poisson(3.0, size=(n_obs, 10)).astype(np.float32))
    adt = AnnData(X=X)
    adt.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adt.var_names = [f"CD{i}" for i in range(10)]
    labels = np.repeat(["CD4 T", "CD8 T", "B", "Mono"], n_obs // 4)
    adt.obs["adt_celltype_manual"] = labels
    return adt


def _make_deg_dict(rna: AnnData) -> dict:
    """Build a minimal cite_deg_dict with synthetic RNA cross-modal results."""
    gene_names = rna.var_names.tolist()
    rna_crossmodal = {}
    for cell_type in ["CD4 T", "CD8 T", "B", "Mono"]:
        rna_crossmodal[cell_type] = pd.DataFrame({
            "gene":     gene_names[:10],
            "score":    np.linspace(5, 1, 10),
            "pval":     np.full(10, 0.001),
            "logfc":    np.linspace(2.0, 0.5, 10),
            "pval_adj": np.full(10, 0.01),
        })
    return {
        "rna_crossmodal": rna_crossmodal,
        "provenance": {
            "module": "cite_deg",
            "groupby": "adt_celltype_manual",
            "input_type": "mudata",
        },
        "input_type": "mudata",
    }


@pytest.fixture
def mdata_fixture():
    rna = _make_rna()
    adt = _make_adt()
    mdata = MuData({"rna": rna, "adt": adt})
    mdata.uns["omicsage_cite_deg"] = {
        "module": "cite_deg",
        "groupby": "adt_celltype_manual",
        "input_type": "mudata",
    }
    return mdata


@pytest.fixture
def deg_dict_fixture(mdata_fixture):
    return _make_deg_dict(mdata_fixture["rna"])


# ---------------------------------------------------------------------------
# Import
# ---------------------------------------------------------------------------

from pipeline.modules.cite.cite_gsea import cite_gsea


# ---------------------------------------------------------------------------
# Validation tests
# ---------------------------------------------------------------------------

class TestValidation:

    def test_raises_when_rna_crossmodal_empty(self, mdata_fixture):
        empty_deg_dict = {
            "rna_crossmodal": {},
            "provenance": {},
            "input_type": "anndata",
        }
        with pytest.raises(ValueError, match="rna_crossmodal"):
            cite_gsea(mdata_fixture, cite_deg_dict=empty_deg_dict, inplace=False)

    def test_warns_when_omicsage_cite_deg_absent(self, mdata_fixture, deg_dict_fixture):
        del mdata_fixture.uns["omicsage_cite_deg"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            try:
                cite_gsea(mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False)
            except Exception:
                pass  # gsea may fail due to network in CI — warning is what we test
        assert any("omicsage_cite_deg" in str(warning.message) for warning in w)


# ---------------------------------------------------------------------------
# Return contract tests
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {"CD4 T": pd.DataFrame()},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            result, gsea_dict = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        assert isinstance(result, MuData)
        assert isinstance(gsea_dict, dict)

    def test_dict_has_expected_keys(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            _, gsea_dict = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        for key in ["results", "summary_df", "provenance", "skipped"]:
            assert key in gsea_dict, f"Missing key: {key}"

    def test_provenance_module_is_cite_gsea(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            result, gsea_dict = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        assert gsea_dict["provenance"]["module"] == "cite_gsea"

    def test_provenance_written_to_uns(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            result, _ = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        assert "omicsage_cite_gsea" in result.uns
        assert result.uns["omicsage_cite_gsea"]["module"] == "cite_gsea"

    def test_skipped_list_present(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [("CD4 T", "up")],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            _, gsea_dict = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        assert isinstance(gsea_dict["skipped"], list)

    def test_source_deg_module_in_provenance(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            _, gsea_dict = cite_gsea(
                mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False
            )
        assert gsea_dict["provenance"]["source_deg_module"] == "cite_deg"


# ---------------------------------------------------------------------------
# Inplace tests
# ---------------------------------------------------------------------------

class TestInplace:

    def test_inplace_false_does_not_modify_original(self, mdata_fixture, deg_dict_fixture):
        original_uns_keys = set(mdata_fixture.uns.keys())
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            cite_gsea(mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=False)
        assert set(mdata_fixture.uns.keys()) == original_uns_keys

    def test_inplace_true_modifies_original(self, mdata_fixture, deg_dict_fixture):
        gsea_mock_result = (
            mdata_fixture["rna"],
            {
                "results": {},
                "summary_df": pd.DataFrame(),
                "provenance": {"module": "gsea", "timestamp": "2026-01-01"},
                "skipped": [],
            },
        )
        with patch(
            "pipeline.modules.cite.cite_gsea._rna_gsea",
            return_value=gsea_mock_result,
        ):
            cite_gsea(mdata_fixture, cite_deg_dict=deg_dict_fixture, inplace=True)
        assert "omicsage_cite_gsea" in mdata_fixture.uns
