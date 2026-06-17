"""
tests/test_cite_epitope.py — Tests for cite_epitope (Phase 4, Step 10)

Coverage
--------
Input validation:
  - raises ValueError on AnnData input

Panel resolution:
  - preset="bmmc" loads built-in panel
  - custom epitope_panels overrides preset
  - unknown preset warns and uses empty panels
  - no preset and no custom warns
  - custom > preset precedence

Epitope scoring (sub-analysis A):
  - obs columns written for each panel
  - obs column names follow "epitope_score_<panel>" convention
  - scores are numeric (float)
  - panels with no matching proteins score all zeros with warning

Groupby fallback (same as cite_deg):
  - falls back to adt_celltype_score when manual absent
  - falls back to leiden when both absent
  - raises ValueError when none found

Cross-modal marker table (sub-analysis D):
  - empty DataFrame when dpe_results=None
  - correct columns when dpe_results provided
  - r column populated when corr_results provided
  - r column is NaN placeholder when corr_results=None

Return contract:
  - returns (MuData, dict)
  - dict has keys: scores_df, panel_used, n_panels, marker_table, provenance
  - provenance written to mdata.uns["omicsage_cite_epitope"]
  - inplace=False does not modify original
  - inplace=True modifies original
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

def _make_rna(n_obs: int = 80, n_genes: int = 30, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32))
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = (
        ["CD3E", "CD4", "CD8A", "CD19", "MS4A1",
         "CD14", "FCGR3A", "NKG7", "GNLY", "CD56"]
        + [f"Gene{i}" for i in range(20)]
    )
    adata.layers["logcounts"] = X.toarray().copy()
    return adata


def _make_adt(n_obs: int = 80, seed: int = 1) -> AnnData:
    rng = np.random.default_rng(seed)
    proteins = [
        "CD3E", "CD4", "CD8a", "CD19", "CD20",
        "CD14", "CD11c", "CD56", "CD71", "CD36",
        "Mouse-IgG1", "Rat-IgG2b",
    ]
    clr = np.abs(rng.normal(0, 1, size=(n_obs, len(proteins)))).astype(np.float32)
    adt = AnnData(X=clr.copy())
    adt.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adt.var_names = proteins
    adt.layers["adt_clr"] = clr.copy()
    labels = np.repeat(["CD4 T", "CD8 T", "B", "Mono"], n_obs // 4)
    adt.obs["adt_celltype_manual"] = labels
    adt.obs["adt_celltype_score"]  = labels
    adt.obs["leiden"] = [str(i % 4) for i in range(n_obs)]
    return adt


def _make_dpe_results() -> dict[str, pd.DataFrame]:
    cell_types = ["CD4 T", "CD8 T", "B", "Mono"]
    proteins = ["CD3E", "CD4", "CD8a", "CD19", "CD14"]
    results = {}
    for i, ct in enumerate(cell_types):
        results[ct] = pd.DataFrame({
            "protein":  proteins,
            "logfc":    np.linspace(2.0, 0.5, len(proteins)),
            "pval":     np.full(len(proteins), 0.001),
            "pval_adj": np.full(len(proteins), 0.01),
            "score":    np.linspace(5, 1, len(proteins)),
        })
    return results


def _make_corr_results() -> pd.DataFrame:
    return pd.DataFrame({
        "protein":    ["CD3E", "CD4", "CD8a", "CD19", "CD14"],
        "gene":       ["CD3E", "CD4", "CD8A", "CD19", "CD14"],
        "r":          [0.55,   0.48,  0.42,   0.61,   0.31],
        "pval":       [0.001] * 5,
        "pval_adj":   [0.01]  * 5,
        "n_cells":    [60] * 5,
        "matched_by": ["exact"] * 5,
    })


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

from pipeline.modules.scripts.cite.cite_epitope import cite_epitope


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_raises_on_anndata_input(self, adt_only_fixture):
        with pytest.raises(ValueError, match="MuData"):
            cite_epitope(adt_only_fixture)


# ---------------------------------------------------------------------------
# Panel resolution
# ---------------------------------------------------------------------------

class TestPanelResolution:

    def test_preset_bmmc_loads(self, mdata_fixture):
        _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert d["n_panels"] > 0
        assert "T_cell" in d["panel_used"] or "CD4_T" in d["panel_used"]

    def test_custom_panels_override_preset(self, mdata_fixture):
        custom = {"MyPanel": ["CD3", "CD4"]}
        _, d = cite_epitope(
            mdata_fixture,
            preset="bmmc",
            epitope_panels=custom,
            inplace=False,
        )
        assert list(d["panel_used"].keys()) == ["MyPanel"]

    def test_unknown_preset_warns(self, mdata_fixture):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, d = cite_epitope(mdata_fixture, preset="unknown_xyz", inplace=False)
        assert any("unknown_xyz" in str(warning.message).lower() for warning in w)
        assert d["n_panels"] == 0

    def test_no_panel_warns(self, mdata_fixture):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            cite_epitope(mdata_fixture, inplace=False)
        assert any("no epitope_panels" in str(warning.message).lower()
                   or "preset" in str(warning.message).lower()
                   for warning in w)

    def test_custom_precedence_over_preset(self, mdata_fixture):
        custom = {"Custom": ["CD3"]}
        _, d = cite_epitope(
            mdata_fixture,
            preset="bmmc",
            epitope_panels=custom,
            inplace=False,
        )
        assert "Custom" in d["panel_used"]
        assert "T_cell" not in d["panel_used"]


# ---------------------------------------------------------------------------
# Epitope scoring
# ---------------------------------------------------------------------------

class TestEpitopeScoring:

    def test_obs_columns_written(self, mdata_fixture):
        result, d = cite_epitope(
            mdata_fixture, preset="bmmc", inplace=False
        )
        for panel_name in d["panel_used"]:
            safe = panel_name.replace(" ", "_").replace("-", "_")
            col = f"epitope_score_{safe}"
            assert col in result["adt"].obs.columns, f"Missing obs column: {col}"

    def test_scores_are_numeric(self, mdata_fixture):
        result, d = cite_epitope(
            mdata_fixture, preset="bmmc", inplace=False
        )
        for panel_name in d["panel_used"]:
            safe = panel_name.replace(" ", "_").replace("-", "_")
            col = f"epitope_score_{safe}"
            assert pd.api.types.is_numeric_dtype(result["adt"].obs[col])

    def test_unmatched_panel_warns_and_scores_zero(self, mdata_fixture):
        custom = {"FakePanel": ["NonExistentProtein999"]}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result, _ = cite_epitope(
                mdata_fixture, epitope_panels=custom, inplace=False
            )
        assert any("FakePanel" in str(warning.message) for warning in w)
        assert (result["adt"].obs["epitope_score_FakePanel"] == 0).all()

    def test_scores_df_has_cell_column(self, mdata_fixture):
        _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert "cell" in d["scores_df"].columns


# ---------------------------------------------------------------------------
# Groupby fallback
# ---------------------------------------------------------------------------

class TestGroupbyFallback:

    def test_falls_back_to_score(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual"], inplace=True
        )
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert d["provenance"]["groupby"] == "adt_celltype_score"

    def test_falls_back_to_leiden(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual", "adt_celltype_score"], inplace=True
        )
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert d["provenance"]["groupby"] == "leiden"

    def test_raises_when_no_groupby(self, mdata_fixture):
        mdata_fixture["adt"].obs.drop(
            columns=["adt_celltype_manual", "adt_celltype_score", "leiden"],
            inplace=True,
        )
        with pytest.raises(ValueError, match="None of"):
            cite_epitope(mdata_fixture, preset="bmmc", inplace=False)


# ---------------------------------------------------------------------------
# Cross-modal marker table
# ---------------------------------------------------------------------------

class TestMarkerTable:

    def test_empty_when_no_dpe(self, mdata_fixture):
        _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert d["marker_table"].empty

    def test_has_columns_with_dpe(self, mdata_fixture):
        _, d = cite_epitope(
            mdata_fixture,
            preset="bmmc",
            dpe_results=_make_dpe_results(),
            inplace=False,
        )
        expected = {"cell_type", "protein", "logfc", "pval_adj", "rna_gene", "r"}
        assert expected.issubset(set(d["marker_table"].columns))

    def test_r_populated_with_corr(self, mdata_fixture):
        _, d = cite_epitope(
            mdata_fixture,
            preset="bmmc",
            dpe_results=_make_dpe_results(),
            corr_results=_make_corr_results(),
            inplace=False,
        )
        mt = d["marker_table"]
        # CD3E has r=0.55 in our fixture
        cd3e_rows = mt[mt["protein"] == "CD3E"]
        if not cd3e_rows.empty:
            assert not np.isnan(cd3e_rows.iloc[0]["r"])

    def test_r_nan_without_corr(self, mdata_fixture):
        _, d = cite_epitope(
            mdata_fixture,
            preset="bmmc",
            dpe_results=_make_dpe_results(),
            corr_results=None,
            inplace=False,
        )
        mt = d["marker_table"]
        if not mt.empty:
            assert mt["r"].isna().all() or (mt["rna_gene"] == "—").all()


# ---------------------------------------------------------------------------
# Return contract
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self, mdata_fixture):
        result, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert isinstance(result, MuData)
        assert isinstance(d, dict)

    def test_dict_keys(self, mdata_fixture):
        _, d = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        for key in ["scores_df", "panel_used", "n_panels", "marker_table", "provenance"]:
            assert key in d

    def test_provenance_in_uns(self, mdata_fixture):
        result, _ = cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert "omicsage_cite_epitope" in result.uns
        prov = result.uns["omicsage_cite_epitope"]
        assert prov["module"] == "cite_epitope"
        assert "timestamp" in prov

    def test_inplace_false_no_modify(self, mdata_fixture):
        original_uns = set(mdata_fixture.uns.keys())
        cite_epitope(mdata_fixture, preset="bmmc", inplace=False)
        assert set(mdata_fixture.uns.keys()) == original_uns

    def test_inplace_true_modifies(self, mdata_fixture):
        cite_epitope(mdata_fixture, preset="bmmc", inplace=True)
        assert "omicsage_cite_epitope" in mdata_fixture.uns
