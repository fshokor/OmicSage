"""
test_gsea.py — Tests for pipeline/modules/downstream/gsea.py

Run with:
    cd ~/OmicSage
    conda activate omicsage
    python -m pytest tests/test_gsea.py -v

Expected: 8 passed

Notes
-----
- All Enrichr API calls are mocked — tests are offline-safe and CI-safe.
- Query gene lists use real HGNC symbols so fixtures are biologically plausible.
- The mock fixture returns a minimal but correctly structured enrichr results
  DataFrame, matching what gseapy.enrichr().results produces.
- gseapy must be installed for these tests to run. In CI environments where
  gseapy is not available, the entire module is skipped gracefully via
  pytest.importorskip — no ERRORs, just SKIPs.
"""

from __future__ import annotations

import os
import warnings
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

# ---------------------------------------------------------------------------
# Module-level skip guard — must come before any gsea import.
# If gseapy is not installed the whole test file is skipped cleanly.
# This also respects the OMICSAGE_CI flag used elsewhere in the test suite.
# ---------------------------------------------------------------------------
gseapy = pytest.importorskip(
    "gseapy",
    reason="gseapy not installed — skipping GSEA tests. "
           "Install with: pip install gseapy",
)

from pipeline.modules.scripts.downstream.gsea import gsea  # noqa: E402  (after importorskip)


# ---------------------------------------------------------------------------
# Mock helpers
# ---------------------------------------------------------------------------

def _mock_enrichr_result(n_terms: int = 10) -> MagicMock:
    """
    Return a mock gseapy enrichr result object whose .results attribute
    is a DataFrame matching the expected gseapy output schema.
    """
    terms = [f"Mock pathway term {i}" for i in range(n_terms)]
    df = pd.DataFrame({
        "Gene_set":          ["GO_Biological_Process_2023"] * n_terms,
        "Term":              terms,
        "Overlap":           ["3/50"] * n_terms,
        "P-value":           np.linspace(0.001, 0.05, n_terms),
        "Adjusted P-value":  np.linspace(0.01,  0.05, n_terms),
        "Old P-value":       [0.0] * n_terms,
        "Old Adjusted P-value": [0.0] * n_terms,
        "Odds Ratio":        [1.5] * n_terms,
        "Combined Score":    [10.0] * n_terms,
        "Genes":             ["CD3D;CD3E;CD4"] * n_terms,
    })
    mock_result = MagicMock()
    mock_result.results = df
    return mock_result


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Real HGNC gene symbols — gives the fixture biological plausibility
_T_CELL_GENES = [
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL7R", "TCF7",
    "CCR7", "SELL", "LEF1", "FOXP3", "TIGIT", "LAG3", "PDCD1",
    "GZMB", "GZMK", "PRF1", "NKG7", "IFNG",
]
_MONOCYTE_GENES = [
    "CD14", "LYZ", "S100A8", "S100A9", "VCAN", "FCN1", "PSAP",
    "CTSS", "AIF1", "CST3", "TYROBP", "LST1", "SERPINA1",
    "CFP", "NEAT1", "FGL2", "WARS", "GBP1", "IFITM3", "IFI30",
]
_B_CELL_GENES = [
    "CD19", "MS4A1", "CD79A", "CD79B", "HLA-DRA", "HLA-DRB1",
    "IGHM", "IGKC", "IGLC2", "BANK1", "PAX5", "CD22",
    "FCER2", "CR2", "TCL1A", "BLK", "FCRL1", "FCRL2",
]

_ALL_MARKER_GENES = list(set(_T_CELL_GENES + _MONOCYTE_GENES + _B_CELL_GENES))

# Pad to a larger realistic var_names list (~200 genes)
_BACKGROUND_GENES = [f"GENE{i:04d}" for i in range(200 - len(_ALL_MARKER_GENES))]
_ALL_GENES = _ALL_MARKER_GENES + _BACKGROUND_GENES


def _make_adata(
    n_cells: int = 300,
    seed: int = 42,
) -> AnnData:
    """
    Minimal AnnData with realistic gene names and logcounts layer.
    Three cell types: T_cell, Monocyte, B_cell.
    """
    rng = np.random.default_rng(seed)
    n_genes = len(_ALL_GENES)
    n_per_group = n_cells // 3

    X = rng.negative_binomial(2, 0.3, size=(n_cells, n_genes)).astype(float)
    logcounts = np.log1p(X / (X.sum(axis=1, keepdims=True) + 1e-9) * 1e4)

    adata = AnnData(X=logcounts)
    adata.layers["logcounts"] = logcounts.copy()
    adata.obs["cell_type_vote"] = pd.Categorical(
        ["T_cell"] * n_per_group +
        ["Monocyte"] * n_per_group +
        ["B_cell"] * (n_cells - 2 * n_per_group)
    )
    adata.var_names = _ALL_GENES
    return adata


def _make_deg_dict(
    groups: list[str] | None = None,
    n_degs_per_group: int = 15,
) -> dict:
    """
    Synthetic deg_dict['results'] with real gene symbols as DEGs.
    All genes have logfc > 0 and pval_adj < 0.05 so they pass ORA filters.
    """
    if groups is None:
        groups = ["T_cell", "Monocyte", "B_cell"]

    gene_pools = {
        "T_cell":   _T_CELL_GENES,
        "Monocyte": _MONOCYTE_GENES,
        "B_cell":   _B_CELL_GENES,
    }

    results = {}
    for group in groups:
        pool = gene_pools.get(group, _T_CELL_GENES)
        genes = pool[:n_degs_per_group]
        df = pd.DataFrame({
            "gene":     genes,
            "score":    np.linspace(10, 2, len(genes)),
            "pval":     np.linspace(0.001, 0.04, len(genes)),
            "logfc":    np.linspace(2.0, 0.3, len(genes)),
            "pval_adj": np.linspace(0.005, 0.04, len(genes)),
        })
        results[group] = df

    return {
        "results": results,
        "provenance": {
            "groupby": "cell_type_vote",
            "method": "wilcoxon",
        },
    }


@pytest.fixture(scope="module")
def adata_base():
    return _make_adata()


@pytest.fixture(scope="module")
def deg_dict_base():
    return _make_deg_dict()


@pytest.fixture(scope="module")
def gsea_results(adata_base, deg_dict_base):
    """Run gsea() once (mocked) and share across structural tests."""
    with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
         patch("gseapy.get_library_name", return_value=_DEFAULT_SETS):
        adata_out, gsea_dict = gsea(
            adata_base,
            deg_dict=deg_dict_base,
            gene_sets=["GO_Biological_Process_2023"],
            inplace=False,
        )
    return adata_out, gsea_dict


_DEFAULT_SETS = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestGseaReturnsCorrectTypes:
    def test_gsea_returns_anndata_and_dict(self, adata_base, deg_dict_base):
        """gsea() must return (AnnData, dict)."""
        with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
             patch("gseapy.get_library_name", return_value=_DEFAULT_SETS):
            result = gsea(
                adata_base,
                deg_dict=deg_dict_base,
                gene_sets=["GO_Biological_Process_2023"],
                inplace=False,
            )
        assert isinstance(result, tuple), "gsea() should return a tuple"
        assert len(result) == 2, "gsea() should return exactly 2 items"
        adata_out, gsea_dict = result
        assert isinstance(adata_out, AnnData)
        assert isinstance(gsea_dict, dict)


class TestGseaProvenance:
    def test_gsea_uns_provenance_keys(self, gsea_results):
        """adata.uns['omicsage_gsea'] must contain all required provenance keys."""
        adata_out, _ = gsea_results
        assert "omicsage_gsea" in adata_out.uns, \
            "adata.uns['omicsage_gsea'] missing"

        required_keys = {
            "groupby", "gene_sets", "organism", "min_logfc",
            "max_pval_adj", "top_n_genes", "min_genes",
            "n_groups_tested", "n_groups_skipped", "groups_skipped",
            "gseapy_version", "omicsage_module", "timestamp",
        }
        present = set(adata_out.uns["omicsage_gsea"].keys())
        missing = required_keys - present
        assert not missing, f"Missing provenance keys: {missing}"

    def test_gsea_provenance_omicsage_module_value(self, gsea_results):
        """omicsage_module provenance value must be 'gsea'."""
        adata_out, _ = gsea_results
        assert adata_out.uns["omicsage_gsea"]["omicsage_module"] == "gsea"


class TestGseaOutputStructure:
    def test_gsea_output_columns(self, gsea_results):
        """Each group DataFrame must contain the 5 required columns."""
        _, gsea_dict = gsea_results
        required_cols = {"Term", "Overlap", "P-value", "Adjusted P-value", "Genes"}
        for group, df in gsea_dict["results"].items():
            assert isinstance(df, pd.DataFrame), \
                f"results['{group}'] should be a DataFrame"
            missing = required_cols - set(df.columns)
            assert not missing, \
                f"results['{group}'] missing columns: {missing}"

    def test_gsea_every_group_has_results_or_is_skipped(self, adata_base, deg_dict_base):
        """Every group must be in results OR in skipped — none silently dropped."""
        with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
             patch("gseapy.get_library_name", return_value=_DEFAULT_SETS):
            _, gsea_dict = gsea(
                adata_base,
                deg_dict=deg_dict_base,
                gene_sets=["GO_Biological_Process_2023"],
                inplace=False,
            )
        all_groups    = set(deg_dict_base["results"].keys())
        result_groups = set(gsea_dict["results"].keys())
        skipped       = set(gsea_dict["skipped"])
        assert all_groups == result_groups | skipped, \
            f"Groups unaccounted for: {all_groups - result_groups - skipped}"

    def test_gsea_pval_range(self, gsea_results):
        """All P-value and Adjusted P-value values must be in [0, 1]."""
        _, gsea_dict = gsea_results
        for group, df in gsea_dict["results"].items():
            if df.empty:
                continue
            pvals = pd.to_numeric(df["P-value"], errors="coerce").dropna()
            padjs = pd.to_numeric(df["Adjusted P-value"], errors="coerce").dropna()
            assert pvals.between(0.0, 1.0).all(), \
                f"P-value out of [0,1] range in group '{group}'"
            assert padjs.between(0.0, 1.0).all(), \
                f"Adjusted P-value out of [0,1] range in group '{group}'"


class TestGseaMinGenesSkip:
    def test_gsea_min_genes_skip_warns_and_skips(self, adata_base):
        """Groups with fewer DEGs than min_genes must be skipped with a UserWarning."""
        # Only 2 DEGs per group — below default min_genes=5
        deg_dict_tiny = _make_deg_dict(n_degs_per_group=2)

        with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
             patch("gseapy.get_library_name", return_value=_DEFAULT_SETS), \
             warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _, gsea_dict = gsea(
                adata_base,
                deg_dict=deg_dict_tiny,
                gene_sets=["GO_Biological_Process_2023"],
                min_genes=5,
                inplace=False,
            )

        warning_messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("min_genes" in msg or "Skipping" in msg for msg in warning_messages), \
            "Expected UserWarning about skipping group due to min_genes"
        assert len(gsea_dict["skipped"]) == 3, \
            "All 3 groups should be in skipped when DEGs < min_genes"
        assert len(gsea_dict["results"]) == 0, \
            "results should be empty when all groups are skipped"


class TestGseaInplace:
    def test_gsea_inplace_false(self, adata_base, deg_dict_base):
        """inplace=False must not modify the caller's AnnData."""
        original_uns_keys = set(adata_base.uns.keys())
        original_X_sum    = float(np.array(adata_base.X).sum())

        with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
             patch("gseapy.get_library_name", return_value=_DEFAULT_SETS):
            gsea(
                adata_base,
                deg_dict=deg_dict_base,
                gene_sets=["GO_Biological_Process_2023"],
                inplace=False,
            )

        assert set(adata_base.uns.keys()) == original_uns_keys, \
            "inplace=False must not add keys to original adata.uns"
        assert float(np.array(adata_base.X).sum()) == pytest.approx(original_X_sum), \
            "inplace=False must not modify original adata.X"


class TestGseaGeneSetsParam:
    def test_gsea_single_gene_set_string_accepted(self, adata_base, deg_dict_base):
        """Passing a single gene set as a string (not list) must work."""
        with patch("gseapy.enrichr", return_value=_mock_enrichr_result()), \
             patch("gseapy.get_library_name", return_value=_DEFAULT_SETS):
            adata_out, gsea_dict = gsea(
                adata_base,
                deg_dict=deg_dict_base,
                gene_sets="GO_Biological_Process_2023",   # ← string, not list
                inplace=False,
            )
        # Provenance should store it as a list
        assert isinstance(adata_out.uns["omicsage_gsea"]["gene_sets"], list), \
            "gene_sets should be stored as a list in provenance even if passed as string"
        assert "GO_Biological_Process_2023" in adata_out.uns["omicsage_gsea"]["gene_sets"]
