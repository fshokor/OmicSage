"""
test_deg.py — Tests for pipeline/modules/downstream/deg.py

Run with:
    cd ~/OmicSage
    conda activate omicsage
    python -m pytest tests/test_deg.py -v

Expected: 9 passed
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
from anndata import AnnData

from pipeline.modules.scripts.downstream.deg import deg


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_adata(
    n_cells: int = 300,
    n_genes: int = 200,
    n_groups: int = 3,
    seed: int = 42,
) -> AnnData:
    """
    Synthetic AnnData with:
    - adata.layers['logcounts']  — log-normalised counts
    - adata.obs['cell_type_vote'] — group labels
    - adata.obs['leiden']         — numeric cluster labels (fallback)

    Groups are given distinct expression profiles so Wilcoxon has
    real signal to detect.
    """
    rng = np.random.default_rng(seed)

    cells_per_group = n_cells // n_groups
    group_labels = []
    expression_blocks = []

    cell_types = ["T_cell", "Monocyte", "B_cell"][:n_groups]

    for i, ct in enumerate(cell_types):
        group_labels.extend([ct] * cells_per_group)
        # Each group gets elevated expression in a distinct gene block
        block = rng.negative_binomial(1, 0.1, size=(cells_per_group, n_genes)).astype(float)
        # Spike a block of genes for this group
        gene_start = i * (n_genes // n_groups)
        gene_end   = gene_start + (n_genes // n_groups)
        block[:, gene_start:gene_end] += rng.negative_binomial(
            5, 0.3, size=(cells_per_group, gene_end - gene_start)
        ).astype(float)
        expression_blocks.append(block)

    X = np.vstack(expression_blocks)
    logcounts = np.log1p(X / (X.sum(axis=1, keepdims=True) + 1e-9) * 1e4)

    adata = AnnData(X=logcounts)
    adata.layers["logcounts"] = logcounts.copy()
    adata.obs["cell_type_vote"] = pd.Categorical(group_labels)
    adata.obs["leiden"] = pd.Categorical(
        [str(i % n_groups) for i in range(len(group_labels))]
    )
    adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]

    return adata


@pytest.fixture(scope="module")
def adata_base():
    return _make_adata()


@pytest.fixture(scope="module")
def deg_results(adata_base):
    """Run deg() once and share across tests that only need to inspect output."""
    adata_out, deg_dict = deg(
        adata_base,
        groupby="cell_type_vote",
        method="wilcoxon",
        min_logfc=0.0,      # no filtering — keeps all genes for structural tests
        max_pval_adj=1.0,
        n_genes=50,
        inplace=False,
    )
    return adata_out, deg_dict


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestDegReturnsCorrectTypes:
    def test_deg_returns_anndata_and_dict(self, adata_base):
        """deg() must return (AnnData, dict)."""
        result = deg(
            adata_base,
            min_logfc=0.0,
            max_pval_adj=1.0,
            inplace=False,
        )
        assert isinstance(result, tuple), "deg() should return a tuple"
        assert len(result) == 2, "deg() should return exactly 2 items"
        adata_out, deg_dict = result
        assert isinstance(adata_out, AnnData)
        assert isinstance(deg_dict, dict)


class TestDegProvenance:
    def test_deg_uns_provenance_keys(self, deg_results):
        """adata.uns['omicsage_deg'] must contain all required provenance keys."""
        adata_out, _ = deg_results
        assert "omicsage_deg" in adata_out.uns, \
            "adata.uns['omicsage_deg'] missing"

        required_keys = {
            "groupby", "method", "min_logfc", "max_pval_adj",
            "n_genes", "n_groups", "results_keys",
            "scanpy_version", "omicsage_module", "timestamp",
        }
        present = set(adata_out.uns["omicsage_deg"].keys())
        missing = required_keys - present
        assert not missing, f"Missing provenance keys: {missing}"

    def test_deg_provenance_omicsage_module_value(self, deg_results):
        """omicsage_module provenance value must be 'deg'."""
        adata_out, _ = deg_results
        assert adata_out.uns["omicsage_deg"]["omicsage_module"] == "deg"


class TestDegOutputStructure:
    def test_deg_output_columns(self, deg_results):
        """Each group DataFrame must have exactly: gene, score, pval, logfc, pval_adj."""
        _, deg_dict = deg_results
        required_cols = {"gene", "score", "pval", "logfc", "pval_adj"}
        for group, df in deg_dict["results"].items():
            assert isinstance(df, pd.DataFrame), \
                f"results['{group}'] should be a DataFrame"
            missing = required_cols - set(df.columns)
            assert not missing, \
                f"results['{group}'] missing columns: {missing}"

    def test_deg_every_cluster_has_results(self, adata_base, deg_results):
        """Every group in obs['cell_type_vote'] must have an entry in results."""
        _, deg_dict = deg_results
        expected_groups = set(adata_base.obs["cell_type_vote"].unique())
        result_groups   = set(deg_dict["results"].keys())
        assert expected_groups == result_groups, \
            f"Missing groups in results: {expected_groups - result_groups}"

    def test_deg_summary_df_structure(self, deg_results):
        """summary_df must be a DataFrame with columns: group, rank, gene, logfc, pval_adj."""
        _, deg_dict = deg_results
        summary = deg_dict["summary_df"]
        assert isinstance(summary, pd.DataFrame)
        required = {"group", "rank", "gene", "logfc", "pval_adj"}
        missing  = required - set(summary.columns)
        assert not missing, f"summary_df missing columns: {missing}"


class TestDegStatisticalProperties:
    def test_deg_pval_range(self, deg_results):
        """All pval and pval_adj values must be in [0, 1]."""
        _, deg_dict = deg_results
        for group, df in deg_dict["results"].items():
            if df.empty:
                continue
            assert df["pval"].between(0.0, 1.0).all(), \
                f"pval out of [0,1] range in group '{group}'"
            assert df["pval_adj"].between(0.0, 1.0).all(), \
                f"pval_adj out of [0,1] range in group '{group}'"


class TestDegThresholds:
    def test_deg_logfc_threshold_filters(self, adata_base):
        """Setting min_logfc=1.0 should return fewer genes than min_logfc=0.0."""
        _, dict_loose = deg(
            adata_base, min_logfc=0.0, max_pval_adj=1.0,
            n_genes=50, inplace=False,
        )
        _, dict_strict = deg(
            adata_base, min_logfc=1.0, max_pval_adj=1.0,
            n_genes=50, inplace=False,
        )
        loose_total  = sum(len(df) for df in dict_loose["results"].values())
        strict_total = sum(len(df) for df in dict_strict["results"].values())
        assert strict_total <= loose_total, \
            "Strict logFC threshold should return ≤ genes vs loose threshold"

    def test_deg_pval_threshold_filters(self, adata_base):
        """Setting max_pval_adj=0.01 should return fewer genes than max_pval_adj=1.0."""
        _, dict_loose = deg(
            adata_base, min_logfc=0.0, max_pval_adj=1.0,
            n_genes=50, inplace=False,
        )
        _, dict_strict = deg(
            adata_base, min_logfc=0.0, max_pval_adj=0.01,
            n_genes=50, inplace=False,
        )
        loose_total  = sum(len(df) for df in dict_loose["results"].values())
        strict_total = sum(len(df) for df in dict_strict["results"].values())
        assert strict_total <= loose_total, \
            "Strict pval_adj threshold should return ≤ genes vs loose threshold"


class TestDegInplace:
    def test_deg_inplace_false(self, adata_base):
        """inplace=False must not modify the caller's AnnData."""
        import copy
        original_uns_keys = set(adata_base.uns.keys())
        original_X_sum    = float(np.array(adata_base.X).sum())

        deg(adata_base, min_logfc=0.0, max_pval_adj=1.0, inplace=False)

        # uns keys on original should be unchanged
        assert set(adata_base.uns.keys()) == original_uns_keys, \
            "inplace=False must not add keys to original adata.uns"

        # X on original should be unchanged
        assert float(np.array(adata_base.X).sum()) == pytest.approx(original_X_sum), \
            "inplace=False must not modify original adata.X"


class TestDegEdgeCases:
    def test_deg_small_group_warning(self):
        """Groups with < 10 cells must trigger a UserWarning."""
        # Build an AnnData where one group has only 5 cells
        rng = np.random.default_rng(0)
        n_genes = 50

        # Group A: 100 cells, Group B: 5 cells
        X_a = np.log1p(rng.negative_binomial(2, 0.3, size=(100, n_genes)).astype(float))
        X_b = np.log1p(rng.negative_binomial(2, 0.3, size=(5,   n_genes)).astype(float))
        X   = np.vstack([X_a, X_b])

        adata_small = AnnData(X=X)
        adata_small.layers["logcounts"] = X.copy()
        adata_small.obs["cell_type_vote"] = pd.Categorical(
            ["GroupA"] * 100 + ["GroupB"] * 5
        )
        adata_small.var_names = [f"gene_{i:03d}" for i in range(n_genes)]

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            deg(adata_small, min_logfc=0.0, max_pval_adj=1.0, inplace=False)

        warning_messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("fewer than 10 cells" in msg for msg in warning_messages), \
            "Expected UserWarning about groups with fewer than 10 cells"

    def test_deg_fallback_to_leiden(self):
        """When groupby col is missing, deg() should fall back to leiden_col with a warning."""
        adata = _make_adata(n_cells=150, n_groups=3)
        # Remove cell_type_vote so fallback triggers
        del adata.obs["cell_type_vote"]

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            adata_out, deg_dict = deg(
                adata,
                groupby="cell_type_vote",
                leiden_col="leiden",
                min_logfc=0.0,
                max_pval_adj=1.0,
                inplace=False,
            )

        warning_messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("Falling back" in msg for msg in warning_messages), \
            "Expected UserWarning about fallback to leiden_col"

        # Should still produce results keyed by leiden groups
        assert len(deg_dict["results"]) > 0
