"""
test_pseudobulk_deg.py — Tests for pipeline/modules/downstream/pseudobulk_deg.py

Run with:
    cd ~/OmicSage
    conda activate omicsage
    python -m pytest tests/test_pseudobulk_deg.py -v

Expected: 14 passed
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from pipeline.modules.downstream.pseudobulk_deg import pseudobulk_deg


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_adata(
    n_cells_per_group: int = 120,
    n_genes: int = 80,
    n_cell_types: int = 3,
    n_donors: int = 6,
    seed: int = 42,
) -> AnnData:
    """
    Synthetic AnnData suitable for pseudobulk DEG testing.

    Design
    ------
    - n_cell_types groups × n_donors donors = n_cell_types × n_donors pseudo-samples
    - Each cell type has a distinct gene block elevated ~5× above background
    - layers['counts'] holds raw integer counts (required by pseudobulk_deg)
    - layers['logcounts'] is also present (pipeline convention)
    - obs['cell_type_vote'] — cell-type labels
    - obs['batch']          — donor labels  (d0 … d{n_donors-1})

    Cell count per (cell_type, donor): n_cells_per_group // n_donors
    All (cell_type, donor) combos are present — no missing combinations.
    """
    rng = np.random.default_rng(seed)

    cell_types = [f"CellType_{i}" for i in range(n_cell_types)]
    donors     = [f"d{j}" for j in range(n_donors)]
    cells_per_combo = max(n_cells_per_group // n_donors, 1)

    group_list  = []
    donor_list  = []
    count_blocks = []

    for ct_idx, ct in enumerate(cell_types):
        for donor in donors:
            # Background
            block = rng.negative_binomial(1, 0.2, size=(cells_per_combo, n_genes))
            # Signal: each cell type elevated in its own gene block
            start = ct_idx * (n_genes // n_cell_types)
            end   = start  + (n_genes // n_cell_types)
            block[:, start:end] += rng.negative_binomial(
                5, 0.3, size=(cells_per_combo, end - start)
            )
            count_blocks.append(block)
            group_list.extend([ct]    * cells_per_combo)
            donor_list.extend([donor] * cells_per_combo)

    X_counts = np.vstack(count_blocks).astype(np.float64)
    totals   = X_counts.sum(axis=1, keepdims=True) + 1e-9
    logcounts = np.log1p(X_counts / totals * 1e4)

    adata = AnnData(X=logcounts)
    adata.layers["counts"]    = X_counts.copy()
    adata.layers["logcounts"] = logcounts.copy()
    adata.obs["cell_type_vote"] = pd.Categorical(group_list)
    adata.obs["batch"]          = pd.Categorical(donor_list)
    adata.var_names = [f"gene_{i:04d}" for i in range(n_genes)]

    return adata


@pytest.fixture(scope="module")
def adata_base():
    return _make_adata()


@pytest.fixture(scope="module")
def pb_results(adata_base):
    """Run pseudobulk_deg() once and share the output across tests."""
    adata_out, pb_dict = pseudobulk_deg(
        adata_base,
        groupby="cell_type_vote",
        donor_key="batch",
        min_cells=1,       # small fixture — relax cell threshold
        min_samples=2,     # minimum 2 donors for Wald test
        min_logfc=0.0,     # no filtering — keeps all genes for structural tests
        max_pval_adj=1.0,
        inplace=False,
    )
    return adata_out, pb_dict


# ---------------------------------------------------------------------------
# Tests: return types
# ---------------------------------------------------------------------------

class TestPseudobulkReturnsCorrectTypes:
    def test_returns_tuple_of_two(self, adata_base):
        result = pseudobulk_deg(
            adata_base,
            min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=1.0, inplace=False,
        )
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self, adata_base):
        adata_out, pb_dict = pseudobulk_deg(
            adata_base,
            min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=1.0, inplace=False,
        )
        assert isinstance(adata_out, AnnData)
        assert isinstance(pb_dict, dict)


# ---------------------------------------------------------------------------
# Tests: provenance
# ---------------------------------------------------------------------------

class TestPseudobulkProvenance:
    def test_uns_key_present(self, pb_results):
        adata_out, _ = pb_results
        assert "omicsage_pseudobulk_deg" in adata_out.uns

    def test_provenance_required_keys(self, pb_results):
        adata_out, _ = pb_results
        prov = adata_out.uns["omicsage_pseudobulk_deg"]
        required = {
            "groupby", "donor_key", "counts_layer", "method",
            "min_cells", "min_samples", "min_logfc", "max_pval_adj",
            "n_groups", "n_skipped", "results_keys", "skipped_keys",
            "omicsage_module", "timestamp",
        }
        missing = required - set(prov.keys())
        assert not missing, f"Missing provenance keys: {missing}"

    def test_provenance_module_name(self, pb_results):
        adata_out, _ = pb_results
        assert adata_out.uns["omicsage_pseudobulk_deg"]["omicsage_module"] == "pseudobulk_deg"


# ---------------------------------------------------------------------------
# Tests: output schema — must match deg.py deg_dict exactly
# ---------------------------------------------------------------------------

class TestPseudobulkOutputSchema:
    def test_pb_dict_top_level_keys(self, pb_results):
        """pb_dict must have the same top-level keys as deg_dict."""
        _, pb_dict = pb_results
        required = {"results", "summary_df", "provenance", "pairwise", "skipped"}
        missing  = required - set(pb_dict.keys())
        assert not missing, f"pb_dict missing keys: {missing}"

    def test_result_dataframe_columns(self, pb_results):
        """Each group DataFrame must have: gene, score, pval, logfc, pval_adj."""
        _, pb_dict = pb_results
        required = {"gene", "score", "pval", "logfc", "pval_adj"}
        for group, df in pb_dict["results"].items():
            assert isinstance(df, pd.DataFrame)
            missing = required - set(df.columns)
            assert not missing, f"results['{group}'] missing columns: {missing}"

    def test_summary_df_columns(self, pb_results):
        """summary_df must have: group, rank, gene, logfc, pval_adj."""
        _, pb_dict = pb_results
        summary = pb_dict["summary_df"]
        assert isinstance(summary, pd.DataFrame)
        required = {"group", "rank", "gene", "logfc", "pval_adj"}
        missing  = required - set(summary.columns)
        assert not missing, f"summary_df missing columns: {missing}"

    def test_pairwise_is_empty_dict(self, pb_results):
        """pairwise must be an empty dict (not implemented for pseudobulk)."""
        _, pb_dict = pb_results
        assert pb_dict["pairwise"] == {}

    def test_all_groups_have_results_or_skipped(self, adata_base, pb_results):
        """Every cell type must appear in either results or skipped."""
        _, pb_dict = pb_results
        expected = set(adata_base.obs["cell_type_vote"].astype(str).unique())
        covered  = set(pb_dict["results"].keys()) | set(pb_dict["skipped"].keys())
        assert expected == covered, f"Cell types unaccounted for: {expected - covered}"


# ---------------------------------------------------------------------------
# Tests: statistical properties
# ---------------------------------------------------------------------------

class TestPseudobulkStatisticalProperties:
    def test_pval_adj_in_unit_interval(self, pb_results):
        """All pval_adj values must lie in [0, 1]."""
        _, pb_dict = pb_results
        for group, df in pb_dict["results"].items():
            if df.empty:
                continue
            assert df["pval_adj"].between(0.0, 1.0).all(), \
                f"pval_adj out of [0,1] in group '{group}'"

    def test_pval_in_unit_interval(self, pb_results):
        """All pval values must lie in [0, 1]."""
        _, pb_dict = pb_results
        for group, df in pb_dict["results"].items():
            if df.empty:
                continue
            assert df["pval"].between(0.0, 1.0).all(), \
                f"pval out of [0,1] in group '{group}'"


# ---------------------------------------------------------------------------
# Tests: thresholds
# ---------------------------------------------------------------------------

class TestPseudobulkThresholds:
    def test_logfc_threshold_reduces_results(self, adata_base):
        """Strict min_logfc should return ≤ genes than loose threshold."""
        _, dict_loose = pseudobulk_deg(
            adata_base, min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=1.0, inplace=False,
        )
        _, dict_strict = pseudobulk_deg(
            adata_base, min_cells=1, min_samples=2,
            min_logfc=2.0, max_pval_adj=1.0, inplace=False,
        )
        loose  = sum(len(df) for df in dict_loose["results"].values())
        strict = sum(len(df) for df in dict_strict["results"].values())
        assert strict <= loose

    def test_pval_threshold_reduces_results(self, adata_base):
        """Strict max_pval_adj should return ≤ genes than loose threshold."""
        _, dict_loose = pseudobulk_deg(
            adata_base, min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=1.0, inplace=False,
        )
        _, dict_strict = pseudobulk_deg(
            adata_base, min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=0.01, inplace=False,
        )
        loose  = sum(len(df) for df in dict_loose["results"].values())
        strict = sum(len(df) for df in dict_strict["results"].values())
        assert strict <= loose


# ---------------------------------------------------------------------------
# Tests: edge cases
# ---------------------------------------------------------------------------

class TestPseudobulkEdgeCases:
    def test_inplace_false_does_not_modify_input(self, adata_base):
        """inplace=False must not add keys to the caller's adata.uns."""
        original_uns_keys = set(adata_base.uns.keys())
        pseudobulk_deg(
            adata_base, min_cells=1, min_samples=2,
            min_logfc=0.0, max_pval_adj=1.0, inplace=False,
        )
        assert set(adata_base.uns.keys()) == original_uns_keys

    def test_cell_type_with_too_few_donors_is_skipped(self):
        """A cell type with < min_samples donors must be in skipped, not results."""
        rng   = np.random.default_rng(7)
        n_g   = 40

        # Design: three cell types.
        # CellType_A (4 donors) and CellType_B (4 donors) both pass min_samples=3.
        # When CellType_A is the target, its rest group = CellType_B (4) + CellType_C (2)
        # = 6 pseudo-samples ≥ 3, so CellType_A gets results.
        # CellType_C (2 donors) fails the target-group check (2 < 3) → skipped.
        blocks, group_l, donor_l = [], [], []
        # CellType_A: 4 donors  (target=4 ≥ 3, rest=CellType_B has 4 ≥ 3 → results)
        for donor in ["d0", "d1", "d2", "d3"]:
            block = rng.negative_binomial(1, 0.2, size=(30, n_g)).astype(float)
            blocks.append(block)
            group_l.extend(["CellType_A"] * 30)
            donor_l.extend([donor] * 30)
        # CellType_B: 4 donors  (rest group for CellType_A — must be ≥ 3)
        for donor in ["d0", "d1", "d2", "d3"]:
            block = rng.negative_binomial(1, 0.2, size=(30, n_g)).astype(float)
            blocks.append(block)
            group_l.extend(["CellType_B"] * 30)
            donor_l.extend([donor] * 30)
        # CellType_C: 2 donors  (< min_samples=3 → skipped)
        for donor in ["d0", "d1"]:
            block = rng.negative_binomial(1, 0.2, size=(30, n_g)).astype(float)
            blocks.append(block)
            group_l.extend(["CellType_C"] * 30)
            donor_l.extend([donor] * 30)

        X = np.vstack(blocks)
        adata_small = AnnData(X=np.log1p(X / (X.sum(1, keepdims=True) + 1e-9) * 1e4))
        adata_small.layers["counts"]    = X.copy()
        adata_small.layers["logcounts"] = adata_small.X.copy()
        adata_small.obs["cell_type_vote"] = pd.Categorical(group_l)
        adata_small.obs["batch"]          = pd.Categorical(donor_l)
        adata_small.var_names = [f"g{i}" for i in range(n_g)]

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            _, pb_dict = pseudobulk_deg(
                adata_small, min_cells=1, min_samples=3,
                min_logfc=0.0, max_pval_adj=1.0, inplace=False,
            )

        assert "CellType_C" in pb_dict["skipped"], \
            "CellType_C should be in skipped (too few donors)"
        assert "CellType_A" in pb_dict["results"], \
            "CellType_A should have results"

        # Warning must have been issued
        msgs = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
        assert any("CellType_C" in m for m in msgs), \
            "Expected UserWarning naming the skipped cell type"

    def test_missing_counts_layer_raises(self, adata_base):
        """Passing a non-existent counts_layer must raise ValueError."""
        with pytest.raises(ValueError, match="layers\\["):
            pseudobulk_deg(
                adata_base,
                counts_layer="does_not_exist",
                min_cells=1, min_samples=2, inplace=False,
            )

    def test_missing_groupby_column_raises(self, adata_base):
        """Passing a non-existent groupby column must raise ValueError."""
        with pytest.raises(ValueError, match="obs column"):
            pseudobulk_deg(
                adata_base,
                groupby="nonexistent_col",
                min_cells=1, min_samples=2, inplace=False,
            )
