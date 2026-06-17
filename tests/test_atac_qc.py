"""
tests/test_atac_qc.py

Unit tests for pipeline.modules.multiome.atac_qc.atac_qc()

Fixture conventions
-------------------
- Synthetic peak count matrices only — no real data loaded in tests
- n_cells=200, n_peaks=500 for most fixtures (fast, CI-safe)
- Fixtures include CellRanger-ARC-style obs columns to mimic real GSE194122 data
- All fixtures use raw integer counts in .X (sparse CSR)
- filter_cells=False is the default; tests that need filtering pass it explicitly
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import anndata as ad

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.modules.scripts.multiome.atac_qc import atac_qc


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_atac_adata(
    n_cells: int = 200,
    n_peaks: int = 500,
    seed: int = 42,
    with_cellranger_obs: bool = True,
    with_cell_type: bool = True,
    sparse: bool = True,
) -> ad.AnnData:
    """
    Build a minimal ATAC AnnData mimicking mdata["atac"] output of qc.py.

    Uses Gaussian clusters to guarantee non-trivial peak distributions.
    CellRanger-ARC obs columns added by default to mimic real data.
    """
    rng = np.random.default_rng(seed)

    # Sparse binary-like integer peak matrix (0/1 with some 2s)
    # ATAC data is very sparse (~1-10% non-zero)
    X = rng.integers(0, 3, size=(n_cells, n_peaks)).astype(float)
    X[X > 1] = 1.0   # mostly binary
    # Zero out ~90% entries to mimic ATAC sparsity
    mask = rng.random(size=X.shape) > 0.10
    X[mask] = 0.0
    if sparse:
        X = sp.csr_matrix(X)

    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"chr1:{i*500}-{i*500+499}" for i in range(n_peaks)]
    adata.var["feature_types"] = "Peaks"

    if with_cell_type:
        cell_types = rng.choice(
            ["CD4 T", "CD8 T", "B cell", "NK", "Monocyte"],
            size=n_cells,
        )
        adata.obs["cell_type"] = cell_types

    if with_cellranger_obs:
        # Mimic CellRanger-ARC pre-computed metrics
        adata.obs["ATAC_nCount_peaks"]        = rng.integers(1500, 50000, size=n_cells).astype(float)
        adata.obs["ATAC_atac_fragments"]      = rng.integers(2000, 60000, size=n_cells).astype(float)
        adata.obs["ATAC_reads_in_peaks_frac"] = rng.uniform(0.4, 0.95, size=n_cells)
        adata.obs["ATAC_blacklist_fraction"]  = rng.uniform(0.0, 0.02, size=n_cells)
        adata.obs["ATAC_nucleosome_signal"]   = rng.uniform(0.5, 1.8, size=n_cells)
        adata.obs["batch"]                    = rng.choice(["s1d1", "s1d2"], size=n_cells)

    return adata


# ---------------------------------------------------------------------------
# 1. Return contract
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self):
        adata = _make_atac_adata()
        result = atac_qc(adata)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_atac_adata()
        adata_out, metrics = atac_qc(adata)
        assert isinstance(adata_out, ad.AnnData)
        assert isinstance(metrics, dict)

    def test_obs_names_preserved(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert list(adata_out.obs_names) == list(adata.obs_names)

    def test_var_names_subset_of_input(self):
        """Feature filter may reduce n_peaks but var_names must be a subset."""
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert set(adata_out.var_names).issubset(set(adata.var_names))


# ---------------------------------------------------------------------------
# 2. QC metrics computed correctly
# ---------------------------------------------------------------------------

class TestAtacQcMetrics:

    def test_n_peaks_by_counts_added(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert "n_peaks_by_counts" in adata_out.obs.columns

    def test_n_peaks_by_counts_non_negative(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert (adata_out.obs["n_peaks_by_counts"] >= 0).all()

    def test_total_peak_counts_added(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert "total_peak_counts" in adata_out.obs.columns

    def test_cellranger_cols_renamed(self):
        """CellRanger-ARC columns should be renamed to OmicSage names."""
        adata = _make_atac_adata(with_cellranger_obs=True)
        adata_out, _ = atac_qc(adata)
        assert "nucleosome_signal"   in adata_out.obs.columns
        assert "reads_in_peaks_frac" in adata_out.obs.columns
        assert "blacklist_fraction"  in adata_out.obs.columns
        assert "total_fragment_counts" in adata_out.obs.columns

    def test_atac_qc_pass_added(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert "atac_qc_pass" in adata_out.obs.columns

    def test_atac_qc_pass_is_boolean(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert adata_out.obs["atac_qc_pass"].dtype == bool or \
               adata_out.obs["atac_qc_pass"].isin([True, False]).all()

    def test_without_cellranger_obs_still_runs(self):
        """Module must work even without CellRanger-ARC obs columns."""
        adata = _make_atac_adata(with_cellranger_obs=False)
        adata_out, metrics = atac_qc(adata)
        assert "n_peaks_by_counts" in adata_out.obs.columns
        assert metrics["cellranger_metrics_available"] is False


# ---------------------------------------------------------------------------
# 3. filter_cells=True removes low-quality cells
# ---------------------------------------------------------------------------

class TestAtacQcFilter:

    def test_filter_cells_false_no_cells_removed(self):
        """Default filter_cells=False must not remove any cells."""
        adata = _make_atac_adata()
        adata_out, metrics = atac_qc(adata, filter_cells=False)
        assert adata_out.n_obs == adata.n_obs

    def test_filter_cells_true_removes_failing_cells(self):
        """filter_cells=True must remove cells that fail QC thresholds."""
        adata = _make_atac_adata(n_cells=200)
        # Force some cells to fail by setting very low peak counts
        if sp.issparse(adata.X):
            X = adata.X.toarray()
        else:
            X = adata.X.copy()
        # First 20 cells get near-zero counts → will fail min_peak_counts
        X[:20, :] = 0.0
        adata.X = sp.csr_matrix(X)
        if "ATAC_nCount_peaks" in adata.obs.columns:
            adata.obs.loc[adata.obs_names[:20], "ATAC_nCount_peaks"] = 10.0

        adata_out, metrics = atac_qc(adata, filter_cells=True)
        assert adata_out.n_obs < adata.n_obs
        assert metrics["filter_cells_applied"] is True

    def test_filter_cells_true_all_pass(self):
        """If all cells pass QC, filter_cells=True must not remove any."""
        adata = _make_atac_adata()
        # Ensure all cells pass by using very permissive thresholds
        adata_out, _ = atac_qc(
            adata,
            min_peaks=0, max_peaks=10_000_000,
            min_peak_counts=0, max_peak_counts=10_000_000,
            max_nucleosome_signal=100.0,
            filter_cells=True,
        )
        assert adata_out.n_obs == adata.n_obs


# ---------------------------------------------------------------------------
# 4. Provenance
# ---------------------------------------------------------------------------

class TestAtacQcProvenance:

    def test_uns_key_present(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert "omicsage_atac_qc" in adata_out.uns

    def test_uns_required_fields(self):
        required = {
            "omicsage_module",
            "omicsage_version",
            "timestamp",
            "n_cells_before",
            "n_cells_after",
            "n_peaks_before",
            "n_peaks_after",
            "fragment_file_available",
            "cellranger_metrics_used",
            "filter_cells",
            "thresholds",
            "note",
        }
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        prov = adata_out.uns["omicsage_atac_qc"]
        missing = required - set(prov.keys())
        assert not missing, f"Provenance uns missing keys: {missing}"

    def test_uns_module_name(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert (
            adata_out.uns["omicsage_atac_qc"]["omicsage_module"]
            == "pipeline.modules.multiome.atac_qc"
        )

    def test_uns_timestamp_is_iso_string(self):
        from datetime import datetime
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        ts = adata_out.uns["omicsage_atac_qc"]["timestamp"]
        datetime.fromisoformat(ts)

    def test_uns_fragment_file_false(self):
        """fragment_file_available must always be False for this dataset."""
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert adata_out.uns["omicsage_atac_qc"]["fragment_file_available"] is False

    def test_uns_cell_counts_consistent(self):
        adata = _make_atac_adata()
        adata_out, metrics = atac_qc(adata)
        prov = adata_out.uns["omicsage_atac_qc"]
        assert prov["n_cells_before"] == adata.n_obs
        assert prov["n_cells_after"]  == metrics["n_cells_after"]


# ---------------------------------------------------------------------------
# 5. Ground truth label preservation
# ---------------------------------------------------------------------------

class TestGroundTruthPreservation:

    def test_cell_type_preserved_as_groundtruth(self):
        adata = _make_atac_adata(with_cell_type=True)
        original_types = adata.obs["cell_type"].values.copy()
        adata_out, _ = atac_qc(adata)
        assert "cell_type_groundtruth" in adata_out.obs.columns
        np.testing.assert_array_equal(
            adata_out.obs["cell_type_groundtruth"].values, original_types
        )

    def test_no_cell_type_col_no_error(self):
        """If cell_type is absent, atac_qc must not raise."""
        adata = _make_atac_adata(with_cell_type=False)
        adata_out, _ = atac_qc(adata)
        assert "cell_type_groundtruth" not in adata_out.obs.columns


# ---------------------------------------------------------------------------
# 6. Fragment file fallback
# ---------------------------------------------------------------------------

class TestFragmentFallback:

    def test_no_cellranger_obs_nucleosome_filter_skipped(self):
        """Without CellRanger obs, nucleosome filter must not be applied."""
        adata = _make_atac_adata(with_cellranger_obs=False)
        adata_out, metrics = atac_qc(adata)
        assert metrics["nucleosome_filter_applied"] is False

    def test_fragment_file_available_always_false(self):
        adata = _make_atac_adata(with_cellranger_obs=True)
        _, metrics = atac_qc(adata)
        assert metrics["fragment_file_available"] is False

    def test_cellranger_metrics_available_true_when_cols_present(self):
        adata = _make_atac_adata(with_cellranger_obs=True)
        _, metrics = atac_qc(adata)
        assert metrics["cellranger_metrics_available"] is True

    def test_cellranger_metrics_available_false_without_cols(self):
        adata = _make_atac_adata(with_cellranger_obs=False)
        _, metrics = atac_qc(adata)
        assert metrics["cellranger_metrics_available"] is False


# ---------------------------------------------------------------------------
# 7. Layers contract
# ---------------------------------------------------------------------------

class TestLayers:

    def test_counts_layer_created(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        assert "counts" in adata_out.layers

    def test_counts_layer_contains_integers(self):
        """Raw counts layer must contain non-negative values."""
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata)
        counts = adata_out.layers["counts"]
        if sp.issparse(counts):
            counts = counts.toarray()
        assert np.all(counts >= 0)

    def test_existing_counts_layer_not_overwritten(self):
        adata = _make_atac_adata()
        sentinel = np.ones(adata.shape) * 999.0
        adata.layers["counts"] = sp.csr_matrix(sentinel)
        adata_out, _ = atac_qc(adata)
        counts = adata_out.layers["counts"]
        if sp.issparse(counts):
            counts = counts.toarray()
        # Feature filter may reduce n_peaks; every surviving value must still be 999
        assert counts.shape[0] == adata_out.n_obs
        assert counts.shape[1] == adata_out.n_vars
        np.testing.assert_array_equal(counts, np.ones(counts.shape) * 999.0)


# ---------------------------------------------------------------------------
# 8. Scrublet doublet detection
# ---------------------------------------------------------------------------

class TestScrublet:

    def test_doublet_score_column_added(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata, run_scrublet=True)
        assert "atac_doublet_score" in adata_out.obs.columns

    def test_doublet_flag_column_added(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata, run_scrublet=True)
        assert "atac_predicted_doublet" in adata_out.obs.columns

    def test_scrublet_false_adds_nan_score(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata, run_scrublet=False)
        assert "atac_doublet_score" in adata_out.obs.columns
        assert adata_out.obs["atac_doublet_score"].isna().all()

    def test_scrublet_false_adds_false_flag(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata, run_scrublet=False)
        assert not adata_out.obs["atac_predicted_doublet"].any()

    def test_cells_never_removed_by_doublet_flag(self):
        """Scrublet flags cells but must never remove them."""
        adata = _make_atac_adata()
        adata_out, _ = atac_qc(adata, run_scrublet=True, filter_cells=False)
        assert adata_out.n_obs == adata.n_obs


# ---------------------------------------------------------------------------
# 9. Metrics dict contract
# ---------------------------------------------------------------------------

class TestMetricsDict:

    def test_required_keys_present(self):
        required = {
            "n_cells_before",
            "n_cells_after",
            "n_cells_removed_qc",
            "n_peaks_before",
            "n_peaks_after",
            "n_peaks_removed_feature_filter",
            "n_qc_pass",
            "n_qc_fail",
            "filter_cells_applied",
            "cellranger_metrics_available",
            "fragment_file_available",
            "nucleosome_filter_applied",
            "n_fail_low_peaks",
            "n_fail_high_peaks",
            "n_fail_low_counts",
            "n_fail_high_counts",
            "n_fail_high_nucleosome",
            "median_n_peaks_by_counts",
            "median_total_peak_counts",
            "thresholds",
        }
        adata = _make_atac_adata()
        _, metrics = atac_qc(adata)
        missing = required - set(metrics.keys())
        assert not missing, f"Metrics dict missing keys: {missing}"

    def test_n_cells_before_correct(self):
        adata = _make_atac_adata(n_cells=150)
        _, metrics = atac_qc(adata)
        assert metrics["n_cells_before"] == 150

    def test_n_qc_pass_plus_fail_equals_before(self):
        adata = _make_atac_adata()
        _, metrics = atac_qc(adata)
        assert metrics["n_qc_pass"] + metrics["n_qc_fail"] == metrics["n_cells_before"]

    def test_filter_cells_applied_false_by_default(self):
        adata = _make_atac_adata()
        _, metrics = atac_qc(adata)
        assert metrics["filter_cells_applied"] is False

    def test_thresholds_recorded(self):
        adata = _make_atac_adata()
        _, metrics = atac_qc(adata, min_peaks=800, min_peak_counts=2000)
        assert metrics["thresholds"]["min_peaks"] == 800
        assert metrics["thresholds"]["min_peak_counts"] == 2000


# ---------------------------------------------------------------------------
# 10. Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_rejects_non_anndata(self):
        with pytest.raises(TypeError, match="AnnData"):
            atac_qc("not_an_anndata")

    def test_rejects_empty_cells(self):
        adata = _make_atac_adata()
        empty = adata[:0].copy()
        with pytest.raises(ValueError, match="0 cells"):
            atac_qc(empty)

    def test_rejects_empty_peaks(self):
        adata = _make_atac_adata()
        empty = adata[:, :0].copy()
        with pytest.raises(ValueError, match="0 peaks"):
            atac_qc(empty)


# ---------------------------------------------------------------------------
# 11. inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:

    def test_default_not_inplace(self):
        adata = _make_atac_adata()
        original_n_obs = adata.n_obs
        original_obs_cols = set(adata.obs.columns)
        atac_qc(adata)
        # Input must not have gained new obs columns
        assert adata.n_obs == original_n_obs
        assert set(adata.obs.columns) == original_obs_cols

    def test_inplace_true_mutates_input(self):
        adata = _make_atac_adata()
        atac_qc(adata, inplace=True)
        assert "n_peaks_by_counts" in adata.obs.columns
        assert "omicsage_atac_qc"  in adata.uns
