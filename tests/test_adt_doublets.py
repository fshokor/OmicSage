"""
tests/test_adt_doublets.py

Tests for pipeline/modules/cite/adt_doublets.py

Test groups
-----------
T01  Fixture sanity
T02  Marker resolution (prefix matching, case-insensitive, missing markers)
T03  Doublet scoring — core arithmetic
T04  Predicted doublet boolean flag
T05  No evaluable pairs → neutral output
T06  filter_doublets=False — cells not removed
T07  filter_doublets=True  — cross-modal filtering
T08  inplace flag
T09  Provenance written to uns
T10  Metrics dict keys and types
T11  Sparse matrix input
T12  Edge cases (all doublets, zero doublets, single pair, single cell)
T13  Threshold sensitivity
T14  Custom marker pairs
T15  Missing adt_clr layer raises KeyError
"""

from __future__ import annotations

import copy

import numpy as np
import pytest
import anndata as ad
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Lazy import so collection works even before install
# ---------------------------------------------------------------------------
def import_module():
    import importlib.util, sys, pathlib
    path = pathlib.Path(__file__).parent.parent / "pipeline/modules/cite/adt_doublets.py"
    spec = importlib.util.spec_from_file_location("adt_doublets", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

class MockMuData:
    """Minimal MuData stand-in that mirrors the parts we use."""

    def __init__(self, adt: ad.AnnData, rna: ad.AnnData | None = None):
        self.mod = {"adt": adt}
        if rna is not None:
            self.mod["rna"] = rna

    def __getitem__(self, key):
        return self.mod[key]

    def __contains__(self, key):
        return key in self.mod

    @property
    def obs_names(self):
        return self.mod["adt"].obs_names


# Marker names that match DEFAULT_MARKER_PAIRS via prefix
ADT_VARS = [
    "CD3-TotalSeqB",   # matches "CD3"
    "CD19-1",          # matches "CD19"
    "CD14-1",          # matches "CD14"
    "CD56-1",
    "CD4-1",
    "CD8a-1",
]

N_CELLS = 50
N_VARS = len(ADT_VARS)
RNG = np.random.default_rng(42)


def _make_clr(n_cells=N_CELLS, rng=RNG) -> np.ndarray:
    """CLR-like float array in roughly [-2, 5] range."""
    return rng.normal(loc=0.5, scale=1.5, size=(n_cells, N_VARS)).astype(np.float32)


def _make_adt(clr: np.ndarray | None = None, n_cells: int = N_CELLS) -> ad.AnnData:
    if clr is None:
        clr = _make_clr(n_cells=n_cells)
    barcodes = [f"cell_{i}" for i in range(n_cells)]
    adata = ad.AnnData(
        X=clr.copy(),
        obs={"barcode": barcodes},
        var={"marker": ADT_VARS},
    )
    adata.obs_names = barcodes
    adata.var_names = ADT_VARS
    adata.layers["adt_clr"] = clr.copy()
    return adata


def _make_rna(n_cells: int = N_CELLS) -> ad.AnnData:
    barcodes = [f"cell_{i}" for i in range(n_cells)]
    X = RNG.integers(0, 100, size=(n_cells, 20)).astype(np.float32)
    rna = ad.AnnData(X=X)
    rna.obs_names = barcodes
    rna.var_names = [f"gene_{i}" for i in range(20)]
    return rna


def _make_mdata(clr=None, n_cells=N_CELLS, with_rna=True):
    adt = _make_adt(clr, n_cells=n_cells)
    rna = _make_rna(n_cells) if with_rna else None
    return MockMuData(adt, rna)


# ---------------------------------------------------------------------------
# T01 — Fixture sanity
# ---------------------------------------------------------------------------

class TestFixtureSanity:
    def test_adt_has_clr_layer(self):
        adt = _make_adt()
        assert "adt_clr" in adt.layers

    def test_adt_var_names(self):
        adt = _make_adt()
        assert "CD3-TotalSeqB" in adt.var_names

    def test_mdata_has_rna(self):
        mdata = _make_mdata()
        assert "rna" in mdata.mod

    def test_mdata_barcodes_match(self):
        mdata = _make_mdata()
        assert list(mdata["adt"].obs_names) == list(mdata["rna"].obs_names)


# ---------------------------------------------------------------------------
# T02 — Marker resolution
# ---------------------------------------------------------------------------

class TestMarkerResolution:
    def setup_method(self):
        self.m = import_module()

    def test_prefix_match_cd3(self):
        # _resolve_marker is a closure — verify prefix matching via the public API:
        # CD3-TotalSeqB in var_names should be matched by marker name "CD3"
        clr = _make_clr()
        mdata = _make_mdata(clr)
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19")], threshold=2.5
        )
        assert ("CD3", "CD19") in [tuple(p) for p in metrics["pairs_evaluated"]]

    def test_case_insensitive_match(self):
        clr = _make_clr()
        mdata = _make_mdata(clr)
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("cd3", "cd19")], threshold=2.5
        )
        assert len(metrics["pairs_evaluated"]) == 1

    def test_missing_marker_skipped(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CXCR5", "CCR7")], threshold=2.5
        )
        assert ("CXCR5", "CCR7") in [tuple(p) for p in metrics["pairs_skipped"]]
        assert len(metrics["pairs_evaluated"]) == 0

    def test_partial_pair_skipped(self):
        """One marker present, one missing → pair skipped."""
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "MISSING_MARKER")], threshold=2.5
        )
        assert len(metrics["pairs_skipped"]) == 1
        assert len(metrics["pairs_evaluated"]) == 0

    def test_default_pairs_evaluated(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(mdata, threshold=2.5)
        # Both default pairs should resolve with our fixture vars
        assert len(metrics["pairs_evaluated"]) == 2

    def test_extra_var_suffix_still_matches(self):
        """CD3-TotalSeqB in var_names should be found by marker name "CD3"."""
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19")], threshold=2.5
        )
        assert ("CD3", "CD19") in [tuple(p) for p in metrics["pairs_evaluated"]]


# ---------------------------------------------------------------------------
# T03 — Doublet scoring arithmetic
# ---------------------------------------------------------------------------

class TestDoubletScoring:
    def setup_method(self):
        self.m = import_module()

    def _score_for(self, clr, **kw):
        mdata = _make_mdata(clr)
        mdata2, _ = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19"), ("CD3", "CD14")],
            threshold=2.5, **kw
        )
        return mdata2["adt"].obs["adt_doublet_score"].values

    def test_score_range_0_to_1(self):
        clr = _make_clr()
        scores = self._score_for(clr)
        assert np.all(scores >= 0.0)
        assert np.all(scores <= 1.0)

    def test_score_zero_when_no_coexpression(self):
        """Set all values below threshold → all scores == 0."""
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        scores = self._score_for(clr)
        np.testing.assert_array_equal(scores, 0.0)

    def test_score_one_when_all_pairs_coexpress(self):
        """Set CD3, CD19, CD14 all above threshold → score == 1.0 (both pairs fire)."""
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        # CD3=idx0, CD19=idx1, CD14=idx2
        clr[:, 0] = 3.0  # CD3
        clr[:, 1] = 3.0  # CD19
        clr[:, 2] = 3.0  # CD14
        scores = self._score_for(clr)
        np.testing.assert_array_almost_equal(scores, 1.0)

    def test_score_half_when_one_of_two_pairs_fires(self):
        """Only CD3+CD19 fires, CD3+CD14 does not → score == 0.5."""
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[:, 0] = 3.0  # CD3
        clr[:, 1] = 3.0  # CD19  — pair 1 fires
        # CD14 stays at 0.0 — pair 2 does not fire
        scores = self._score_for(clr)
        np.testing.assert_array_almost_equal(scores, 0.5)

    def test_score_dtype_float64(self):
        scores = self._score_for(_make_clr())
        assert scores.dtype == np.float64

    def test_per_cell_independence(self):
        """First cell is doublet, rest are not."""
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[0, 0] = 3.0  # CD3
        clr[0, 1] = 3.0  # CD19
        scores = self._score_for(clr)
        assert scores[0] > 0.0
        np.testing.assert_array_equal(scores[1:], 0.0)


# ---------------------------------------------------------------------------
# T04 — Predicted doublet boolean flag
# ---------------------------------------------------------------------------

class TestPredictedDoublet:
    def setup_method(self):
        self.m = import_module()

    def test_predicted_doublet_dtype_bool(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        col = mdata2["adt"].obs["adt_predicted_doublet"]
        assert col.dtype == bool or col.dtype == np.bool_

    def test_predicted_doublet_true_when_score_gt_zero(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[5, 0] = 3.0; clr[5, 1] = 3.0  # cell 5 is a doublet
        mdata = _make_mdata(clr)
        mdata2, _ = self.m.detect_adt_doublets(mdata, threshold=2.5)
        flags = mdata2["adt"].obs["adt_predicted_doublet"].values
        assert flags[5] == True
        # all others should be False
        other_idx = [i for i in range(N_CELLS) if i != 5]
        assert not np.any(flags[other_idx])

    def test_predicted_false_when_score_zero(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        mdata = _make_mdata(clr)
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        flags = mdata2["adt"].obs["adt_predicted_doublet"].values
        assert not np.any(flags)

    def test_predicted_columns_present(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        assert "adt_doublet_score" in mdata2["adt"].obs.columns
        assert "adt_predicted_doublet" in mdata2["adt"].obs.columns


# ---------------------------------------------------------------------------
# T05 — No evaluable pairs → neutral output
# ---------------------------------------------------------------------------

class TestNoEvaluablePairs:
    def setup_method(self):
        self.m = import_module()

    def test_score_zero_when_no_pairs(self):
        mdata = _make_mdata()
        mdata2, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("MISSING_A", "MISSING_B")]
        )
        scores = mdata2["adt"].obs["adt_doublet_score"].values
        np.testing.assert_array_equal(scores, 0.0)

    def test_no_doublets_flagged_when_no_pairs(self):
        mdata = _make_mdata()
        mdata2, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("MISSING_A", "MISSING_B")]
        )
        flags = mdata2["adt"].obs["adt_predicted_doublet"].values
        assert not np.any(flags)

    def test_metrics_reflect_no_pairs(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("MISSING_A", "MISSING_B")]
        )
        assert metrics["n_doublets_detected"] == 0
        assert len(metrics["pairs_evaluated"]) == 0
        assert len(metrics["pairs_skipped"]) == 1

    def test_cell_count_unchanged_when_no_pairs(self):
        mdata = _make_mdata()
        mdata2, metrics = self.m.detect_adt_doublets(
            mdata,
            marker_pairs=[("MISSING_A", "MISSING_B")],
            filter_doublets=True,
        )
        assert metrics["n_cells_after"] == N_CELLS


# ---------------------------------------------------------------------------
# T06 — filter_doublets=False
# ---------------------------------------------------------------------------

class TestFilterDoubletsFalse:
    def setup_method(self):
        self.m = import_module()

    def _mdata_with_known_doublets(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[:5, 0] = 3.0; clr[:5, 1] = 3.0  # cells 0-4 are doublets
        return _make_mdata(clr)

    def test_cell_count_unchanged(self):
        mdata = self._mdata_with_known_doublets()
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=False)
        assert mdata2["adt"].n_obs == N_CELLS

    def test_rna_count_unchanged(self):
        mdata = self._mdata_with_known_doublets()
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=False)
        assert mdata2["rna"].n_obs == N_CELLS

    def test_flags_written_even_without_filter(self):
        mdata = self._mdata_with_known_doublets()
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=False)
        n_flagged = mdata2["adt"].obs["adt_predicted_doublet"].sum()
        assert n_flagged == 5

    def test_metrics_n_cells_after_equals_before(self):
        mdata = self._mdata_with_known_doublets()
        _, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=False)
        assert metrics["n_cells_after"] == metrics["n_cells_before"]


# ---------------------------------------------------------------------------
# T07 — filter_doublets=True (cross-modal filtering)
# ---------------------------------------------------------------------------

class TestFilterDoubletsTrue:
    def setup_method(self):
        self.m = import_module()

    def _mdata_with_n_doublets(self, n=5):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[:n, 0] = 3.0
        clr[:n, 1] = 3.0
        return _make_mdata(clr), n

    def test_adt_cells_removed(self):
        mdata, n_d = self._mdata_with_n_doublets(5)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert mdata2["adt"].n_obs == N_CELLS - n_d

    def test_rna_cells_removed(self):
        mdata, n_d = self._mdata_with_n_doublets(5)
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert mdata2["rna"].n_obs == N_CELLS - n_d

    def test_rna_adt_obs_names_in_sync(self):
        mdata, _ = self._mdata_with_n_doublets(5)
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert list(mdata2["adt"].obs_names) == list(mdata2["rna"].obs_names)

    def test_correct_barcodes_removed(self):
        mdata, n_d = self._mdata_with_n_doublets(5)
        mdata2, _ = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        remaining = list(mdata2["adt"].obs_names)
        for i in range(n_d):
            assert f"cell_{i}" not in remaining
        for i in range(n_d, N_CELLS):
            assert f"cell_{i}" in remaining

    def test_metrics_n_cells_after_correct(self):
        mdata, n_d = self._mdata_with_n_doublets(5)
        _, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert metrics["n_cells_after"] == N_CELLS - n_d

    def test_filter_without_rna_modality(self):
        """filter_doublets=True should work even if there's no RNA modality."""
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[:3, 0] = 3.0; clr[:3, 1] = 3.0
        mdata = _make_mdata(clr, with_rna=False)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert mdata2["adt"].n_obs == N_CELLS - 3

    def test_no_doublets_filter_is_noop(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        mdata = _make_mdata(clr)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert mdata2["adt"].n_obs == N_CELLS
        assert metrics["n_cells_after"] == N_CELLS


# ---------------------------------------------------------------------------
# T08 — inplace flag
# ---------------------------------------------------------------------------

class TestInplaceFlag:
    def setup_method(self):
        self.m = import_module()

    def test_inplace_false_does_not_modify_original(self):
        mdata = _make_mdata()
        original_obs_cols = set(mdata["adt"].obs.columns)
        self.m.detect_adt_doublets(mdata, inplace=False)
        assert "adt_doublet_score" not in original_obs_cols
        assert "adt_doublet_score" not in mdata["adt"].obs.columns

    def test_inplace_true_modifies_original(self):
        mdata = _make_mdata()
        self.m.detect_adt_doublets(mdata, inplace=True)
        assert "adt_doublet_score" in mdata["adt"].obs.columns

    def test_inplace_false_returns_copy(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata, inplace=False)
        assert mdata2 is not mdata


# ---------------------------------------------------------------------------
# T09 — Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    def setup_method(self):
        self.m = import_module()

    def test_uns_key_exists(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        assert "omicsage_adt_doublets" in mdata2["adt"].uns

    def test_uns_has_timestamp(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        prov = mdata2["adt"].uns["omicsage_adt_doublets"]
        assert "timestamp" in prov

    def test_uns_has_module_name(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata)
        prov = mdata2["adt"].uns["omicsage_adt_doublets"]
        assert prov["omicsage_module"] == "adt_doublets"

    def test_uns_threshold_matches_input(self):
        mdata = _make_mdata()
        mdata2, _ = self.m.detect_adt_doublets(mdata, threshold=3.0)
        prov = mdata2["adt"].uns["omicsage_adt_doublets"]
        assert prov["threshold"] == 3.0


# ---------------------------------------------------------------------------
# T10 — Metrics dict
# ---------------------------------------------------------------------------

class TestMetrics:
    def setup_method(self):
        self.m = import_module()

    def _run(self, **kw):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(mdata, **kw)
        return metrics

    def test_all_required_keys_present(self):
        metrics = self._run()
        required = {
            "n_cells_before", "n_doublets_detected", "pct_doublets",
            "pairs_evaluated", "pairs_skipped", "n_cells_after",
            "threshold", "filter_doublets",
        }
        assert required.issubset(metrics.keys())

    def test_n_cells_before_correct(self):
        metrics = self._run()
        assert metrics["n_cells_before"] == N_CELLS

    def test_pct_doublets_is_float(self):
        metrics = self._run()
        assert isinstance(metrics["pct_doublets"], float)

    def test_pct_doublets_in_range(self):
        metrics = self._run()
        assert 0.0 <= metrics["pct_doublets"] <= 100.0

    def test_n_doublets_plus_remaining_equals_before_when_filtered(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        clr[:7, 0] = 3.0; clr[:7, 1] = 3.0
        mdata = _make_mdata(clr)
        _, metrics = self.m.detect_adt_doublets(mdata, filter_doublets=True)
        assert metrics["n_doublets_detected"] + metrics["n_cells_after"] == metrics["n_cells_before"]

    def test_pairs_evaluated_is_list(self):
        metrics = self._run()
        assert isinstance(metrics["pairs_evaluated"], list)

    def test_filter_doublets_flag_in_metrics(self):
        m_false = self._run(filter_doublets=False)
        m_true_mdata = _make_mdata()
        _, m_true = self.m.detect_adt_doublets(m_true_mdata, filter_doublets=True)
        assert m_false["filter_doublets"] is False
        assert m_true["filter_doublets"] is True


# ---------------------------------------------------------------------------
# T11 — Sparse matrix input
# ---------------------------------------------------------------------------

class TestSparseInput:
    def setup_method(self):
        self.m = import_module()

    def test_sparse_clr_layer(self):
        clr = _make_clr()
        clr_sparse = sp.csr_matrix(clr)
        barcodes = [f"cell_{i}" for i in range(N_CELLS)]
        adata = ad.AnnData(X=clr_sparse.copy())
        adata.obs_names = barcodes
        adata.var_names = ADT_VARS
        adata.layers["adt_clr"] = clr_sparse.copy()
        mdata = MockMuData(adata)
        mdata2, metrics = self.m.detect_adt_doublets(mdata)
        assert "adt_doublet_score" in mdata2["adt"].obs.columns

    def test_sparse_scores_match_dense(self):
        """Sparse and dense CLR should give identical scores."""
        clr = _make_clr()
        clr_sparse = sp.csr_matrix(clr)

        # dense version
        mdata_d = _make_mdata(clr)
        mdata_d2, _ = self.m.detect_adt_doublets(mdata_d)
        scores_dense = mdata_d2["adt"].obs["adt_doublet_score"].values

        # sparse version
        barcodes = [f"cell_{i}" for i in range(N_CELLS)]
        adata_s = ad.AnnData(X=clr_sparse.copy())
        adata_s.obs_names = barcodes
        adata_s.var_names = ADT_VARS
        adata_s.layers["adt_clr"] = clr_sparse.copy()
        mdata_s = MockMuData(adata_s)
        mdata_s2, _ = self.m.detect_adt_doublets(mdata_s)
        scores_sparse = mdata_s2["adt"].obs["adt_doublet_score"].values

        np.testing.assert_array_almost_equal(scores_dense, scores_sparse)


# ---------------------------------------------------------------------------
# T12 — Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def setup_method(self):
        self.m = import_module()

    def test_all_cells_are_doublets(self):
        clr = np.full((N_CELLS, N_VARS), 3.0, dtype=np.float32)
        mdata = _make_mdata(clr)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, threshold=2.5)
        assert metrics["n_doublets_detected"] == N_CELLS
        assert np.all(mdata2["adt"].obs["adt_predicted_doublet"].values)

    def test_zero_doublets(self):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        mdata = _make_mdata(clr)
        _, metrics = self.m.detect_adt_doublets(mdata, threshold=2.5)
        assert metrics["n_doublets_detected"] == 0

    def test_single_pair(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19")]
        )
        assert len(metrics["pairs_evaluated"]) == 1

    def test_single_cell(self):
        clr = np.array([[3.0, 3.0, 0.0, 0.0, 0.0, 0.0]], dtype=np.float32)
        mdata = _make_mdata(clr, n_cells=1)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, threshold=2.5)
        assert metrics["n_doublets_detected"] == 1
        assert mdata2["adt"].obs["adt_predicted_doublet"].values[0] == True

    def test_single_cell_not_doublet(self):
        clr = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]], dtype=np.float32)
        mdata = _make_mdata(clr, n_cells=1)
        mdata2, metrics = self.m.detect_adt_doublets(mdata, threshold=2.5)
        assert metrics["n_doublets_detected"] == 0


# ---------------------------------------------------------------------------
# T13 — Threshold sensitivity
# ---------------------------------------------------------------------------

class TestThresholdSensitivity:
    def setup_method(self):
        self.m = import_module()

    def _n_doublets(self, threshold):
        clr = np.full((N_CELLS, N_VARS), 0.0, dtype=np.float32)
        # Set CD3 and CD19 to exactly 2.5
        clr[:, 0] = 2.5
        clr[:, 1] = 2.5
        mdata = _make_mdata(clr)
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19")], threshold=threshold
        )
        return metrics["n_doublets_detected"]

    def test_at_threshold_not_flagged(self):
        """Strict greater-than: value exactly at threshold is NOT a doublet."""
        assert self._n_doublets(2.5) == 0

    def test_below_threshold_not_flagged(self):
        assert self._n_doublets(3.0) == 0

    def test_above_threshold_flagged(self):
        assert self._n_doublets(2.0) == N_CELLS

    def test_very_high_threshold_no_doublets(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(mdata, threshold=999.0)
        assert metrics["n_doublets_detected"] == 0


# ---------------------------------------------------------------------------
# T14 — Custom marker pairs
# ---------------------------------------------------------------------------

class TestCustomMarkerPairs:
    def setup_method(self):
        self.m = import_module()

    def test_custom_pair_evaluated(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD4", "CD8a")]
        )
        assert len(metrics["pairs_evaluated"]) == 1

    def test_multiple_custom_pairs(self):
        mdata = _make_mdata()
        _, metrics = self.m.detect_adt_doublets(
            mdata, marker_pairs=[("CD3", "CD19"), ("CD4", "CD8a"), ("CD3", "CD14")]
        )
        assert len(metrics["pairs_evaluated"]) == 3

    def test_empty_marker_pairs_list(self):
        mdata = _make_mdata()
        mdata2, metrics = self.m.detect_adt_doublets(mdata, marker_pairs=[])
        assert metrics["n_doublets_detected"] == 0
        assert not np.any(mdata2["adt"].obs["adt_predicted_doublet"].values)


# ---------------------------------------------------------------------------
# T15 — Missing adt_clr layer
# ---------------------------------------------------------------------------

class TestMissingLayer:
    def setup_method(self):
        self.m = import_module()

    def test_raises_key_error_without_adt_clr(self):
        adata = ad.AnnData(
            X=np.zeros((N_CELLS, N_VARS)),
        )
        adata.obs_names = [f"cell_{i}" for i in range(N_CELLS)]
        adata.var_names = ADT_VARS
        # deliberately omit adt_clr layer
        mdata = MockMuData(adata)
        with pytest.raises(KeyError, match="adt_clr"):
            self.m.detect_adt_doublets(mdata)

    def test_error_message_mentions_normalize(self):
        adata = ad.AnnData(X=np.zeros((N_CELLS, N_VARS)))
        adata.obs_names = [f"cell_{i}" for i in range(N_CELLS)]
        adata.var_names = ADT_VARS
        mdata = MockMuData(adata)
        with pytest.raises(KeyError, match="normalize_adt"):
            self.m.detect_adt_doublets(mdata)
