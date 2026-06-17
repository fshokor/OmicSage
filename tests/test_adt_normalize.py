"""
tests/test_adt_normalize.py

Unit tests for pipeline.modules.cite.adt_normalize.normalize_adt()

Fixture conventions
-------------------
- Small fixtures only (integer ADT counts, sparse CSR) — no real dataset
- flavor / axis tests kept isolated so failures are unambiguous
- All fixtures use raw integer counts in .X to match qc.py output contract
- muon CLR operates in-place; tests verify both .X and layers["adt_clr"]
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import anndata as ad

# ---------------------------------------------------------------------------
# Import the module under test
# (mirrors the path pipeline/modules/scripts/cite/adt_normalize.py)
# ---------------------------------------------------------------------------
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.modules.scripts.cite.adt_normalize import normalize_adt


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_adt_adata(
    n_cells: int = 20,
    n_proteins: int = 10,
    seed: int = 42,
    sparse: bool = True,
) -> ad.AnnData:
    """
    Build a minimal AnnData with raw integer ADT counts.
    Mimics mdata["adt"] output of qc.py.
    """
    rng = np.random.default_rng(seed)
    # ADT counts are typically low integers (0-500); use Poisson-like distribution
    X = rng.integers(0, 50, size=(n_cells, n_proteins)).astype(float)
    # Introduce some zeros (ADT is sparser than RNA for negative markers)
    X[X < 5] = 0.0
    if sparse:
        X = sp.csr_matrix(X)
    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"protein_{i}" for i in range(n_proteins)]
    return adata


# ---------------------------------------------------------------------------
# 1. Return type and shape
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self):
        adata = _make_adt_adata()
        result = normalize_adt(adata)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_adt_adata()
        adata_out, metrics = normalize_adt(adata)
        assert isinstance(adata_out, ad.AnnData)
        assert isinstance(metrics, dict)

    def test_shape_preserved(self):
        adata = _make_adt_adata(n_cells=30, n_proteins=15)
        adata_out, _ = normalize_adt(adata)
        assert adata_out.shape == (30, 15)

    def test_obs_names_preserved(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert list(adata_out.obs_names) == list(adata.obs_names)

    def test_var_names_preserved(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert list(adata_out.var_names) == list(adata.var_names)


# ---------------------------------------------------------------------------
# 2. Layer contract
# ---------------------------------------------------------------------------

class TestLayers:

    def test_counts_layer_created(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert "counts" in adata_out.layers

    def test_adt_clr_layer_created(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert "adt_clr" in adata_out.layers

    def test_counts_layer_contains_raw_integers(self):
        adata = _make_adt_adata()
        raw_X = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        )
        adata_out, _ = normalize_adt(adata)
        counts = (
            adata_out.layers["counts"].toarray()
            if sp.issparse(adata_out.layers["counts"])
            else np.asarray(adata_out.layers["counts"])
        )
        np.testing.assert_array_equal(counts, raw_X)

    def test_adt_clr_matches_X(self):
        """layers['adt_clr'] must be an exact copy of .X post-normalization."""
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        X = (
            adata_out.X.toarray()
            if sp.issparse(adata_out.X)
            else np.asarray(adata_out.X)
        )
        clr = (
            adata_out.layers["adt_clr"].toarray()
            if sp.issparse(adata_out.layers["adt_clr"])
            else np.asarray(adata_out.layers["adt_clr"])
        )
        np.testing.assert_array_almost_equal(X, clr)

    def test_clr_values_differ_from_raw(self):
        """After CLR, .X must differ from raw counts."""
        adata = _make_adt_adata()
        raw_X = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        )
        adata_out, _ = normalize_adt(adata)
        clr = (
            adata_out.X.toarray()
            if sp.issparse(adata_out.X)
            else np.asarray(adata_out.X)
        )
        # CLR values and raw counts must not be identical
        assert not np.allclose(clr, raw_X), "CLR output is identical to raw counts"

    def test_existing_counts_layer_not_overwritten(self):
        """If layers['counts'] already exists, normalize_adt must not overwrite it."""
        adata = _make_adt_adata()
        sentinel = np.ones(adata.shape) * 999.0
        adata.layers["counts"] = sp.csr_matrix(sentinel)
        adata_out, _ = normalize_adt(adata)
        counts = (
            adata_out.layers["counts"].toarray()
            if sp.issparse(adata_out.layers["counts"])
            else np.asarray(adata_out.layers["counts"])
        )
        np.testing.assert_array_equal(counts, sentinel)


# ---------------------------------------------------------------------------
# 3. CLR correctness
# ---------------------------------------------------------------------------

class TestCLRCorrectness:

    def test_clr_zeros_remain_zero(self):
        """Zero counts must remain zero after CLR (log1p(0/geomean) = 0)."""
        adata = _make_adt_adata()
        raw = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        )
        adata_out, _ = normalize_adt(adata)
        clr = (
            adata_out.X.toarray()
            if sp.issparse(adata_out.X)
            else np.asarray(adata_out.X)
        )
        # Wherever raw count is 0, CLR must also be 0
        assert np.all(clr[raw == 0] == 0.0)

    def test_clr_non_negative(self):
        """CLR values on raw counts must be non-negative."""
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        clr = (
            adata_out.X.toarray()
            if sp.issparse(adata_out.X)
            else np.asarray(adata_out.X)
        )
        assert np.all(clr >= 0.0), f"CLR produced negative values: min={clr.min()}"

    def test_clr_axis0_default(self):
        """Default axis=0: normalization is per-protein across cells."""
        adata = _make_adt_adata()
        adata_out_default, _ = normalize_adt(adata)
        adata_out_axis0, _ = normalize_adt(adata, clr_axis=0)
        X_default = (
            adata_out_default.X.toarray()
            if sp.issparse(adata_out_default.X)
            else np.asarray(adata_out_default.X)
        )
        X_axis0 = (
            adata_out_axis0.X.toarray()
            if sp.issparse(adata_out_axis0.X)
            else np.asarray(adata_out_axis0.X)
        )
        np.testing.assert_array_almost_equal(X_default, X_axis0)

    def test_clr_axis0_and_axis1_differ(self):
        """axis=0 and axis=1 must produce different results."""
        adata = _make_adt_adata()
        adata_out_0, _ = normalize_adt(adata, clr_axis=0)
        adata_out_1, _ = normalize_adt(adata, clr_axis=1)
        X0 = (
            adata_out_0.X.toarray()
            if sp.issparse(adata_out_0.X)
            else np.asarray(adata_out_0.X)
        )
        X1 = (
            adata_out_1.X.toarray()
            if sp.issparse(adata_out_1.X)
            else np.asarray(adata_out_1.X)
        )
        assert not np.allclose(X0, X1), "axis=0 and axis=1 produced identical results"

    def test_clr_dense_input(self):
        """CLR must work on dense (non-sparse) input."""
        adata = _make_adt_adata(sparse=False)
        adata_out, metrics = normalize_adt(adata)
        assert adata_out.n_obs == adata.n_obs
        assert metrics["n_proteins"] == adata.n_vars


# ---------------------------------------------------------------------------
# 4. inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:

    def test_default_not_inplace(self):
        """Default inplace=False must not mutate the input."""
        adata = _make_adt_adata()
        raw_X = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        ).copy()
        normalize_adt(adata)
        X_after = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        )
        np.testing.assert_array_equal(X_after, raw_X)

    def test_inplace_true_mutates_input(self):
        """inplace=True must modify adata.X directly."""
        adata = _make_adt_adata()
        raw_X = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        ).copy()
        normalize_adt(adata, inplace=True)
        X_after = (
            adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        )
        assert not np.allclose(X_after, raw_X), "inplace=True did not mutate adata.X"


# ---------------------------------------------------------------------------
# 5. Metrics dict contract
# ---------------------------------------------------------------------------

class TestMetrics:

    def test_required_keys_present(self):
        required = {
            "n_cells",
            "n_proteins",
            "clr_axis",
            "dsb_applied",
            "raw_counts_in_layer",
            "clr_in_layer",
            "clr_mean_per_cell",
            "clr_max",
            "clr_min",
            "raw_median_total_counts_per_cell",
            "raw_max_count",
        }
        adata = _make_adt_adata()
        _, metrics = normalize_adt(adata)
        missing = required - set(metrics.keys())
        assert not missing, f"Metrics dict missing keys: {missing}"

    def test_n_cells_correct(self):
        adata = _make_adt_adata(n_cells=25)
        _, metrics = normalize_adt(adata)
        assert metrics["n_cells"] == 25

    def test_n_proteins_correct(self):
        adata = _make_adt_adata(n_proteins=12)
        _, metrics = normalize_adt(adata)
        assert metrics["n_proteins"] == 12

    def test_dsb_applied_false(self):
        adata = _make_adt_adata()
        _, metrics = normalize_adt(adata)
        assert metrics["dsb_applied"] is False

    def test_layer_name_values(self):
        adata = _make_adt_adata()
        _, metrics = normalize_adt(adata)
        assert metrics["raw_counts_in_layer"] == "counts"
        assert metrics["clr_in_layer"] == "adt_clr"

    def test_clr_axis_recorded(self):
        adata = _make_adt_adata()
        _, metrics = normalize_adt(adata, clr_axis=1)
        assert metrics["clr_axis"] == 1


# ---------------------------------------------------------------------------
# 6. Provenance (uns)
# ---------------------------------------------------------------------------

class TestProvenance:

    def test_uns_key_present(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert "omicsage_adt_normalize" in adata_out.uns

    def test_uns_required_fields(self):
        required = {
            "clr_axis",
            "clr_axis_description",
            "dsb_applied",
            "normalized_in_layer",
            "raw_in_layer",
            "muon_version",
            "omicsage_module",
            "omicsage_version",
            "timestamp",
        }
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        prov = adata_out.uns["omicsage_adt_normalize"]
        missing = required - set(prov.keys())
        assert not missing, f"Provenance uns missing keys: {missing}"

    def test_uns_module_name(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        assert (
            adata_out.uns["omicsage_adt_normalize"]["omicsage_module"]
            == "pipeline.modules.cite.adt_normalize"
        )

    def test_uns_timestamp_is_iso_string(self):
        from datetime import datetime
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata)
        ts = adata_out.uns["omicsage_adt_normalize"]["timestamp"]
        # Should parse without error
        datetime.fromisoformat(ts)

    def test_uns_axis_description_axis0(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata, clr_axis=0)
        desc = adata_out.uns["omicsage_adt_normalize"]["clr_axis_description"]
        assert "per-protein" in desc

    def test_uns_axis_description_axis1(self):
        adata = _make_adt_adata()
        adata_out, _ = normalize_adt(adata, clr_axis=1)
        desc = adata_out.uns["omicsage_adt_normalize"]["clr_axis_description"]
        assert "per-cell" in desc


# ---------------------------------------------------------------------------
# 7. Input validation
# ---------------------------------------------------------------------------

class TestValidation:

    def test_rejects_non_anndata(self):
        with pytest.raises(TypeError, match="AnnData"):
            normalize_adt("not_an_anndata")

    def test_rejects_empty_cells(self):
        adata = _make_adt_adata()
        empty = adata[:0].copy()
        with pytest.raises(ValueError, match="0 cells"):
            normalize_adt(empty)

    def test_rejects_empty_proteins(self):
        adata = _make_adt_adata()
        empty = adata[:, :0].copy()
        with pytest.raises(ValueError, match="0 proteins"):
            normalize_adt(empty)

    def test_rejects_invalid_axis(self):
        adata = _make_adt_adata()
        with pytest.raises(ValueError, match="clr_axis must be 0 or 1"):
            normalize_adt(adata, clr_axis=2)

    def test_rejects_axis_minus_one(self):
        adata = _make_adt_adata()
        with pytest.raises(ValueError, match="clr_axis must be 0 or 1"):
            normalize_adt(adata, clr_axis=-1)


    def test_non_integer_input_warns(self, caplog):
        """Non-integer .X should trigger a warning, not an error."""
        import logging
        adata = _make_adt_adata()
        # Replace X with non-integer (already-normalized) data
        X_float = np.random.default_rng(0).random(adata.shape)
        adata.X = sp.csr_matrix(X_float)
        with caplog.at_level(logging.WARNING, logger="pipeline.modules.cite.adt_normalize"):
            normalize_adt(adata)  # should not raise
        assert any("non-integer" in r.message.lower() for r in caplog.records)
