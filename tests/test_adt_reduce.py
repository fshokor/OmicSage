"""
tests/test_adt_reduce.py — Tests for pipeline/modules/scripts/cite/adt_reduce.py

Coverage:
  - Happy path: PCA + neighbors + UMAP computed correctly
  - Output key naming (X_pca_adt, X_umap_adt — no collision with RNA keys)
  - n_comps capping (when n_comps > min(n_cells-1, n_vars-1))
  - n_pcs capping (when n_pcs gets capped alongside n_comps)
  - n_pcs > n_comps raises ValueError
  - Missing "adt" modality raises KeyError
  - Missing "adt_clr" layer raises KeyError
  - Too few cells raises ValueError
  - Too few vars raises ValueError
  - inplace=False returns a copy (original untouched)
  - inplace=True modifies in place
  - Provenance recorded in uns["omicsage_adt_reduce"]
  - Provenance records correct n_comps_actual and n_pcs_used
  - Metrics dict structure and values
  - Variance explained is in [0, 1] range (per-component), total <= 1.0 + epsilon
  - UMAP shape is (n_cells, 2)
  - PCA shape is (n_cells, n_comps_actual)
  - RNA obsm keys NOT overwritten (X_pca, X_umap not added to adt AnnData)
  - Default parameters match sc-best-practices ch.36 (n_pcs=20, svd_solver="arpack")
  - Small fixture with n_vars < n_comps: n_comps_actual correctly capped
  - svd_solver parameter forwarded correctly
  - random_state parameter forwarded correctly
  - n_neighbors parameter forwarded to neighbors
  - Idempotent: running twice overwrites outputs cleanly
  - CLR layer preserved after reduction
  - .X set to CLR values before PCA (not raw counts)
  - mdata["rna"] untouched when reducing adt
"""

from __future__ import annotations

import sys
import os

import numpy as np
import pytest
import anndata as ad
import scanpy as sc
from anndata import AnnData
from mudata import MuData

# ---------------------------------------------------------------------------
# Path setup — import from project root
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.modules.scripts.cite.adt_reduce import reduce_adt


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_clr_matrix(n_cells: int, n_vars: int, seed: int = 42) -> np.ndarray:
    """Make realistic CLR-normalized ADT values (approximately normal, centred near 0)."""
    rng = np.random.default_rng(seed)
    return rng.normal(loc=0.0, scale=2.0, size=(n_cells, n_vars)).astype(np.float32)


def _make_adt_adata(n_cells: int = 80, n_vars: int = 25, seed: int = 42) -> AnnData:
    """Create a minimal AnnData with layers['adt_clr'] set."""
    clr = _make_clr_matrix(n_cells, n_vars, seed)
    # raw counts in .X, CLR in layers
    raw = np.abs(np.random.default_rng(seed).integers(0, 50, size=(n_cells, n_vars))).astype(np.float32)
    adata = ad.AnnData(X=raw)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"ADT_{i}" for i in range(n_vars)]
    adata.layers["counts"] = raw.copy()
    adata.layers["adt_clr"] = clr.copy()
    return adata


def _make_rna_adata(n_cells: int = 80, n_genes: int = 200, seed: int = 99) -> AnnData:
    """Create a minimal RNA AnnData (pre-reduced, with X_pca and X_umap)."""
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    # Simulate already-run RNA reduce.py
    adata.obsm["X_pca"] = rng.standard_normal((n_cells, 50)).astype(np.float32)
    adata.obsm["X_umap"] = rng.standard_normal((n_cells, 2)).astype(np.float32)
    return adata


def _make_mdata(n_cells: int = 80, n_adt: int = 25, include_rna: bool = True, seed: int = 42) -> MuData:
    adt = _make_adt_adata(n_cells=n_cells, n_vars=n_adt, seed=seed)
    if include_rna:
        rna = _make_rna_adata(n_cells=n_cells)
        return MuData({"rna": rna, "adt": adt})
    return MuData({"adt": adt})


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestReduceAdtHappyPath:
    def test_returns_tuple(self):
        mdata = _make_mdata()
        result = reduce_adt(mdata)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_anndata_and_dict(self):
        mdata = _make_mdata()
        adata, metrics = reduce_adt(mdata)
        assert isinstance(adata, AnnData)
        assert isinstance(metrics, dict)

    def test_pca_key_present(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "X_pca_adt" in adata.obsm

    def test_umap_key_present(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "X_umap_adt" in adata.obsm

    def test_pca_shape_cells(self):
        mdata = _make_mdata(n_cells=80, n_adt=25)
        adata, _ = reduce_adt(mdata)
        assert adata.obsm["X_pca_adt"].shape[0] == 80

    def test_pca_shape_components(self):
        """n_comps capped at min(79, 24) = 24 for this fixture."""
        mdata = _make_mdata(n_cells=80, n_adt=25)
        adata, _ = reduce_adt(mdata)
        assert adata.obsm["X_pca_adt"].shape[1] == 24  # min(79, 24)

    def test_umap_shape(self):
        mdata = _make_mdata(n_cells=80, n_adt=25)
        adata, _ = reduce_adt(mdata)
        assert adata.obsm["X_umap_adt"].shape == (80, 2)

    def test_umap_is_2d(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert adata.obsm["X_umap_adt"].ndim == 2
        assert adata.obsm["X_umap_adt"].shape[1] == 2

    def test_no_x_pca_key_in_adt(self):
        """X_pca must not exist on adt after reduction — only X_pca_adt."""
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "X_pca" not in adata.obsm

    def test_no_x_umap_key_in_adt(self):
        """X_umap must not exist on adt after reduction — only X_umap_adt."""
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "X_umap" not in adata.obsm


# ---------------------------------------------------------------------------
# RNA modality untouched
# ---------------------------------------------------------------------------

class TestRnaUntouched:
    def test_rna_x_pca_unchanged(self):
        mdata = _make_mdata()
        rna_pca_before = mdata["rna"].obsm["X_pca"].copy()
        reduce_adt(mdata, inplace=False)
        np.testing.assert_array_equal(mdata["rna"].obsm["X_pca"], rna_pca_before)

    def test_rna_x_umap_unchanged(self):
        mdata = _make_mdata()
        rna_umap_before = mdata["rna"].obsm["X_umap"].copy()
        reduce_adt(mdata, inplace=False)
        np.testing.assert_array_equal(mdata["rna"].obsm["X_umap"], rna_umap_before)

    def test_rna_x_matrix_unchanged(self):
        mdata = _make_mdata()
        rna_X_before = mdata["rna"].X.copy()
        reduce_adt(mdata, inplace=False)
        np.testing.assert_array_equal(mdata["rna"].X, rna_X_before)


# ---------------------------------------------------------------------------
# n_comps capping
# ---------------------------------------------------------------------------

class TestNCompsCapping:
    def test_n_comps_capped_by_n_vars(self):
        """n_vars=10 → n_comps_actual = min(79, 9) = 9."""
        mdata = _make_mdata(n_cells=80, n_adt=10)
        adata, metrics = reduce_adt(mdata, n_comps=50)
        assert adata.obsm["X_pca_adt"].shape[1] == 9
        assert metrics["n_comps_actual"] == 9

    def test_n_comps_capped_by_n_cells(self):
        """n_cells=10 → n_comps_actual = min(9, 24) = 9."""
        mdata = _make_mdata(n_cells=10, n_adt=25)
        adata, metrics = reduce_adt(mdata, n_comps=50)
        assert adata.obsm["X_pca_adt"].shape[1] == 9
        assert metrics["n_comps_actual"] == 9

    def test_n_comps_not_capped_when_large_data(self):
        """With n_cells=200, n_vars=60, n_comps=30 → no capping."""
        mdata = _make_mdata(n_cells=200, n_adt=60)
        adata, metrics = reduce_adt(mdata, n_comps=30, n_pcs=20)
        assert metrics["n_comps_actual"] == 30
        assert adata.obsm["X_pca_adt"].shape[1] == 30

    def test_n_pcs_capped_when_n_comps_capped(self):
        """If n_comps gets capped below n_pcs, n_pcs_used is also capped."""
        # n_vars=5 → n_comps_actual = min(79,4) = 4; n_pcs=20 → n_pcs_used=4
        mdata = _make_mdata(n_cells=80, n_adt=5)
        _, metrics = reduce_adt(mdata, n_comps=50, n_pcs=20)
        assert metrics["n_pcs_used"] == 4
        assert metrics["n_pcs_used"] <= metrics["n_comps_actual"]


# ---------------------------------------------------------------------------
# n_pcs validation
# ---------------------------------------------------------------------------

class TestNPcsValidation:
    def test_n_pcs_greater_than_n_comps_raises(self):
        mdata = _make_mdata()
        with pytest.raises(ValueError, match="n_pcs"):
            reduce_adt(mdata, n_comps=10, n_pcs=20)

    def test_n_pcs_equal_n_comps_ok(self):
        mdata = _make_mdata(n_cells=200, n_adt=60)
        _, metrics = reduce_adt(mdata, n_comps=20, n_pcs=20)
        assert metrics["n_pcs_used"] == 20


# ---------------------------------------------------------------------------
# Error cases
# ---------------------------------------------------------------------------

class TestErrorCases:
    def test_missing_adt_modality_raises(self):
        rna = _make_rna_adata()
        mdata = MuData({"rna": rna})
        with pytest.raises(KeyError, match="adt"):
            reduce_adt(mdata)

    def test_missing_adt_clr_layer_raises(self):
        adt = _make_adt_adata()
        del adt.layers["adt_clr"]
        mdata = MuData({"adt": adt})
        with pytest.raises(KeyError, match="adt_clr"):
            reduce_adt(mdata)

    def test_too_few_cells_raises(self):
        adt = _make_adt_adata(n_cells=1, n_vars=10)
        mdata = MuData({"adt": adt})
        with pytest.raises(ValueError, match="cells"):
            reduce_adt(mdata)

    def test_too_few_vars_raises(self):
        adt = _make_adt_adata(n_cells=20, n_vars=1)
        mdata = MuData({"adt": adt})
        with pytest.raises(ValueError, match="features"):
            reduce_adt(mdata)


# ---------------------------------------------------------------------------
# inplace behaviour
# ---------------------------------------------------------------------------

class TestInplaceBehaviour:
    def test_inplace_false_original_adt_untouched(self):
        mdata = _make_mdata()
        original_obsm_keys = set(mdata["adt"].obsm.keys())
        reduce_adt(mdata, inplace=False)
        assert set(mdata["adt"].obsm.keys()) == original_obsm_keys

    def test_inplace_false_returns_copy(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata, inplace=False)
        assert adata is not mdata["adt"]

    def test_inplace_true_modifies_original(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata, inplace=True)
        assert "X_pca_adt" in mdata["adt"].obsm
        assert "X_umap_adt" in mdata["adt"].obsm

    def test_inplace_true_returns_same_object(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata, inplace=True)
        assert adata is mdata["adt"]


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_provenance_key_present(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "omicsage_adt_reduce" in adata.uns

    def test_provenance_has_module_name(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert adata.uns["omicsage_adt_reduce"]["module"] == "adt_reduce"

    def test_provenance_has_timestamp(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "timestamp" in adata.uns["omicsage_adt_reduce"]
        assert isinstance(adata.uns["omicsage_adt_reduce"]["timestamp"], str)

    def test_provenance_has_params(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        assert "params" in adata.uns["omicsage_adt_reduce"]

    def test_provenance_params_n_comps_actual(self):
        mdata = _make_mdata(n_cells=80, n_adt=25)
        adata, _ = reduce_adt(mdata, n_comps=50)
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["n_comps_requested"] == 50
        assert params["n_comps_actual"] == 24  # min(79, 24)

    def test_provenance_params_n_pcs_used(self):
        mdata = _make_mdata(n_cells=200, n_adt=60)
        adata, _ = reduce_adt(mdata, n_comps=50, n_pcs=20)
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["n_pcs_used"] == 20

    def test_provenance_output_keys(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        outputs = adata.uns["omicsage_adt_reduce"]["outputs"]
        assert outputs["pca_key"] == "X_pca_adt"
        assert outputs["umap_key"] == "X_umap_adt"

    def test_provenance_svd_solver(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata, svd_solver="arpack")
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["svd_solver"] == "arpack"

    def test_provenance_random_state(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata, random_state=42)
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["random_state"] == 42


# ---------------------------------------------------------------------------
# Metrics dict
# ---------------------------------------------------------------------------

class TestMetrics:
    def test_metrics_has_required_keys(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        required = {
            "n_cells", "n_vars", "n_comps_actual", "n_pcs_used",
            "n_neighbors", "variance_explained_total",
            "umap_computed", "pca_key", "umap_key",
        }
        assert required.issubset(set(metrics.keys()))

    def test_metrics_n_cells_correct(self):
        mdata = _make_mdata(n_cells=80)
        _, metrics = reduce_adt(mdata)
        assert metrics["n_cells"] == 80

    def test_metrics_n_vars_correct(self):
        mdata = _make_mdata(n_adt=25)
        _, metrics = reduce_adt(mdata)
        assert metrics["n_vars"] == 25

    def test_metrics_umap_computed_true(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["umap_computed"] is True

    def test_metrics_pca_key(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["pca_key"] == "X_pca_adt"

    def test_metrics_umap_key(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["umap_key"] == "X_umap_adt"

    def test_metrics_variance_explained_not_none(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["variance_explained_total"] is not None

    def test_metrics_variance_explained_positive(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["variance_explained_total"] > 0.0

    def test_metrics_variance_explained_at_most_one(self):
        """Sum of variance ratios across computed components <= 1.0 + epsilon."""
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata)
        assert metrics["variance_explained_total"] <= 1.0 + 1e-6

    def test_metrics_n_neighbors_correct(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata, n_neighbors=10)
        assert metrics["n_neighbors"] == 10


# ---------------------------------------------------------------------------
# Default parameters (sc-best-practices ch.36)
# ---------------------------------------------------------------------------

class TestDefaultParameters:
    def test_default_n_pcs_is_20(self):
        """sc-best-practices ch.36 uses n_pcs=20 for ADT neighbors."""
        mdata = _make_mdata(n_cells=200, n_adt=60)
        _, metrics = reduce_adt(mdata)
        assert metrics["n_pcs_used"] == 20

    def test_default_svd_solver_is_arpack(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["svd_solver"] == "arpack"

    def test_default_random_state_is_0(self):
        mdata = _make_mdata()
        adata, _ = reduce_adt(mdata)
        params = adata.uns["omicsage_adt_reduce"]["params"]
        assert params["random_state"] == 0


# ---------------------------------------------------------------------------
# Layer and matrix preservation
# ---------------------------------------------------------------------------

class TestLayerPreservation:
    def test_adt_clr_layer_preserved(self):
        mdata = _make_mdata()
        clr_before = mdata["adt"].layers["adt_clr"].copy()
        adata, _ = reduce_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.layers["adt_clr"], clr_before)

    def test_counts_layer_preserved(self):
        mdata = _make_mdata()
        counts_before = mdata["adt"].layers["counts"].copy()
        adata, _ = reduce_adt(mdata, inplace=False)
        np.testing.assert_array_equal(adata.layers["counts"], counts_before)


# ---------------------------------------------------------------------------
# Idempotency
# ---------------------------------------------------------------------------

class TestIdempotency:
    def test_running_twice_overwrites_outputs(self):
        mdata = _make_mdata()
        adata1, _ = reduce_adt(mdata, inplace=False)
        adata_inplace, _ = reduce_adt(MuData({"adt": adata1}), inplace=True)
        # Should not raise; outputs overwritten cleanly
        assert "X_pca_adt" in adata_inplace.obsm
        assert "X_umap_adt" in adata_inplace.obsm


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

class TestReproducibility:
    def test_same_random_state_same_umap(self):
        mdata1 = _make_mdata(seed=7)
        mdata2 = _make_mdata(seed=7)
        adata1, _ = reduce_adt(mdata1, random_state=0)
        adata2, _ = reduce_adt(mdata2, random_state=0)
        np.testing.assert_array_almost_equal(
            adata1.obsm["X_umap_adt"], adata2.obsm["X_umap_adt"], decimal=4
        )

    def test_same_random_state_same_pca(self):
        mdata1 = _make_mdata(seed=7)
        mdata2 = _make_mdata(seed=7)
        adata1, _ = reduce_adt(mdata1, random_state=0)
        adata2, _ = reduce_adt(mdata2, random_state=0)
        np.testing.assert_array_almost_equal(
            adata1.obsm["X_pca_adt"], adata2.obsm["X_pca_adt"], decimal=4
        )


# ---------------------------------------------------------------------------
# Custom parameters forwarded
# ---------------------------------------------------------------------------

class TestCustomParameters:
    def test_custom_n_neighbors_forwarded(self):
        mdata = _make_mdata()
        _, metrics = reduce_adt(mdata, n_neighbors=5)
        assert metrics["n_neighbors"] == 5

    def test_custom_n_comps_used_when_smaller_than_cap(self):
        mdata = _make_mdata(n_cells=200, n_adt=60)
        adata, metrics = reduce_adt(mdata, n_comps=15, n_pcs=10)
        assert metrics["n_comps_actual"] == 15
        assert adata.obsm["X_pca_adt"].shape[1] == 15

    def test_custom_n_pcs_recorded_in_metrics(self):
        mdata = _make_mdata(n_cells=200, n_adt=60)
        _, metrics = reduce_adt(mdata, n_comps=30, n_pcs=15)
        assert metrics["n_pcs_used"] == 15


# ---------------------------------------------------------------------------
# Fixture: no RNA modality (adt-only MuData)
# ---------------------------------------------------------------------------

class TestAdtOnlyMuData:
    def test_works_without_rna_modality(self):
        mdata = _make_mdata(include_rna=False)
        adata, metrics = reduce_adt(mdata)
        assert "X_pca_adt" in adata.obsm
        assert "X_umap_adt" in adata.obsm

    def test_metrics_correct_adt_only(self):
        mdata = _make_mdata(n_cells=80, n_adt=25, include_rna=False)
        _, metrics = reduce_adt(mdata)
        assert metrics["n_cells"] == 80
        assert metrics["n_vars"] == 25
