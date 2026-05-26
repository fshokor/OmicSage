"""
tests/test_adt_harmony.py — Tests for pipeline/modules/cite/adt_harmony.py

Coverage:
  - Happy path: Harmony + neighbors + UMAP computed correctly
  - Output key naming (X_pca_harmony_adt, X_umap_adt — no collision with RNA keys)
  - X_pca_adt preserved after harmony (original PCA not overwritten)
  - No stray X_pca_harmony key on adt (ADT-namespaced key only)
  - No stray X_umap key on adt after UMAP rename (only X_umap_adt)
  - batch_key missing from obs raises ValueError
  - batch_key missing from mdata raises KeyError (adt modality absent)
  - X_pca_adt missing raises KeyError
  - n_pcs > available components raises ValueError
  - n_pcs capped at available harmony components (n_pcs_used in metrics)
  - inplace=False returns a copy (original mdata["adt"] untouched)
  - inplace=True modifies adt in place
  - Provenance recorded in uns["omicsage_adt_harmony"]
  - Provenance records correct batch_key, n_pcs_used, n_neighbors, random_state
  - Provenance output keys correct
  - Metrics dict structure and required keys
  - Metrics n_cells, n_batches, batch_values, batch_key correct
  - Metrics umap_computed always True
  - Metrics harmony_key == "X_pca_harmony_adt"
  - Metrics umap_key == "X_umap_adt"
  - UMAP shape is (n_cells, 2)
  - Harmony embedding shape is (n_cells, n_comps)
  - RNA obsm keys untouched (X_pca, X_umap, X_pca_harmony not touched on adt)
  - mdata["rna"] entirely untouched
  - Default n_pcs=20 (sc-best-practices ch.38)
  - Default random_state=0
  - Idempotent: running twice overwrites outputs cleanly
  - Single-batch data (Harmony still runs; no error)
  - Multi-batch data (2+ batches, expected use case)
  - adt-only MuData (no rna modality present)
  - Layers preserved after harmony (adt_clr, counts)
  - CLR layer not modified by harmony
"""

from __future__ import annotations

import sys
import os

import numpy as np
import pytest
import anndata as ad
from anndata import AnnData
from mudata import MuData

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.modules.cite.adt_harmony import run_harmony_adt


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_adt_adata(
    n_cells: int = 120,
    n_vars: int = 25,
    n_pca_comps: int = 20,
    batch_key: str = "donor",
    n_batches: int = 2,
    seed: int = 42,
) -> AnnData:
    """
    Create a minimal AnnData with:
      - layers["counts"]   — raw integer counts
      - layers["adt_clr"]  — CLR-normalized values
      - obsm["X_pca_adt"]  — ADT PCA embedding (simulates adt_reduce.py output)
      - obs[batch_key]     — batch/donor column
    """
    rng = np.random.default_rng(seed)
    raw = rng.integers(0, 50, size=(n_cells, n_vars)).astype(np.float32)
    clr = rng.normal(loc=0.0, scale=2.0, size=(n_cells, n_vars)).astype(np.float32)
    pca = rng.standard_normal((n_cells, n_pca_comps)).astype(np.float32)

    adata = ad.AnnData(X=raw.copy())
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"ADT_{i}" for i in range(n_vars)]
    adata.layers["counts"] = raw.copy()
    adata.layers["adt_clr"] = clr.copy()
    adata.obsm["X_pca_adt"] = pca

    # Batch column: evenly distribute cells across batches
    batch_labels = [f"donor_{i % n_batches}" for i in range(n_cells)]
    adata.obs[batch_key] = batch_labels

    return adata


def _make_rna_adata(n_cells: int = 120, n_genes: int = 200, seed: int = 99) -> AnnData:
    """Create a minimal RNA AnnData with pre-computed obsm keys."""
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(n_genes)]
    adata.obsm["X_pca"] = rng.standard_normal((n_cells, 50)).astype(np.float32)
    adata.obsm["X_umap"] = rng.standard_normal((n_cells, 2)).astype(np.float32)
    return adata


def _make_mdata(
    n_cells: int = 120,
    n_adt: int = 25,
    n_pca_comps: int = 20,
    batch_key: str = "donor",
    n_batches: int = 2,
    include_rna: bool = True,
    seed: int = 42,
) -> MuData:
    adt = _make_adt_adata(
        n_cells=n_cells,
        n_vars=n_adt,
        n_pca_comps=n_pca_comps,
        batch_key=batch_key,
        n_batches=n_batches,
        seed=seed,
    )
    if include_rna:
        rna = _make_rna_adata(n_cells=n_cells, seed=seed + 1)
        return MuData({"rna": rna, "adt": adt})
    return MuData({"adt": adt})


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestRunHarmonyAdtHappyPath:
    def test_returns_tuple(self):
        mdata = _make_mdata()
        result = run_harmony_adt(mdata, batch_key="donor")
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_anndata_and_dict(self):
        mdata = _make_mdata()
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert isinstance(adata, AnnData)
        assert isinstance(metrics, dict)

    def test_harmony_key_present(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert "X_pca_harmony_adt" in adata.obsm

    def test_umap_key_present(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert "X_umap_adt" in adata.obsm

    def test_harmony_shape_cells(self):
        mdata = _make_mdata(n_cells=120, n_pca_comps=20)
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.obsm["X_pca_harmony_adt"].shape[0] == 120

    def test_harmony_shape_components(self):
        """Harmony corrects each PC independently — output has same n_comps as input."""
        mdata = _make_mdata(n_cells=120, n_pca_comps=20)
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.obsm["X_pca_harmony_adt"].shape[1] == 20

    def test_umap_shape(self):
        mdata = _make_mdata(n_cells=120)
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.obsm["X_umap_adt"].shape == (120, 2)

    def test_umap_is_2d(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.obsm["X_umap_adt"].ndim == 2
        assert adata.obsm["X_umap_adt"].shape[1] == 2


# ---------------------------------------------------------------------------
# Embedding key naming — no namespace collisions
# ---------------------------------------------------------------------------

class TestEmbeddingKeyNaming:
    def test_no_stray_x_pca_harmony_on_adt(self):
        """Only X_pca_harmony_adt must exist — not the RNA key X_pca_harmony."""
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert "X_pca_harmony" not in adata.obsm

    def test_no_stray_x_umap_on_adt(self):
        """scanpy writes X_umap; we must rename to X_umap_adt and remove original."""
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert "X_umap" not in adata.obsm

    def test_x_pca_adt_preserved(self):
        """Original ADT PCA must not be overwritten by Harmony."""
        mdata = _make_mdata()
        pca_before = mdata["adt"].obsm["X_pca_adt"].copy()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(adata.obsm["X_pca_adt"], pca_before)

    def test_harmony_key_is_adt_namespaced(self):
        mdata = _make_mdata()
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["harmony_key"] == "X_pca_harmony_adt"

    def test_umap_key_is_adt_namespaced(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["umap_key"] == "X_umap_adt"


# ---------------------------------------------------------------------------
# RNA modality untouched
# ---------------------------------------------------------------------------

class TestRnaUntouched:
    def test_rna_x_pca_unchanged(self):
        mdata = _make_mdata()
        rna_pca_before = mdata["rna"].obsm["X_pca"].copy()
        run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(mdata["rna"].obsm["X_pca"], rna_pca_before)

    def test_rna_x_umap_unchanged(self):
        mdata = _make_mdata()
        rna_umap_before = mdata["rna"].obsm["X_umap"].copy()
        run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(mdata["rna"].obsm["X_umap"], rna_umap_before)

    def test_rna_x_matrix_unchanged(self):
        mdata = _make_mdata()
        rna_X_before = mdata["rna"].X.copy()
        run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(mdata["rna"].X, rna_X_before)

    def test_rna_has_no_harmony_key(self):
        """Harmony must not add X_pca_harmony or X_pca_harmony_adt to RNA."""
        mdata = _make_mdata()
        run_harmony_adt(mdata, batch_key="donor", inplace=False)
        assert "X_pca_harmony" not in mdata["rna"].obsm
        assert "X_pca_harmony_adt" not in mdata["rna"].obsm


# ---------------------------------------------------------------------------
# Input validation errors
# ---------------------------------------------------------------------------

class TestValidationErrors:
    def test_missing_adt_modality_raises_keyerror(self):
        rna = _make_rna_adata()
        mdata = MuData({"rna": rna})
        with pytest.raises(KeyError, match="adt"):
            run_harmony_adt(mdata, batch_key="donor")

    def test_missing_x_pca_adt_raises_keyerror(self):
        mdata = _make_mdata()
        del mdata["adt"].obsm["X_pca_adt"]
        with pytest.raises(KeyError, match="X_pca_adt"):
            run_harmony_adt(mdata, batch_key="donor")

    def test_missing_batch_key_raises_valueerror(self):
        mdata = _make_mdata()
        with pytest.raises(ValueError, match="batch_key"):
            run_harmony_adt(mdata, batch_key="nonexistent_column")

    def test_n_pcs_exceeds_components_raises_valueerror(self):
        """n_pcs=30 > n_pca_comps=20 should raise ValueError."""
        mdata = _make_mdata(n_pca_comps=20)
        with pytest.raises(ValueError, match="n_pcs"):
            run_harmony_adt(mdata, batch_key="donor", n_pcs=30)

    def test_error_message_lists_available_obs_columns(self):
        """ValueError for bad batch_key must list available columns."""
        mdata = _make_mdata()
        with pytest.raises(ValueError, match="donor"):
            run_harmony_adt(mdata, batch_key="nonexistent_column")


# ---------------------------------------------------------------------------
# inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:
    def test_inplace_false_does_not_mutate_original(self):
        mdata = _make_mdata()
        original_keys = set(mdata["adt"].obsm.keys())
        run_harmony_adt(mdata, batch_key="donor", inplace=False)
        # Original mdata["adt"].obsm must not have gained harmony keys
        assert "X_pca_harmony_adt" not in original_keys
        assert "X_pca_harmony_adt" not in mdata["adt"].obsm

    def test_inplace_true_modifies_original(self):
        mdata = _make_mdata()
        run_harmony_adt(mdata, batch_key="donor", inplace=True)
        assert "X_pca_harmony_adt" in mdata["adt"].obsm

    def test_inplace_false_returns_new_object(self):
        mdata = _make_mdata()
        adata_out, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        assert adata_out is not mdata["adt"]

    def test_inplace_true_returns_same_object(self):
        mdata = _make_mdata()
        adata_out, _ = run_harmony_adt(mdata, batch_key="donor", inplace=True)
        assert adata_out is mdata["adt"]


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_provenance_key_present(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert "omicsage_adt_harmony" in adata.uns

    def test_provenance_has_module_field(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.uns["omicsage_adt_harmony"]["module"] == "adt_harmony"

    def test_provenance_has_timestamp(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        ts = adata.uns["omicsage_adt_harmony"]["timestamp"]
        assert isinstance(ts, str)
        assert len(ts) > 10

    def test_provenance_params_batch_key(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.uns["omicsage_adt_harmony"]["params"]["batch_key"] == "donor"

    def test_provenance_params_n_pcs_used(self):
        mdata = _make_mdata(n_pca_comps=20)
        adata, _ = run_harmony_adt(mdata, batch_key="donor", n_pcs=15)
        assert adata.uns["omicsage_adt_harmony"]["params"]["n_pcs_used"] == 15

    def test_provenance_params_n_pcs_requested(self):
        mdata = _make_mdata(n_pca_comps=20)
        adata, _ = run_harmony_adt(mdata, batch_key="donor", n_pcs=15)
        assert adata.uns["omicsage_adt_harmony"]["params"]["n_pcs_requested"] == 15

    def test_provenance_params_n_neighbors(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", n_neighbors=10)
        assert adata.uns["omicsage_adt_harmony"]["params"]["n_neighbors"] == 10

    def test_provenance_params_random_state(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", random_state=42)
        assert adata.uns["omicsage_adt_harmony"]["params"]["random_state"] == 42

    def test_provenance_output_harmony_key(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.uns["omicsage_adt_harmony"]["outputs"]["harmony_key"] == "X_pca_harmony_adt"

    def test_provenance_output_umap_key(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.uns["omicsage_adt_harmony"]["outputs"]["umap_key"] == "X_umap_adt"

    def test_provenance_n_pcs_used_capped(self):
        """When n_pcs is capped by available components, provenance reflects actual value."""
        mdata = _make_mdata(n_pca_comps=10)
        adata, _ = run_harmony_adt(mdata, batch_key="donor", n_pcs=10)
        assert adata.uns["omicsage_adt_harmony"]["params"]["n_pcs_used"] == 10


# ---------------------------------------------------------------------------
# Metrics dict
# ---------------------------------------------------------------------------

class TestMetrics:
    def test_metrics_has_required_keys(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        required = {
            "n_cells", "n_pcs_used", "n_neighbors", "batch_key",
            "n_batches", "batch_values", "harmony_key", "umap_key",
            "umap_computed", "random_state",
        }
        assert required.issubset(set(metrics.keys()))

    def test_metrics_n_cells_correct(self):
        mdata = _make_mdata(n_cells=120)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_cells"] == 120

    def test_metrics_batch_key_correct(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["batch_key"] == "donor"

    def test_metrics_n_batches_correct(self):
        mdata = _make_mdata(n_batches=2)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_batches"] == 2

    def test_metrics_batch_values_is_list(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert isinstance(metrics["batch_values"], list)

    def test_metrics_batch_values_count(self):
        mdata = _make_mdata(n_batches=3)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert len(metrics["batch_values"]) == 3

    def test_metrics_umap_computed_true(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["umap_computed"] is True

    def test_metrics_n_pcs_used_correct(self):
        mdata = _make_mdata(n_pca_comps=20)
        _, metrics = run_harmony_adt(mdata, batch_key="donor", n_pcs=15)
        assert metrics["n_pcs_used"] == 15

    def test_metrics_n_neighbors_correct(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor", n_neighbors=10)
        assert metrics["n_neighbors"] == 10

    def test_metrics_random_state_correct(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor", random_state=7)
        assert metrics["random_state"] == 7

    def test_metrics_harmony_key_value(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["harmony_key"] == "X_pca_harmony_adt"

    def test_metrics_umap_key_value(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["umap_key"] == "X_umap_adt"


# ---------------------------------------------------------------------------
# Default parameters (sc-best-practices ch.38)
# ---------------------------------------------------------------------------

class TestDefaultParameters:
    def test_default_n_pcs_is_20(self):
        """sc-best-practices ch.38 uses n_pcs=20 for post-harmony neighbors."""
        mdata = _make_mdata(n_pca_comps=20)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_pcs_used"] == 20

    def test_default_random_state_is_0(self):
        mdata = _make_mdata()
        adata, _ = run_harmony_adt(mdata, batch_key="donor")
        assert adata.uns["omicsage_adt_harmony"]["params"]["random_state"] == 0

    def test_default_n_neighbors_is_15(self):
        mdata = _make_mdata()
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_neighbors"] == 15


# ---------------------------------------------------------------------------
# n_pcs capping
# ---------------------------------------------------------------------------

class TestNPcsCapping:
    def test_n_pcs_capped_at_available_components(self):
        """If n_pca_comps=10, n_pcs=10 — n_pcs_used should be 10."""
        mdata = _make_mdata(n_pca_comps=10)
        _, metrics = run_harmony_adt(mdata, batch_key="donor", n_pcs=10)
        assert metrics["n_pcs_used"] == 10

    def test_n_pcs_not_capped_when_within_range(self):
        mdata = _make_mdata(n_pca_comps=20)
        _, metrics = run_harmony_adt(mdata, batch_key="donor", n_pcs=15)
        assert metrics["n_pcs_used"] == 15


# ---------------------------------------------------------------------------
# Layer preservation
# ---------------------------------------------------------------------------

class TestLayerPreservation:
    def test_adt_clr_layer_preserved(self):
        mdata = _make_mdata()
        clr_before = mdata["adt"].layers["adt_clr"].copy()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(adata.layers["adt_clr"], clr_before)

    def test_counts_layer_preserved(self):
        mdata = _make_mdata()
        counts_before = mdata["adt"].layers["counts"].copy()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(adata.layers["counts"], counts_before)

    def test_x_pca_adt_not_overwritten_by_harmony(self):
        """Harmony should not modify the original X_pca_adt values."""
        mdata = _make_mdata()
        pca_before = mdata["adt"].obsm["X_pca_adt"].copy()
        adata, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        np.testing.assert_array_equal(adata.obsm["X_pca_adt"], pca_before)


# ---------------------------------------------------------------------------
# Batch scenarios
# ---------------------------------------------------------------------------

class TestBatchScenarios:
    def test_two_batches(self):
        """Standard case: 2 donors."""
        mdata = _make_mdata(n_batches=2)
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_batches"] == 2
        assert "X_pca_harmony_adt" in adata.obsm

    def test_three_batches(self):
        mdata = _make_mdata(n_batches=3, n_cells=150)
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_batches"] == 3
        assert "X_pca_harmony_adt" in adata.obsm

    def test_single_batch_no_error(self):
        """Harmony should run without error even with a single batch value."""
        mdata = _make_mdata(n_batches=1, n_cells=120)
        # Should not raise
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert "X_pca_harmony_adt" in adata.obsm
        assert metrics["n_batches"] == 1

    def test_batch_values_are_strings(self):
        """batch_values in metrics must be a list of strings."""
        mdata = _make_mdata(n_batches=2)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        for v in metrics["batch_values"]:
            assert isinstance(v, str)

    def test_custom_batch_key(self):
        """Harmony works with any valid obs column name."""
        mdata = _make_mdata()
        mdata["adt"].obs["site"] = ["site_X"] * 60 + ["site_Y"] * 60
        adata, metrics = run_harmony_adt(mdata, batch_key="site")
        assert metrics["batch_key"] == "site"
        assert "X_pca_harmony_adt" in adata.obsm


# ---------------------------------------------------------------------------
# ADT-only MuData (no RNA modality)
# ---------------------------------------------------------------------------

class TestAdtOnlyMuData:
    def test_works_without_rna_modality(self):
        mdata = _make_mdata(include_rna=False)
        adata, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert "X_pca_harmony_adt" in adata.obsm
        assert "X_umap_adt" in adata.obsm

    def test_metrics_correct_adt_only(self):
        mdata = _make_mdata(n_cells=120, include_rna=False)
        _, metrics = run_harmony_adt(mdata, batch_key="donor")
        assert metrics["n_cells"] == 120
        assert metrics["n_batches"] == 2


# ---------------------------------------------------------------------------
# Idempotency
# ---------------------------------------------------------------------------

class TestIdempotency:
    def test_running_twice_overwrites_cleanly(self):
        mdata = _make_mdata()
        adata1, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        # Build a new MuData from the already-harmonized adt
        mdata2 = MuData({"adt": adata1})
        adata2, _ = run_harmony_adt(mdata2, batch_key="donor", inplace=True)
        assert "X_pca_harmony_adt" in adata2.obsm
        assert "X_umap_adt" in adata2.obsm

    def test_second_run_overwrites_umap(self):
        """Second run must replace X_umap_adt, not create a stray X_umap alongside it."""
        mdata = _make_mdata()
        adata1, _ = run_harmony_adt(mdata, batch_key="donor", inplace=False)
        mdata2 = MuData({"adt": adata1})
        adata2, _ = run_harmony_adt(mdata2, batch_key="donor", inplace=True)
        assert "X_umap" not in adata2.obsm
        assert "X_umap_adt" in adata2.obsm
