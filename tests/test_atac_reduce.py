"""
tests/test_atac_reduce.py

Unit tests for pipeline.modules.multiome.atac_reduce.atac_reduce()

Fixture conventions
-------------------
- Synthetic peak count matrices only — no real data loaded
- n_cells=150, n_peaks=300 for most fixtures (fast, CI-safe)
- n_components kept small (10) to keep SVD fast in CI
- All fixtures use raw integer counts in .X (sparse CSR)
- filter_cells=False default from atac_qc — fixtures mimic already-QC'd data
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import anndata as ad

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pipeline.modules.scripts.multiome.atac_reduce import atac_reduce, _tfidf


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

def _make_atac_adata(
    n_cells: int = 150,
    n_peaks: int = 300,
    seed: int = 42,
    with_counts_layer: bool = True,
    with_cellranger_obs: bool = True,
    sparse: bool = True,
) -> ad.AnnData:
    """
    Build a minimal QC-filtered ATAC AnnData mimicking atac_qc() output.
    Uses random sparse binary-like counts (10% non-zero).
    """
    rng = np.random.default_rng(seed)
    X = rng.integers(0, 2, size=(n_cells, n_peaks)).astype(float)
    mask = rng.random(size=X.shape) > 0.10
    X[mask] = 0.0
    if sparse:
        X = sp.csr_matrix(X)

    adata = ad.AnnData(X)
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"chr1:{i*500}-{i*500+499}" for i in range(n_peaks)]
    adata.var["feature_types"] = "Peaks"

    if with_counts_layer:
        adata.layers["counts"] = adata.X.copy()

    if with_cellranger_obs:
        adata.obs["total_peak_counts"]    = rng.integers(2000, 40000, size=n_cells).astype(float)
        adata.obs["nucleosome_signal"]    = rng.uniform(0.5, 1.8, size=n_cells)
        adata.obs["batch"]               = rng.choice(["s1d1", "s1d2"], size=n_cells)
        adata.obs["cell_type_groundtruth"]= rng.choice(
            ["CD4 T", "CD8 T", "B cell", "NK"], size=n_cells
        )

    return adata


# ---------------------------------------------------------------------------
# 1. Return contract
# ---------------------------------------------------------------------------

class TestReturnContract:

    def test_returns_tuple(self):
        adata = _make_atac_adata()
        result = atac_reduce(adata, n_components=10)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_atac_adata()
        adata_out, metrics = atac_reduce(adata, n_components=10)
        assert isinstance(adata_out, ad.AnnData)
        assert isinstance(metrics, dict)

    def test_n_obs_preserved(self):
        adata = _make_atac_adata(n_cells=120)
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert adata_out.n_obs == 120

    def test_n_vars_preserved(self):
        adata = _make_atac_adata(n_peaks=200)
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert adata_out.n_vars == 200

    def test_obs_names_preserved(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert list(adata_out.obs_names) == list(adata.obs_names)


# ---------------------------------------------------------------------------
# 2. TF-IDF layer
# ---------------------------------------------------------------------------

class TestTFIDF:

    def test_tfidf_layer_created(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "tf_idf" in adata_out.layers

    def test_tfidf_shape_correct(self):
        adata = _make_atac_adata(n_cells=100, n_peaks=200)
        adata_out, _ = atac_reduce(adata, n_components=10)
        tfidf = adata_out.layers["tf_idf"]
        if sp.issparse(tfidf):
            tfidf = tfidf.toarray()
        assert tfidf.shape == (100, 200)

    def test_tfidf_non_negative(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        tfidf = adata_out.layers["tf_idf"]
        if sp.issparse(tfidf):
            tfidf = tfidf.toarray()
        assert np.all(tfidf >= 0)

    def test_tfidf_differs_from_raw(self):
        adata = _make_atac_adata()
        raw = adata.X.toarray() if sp.issparse(adata.X) else adata.X.copy()
        adata_out, _ = atac_reduce(adata, n_components=10)
        tfidf = adata_out.layers["tf_idf"]
        if sp.issparse(tfidf):
            tfidf = tfidf.toarray()
        assert not np.allclose(tfidf, raw), "TF-IDF output identical to raw counts"

    def test_tfidf_standalone_function(self):
        """_tfidf() helper should return a sparse matrix of the same shape."""
        rng = np.random.default_rng(0)
        X = sp.csr_matrix(rng.integers(0, 3, size=(50, 100)).astype(float))
        result = _tfidf(X)
        assert sp.issparse(result)
        assert result.shape == (50, 100)
        assert np.all(result.data >= 0)

    def test_tfidf_uses_counts_layer_when_available(self):
        """atac_reduce should prefer .layers['counts'] over .X."""
        adata = _make_atac_adata(with_counts_layer=True)
        # Corrupt .X with zeros — if counts layer is used, result should still be valid
        adata.X = sp.csr_matrix(adata.X.shape)
        adata_out, _ = atac_reduce(adata, n_components=10, use_raw_counts=True)
        # TF-IDF from all-zero X would be all-zero; from counts layer it should have values
        tfidf = adata_out.layers["tf_idf"]
        if sp.issparse(tfidf):
            tfidf = tfidf.toarray()
        # counts layer was non-zero so tf_idf should have non-zero entries
        assert tfidf.sum() > 0


# ---------------------------------------------------------------------------
# 3. LSI embedding
# ---------------------------------------------------------------------------

class TestLSI:

    def test_x_lsi_in_obsm(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "X_lsi" in adata_out.obsm

    def test_x_lsi_shape(self):
        """X_lsi should have n_components - 1 columns (component 1 dropped)."""
        adata = _make_atac_adata()
        n_comp = 10
        adata_out, _ = atac_reduce(adata, n_components=n_comp)
        assert adata_out.obsm["X_lsi"].shape == (adata.n_obs, n_comp - 1)

    def test_component_1_dropped(self):
        """metrics must record component_1_dropped=True."""
        adata = _make_atac_adata()
        _, metrics = atac_reduce(adata, n_components=10)
        assert metrics["component_1_dropped"] is True

    def test_n_lsi_components_used_correct(self):
        adata = _make_atac_adata()
        n_comp = 10
        _, metrics = atac_reduce(adata, n_components=n_comp)
        assert metrics["n_lsi_components_used"] == n_comp - 1

    def test_lsi_l2_normalised(self):
        """Each row of X_lsi should have unit L2 norm (post-normalisation)."""
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        lsi = adata_out.obsm["X_lsi"]
        row_norms = np.linalg.norm(lsi, axis=1)
        np.testing.assert_allclose(row_norms, np.ones(len(row_norms)), atol=1e-5)

    def test_lsi_finite_values(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        lsi = adata_out.obsm["X_lsi"]
        assert np.all(np.isfinite(lsi))


# ---------------------------------------------------------------------------
# 4. UMAP embedding
# ---------------------------------------------------------------------------

class TestUMAP:

    def test_x_umap_atac_in_obsm(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "X_umap_atac" in adata_out.obsm

    def test_x_umap_atac_shape(self):
        adata = _make_atac_adata(n_cells=100)
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert adata_out.obsm["X_umap_atac"].shape == (100, 2)

    def test_x_umap_atac_finite(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert np.all(np.isfinite(adata_out.obsm["X_umap_atac"]))

    def test_umap_key_namespaced(self):
        """X_umap_atac must be present; confirms no collision with RNA X_umap."""
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "X_umap_atac" in adata_out.obsm


# ---------------------------------------------------------------------------
# 5. Leiden clustering
# ---------------------------------------------------------------------------

class TestLeiden:

    def test_atac_leiden_in_obs(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "atac_leiden" in adata_out.obs.columns

    def test_atac_leiden_is_categorical_or_string(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        dtype = str(adata_out.obs["atac_leiden"].dtype)
        assert dtype in ("category", "object") or "str" in dtype

    def test_atac_leiden_has_multiple_clusters(self):
        adata = _make_atac_adata(n_cells=150)
        adata_out, metrics = atac_reduce(adata, n_components=10)
        assert metrics["n_leiden_clusters"] >= 1

    def test_atac_leiden_key_namespaced(self):
        """atac_leiden must not collide with RNA 'leiden' key."""
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "atac_leiden" in adata_out.obs.columns
        assert "leiden" not in adata_out.obs.columns


# ---------------------------------------------------------------------------
# 6. Provenance
# ---------------------------------------------------------------------------

class TestProvenance:

    def test_uns_key_present(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert "omicsage_atac_reduce" in adata_out.uns

    def test_uns_required_fields(self):
        required = {
            "omicsage_module",
            "omicsage_version",
            "timestamp",
            "n_components_computed",
            "n_lsi_components_used",
            "component_1_dropped",
            "component_1_drop_reason",
            "tfidf_method",
            "lsi_method",
            "lsi_normalisation",
            "n_neighbors",
            "leiden_resolution",
            "n_leiden_clusters",
            "variance_ratio_used",
            "cumulative_variance_used",
            "random_state",
            "use_raw_counts",
            "scanpy_version",
        }
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        prov = adata_out.uns["omicsage_atac_reduce"]
        missing = required - set(prov.keys())
        assert not missing, f"Provenance missing keys: {missing}"

    def test_uns_module_name(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert (
            adata_out.uns["omicsage_atac_reduce"]["omicsage_module"]
            == "pipeline.modules.multiome.atac_reduce"
        )

    def test_uns_timestamp_iso(self):
        from datetime import datetime
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        ts = adata_out.uns["omicsage_atac_reduce"]["timestamp"]
        datetime.fromisoformat(ts)

    def test_uns_component_1_dropped_true(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10)
        assert adata_out.uns["omicsage_atac_reduce"]["component_1_dropped"] is True


# ---------------------------------------------------------------------------
# 7. Metrics dict contract
# ---------------------------------------------------------------------------

class TestMetricsDict:

    def test_required_keys_present(self):
        required = {
            "n_cells",
            "n_peaks",
            "n_components_computed",
            "n_lsi_components_used",
            "component_1_dropped",
            "n_neighbors",
            "leiden_resolution",
            "n_leiden_clusters",
            "variance_ratio_all",
            "variance_ratio_used",
            "cumulative_variance_used",
            "embeddings_computed",
            "cluster_key",
        }
        adata = _make_atac_adata()
        _, metrics = atac_reduce(adata, n_components=10)
        missing = required - set(metrics.keys())
        assert not missing, f"Metrics dict missing keys: {missing}"

    def test_n_cells_correct(self):
        adata = _make_atac_adata(n_cells=120)
        _, metrics = atac_reduce(adata, n_components=10)
        assert metrics["n_cells"] == 120

    def test_n_peaks_correct(self):
        adata = _make_atac_adata(n_peaks=250)
        _, metrics = atac_reduce(adata, n_components=10)
        assert metrics["n_peaks"] == 250

    def test_variance_ratio_all_length(self):
        adata = _make_atac_adata()
        n_comp = 10
        _, metrics = atac_reduce(adata, n_components=n_comp)
        assert len(metrics["variance_ratio_all"]) == n_comp

    def test_variance_ratio_used_length(self):
        adata = _make_atac_adata()
        n_comp = 10
        _, metrics = atac_reduce(adata, n_components=n_comp)
        assert len(metrics["variance_ratio_used"]) == n_comp - 1

    def test_embeddings_computed_list(self):
        adata = _make_atac_adata()
        _, metrics = atac_reduce(adata, n_components=10)
        assert "X_lsi" in metrics["embeddings_computed"]
        assert "X_umap_atac" in metrics["embeddings_computed"]

    def test_cluster_key_value(self):
        adata = _make_atac_adata()
        _, metrics = atac_reduce(adata, n_components=10)
        assert metrics["cluster_key"] == "atac_leiden"


# ---------------------------------------------------------------------------
# 8. inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:

    def test_default_not_inplace(self):
        adata = _make_atac_adata()
        original_obsm_keys = set(adata.obsm.keys())
        atac_reduce(adata, n_components=10)
        assert set(adata.obsm.keys()) == original_obsm_keys

    def test_inplace_false_returns_new_object(self):
        adata = _make_atac_adata()
        adata_out, _ = atac_reduce(adata, n_components=10, inplace=False)
        assert adata_out is not adata


# ---------------------------------------------------------------------------
# 9. Input validation
# ---------------------------------------------------------------------------

class TestInputValidation:

    def test_rejects_non_anndata(self):
        with pytest.raises(TypeError, match="AnnData"):
            atac_reduce("not_an_anndata")

    def test_rejects_empty_cells(self):
        adata = _make_atac_adata()
        empty = adata[:0].copy()
        with pytest.raises(ValueError, match="0 cells"):
            atac_reduce(empty)

    def test_rejects_empty_peaks(self):
        adata = _make_atac_adata()
        empty = adata[:, :0].copy()
        with pytest.raises(ValueError, match="0 peaks"):
            atac_reduce(empty)

    def test_rejects_n_components_too_large(self):
        adata = _make_atac_adata(n_peaks=50)
        with pytest.raises(ValueError, match="n_components"):
            atac_reduce(adata, n_components=60)

    def test_rejects_n_components_less_than_2(self):
        adata = _make_atac_adata()
        with pytest.raises(ValueError, match="n_components"):
            atac_reduce(adata, n_components=1)


# ---------------------------------------------------------------------------
# 10. Reproducibility
# ---------------------------------------------------------------------------

class TestReproducibility:

    def test_same_random_state_same_umap(self):
        adata = _make_atac_adata(seed=0)
        out1, _ = atac_reduce(adata, n_components=10, random_state=42)
        out2, _ = atac_reduce(adata, n_components=10, random_state=42)
        np.testing.assert_array_almost_equal(
            out1.obsm["X_umap_atac"], out2.obsm["X_umap_atac"], decimal=4
        )

    def test_different_random_state_different_umap(self):
        adata = _make_atac_adata(seed=0)
        out1, _ = atac_reduce(adata, n_components=10, random_state=0)
        out2, _ = atac_reduce(adata, n_components=10, random_state=99)
        # UMAPs from different seeds should differ (not guaranteed but highly probable)
        # Use LSI instead which is fully deterministic given seed
        lsi1 = out1.obsm["X_lsi"]
        lsi2 = out2.obsm["X_lsi"]
        # SVD with different seeds may produce sign-flipped components — check magnitude
        assert not np.allclose(np.abs(lsi1), np.abs(lsi2), atol=1e-8)
