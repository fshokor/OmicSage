"""
test_harmony.py — Tests for harmony_correct.py

Run with:
    cd ~/OmicSage
    conda activate omicsage
    python -m pytest tests/test_harmony.py -v
"""

import numpy as np
import pytest
import scanpy as sc
import anndata as ad


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_adata(n_cells: int = 300, n_genes: int = 200,
                n_pcs: int = 30, n_batches: int = 3,
                batch_key: str = "batch") -> "AnnData":
    """Minimal AnnData with X_pca and a batch column."""
    rng = np.random.default_rng(42)
    X = rng.poisson(1.5, size=(n_cells, n_genes)).astype(np.float32)
    adata = ad.AnnData(X)

    # Fake PCA embedding (each batch shifted to create a real batch effect)
    pca = rng.standard_normal((n_cells, n_pcs)).astype(np.float32)
    batch_labels = np.repeat(
        [f"batch{i}" for i in range(n_batches)],
        n_cells // n_batches,
    )
    # pad if n_cells not divisible
    if len(batch_labels) < n_cells:
        batch_labels = np.append(
            batch_labels,
            [f"batch{n_batches - 1}"] * (n_cells - len(batch_labels)),
        )

    for i, label in enumerate(np.unique(batch_labels)):
        mask = batch_labels == label
        pca[mask] += i * 5.0   # artificial batch offset

    adata.obsm["X_pca"] = pca
    adata.obs[batch_key] = batch_labels
    adata.obs[batch_key] = adata.obs[batch_key].astype("category")
    return adata


@pytest.fixture
def adata_basic():
    return _make_adata()


@pytest.fixture
def adata_two_batch():
    return _make_adata(n_batches=2)


@pytest.fixture
def adata_with_umap():
    """AnnData that already has X_umap so preservation logic is exercised."""
    adata = _make_adata()
    rng = np.random.default_rng(99)
    adata.obsm["X_umap"] = rng.standard_normal((adata.n_obs, 2)).astype(np.float32)
    return adata


# ---------------------------------------------------------------------------
# Import guard
# ---------------------------------------------------------------------------

harmony = pytest.importorskip(
    "harmonypy",
    reason="harmonypy not installed — skipping harmony tests",
)

from pipeline.modules.integration.harmony_correct import (  # noqa: E402
    harmony_correct,
    _validate_inputs,
)


# ---------------------------------------------------------------------------
# Test 1: corrected embedding is created
# ---------------------------------------------------------------------------

def test_harmony_embedding_created(adata_basic):
    out = harmony_correct(adata_basic, batch_key="batch", copy=True)
    assert "X_pca_harmony" in out.obsm, "X_pca_harmony not stored in obsm"


# ---------------------------------------------------------------------------
# Test 2: corrected embedding shape matches input PCA shape
# ---------------------------------------------------------------------------

def test_harmony_embedding_shape(adata_basic):
    out = harmony_correct(adata_basic, batch_key="batch", copy=True)
    n_cells = adata_basic.n_obs
    n_pcs   = min(30, adata_basic.obsm["X_pca"].shape[1])
    expected = (n_cells, n_pcs)
    got      = out.obsm["X_pca_harmony"].shape
    assert got == expected, f"X_pca_harmony shape mismatch: got {got}, expected {expected}"


# ---------------------------------------------------------------------------
# Test 3: UMAP is recomputed
# ---------------------------------------------------------------------------

def test_umap_recomputed(adata_with_umap):
    """X_umap_harmony is created and the original X_umap is preserved."""
    original_umap = adata_with_umap.obsm["X_umap"].copy()
    out = harmony_correct(adata_with_umap, batch_key="batch", copy=True)

    # New harmony UMAP must exist and be 2-D
    assert "X_umap_harmony" in out.obsm, "X_umap_harmony not found after harmony correction"
    assert out.obsm["X_umap_harmony"].shape == (adata_with_umap.n_obs, 2)

    # Pre-correction UMAP must be preserved and unchanged
    assert "X_umap_precorrection" in out.obsm, "X_umap_precorrection not preserved"
    assert out.obsm["X_umap_precorrection"].shape == (adata_with_umap.n_obs, 2)
    np.testing.assert_array_equal(
        out.obsm["X_umap_precorrection"], original_umap,
        err_msg="X_umap_precorrection was mutated — should be an exact copy of the original X_umap",
    )


# ---------------------------------------------------------------------------
# Test 4: neighbor graph key stored
# ---------------------------------------------------------------------------

def test_neighbors_key_stored(adata_basic):
    out = harmony_correct(adata_basic, batch_key="batch", copy=True)
    assert "neighbors_harmony" in out.uns, \
        "neighbors_harmony not found in uns"


# ---------------------------------------------------------------------------
# Test 5: provenance written correctly
# ---------------------------------------------------------------------------

def test_provenance_written(adata_basic):
    out = harmony_correct(adata_basic, batch_key="batch", copy=True)
    prov = out.uns.get("omicsage_harmony", {})
    for key in ["batch_key", "n_batches", "batch_sizes", "n_pcs",
                "n_neighbors", "elapsed_seconds", "embedding_key",
                "umap_key", "umap_precorrection_key"]:
        assert key in prov, f"Provenance missing key: {key}"
    assert prov["batch_key"] == "batch"
    assert prov["n_batches"] == 3
    assert prov["embedding_key"] == "X_pca_harmony"
    assert prov["umap_key"] == "X_umap_harmony"
    assert prov["umap_precorrection_key"] == "X_umap_precorrection"


# ---------------------------------------------------------------------------
# Test 6: in-place modification (copy=False)
# ---------------------------------------------------------------------------

def test_inplace_modification(adata_basic):
    original_id = id(adata_basic)
    out = harmony_correct(adata_basic, batch_key="batch", copy=False)
    assert id(out) == original_id, "copy=False should modify adata in place"
    assert "X_pca_harmony" in adata_basic.obsm


# ---------------------------------------------------------------------------
# Test 7: n_pcs capping when fewer PCs exist
# ---------------------------------------------------------------------------

def test_n_pcs_capping():
    adata = _make_adata(n_pcs=10)  # only 10 PCs in X_pca
    out = harmony_correct(adata, batch_key="batch", n_pcs=50, copy=True)
    # Harmony caps at available PCs (10), not requested (50)
    assert out.obsm["X_pca_harmony"].shape == (adata.n_obs, 10), \
        f"Expected ({adata.n_obs}, 10), got {out.obsm['X_pca_harmony'].shape}"


# ---------------------------------------------------------------------------
# Test 8: missing X_pca raises ValueError
# ---------------------------------------------------------------------------

def test_missing_pca_raises():
    adata = _make_adata()
    del adata.obsm["X_pca"]
    with pytest.raises(ValueError, match="X_pca"):
        harmony_correct(adata, batch_key="batch")


# ---------------------------------------------------------------------------
# Test 9: missing batch_key raises KeyError
# ---------------------------------------------------------------------------

def test_missing_batch_key_raises(adata_basic):
    with pytest.raises(KeyError, match="nonexistent"):
        harmony_correct(adata_basic, batch_key="nonexistent")


# ---------------------------------------------------------------------------
# Test 10: single batch raises ValueError
# ---------------------------------------------------------------------------

def test_single_batch_raises():
    adata = _make_adata(n_batches=1)
    # Override to have only one batch
    adata.obs["batch"] = "batch0"
    adata.obs["batch"] = adata.obs["batch"].astype("category")
    with pytest.raises(ValueError, match="at least 2 batches"):
        harmony_correct(adata, batch_key="batch")


# ---------------------------------------------------------------------------
# Test 11: custom batch_key column works
# ---------------------------------------------------------------------------

def test_custom_batch_key():
    adata = _make_adata(batch_key="donor")
    out = harmony_correct(adata, batch_key="donor", copy=True)
    assert "X_pca_harmony" in out.obsm
    assert out.uns["omicsage_harmony"]["batch_key"] == "donor"


# ---------------------------------------------------------------------------
# Test 12: two-batch correction works
# ---------------------------------------------------------------------------

def test_two_batch_correction(adata_two_batch):
    out = harmony_correct(adata_two_batch, batch_key="batch", copy=True)
    assert "X_pca_harmony" in out.obsm
    assert out.uns["omicsage_harmony"]["n_batches"] == 2
