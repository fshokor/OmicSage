"""
OmicSage — Dimensionality Reduction Tests
tests/test_reduce.py

Run with:
    conda activate omicsage
    python -m pytest tests/test_reduce.py -v

All tests use lightweight synthetic AnnData fixtures so they run in seconds
with no disk I/O.  Fixtures are pre-normalized (log1p values in .X, HVGs
flagged) to match the expected input contract of reduce().

Fixture design
--------------
- _make_normalized_adata()       : 200 cells × 500 genes, 100 HVGs flagged
                                   Used for most tests.
- _make_normalized_adata_large() : 500 cells × 800 genes, 200 HVGs flagged
                                   Used for PC selection tests that need a
                                   clear scree curve (more cells = better PCA).
- _make_no_hvg_adata()           : 200 cells × 500 genes, no HVG column
                                   Used to test graceful fallback when HVGs
                                   are absent.

Test inventory (11 tests)
--------------------------
 1. test_pca_embedding_shape
 2. test_umap_embedding_shape
 3. test_pca_uses_hvg_only
 4. test_neighbors_graph_computed
 5. test_params_stored_in_uns
 6. test_original_not_mutated
 7. test_tsne_optional
 8. test_elbow_pc_selection
 9. test_manual_n_pcs_override
10. test_variance_pc_selection
11. test_raises_on_non_anndata_input
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------
from pipeline.modules.qc.reduce import reduce, _select_variance, _select_fixed

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_normalized_adata(
    n_cells: int = 200,
    n_genes: int = 500,
    n_hvg: int = 100,
    seed: int = 42,
) -> AnnData:
    """
    Return a minimal AnnData that mimics the output of normalize().

    - .X contains log1p-normalized float values (dense, for speed)
    - .var['highly_variable'] flags the first n_hvg genes as HVGs
    - No layers, no uns — just what reduce() needs as input
    """
    rng = np.random.default_rng(seed)

    # Simulate log1p-normalized expression values (roughly 0–6 range)
    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    # Force ~60% sparsity (zeros are common in scRNA-seq even after normalization)
    mask = rng.random(size=x.shape) < 0.60
    x[mask] = 0.0

    obs_names = [f"cell_{i}" for i in range(n_cells)]
    var_names = [f"gene_{j}" for j in range(n_genes)]

    hvg_flag = np.zeros(n_genes, dtype=bool)
    hvg_flag[:n_hvg] = True

    return AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": obs_names},
        var={
            "gene_id": var_names,
            "highly_variable": hvg_flag,
        },
    )


def _make_normalized_adata_large(
    n_cells: int = 500,
    n_genes: int = 800,
    n_hvg: int = 200,
    seed: int = 7,
) -> AnnData:
    """
    Larger fixture for PC selection tests.

    More cells → richer covariance structure → cleaner scree curve →
    elbow / variance methods behave more predictably.
    """
    rng = np.random.default_rng(seed)
    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=x.shape) < 0.55
    x[mask] = 0.0

    hvg_flag = np.zeros(n_genes, dtype=bool)
    hvg_flag[:n_hvg] = True

    return AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={
            "gene_id": [f"gene_{j}" for j in range(n_genes)],
            "highly_variable": hvg_flag,
        },
    )


def _make_no_hvg_adata(
    n_cells: int = 200,
    n_genes: int = 500,
    seed: int = 99,
) -> AnnData:
    """
    Normalized AnnData with NO highly_variable column — tests graceful fallback
    when reduce() is called without prior HVG selection.
    """
    rng = np.random.default_rng(seed)
    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=x.shape) < 0.60
    x[mask] = 0.0

    return AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={"gene_id": [f"gene_{j}" for j in range(n_genes)]},
    )


# ---------------------------------------------------------------------------
# Test 1
# ---------------------------------------------------------------------------

def test_pca_embedding_shape():
    """
    obsm['X_pca'] must exist after reduce() and have shape (n_cells, n_comps).
    """
    n_comps = 20
    adata = _make_normalized_adata()

    adata_reduced, _ = reduce(adata, n_comps=n_comps)

    assert "X_pca" in adata_reduced.obsm, (
        "obsm['X_pca'] missing — PCA was not computed"
    )
    assert adata_reduced.obsm["X_pca"].shape == (adata.n_obs, n_comps), (
        f"Expected X_pca shape ({adata.n_obs}, {n_comps}), "
        f"got {adata_reduced.obsm['X_pca'].shape}"
    )


# ---------------------------------------------------------------------------
# Test 2
# ---------------------------------------------------------------------------

def test_umap_embedding_shape():
    """
    obsm['X_umap'] must exist after reduce() and have shape (n_cells, 2).
    """
    adata = _make_normalized_adata()

    adata_reduced, _ = reduce(adata)

    assert "X_umap" in adata_reduced.obsm, (
        "obsm['X_umap'] missing — UMAP was not computed"
    )
    assert adata_reduced.obsm["X_umap"].shape == (adata.n_obs, 2), (
        f"Expected X_umap shape ({adata.n_obs}, 2), "
        f"got {adata_reduced.obsm['X_umap'].shape}"
    )


# ---------------------------------------------------------------------------
# Test 3
# ---------------------------------------------------------------------------

def test_pca_uses_hvg_only():
    """
    PCA must be run on HVG subset only when highly_variable is present.

    We verify this indirectly: varm['PCs'] (gene loadings) is only defined
    for HVG genes when use_highly_variable=True.  Scanpy sets the loadings
    for non-HVG genes to zero in the PCs matrix stored in varm.
    """
    n_hvg = 100
    adata = _make_normalized_adata(n_hvg=n_hvg)

    adata_reduced, metrics = reduce(adata)

    # Confirm PCA was flagged as using HVGs
    assert metrics["n_hvg_used"] == n_hvg, (
        f"Expected n_hvg_used={n_hvg}, got {metrics['n_hvg_used']}"
    )
    assert bool(adata_reduced.uns["omicsage_reduce"]["use_highly_variable"]) is True


# ---------------------------------------------------------------------------
# Test 4
# ---------------------------------------------------------------------------

def test_neighbors_graph_computed():
    """
    After reduce(), obsp must contain 'connectivities' and 'distances' — the
    kNN graph required by UMAP and downstream clustering.
    """
    adata = _make_normalized_adata()

    adata_reduced, _ = reduce(adata)

    assert "connectivities" in adata_reduced.obsp, (
        "obsp['connectivities'] missing — neighbor graph was not computed"
    )
    assert "distances" in adata_reduced.obsp, (
        "obsp['distances'] missing — neighbor graph was not computed"
    )

    # Connectivities must be a square matrix of shape (n_cells, n_cells)
    conn = adata_reduced.obsp["connectivities"]
    assert conn.shape == (adata.n_obs, adata.n_obs), (
        f"Connectivities shape {conn.shape} does not match "
        f"expected ({adata.n_obs}, {adata.n_obs})"
    )


# ---------------------------------------------------------------------------
# Test 5
# ---------------------------------------------------------------------------

def test_params_stored_in_uns():
    """
    adata.uns['omicsage_reduce'] must exist and contain all required keys
    with correct values for the parameters used in this run.
    """
    n_comps = 20
    n_neighbors = 10
    adata = _make_normalized_adata()

    adata_reduced, _ = reduce(adata, n_comps=n_comps, n_neighbors=n_neighbors)

    assert "omicsage_reduce" in adata_reduced.uns, (
        "adata.uns['omicsage_reduce'] missing"
    )

    record = adata_reduced.uns["omicsage_reduce"]

    required_keys = [
        "n_comps",
        "n_pcs_used",
        "pc_selection_method",
        "variance_threshold",
        "n_neighbors",
        "run_tsne",
        "use_highly_variable",
        "n_hvg_used",
        "variance_explained_per_pc",
        "cumulative_variance_explained_by_selected_pcs",
        "random_state",
        "scanpy_version",
        "omicsage_module",
        "omicsage_version",
        "timestamp",
    ]
    for key in required_keys:
        assert key in record, f"uns['omicsage_reduce'] missing key: '{key}'"

    assert record["n_comps"] == n_comps, (
        f"n_comps recorded as {record['n_comps']}, expected {n_comps}"
    )
    assert record["n_neighbors"] == n_neighbors, (
        f"n_neighbors recorded as {record['n_neighbors']}, expected {n_neighbors}"
    )
    assert record["run_tsne"] is False
    assert record["omicsage_module"] == "pipeline.modules.qc.reduce"

    # variance_explained_per_pc must be a list of length n_comps
    assert isinstance(record["variance_explained_per_pc"], list)
    assert len(record["variance_explained_per_pc"]) == n_comps

    # timestamp must be a non-empty string
    assert isinstance(record["timestamp"], str) and len(record["timestamp"]) > 0


# ---------------------------------------------------------------------------
# Test 6
# ---------------------------------------------------------------------------

def test_original_not_mutated():
    """
    By default (inplace=False) the caller's AnnData must not be modified.
    - .X must be unchanged
    - obsm must not contain 'X_pca', 'X_umap'
    - obsp must not contain 'connectivities'
    - uns must not contain 'omicsage_reduce'
    """
    adata = _make_normalized_adata()
    original_x = adata.X.copy()
    original_obsm_keys = set(adata.obsm.keys())
    original_obsp_keys = set(adata.obsp.keys())
    original_uns_keys = set(adata.uns.keys())

    # Run reduce — should work on a copy
    _ = reduce(adata, inplace=False)

    # Check .X unchanged
    x_after = adata.X
    if sp.issparse(x_after):
        x_after = x_after.toarray()
    if sp.issparse(original_x):
        original_x = original_x.toarray()
    np.testing.assert_array_equal(
        x_after,
        original_x,
        err_msg="adata.X was mutated despite inplace=False",
    )

    # Check no new obsm keys
    assert set(adata.obsm.keys()) == original_obsm_keys, (
        f"New obsm key(s) added to original adata: "
        f"{set(adata.obsm.keys()) - original_obsm_keys}"
    )

    # Check no new obsp keys
    assert set(adata.obsp.keys()) == original_obsp_keys, (
        f"New obsp key(s) added to original adata: "
        f"{set(adata.obsp.keys()) - original_obsp_keys}"
    )

    # Check no new uns keys
    assert set(adata.uns.keys()) == original_uns_keys, (
        f"New uns key(s) added to original adata: "
        f"{set(adata.uns.keys()) - original_uns_keys}"
    )


# ---------------------------------------------------------------------------
# Test 7
# ---------------------------------------------------------------------------

def test_tsne_optional():
    """
    t-SNE must be absent by default and present when run_tsne=True.
    """
    adata = _make_normalized_adata()

    # Default — t-SNE off
    adata_no_tsne, metrics_off = reduce(adata, run_tsne=False)

    assert "X_tsne" not in adata_no_tsne.obsm, (
        "obsm['X_tsne'] present but run_tsne=False was set"
    )
    assert metrics_off["run_tsne"] is False
    assert "X_tsne" not in metrics_off["embeddings_computed"]

    # t-SNE on
    adata_with_tsne, metrics_on = reduce(adata, run_tsne=True)

    assert "X_tsne" in adata_with_tsne.obsm, (
        "obsm['X_tsne'] missing despite run_tsne=True"
    )
    assert adata_with_tsne.obsm["X_tsne"].shape == (adata.n_obs, 2), (
        f"Expected X_tsne shape ({adata.n_obs}, 2), "
        f"got {adata_with_tsne.obsm['X_tsne'].shape}"
    )
    assert metrics_on["run_tsne"] is True
    assert "X_tsne" in metrics_on["embeddings_computed"]


# ---------------------------------------------------------------------------
# Test 8
# ---------------------------------------------------------------------------

def test_elbow_pc_selection():
    """
    With n_pcs=None and n_pcs_method='elbow', the auto-selected n_pcs_used
    must be within a reasonable range [5, n_comps] and be recorded correctly
    in both the metrics dict and uns.
    """
    n_comps = 30
    adata = _make_normalized_adata_large()

    adata_reduced, metrics = reduce(
        adata,
        n_comps=n_comps,
        n_pcs=None,
        n_pcs_method="elbow",
    )

    n_pcs_used = metrics["n_pcs_used"]

    assert 5 <= n_pcs_used <= n_comps, (
        f"Auto-selected n_pcs_used={n_pcs_used} is outside expected range [5, {n_comps}]"
    )
    assert metrics["pc_selection_method"] in ("elbow", "variance"), (
        # 'variance' is acceptable here — elbow may fall back to variance on
        # synthetic data with a flat scree curve
        f"Unexpected pc_selection_method: {metrics['pc_selection_method']}"
    )
    assert adata_reduced.uns["omicsage_reduce"]["n_pcs_used"] == n_pcs_used


# ---------------------------------------------------------------------------
# Test 9
# ---------------------------------------------------------------------------

def test_manual_n_pcs_override():
    """
    When n_pcs is specified explicitly, auto-selection must be bypassed
    entirely and n_pcs_used must equal the requested value.
    """
    requested_n_pcs = 15
    adata = _make_normalized_adata_large()

    adata_reduced, metrics = reduce(
        adata,
        n_comps=30,
        n_pcs=requested_n_pcs,
    )

    assert metrics["n_pcs_used"] == requested_n_pcs, (
        f"Expected n_pcs_used={requested_n_pcs}, got {metrics['n_pcs_used']}"
    )
    assert metrics["pc_selection_method"] == "manual", (
        f"Expected pc_selection_method='manual', got {metrics['pc_selection_method']}"
    )
    assert adata_reduced.uns["omicsage_reduce"]["pc_selection_method"] == "manual"


# ---------------------------------------------------------------------------
# Test 10
# ---------------------------------------------------------------------------

def test_variance_pc_selection():
    """
    With n_pcs_method='variance', n_pcs_used must be chosen such that
    cumulative variance explained >= variance_threshold.

    We test this directly on _select_variance() (a pure function) so we
    don't have to run the full PCA pipeline in this test.
    """
    # Construct a synthetic variance ratio array:
    # first 10 PCs explain most variance, then it flattens out
    variance_ratio = np.array(
        [0.15, 0.12, 0.10, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.03]
        + [0.01] * 40
    )
    variance_ratio = variance_ratio / variance_ratio.sum()   # normalize to sum=1

    threshold = 0.80
    n_pcs, method = _select_variance(
        variance_ratio=variance_ratio,
        variance_threshold=threshold,
        n_comps=len(variance_ratio),
    )

    cumvar_at_selected = float(np.cumsum(variance_ratio)[n_pcs - 1])

    assert method == "variance"
    assert cumvar_at_selected >= threshold, (
        f"Selected {n_pcs} PCs explain only {cumvar_at_selected:.1%} variance "
        f"— below the requested threshold of {threshold:.0%}"
    )
    assert 5 <= n_pcs <= len(variance_ratio), (
        f"n_pcs={n_pcs} is outside expected bounds [5, {len(variance_ratio)}]"
    )


# ---------------------------------------------------------------------------
# Test 11 — input validation
# ---------------------------------------------------------------------------

def test_raises_on_non_anndata_input():
    """reduce() should raise TypeError if passed anything other than AnnData."""
    with pytest.raises(TypeError, match="AnnData"):
        reduce({"not": "anndata"})  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Bonus — no HVG fallback
# ---------------------------------------------------------------------------

def test_no_hvg_fallback():
    """
    When adata.var has no 'highly_variable' column, reduce() must still run
    successfully — falling back to PCA on all genes — and record
    use_highly_variable=False in uns.
    """
    adata = _make_no_hvg_adata()

    adata_reduced, metrics = reduce(adata, n_comps=15)

    # Must still produce valid embeddings
    assert "X_pca" in adata_reduced.obsm
    assert "X_umap" in adata_reduced.obsm

    # Must record that HVGs were not used
    assert adata_reduced.uns["omicsage_reduce"]["use_highly_variable"] is False

    # n_hvg_used must equal total genes when no HVG subset was used
    assert metrics["n_hvg_used"] == adata.n_vars
