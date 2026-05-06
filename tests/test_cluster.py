"""
OmicSage — Leiden Clustering Tests
tests/test_cluster.py

Run with:
    conda activate omicsage
    python -m pytest tests/test_cluster.py -v

All tests use lightweight synthetic AnnData fixtures so they run in seconds
with no disk I/O.  Fixtures are pre-reduced (X_pca in obsm, connectivities
in obsp) to match the expected input contract of cluster().

Fixture design
--------------
- _make_reduced_adata()       : 200 cells × 500 genes, PCA + neighbor graph
                                pre-computed.  Used for most tests.
- _make_reduced_adata_large() : 600 cells × 500 genes, PCA + neighbor graph
                                pre-computed.  Used for silhouette subsample
                                test (n_obs > _SILHOUETTE_MAX_CELLS not
                                practical here, so large fixture validates
                                silhouette is still computed correctly at
                                scale within reason).

Test inventory (8 tests)
--------------------------
 1. test_leiden_clusters_in_obs
 2. test_all_resolutions_computed
 3. test_cluster_labels_are_strings
 4. test_n_clusters_reasonable
 5. test_silhouette_scores_computed
 6. test_best_resolution_selected
 7. test_params_stored_in_uns
 8. test_original_not_mutated
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scanpy as sc
from anndata import AnnData

# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------
from pipeline.modules.clustering.cluster import cluster, _res_key

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_reduced_adata(
    n_cells: int = 200,
    n_genes: int = 500,
    n_hvg: int = 100,
    n_comps: int = 20,
    n_neighbors: int = 10,
    seed: int = 42,
) -> AnnData:
    """
    Return a minimal AnnData that mimics the output of reduce().

    - .X contains log1p-normalized float values
    - .var['highly_variable'] flags the first n_hvg genes
    - obsm['X_pca'] is pre-computed via scanpy
    - obsp['connectivities'] and obsp['distances'] are pre-computed via scanpy
    """
    rng = np.random.default_rng(seed)

    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=x.shape) < 0.60
    x[mask] = 0.0

    hvg_flag = np.zeros(n_genes, dtype=bool)
    hvg_flag[:n_hvg] = True

    adata = AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={
            "gene_id": [f"gene_{j}" for j in range(n_genes)],
            "highly_variable": hvg_flag,
        },
    )

    # Compute PCA + neighbor graph so the fixture matches reduce() output
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=seed)

    return adata


def _make_reduced_adata_large(
    n_cells: int = 600,
    n_genes: int = 500,
    n_hvg: int = 150,
    n_comps: int = 20,
    n_neighbors: int = 10,
    seed: int = 7,
) -> AnnData:
    """
    Larger fixture for tests that need more cells (e.g. silhouette subsample
    path validation, n_clusters reasonableness at higher cell count).
    """
    rng = np.random.default_rng(seed)

    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=x.shape) < 0.55
    x[mask] = 0.0

    hvg_flag = np.zeros(n_genes, dtype=bool)
    hvg_flag[:n_hvg] = True

    adata = AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={
            "gene_id": [f"gene_{j}" for j in range(n_genes)],
            "highly_variable": hvg_flag,
        },
    )

    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=seed)

    return adata


# ---------------------------------------------------------------------------
# Test 1
# ---------------------------------------------------------------------------

def test_leiden_clusters_in_obs():
    """
    obs must contain a 'leiden_{res}' column for each requested resolution,
    and obs['leiden'] (convenience key) must also be present.
    """
    resolutions = [0.3, 0.5]
    adata = _make_reduced_adata()

    adata_clustered, _ = cluster(adata, resolution_range=resolutions)

    for res in resolutions:
        key = _res_key(res)
        assert key in adata_clustered.obs.columns, (
            f"obs['{key}'] missing after cluster() with resolution={res}"
        )

    assert "leiden" in adata_clustered.obs.columns, (
        "obs['leiden'] convenience key missing — best resolution was not copied"
    )


# ---------------------------------------------------------------------------
# Test 2
# ---------------------------------------------------------------------------

def test_all_resolutions_computed():
    """
    Every resolution in resolution_range must produce a corresponding obs
    column.  Duplicates must be deduplicated (only one column per unique res).
    """
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
    adata = _make_reduced_adata()

    adata_clustered, metrics = cluster(adata, resolution_range=resolutions)

    assert set(metrics["resolutions"]) == set(resolutions), (
        f"metrics['resolutions'] = {metrics['resolutions']}, "
        f"expected {resolutions}"
    )

    for res in resolutions:
        key = _res_key(res)
        assert key in adata_clustered.obs.columns, (
            f"obs['{key}'] missing — resolution {res} was not clustered"
        )

    # Duplicate resolution must not produce two columns
    adata_dup, metrics_dup = cluster(
        adata, resolution_range=[0.5, 0.5, 0.5]
    )
    assert metrics_dup["resolutions"] == [0.5], (
        "Duplicate resolutions were not deduplicated"
    )


# ---------------------------------------------------------------------------
# Test 3
# ---------------------------------------------------------------------------

def test_cluster_labels_are_strings():
    """
    Leiden cluster labels must be strings — this is the Scanpy convention
    and downstream steps (annotation, plotting) rely on it.
    """
    adata = _make_reduced_adata()

    adata_clustered, _ = cluster(adata, resolution_range=[0.5])

    labels = adata_clustered.obs[_res_key(0.5)].values

    assert all(isinstance(lbl, str) for lbl in labels), (
        f"Cluster labels are not all strings — "
        f"got types: {set(type(l).__name__ for l in labels)}"
    )


# ---------------------------------------------------------------------------
# Test 4
# ---------------------------------------------------------------------------

def test_n_clusters_reasonable():
    """
    The number of clusters at each resolution must be between 2 and n_cells.
    At higher resolutions, n_clusters should be >= n_clusters at lower
    resolutions (monotonicity is expected but not strict — we just check
    the bounds here).
    """
    # Use resolutions high enough to produce multiple clusters on synthetic data.
    # Low resolutions (e.g. 0.2) can yield 1 cluster on random data with no
    # biological structure — that is valid behaviour, not a bug.
    adata = _make_reduced_adata()
    resolutions = [0.5, 1.0, 1.5]

    adata_clustered, metrics = cluster(adata, resolution_range=resolutions)

    for res in resolutions:
        n = metrics["n_clusters"][res]
        assert 2 <= n <= adata.n_obs, (
            f"n_clusters at resolution={res} is {n}, "
            f"expected in range [2, {adata.n_obs}]"
        )


# ---------------------------------------------------------------------------
# Test 5
# ---------------------------------------------------------------------------

def test_silhouette_scores_computed():
    """
    metrics['silhouette_scores'] must contain a float entry for every
    requested resolution, and each score must be in [-1.0, 1.0].
    """
    resolutions = [0.2, 0.4, 0.6]
    adata = _make_reduced_adata()

    _, metrics = cluster(adata, resolution_range=resolutions)

    assert "silhouette_scores" in metrics, (
        "metrics['silhouette_scores'] key missing"
    )

    for res in resolutions:
        assert res in metrics["silhouette_scores"], (
            f"silhouette_scores missing entry for resolution={res}"
        )
        score = metrics["silhouette_scores"][res]
        assert isinstance(score, float), (
            f"silhouette_scores[{res}] is {type(score).__name__}, expected float"
        )
        assert -1.0 <= score <= 1.0, (
            f"silhouette_scores[{res}] = {score:.4f} is outside [-1, 1]"
        )


# ---------------------------------------------------------------------------
# Test 6
# ---------------------------------------------------------------------------

def test_best_resolution_selected():
    """
    uns['omicsage_cluster']['best_resolution'] must be set to the resolution
    with the highest silhouette score, and obs['leiden'] must match the
    labels at that resolution.
    """
    resolutions = [0.2, 0.4, 0.6, 0.8]
    adata = _make_reduced_adata()

    adata_clustered, metrics = cluster(adata, resolution_range=resolutions)

    best_res = metrics["best_resolution"]

    assert best_res in resolutions, (
        f"best_resolution={best_res} is not one of the tested resolutions {resolutions}"
    )

    # best_resolution must be the argmax of silhouette scores
    expected_best = max(
        metrics["silhouette_scores"], key=lambda r: metrics["silhouette_scores"][r]
    )
    assert best_res == expected_best, (
        f"best_resolution={best_res} does not match the resolution with highest "
        f"silhouette score ({expected_best})"
    )

    # obs['leiden'] must be identical to obs[_res_key(best_res)]
    leiden_best = adata_clustered.obs[_res_key(best_res)].values
    leiden_convenience = adata_clustered.obs["leiden"].values
    np.testing.assert_array_equal(
        leiden_convenience,
        leiden_best,
        err_msg=(
            "obs['leiden'] does not match obs[_res_key(best_resolution)] — "
            "convenience key was not copied correctly"
        ),
    )

    # uns must also record best_resolution
    assert adata_clustered.uns["omicsage_cluster"]["best_resolution"] == best_res


# ---------------------------------------------------------------------------
# Test 7
# ---------------------------------------------------------------------------

def test_params_stored_in_uns():
    """
    uns['omicsage_cluster'] must contain all required provenance keys with
    correct types and values.
    """
    resolutions = [0.2, 0.5, 1.0]
    adata = _make_reduced_adata()

    adata_clustered, metrics = cluster(
        adata,
        resolution_range=resolutions,
        random_state=0,
    )

    record = adata_clustered.uns["omicsage_cluster"]

    required_keys = [
        "resolutions",
        "n_clusters",
        "silhouette_scores",
        "best_resolution",
        "best_silhouette",
        "best_n_clusters",
        "pca_key",
        "connectivities_key",
        "random_state",
        "scanpy_version",
        "omicsage_module",
        "omicsage_version",
        "timestamp",
    ]

    for key in required_keys:
        assert key in record, f"uns['omicsage_cluster'] missing key: '{key}'"

    assert record["random_state"] == 0
    assert record["pca_key"] == "X_pca"
    assert record["connectivities_key"] == "connectivities"
    assert record["omicsage_module"] == "pipeline.modules.qc.cluster"
    assert record["omicsage_version"] == "0.1.0"

    # timestamp must be a non-empty string
    assert isinstance(record["timestamp"], str) and len(record["timestamp"]) > 0

    # scanpy_version must be a non-empty string
    assert isinstance(record["scanpy_version"], str) and len(record["scanpy_version"]) > 0

    # n_clusters must be a dict with an entry per resolution.
    # uns stores string keys (e.g. '0.2') for h5ad/HDF5 compatibility —
    # the metrics dict returned to the caller still uses float keys.
    assert isinstance(record["n_clusters"], dict)
    for res in resolutions:
        assert str(res) in record["n_clusters"], (
            f"n_clusters missing entry for resolution={res} "
            f"(uns uses string keys, got: {list(record['n_clusters'].keys())})"
        )


# ---------------------------------------------------------------------------
# Test 8
# ---------------------------------------------------------------------------

def test_original_not_mutated():
    """
    By default (inplace=False) the caller's AnnData must not be modified.
    - obs must not contain any new 'leiden_*' columns
    - obs must not contain 'leiden'
    - uns must not contain 'omicsage_cluster'
    """
    adata = _make_reduced_adata()
    original_x = adata.X.copy()
    original_obs_cols = set(adata.obs.columns)
    original_uns_keys = set(adata.uns.keys())

    # Run cluster — should work on a copy
    _ = cluster(adata, inplace=False)

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

    # Check no new obs columns
    new_obs_cols = set(adata.obs.columns) - original_obs_cols
    assert len(new_obs_cols) == 0, (
        f"New obs column(s) added to original adata: {new_obs_cols}"
    )

    # Check no new uns keys
    new_uns_keys = set(adata.uns.keys()) - original_uns_keys
    assert len(new_uns_keys) == 0, (
        f"New uns key(s) added to original adata: {new_uns_keys}"
    )
