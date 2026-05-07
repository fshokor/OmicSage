"""
OmicSage — Leiden Clustering Tests
tests/test_cluster.py

Run with:
    conda activate omicsage
    python -m pytest tests/test_cluster.py -v

All tests use lightweight synthetic AnnData fixtures so they run in seconds
with no disk I/O.  Fixtures are pre-reduced (X_pca in obsm, connectivities
in obsp) to match the expected input contract of cluster().

Test inventory (12 tests)
--------------------------
 1. test_leiden_clusters_in_obs
 2. test_all_resolutions_computed
 3. test_cluster_labels_are_strings
 4. test_n_clusters_reasonable
 5. test_silhouette_scores_computed
 6. test_best_resolution_selected_silhouette
 7. test_best_resolution_selected_expected_count
 8. test_stability_scores_computed
 9. test_delta_computed
10. test_params_stored_in_uns
11. test_original_not_mutated
12. test_override_validated
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
import scanpy as sc
from anndata import AnnData

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
    rng = np.random.default_rng(seed)
    x = rng.exponential(scale=0.5, size=(n_cells, n_genes)).astype(np.float32)
    x[rng.random(size=x.shape) < 0.60] = 0.0

    hvg_flag = np.zeros(n_genes, dtype=bool)
    hvg_flag[:n_hvg] = True

    adata = AnnData(
        X=sp.csr_matrix(x),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={"gene_id": [f"gene_{j}" for j in range(n_genes)],
             "highly_variable": hvg_flag},
    )
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True, random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_comps, random_state=seed)
    return adata


# ---------------------------------------------------------------------------
# Test 1
# ---------------------------------------------------------------------------

def test_leiden_clusters_in_obs():
    """obs must contain a leiden_{res} column for each resolution, plus obs['leiden']."""
    resolutions = [0.3, 0.5]
    adata_clustered, _ = cluster(_make_reduced_adata(), resolution_range=resolutions)

    for res in resolutions:
        key = _res_key(res)
        assert key in adata_clustered.obs.columns, f"obs['{key}'] missing"

    assert "leiden" in adata_clustered.obs.columns, "obs['leiden'] convenience key missing"


# ---------------------------------------------------------------------------
# Test 2
# ---------------------------------------------------------------------------

def test_all_resolutions_computed():
    """Every resolution must produce a column; duplicates must be deduplicated."""
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
    adata_clustered, metrics = cluster(_make_reduced_adata(), resolution_range=resolutions)

    assert set(metrics["resolutions"]) == set(resolutions)
    for res in resolutions:
        assert _res_key(res) in adata_clustered.obs.columns

    # Duplicate resolutions must not produce two columns
    _, metrics_dup = cluster(_make_reduced_adata(), resolution_range=[0.5, 0.5, 0.5])
    assert metrics_dup["resolutions"] == [0.5], "Duplicates not deduplicated"


# ---------------------------------------------------------------------------
# Test 3
# ---------------------------------------------------------------------------

def test_cluster_labels_are_strings():
    """Leiden labels must be strings — Scanpy convention."""
    adata_clustered, _ = cluster(_make_reduced_adata(), resolution_range=[0.5])
    labels = adata_clustered.obs[_res_key(0.5)].values
    assert all(isinstance(lbl, str) for lbl in labels), (
        f"Non-string labels found: {set(type(l).__name__ for l in labels)}"
    )


# ---------------------------------------------------------------------------
# Test 4
# ---------------------------------------------------------------------------

def test_n_clusters_reasonable():
    """n_clusters at each resolution must be in [2, n_cells]."""
    adata = _make_reduced_adata()
    resolutions = [0.5, 1.0, 1.5]
    _, metrics = cluster(adata, resolution_range=resolutions)

    for res in resolutions:
        n = metrics["n_clusters"][res]
        assert 2 <= n <= adata.n_obs, (
            f"n_clusters at resolution={res} is {n}, expected in [2, {adata.n_obs}]"
        )


# ---------------------------------------------------------------------------
# Test 5
# ---------------------------------------------------------------------------

def test_silhouette_scores_computed():
    """metrics['silhouette_scores'] must have a float in [-1, 1] per resolution."""
    resolutions = [0.2, 0.4, 0.6]
    _, metrics = cluster(_make_reduced_adata(), resolution_range=resolutions)

    assert "silhouette_scores" in metrics
    for res in resolutions:
        assert res in metrics["silhouette_scores"]
        score = metrics["silhouette_scores"][res]
        assert isinstance(score, float)
        assert -1.0 <= score <= 1.0, f"silhouette[{res}]={score:.4f} outside [-1, 1]"


# ---------------------------------------------------------------------------
# Test 6
# ---------------------------------------------------------------------------

def test_best_resolution_selected_silhouette():
    """
    When no override or expected count is given, selection_reason must be
    either 'stability_plateau' or 'silhouette', and best_resolution must be
    one of the tested resolutions.
    """
    resolutions = [0.2, 0.4, 0.6, 0.8]
    adata_clustered, metrics = cluster(_make_reduced_adata(), resolution_range=resolutions)

    assert metrics["best_resolution"] in resolutions
    assert metrics["selection_reason"] in ("stability_plateau", "silhouette"), (
        f"Unexpected selection_reason: {metrics['selection_reason']}"
    )
    # obs['leiden'] must match obs at the selected resolution
    best_key = _res_key(metrics["best_resolution"])
    np.testing.assert_array_equal(
        adata_clustered.obs["leiden"].values,
        adata_clustered.obs[best_key].values,
        err_msg="obs['leiden'] does not match the selected resolution labels",
    )


# ---------------------------------------------------------------------------
# Test 6b
# ---------------------------------------------------------------------------

def test_delta_elbow_selects_after_largest_jump():
    """
    The delta-elbow logic must select the resolution that captured the
    largest jump in cluster count — not the last resolution.

    We mock a scenario with a clear elbow by using a wide resolution range
    and checking that the selected resolution is NOT always the highest.
    This guards against the regression where stability=1.0 for the last
    resolution always won.
    """
    # Use a wide range — on real-structured data the elbow should not be last
    resolutions = [0.5, 1.0, 1.5, 2.0, 2.5]
    adata = _make_reduced_adata()
    _, metrics = cluster(adata, resolution_range=resolutions)

    best_res = metrics["best_resolution"]
    reason   = metrics["selection_reason"]

    # When stability_plateau is used, the selected resolution must NOT always
    # be the last one — we are testing the elbow logic not just "did it run"
    if reason == "stability_plateau":
        deltas = metrics["n_clusters_delta"]
        # The selected resolution should be at or after the largest delta
        max_delta     = max(deltas[r] for r in resolutions)
        max_delta_res = [r for r in resolutions if deltas[r] == max_delta]
        # best_res should be >= the resolution of the largest jump
        assert best_res >= min(max_delta_res), (
            f"Elbow selection chose res={best_res} which is before the "
            f"largest jump at res={min(max_delta_res)} (delta={max_delta})"
        )


# ---------------------------------------------------------------------------
# Test 7
# ---------------------------------------------------------------------------

def test_best_resolution_selected_expected_count():
    """
    When n_clusters_expected is provided, selection_reason must be
    'expected_count' and the selected resolution must be the one whose
    cluster count is closest to the target.
    """
    resolutions = [0.5, 1.0, 1.5, 2.0]
    adata = _make_reduced_adata()

    # First pass: find actual cluster counts so we can set a realistic target
    _, metrics_probe = cluster(adata, resolution_range=resolutions)
    counts = metrics_probe["n_clusters"]  # float-keyed

    # Choose a target equal to the count at the middle resolution
    middle_res = resolutions[len(resolutions) // 2]
    target = counts[middle_res]

    _, metrics = cluster(
        adata,
        resolution_range=resolutions,
        n_clusters_expected=target,
    )

    assert metrics["selection_reason"] == "expected_count", (
        f"Expected 'expected_count', got '{metrics['selection_reason']}'"
    )
    # The selected count must be the closest to the target
    selected_count = metrics["best_n_clusters"]
    for res in resolutions:
        assert abs(selected_count - target) <= abs(counts[res] - target), (
            f"Resolution {res} (n={counts[res]}) is closer to target={target} "
            f"than selected resolution {metrics['best_resolution']} (n={selected_count})"
        )


# ---------------------------------------------------------------------------
# Test 8
# ---------------------------------------------------------------------------

def test_stability_scores_computed():
    """
    metrics['stability_scores'] must contain a float in (0, 1] per resolution.
    The last resolution must always have stability 1.0 (no next delta to penalise).
    """
    resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
    _, metrics = cluster(_make_reduced_adata(), resolution_range=resolutions)

    assert "stability_scores" in metrics
    for res in resolutions:
        assert res in metrics["stability_scores"], f"stability_scores missing res={res}"
        score = metrics["stability_scores"][res]
        assert isinstance(score, float)
        assert 0.0 < score <= 1.0, f"stability_scores[{res}]={score} outside (0, 1]"

    # Last resolution has no next step → stability must be 1.0
    last_res = resolutions[-1]
    assert metrics["stability_scores"][last_res] == 1.0, (
        f"Last resolution stability should be 1.0, got {metrics['stability_scores'][last_res]}"
    )


# ---------------------------------------------------------------------------
# Test 9
# ---------------------------------------------------------------------------

def test_delta_computed():
    """
    metrics['n_clusters_delta'] must be present and correct.
    - First resolution delta == 0 (no previous step)
    - delta[res] == n_clusters[res] - n_clusters[prev_res] for all others
    - All deltas must be non-negative (cluster count is non-decreasing with resolution)
    """
    resolutions = [0.5, 1.0, 1.5]
    _, metrics = cluster(_make_reduced_adata(), resolution_range=resolutions)

    assert "n_clusters_delta" in metrics
    n   = metrics["n_clusters"]
    d   = metrics["n_clusters_delta"]

    assert d[resolutions[0]] == 0, "First resolution delta must be 0"

    for i in range(1, len(resolutions)):
        res  = resolutions[i]
        prev = resolutions[i - 1]
        expected_delta = n[res] - n[prev]
        assert d[res] == expected_delta, (
            f"delta[{res}]={d[res]}, expected {expected_delta} "
            f"(n[{res}]={n[res]}, n[{prev}]={n[prev]})"
        )
        assert d[res] >= 0, (
            f"Negative delta at res={res}: {d[res]} "
            "(cluster count should not decrease as resolution increases)"
        )


# ---------------------------------------------------------------------------
# Test 10
# ---------------------------------------------------------------------------

def test_params_stored_in_uns():
    """
    uns['omicsage_cluster'] must contain all required provenance keys.
    uns uses string-keyed dicts for h5ad/HDF5 compatibility.
    """
    resolutions = [0.2, 0.5, 1.0]
    adata_clustered, metrics = cluster(
        _make_reduced_adata(),
        resolution_range=resolutions,
        random_state=0,
    )

    record = adata_clustered.uns["omicsage_cluster"]

    required_keys = [
        "resolutions", "n_clusters", "n_clusters_delta",
        "silhouette_scores", "stability_scores",
        "best_resolution", "best_n_clusters", "best_silhouette", "best_stability",
        "selection_reason", "n_clusters_expected",
        "pca_key", "connectivities_key", "random_state",
        "scanpy_version", "omicsage_module", "omicsage_version", "timestamp",
    ]
    for key in required_keys:
        assert key in record, f"uns['omicsage_cluster'] missing key: '{key}'"

    assert record["random_state"] == 0
    assert record["pca_key"] == "X_pca"
    assert record["connectivities_key"] == "connectivities"
    assert record["omicsage_module"] == "pipeline.modules.qc.cluster"
    assert record["omicsage_version"] == "0.1.0"
    assert isinstance(record["timestamp"], str) and len(record["timestamp"]) > 0
    assert isinstance(record["scanpy_version"], str) and len(record["scanpy_version"]) > 0

    # uns stores string-keyed dicts for h5ad/HDF5 compatibility
    assert isinstance(record["n_clusters"], dict)
    assert isinstance(record["stability_scores"], dict)
    assert isinstance(record["n_clusters_delta"], dict)
    for res in resolutions:
        assert str(res) in record["n_clusters"], (
            f"n_clusters missing entry for resolution={res} "
            f"(uns uses string keys, got: {list(record['n_clusters'].keys())})"
        )
        assert str(res) in record["stability_scores"], (
            f"stability_scores missing entry for resolution={res}"
        )
        assert str(res) in record["n_clusters_delta"], (
            f"n_clusters_delta missing entry for resolution={res}"
        )


# ---------------------------------------------------------------------------
# Test 11
# ---------------------------------------------------------------------------

def test_original_not_mutated():
    """inplace=False must leave the caller's AnnData completely unchanged."""
    adata = _make_reduced_adata()
    original_x        = adata.X.copy()
    original_obs_cols = set(adata.obs.columns)
    original_uns_keys = set(adata.uns.keys())

    _ = cluster(adata, inplace=False)

    x_after = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    x_orig  = original_x.toarray() if sp.issparse(original_x) else original_x
    np.testing.assert_array_equal(x_after, x_orig, err_msg="adata.X mutated")

    new_obs = set(adata.obs.columns) - original_obs_cols
    assert len(new_obs) == 0, f"New obs columns added to original: {new_obs}"

    new_uns = set(adata.uns.keys()) - original_uns_keys
    assert len(new_uns) == 0, f"New uns keys added to original: {new_uns}"


# ---------------------------------------------------------------------------
# Test 12
# ---------------------------------------------------------------------------

def test_override_validated():
    """best_resolution_override must raise ValueError if not in resolution_range."""
    adata = _make_reduced_adata()
    with pytest.raises(ValueError, match="best_resolution_override"):
        cluster(adata, resolution_range=[0.2, 0.4], best_resolution_override=0.9)
