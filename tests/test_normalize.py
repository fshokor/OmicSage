"""
OmicSage — Normalization Tests
tests/test_normalize.py

Run with:
    conda activate omicsage
    python -m pytest tests/test_normalize.py -v

All tests use lightweight synthetic AnnData fixtures so they run in seconds
with no disk I/O.

Fixture design
--------------
- _make_raw_adata()       : 80 cells × 300 genes — used for most tests
- _make_large_adata()     : 150 cells × 400 genes — used for HVG count tests
- _make_seurat_v3_adata() : 500 cells × 600 genes — minimum for seurat_v3
                            stability (its variance model needs enough cells)

HVG flavor strategy
-------------------
- Tests 1-7 use flavor='seurat' (runs post log1p, stable on small fixtures)
- test_seurat_v3_runs_on_real_data uses flavor='seurat_v3' with a large fixture
  This is the correct production flavor for normalize() — the small fixture
  tests just don't have enough cells for its variance model to be stable.

Test inventory (9 tests)
------------------------
1.  test_raw_counts_preserved_in_layer
2.  test_x_is_normalized_after_run
3.  test_log1p_applied
4.  test_hvg_selected
5.  test_hvg_count_correct
6.  test_normalization_params_in_uns
7.  test_original_adata_not_mutated
8.  test_raises_on_non_anndata_input        (bonus)
9.  test_raises_on_already_normalized_input (bonus)
10. test_seurat_v3_runs_on_real_data        (flavor-specific)
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp
from anndata import AnnData

# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------
from pipeline.modules.qc.normalize import normalize


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _make_raw_adata(
    n_cells: int = 80,
    n_genes: int = 300,
    seed: int = 42,
) -> AnnData:
    """
    Return a minimal AnnData with:
    - integer counts in .X  (sparse CSR)
    - cell / gene names
    - no layers, no uns
    """
    rng = np.random.default_rng(seed)
    # Negative binomial counts — realistic sparsity
    counts = rng.negative_binomial(n=2, p=0.9, size=(n_cells, n_genes)).astype(np.float32)
    # Force ~70 % sparsity
    mask = rng.random(size=counts.shape) < 0.70
    counts[mask] = 0.0

    obs_names = [f"cell_{i}" for i in range(n_cells)]
    var_names = [f"gene_{j}" for j in range(n_genes)]

    return AnnData(
        X=sp.csr_matrix(counts),
        obs={"cell_id": obs_names},
        var={"gene_id": var_names},
    )


def _make_large_adata(n_cells: int = 150, n_genes: int = 400, seed: int = 7) -> AnnData:
    """Larger fixture — needed for HVG count accuracy test."""
    rng = np.random.default_rng(seed)
    counts = rng.negative_binomial(n=3, p=0.85, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=counts.shape) < 0.65
    counts[mask] = 0.0

    return AnnData(
        X=sp.csr_matrix(counts),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={"gene_id": [f"gene_{j}" for j in range(n_genes)]},
    )


def _make_seurat_v3_adata(n_cells: int = 500, n_genes: int = 600, seed: int = 99) -> AnnData:
    """
    Large fixture specifically for seurat_v3 HVG flavor.

    seurat_v3 fits a regularized negative binomial variance model per gene.
    With too few cells the model hits near-singularities and raises a ValueError.
    500 cells × 600 genes is safely above the practical stability threshold.
    """
    rng = np.random.default_rng(seed)
    counts = rng.negative_binomial(n=5, p=0.8, size=(n_cells, n_genes)).astype(np.float32)
    mask = rng.random(size=counts.shape) < 0.60
    counts[mask] = 0.0

    return AnnData(
        X=sp.csr_matrix(counts),
        obs={"cell_id": [f"cell_{i}" for i in range(n_cells)]},
        var={"gene_id": [f"gene_{j}" for j in range(n_genes)]},
    )


# ---------------------------------------------------------------------------
# Test 1
# ---------------------------------------------------------------------------

def test_raw_counts_preserved_in_layer():
    """
    Raw integer counts must be saved to adata.layers['counts']
    and must be identical to the original .X before normalization.
    """
    adata = _make_raw_adata()
    original_x = adata.X.copy()

    adata_norm, _ = normalize(adata, hvg_flavor="seurat")

    assert "counts" in adata_norm.layers, (
        "adata.layers['counts'] missing — raw counts were not preserved"
    )

    saved = adata_norm.layers["counts"]
    if sp.issparse(saved):
        saved = saved.toarray()
    if sp.issparse(original_x):
        original_x = original_x.toarray()

    np.testing.assert_array_equal(
        saved,
        original_x,
        err_msg="layers['counts'] does not match the original raw counts",
    )


# ---------------------------------------------------------------------------
# Test 2
# ---------------------------------------------------------------------------

def test_x_is_normalized_after_run():
    """
    After normalization each cell's total count in .X must equal target_sum
    (before log1p, so we check on the layers['counts'] ratio — but we can
    verify indirectly: after log1p the row sums are no longer target_sum,
    so instead we verify .X values differ from raw counts, and that the
    pre-log1p per-cell sums were uniform).

    Approach: run normalize, recover pre-log1p values via expm1, check row sums.
    """
    target_sum = 1e4
    adata = _make_raw_adata()
    adata_norm, _ = normalize(adata, target_sum=target_sum, hvg_flavor="seurat")

    x = adata_norm.X
    if sp.issparse(x):
        x = x.toarray()

    # expm1 reverses log1p → should give per-cell normalized counts
    pre_log = np.expm1(x)
    row_sums = pre_log.sum(axis=1)

    # Every row sum should be close to target_sum
    # (small deviations OK due to float32 precision)
    np.testing.assert_allclose(
        row_sums,
        target_sum,
        rtol=0.01,
        err_msg=(
            "After reversing log1p, per-cell sums should equal target_sum. "
            "Normalization may not have run correctly."
        ),
    )


# ---------------------------------------------------------------------------
# Test 3
# ---------------------------------------------------------------------------

def test_log1p_applied():
    """
    Log1p must have been applied: all values in .X must be ≥ 0,
    and the maximum value must be consistent with log1p of normalized counts
    (i.e., substantially smaller than the raw maximum).
    """
    adata = _make_raw_adata()
    raw_max = float(
        adata.X.toarray().max() if sp.issparse(adata.X) else np.asarray(adata.X).max()
    )

    adata_norm, _ = normalize(adata, hvg_flavor="seurat")
    x = adata_norm.X
    if sp.issparse(x):
        x = x.toarray()

    assert x.min() >= 0, "Negative values found after log1p — unexpected."

    # log1p(target_sum) ≈ 9.21 for target_sum=1e4
    # Raw counts times 1e4 / cell_sum can be large; log1p caps it
    assert x.max() < 20, (
        f"Max value after log1p is {x.max():.2f} — suspiciously large, "
        "log1p may not have been applied."
    )

    # Also verify .uns records log1p_applied = True
    assert adata_norm.uns.get("omicsage_normalization", {}).get("log1p_applied") is True


# ---------------------------------------------------------------------------
# Test 4
# ---------------------------------------------------------------------------

def test_hvg_selected():
    """
    adata.var must contain a boolean 'highly_variable' column after normalization.
    """
    adata = _make_raw_adata()
    adata_norm, _ = normalize(adata, hvg_flavor="seurat")

    assert "highly_variable" in adata_norm.var.columns, (
        "adata.var['highly_variable'] missing — HVG selection did not run"
    )
    assert adata_norm.var["highly_variable"].dtype == bool or \
           adata_norm.var["highly_variable"].dtype == np.bool_, (
        "highly_variable column should be boolean"
    )


# ---------------------------------------------------------------------------
# Test 5
# ---------------------------------------------------------------------------

def test_hvg_count_correct():
    """
    The number of True values in var['highly_variable'] must equal n_top_genes,
    capped at n_vars (cannot select more genes than exist).
    """
    n_top = 200          # deliberately less than n_genes=400
    adata = _make_large_adata()
    adata_norm, metrics = normalize(adata, n_top_genes=n_top, hvg_flavor="seurat")

    n_hvg = int(adata_norm.var["highly_variable"].sum())
    expected = min(n_top, adata_norm.n_vars)

    assert n_hvg == expected, (
        f"Expected {expected} HVGs but got {n_hvg}. "
        "Check n_top_genes handling in normalize()."
    )
    assert metrics["n_hvg_selected"] == n_hvg, (
        "metrics['n_hvg_selected'] does not match actual HVG count in var"
    )


# ---------------------------------------------------------------------------
# Test 6
# ---------------------------------------------------------------------------

def test_normalization_params_in_uns():
    """
    adata.uns['omicsage_normalization'] must exist and contain the key
    parameters used during this run.
    """
    target_sum = 5000.0
    n_top = 150
    adata = _make_large_adata()
    adata_norm, _ = normalize(adata, target_sum=target_sum, n_top_genes=n_top, hvg_flavor="seurat")

    assert "omicsage_normalization" in adata_norm.uns, (
        "adata.uns['omicsage_normalization'] missing"
    )

    record = adata_norm.uns["omicsage_normalization"]

    required_keys = [
        "target_sum",
        "n_top_genes",
        "hvg_flavor",
        "n_hvg_selected",
        "log1p_applied",
        "scanpy_version",
        "omicsage_module",
    ]
    for key in required_keys:
        assert key in record, f"uns['omicsage_normalization'] missing key: '{key}'"

    assert record["target_sum"] == target_sum, (
        f"target_sum recorded as {record['target_sum']}, expected {target_sum}"
    )
    assert record["n_top_genes"] == n_top, (
        f"n_top_genes recorded as {record['n_top_genes']}, expected {n_top}"
    )
    assert record["log1p_applied"] is True


# ---------------------------------------------------------------------------
# Test 7
# ---------------------------------------------------------------------------

def test_original_adata_not_mutated():
    """
    By default (inplace=False) the caller's AnnData must not be modified.
    - .X must still contain raw counts
    - .layers must not contain 'counts'
    - .var must not contain 'highly_variable'
    - .uns must not contain 'omicsage_normalization'
    """
    adata = _make_raw_adata()
    original_x = adata.X.copy()
    original_layers = set(adata.layers.keys())
    original_var_cols = set(adata.var.columns)
    original_uns_keys = set(adata.uns.keys())

    # Run normalize — should work on a copy
    _ = normalize(adata, inplace=False, hvg_flavor="seurat")

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

    # Check no new layers added to original
    assert set(adata.layers.keys()) == original_layers, (
        f"New layer(s) added to original adata: {set(adata.layers.keys()) - original_layers}"
    )

    # Check no new var columns
    assert set(adata.var.columns) == original_var_cols, (
        f"New var column(s) added to original adata: {set(adata.var.columns) - original_var_cols}"
    )

    # Check no new uns keys
    assert set(adata.uns.keys()) == original_uns_keys, (
        f"New uns key(s) added to original adata: {set(adata.uns.keys()) - original_uns_keys}"
    )


# ---------------------------------------------------------------------------
# Bonus — input validation tests
# ---------------------------------------------------------------------------

def test_raises_on_non_anndata_input():
    """normalize() should raise TypeError if passed a MuData or dict."""
    with pytest.raises(TypeError, match="AnnData"):
        normalize({"not": "anndata"})  # type: ignore[arg-type]


def test_raises_on_already_normalized_input():
    """
    normalize() should raise ValueError if .X looks like it was already
    log-normalized (non-integer values).
    """
    adata = _make_raw_adata()
    # Simulate already-normalized data by running log1p manually
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    with pytest.raises(ValueError, match="non-integer"):
        normalize(adata)


# ---------------------------------------------------------------------------
# Test 10 — seurat_v3 flavor specifically (needs large fixture)
# ---------------------------------------------------------------------------

def test_seurat_v3_runs_on_real_data():
    """
    seurat_v3 HVG selection requires enough cells for its variance model
    to be numerically stable (near-singularity errors on <200 cells).
    This test confirms the production default flavor works correctly
    when given an appropriately sized dataset (500 cells × 600 genes).
    """
    adata = _make_seurat_v3_adata()          # 500 × 600
    n_top = 300

    adata_norm, metrics = normalize(
        adata,
        n_top_genes=n_top,
        hvg_flavor="seurat_v3",              # production default
    )

    # HVG count correct
    n_hvg = int(adata_norm.var["highly_variable"].sum())
    assert n_hvg == n_top, f"Expected {n_top} HVGs, got {n_hvg}"

    # Raw counts preserved
    assert "counts" in adata_norm.layers

    # log1p applied
    x = adata_norm.X
    if sp.issparse(x):
        x = x.toarray()
    assert x.max() < 20, "log1p may not have been applied"

    # uns record present
    assert "omicsage_normalization" in adata_norm.uns
    assert adata_norm.uns["omicsage_normalization"]["hvg_flavor"] == "seurat_v3"



# ---------------------------------------------------------------------------
# Test — normalized data saved to layers['normalized']
# ---------------------------------------------------------------------------

def test_normalized_data_in_layer():
    """
    adata.layers['normalized'] must contain the log1p-normalized values
    and must be identical to adata.X after normalization.
    """
    adata = _make_raw_adata()
    adata_norm, metrics = normalize(adata, hvg_flavor="seurat")

    assert "logcounts" in adata_norm.layers, (
        "adata.layers['normalized'] missing — log1p values were not saved to layers"
    )

    x = adata_norm.X
    norm_layer = adata_norm.layers["logcounts"]

    if sp.issparse(x):
        x = x.toarray()
    if sp.issparse(norm_layer):
        norm_layer = norm_layer.toarray()

    np.testing.assert_array_equal(
        norm_layer, x,
        err_msg="layers['normalized'] does not match adata.X"
    )

    # Also check metrics records it
    assert metrics.get("normalized_in_layer") == "logcounts"

# ---------------------------------------------------------------------------
# Test 11 — batch_key is respected
# ---------------------------------------------------------------------------

def test_batch_key_hvg_selection():
    """
    When batch_key is provided, highly_variable_genes runs per batch.
    - var['highly_variable'] must still exist and be boolean
    - var['highly_variable_nbatches'] must exist (scanpy adds this per-batch column)
    - metrics and uns must record the batch_key that was used
    - HVG count should still be > 0
    """
    adata = _make_large_adata(n_cells=150, n_genes=400)

    # Add a simple 3-batch label to obs
    import numpy as np
    batch_labels = np.array(["batchA", "batchB", "batchC"])
    adata.obs["batch"] = batch_labels[np.arange(adata.n_obs) % 3]

    n_top = 150
    adata_norm, metrics = normalize(
        adata,
        n_top_genes=n_top,
        hvg_flavor="seurat",        # seurat flavor — stable on small fixtures
        batch_key="batch",
    )

    # HVG column present and boolean
    assert "highly_variable" in adata_norm.var.columns
    assert adata_norm.var["highly_variable"].sum() > 0

    # scanpy adds highly_variable_nbatches when batch_key is used
    assert "highly_variable_nbatches" in adata_norm.var.columns, (
        "var['highly_variable_nbatches'] missing — batch_key may not have been passed"
    )

    # batch_key recorded in metrics and uns
    assert metrics["batch_key"] == "batch"
    assert adata_norm.uns["omicsage_normalization"]["batch_key"] == "batch"
