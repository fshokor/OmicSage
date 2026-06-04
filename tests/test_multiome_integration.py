"""
tests/test_multiome_integration.py
Tests for pipeline/modules/multiome/multiome_integration.py

Coverage (active — MOFA+):
  Happy path
  - Returns (MuData, dict) tuple
  - obsm["X_mofa"] written
  - obsm["X_umap_mofa"] written
  - obsm["X_umap"] NOT left on mdata after MOFA+
  - uns["omicsage_mofa"] written

  Metrics dict — MOFA+
  - Required keys present
  - n_cells correct
  - batch_key recorded
  - n_batches correct
  - mofa_key == "X_mofa"
  - umap_key == "X_umap_mofa"
  - umap_computed True
  - n_factors matches obsm shape

  Provenance — MOFA+
  - module, timestamp, params (batch_key, n_factors, use_layer, random_state)
  - outputs (mofa_key, umap_key, n_factors_actual)
  - n_factors_actual matches obsm["X_mofa"].shape[1]

  inplace — MOFA+
  - inplace=False: original untouched
  - inplace=True: mdata modified in place

  Validation — MOFA+
  - Missing "rna" modality raises KeyError
  - Missing "atac" modality raises KeyError
  - Missing batch_key raises KeyError

Coverage (commented out — MultiVI):
  All MultiVI tests are commented out because scvi-tools training is too slow
  for CI and the MultiVI.setup_mudata API requires a real GPU-ready environment.
  Uncomment and mark @pytest.mark.slow when running integration tests locally.
  Pattern mirrors test_cite_integration.py totalVI test handling.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData

from pipeline.modules.multiome.multiome_integration import run_mofa, run_multivi

try:
    import scvi as _scvi_check
    _SCVI_INSTALLED = True
except Exception:
    _SCVI_INSTALLED = False

skipif_no_scvi = pytest.mark.skipif(
    not _SCVI_INSTALLED,
    reason="scvi-tools not installed or not importable in this environment",
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_CELLS   = 80
N_GENES   = 60
N_PEAKS   = 100
N_PCS     = 15
BATCH_KEY = "batch"


def _make_rna(n_cells: int = N_CELLS, seed: int = 0) -> AnnData:
    """Minimal RNA AnnData resembling the output of normalize.py."""
    rng  = np.random.default_rng(seed)
    raw  = rng.integers(1, 200, size=(n_cells, N_GENES)).astype(np.float32)
    logX = np.log1p(raw)

    adata = ad.AnnData(X=sp.csr_matrix(logX))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(N_GENES)]
    adata.obs[BATCH_KEY] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"] = sp.csr_matrix(raw)
    adata.obsm["X_pca"]   = rng.standard_normal((n_cells, N_PCS)).astype(np.float32)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=5, random_state=0)
    return adata


def _make_atac(n_cells: int = N_CELLS, seed: int = 1) -> AnnData:
    """
    Minimal ATAC AnnData resembling the output of atac_reduce.py.

    Keys present: layers["counts"], layers["tf_idf"], obsm["X_lsi"],
    obsm["X_umap_atac"], obs["atac_leiden"], obs[BATCH_KEY].
    """
    rng    = np.random.default_rng(seed)
    counts = rng.integers(0, 20, size=(n_cells, N_PEAKS)).astype(np.float32)
    tfidf  = rng.random(size=(n_cells, N_PEAKS)).astype(np.float32)

    # Peak var_names in chr:start-end format
    peak_names = [f"chr1:{i*1000}-{i*1000+500}" for i in range(N_PEAKS)]

    adata = ad.AnnData(X=sp.csr_matrix(tfidf))
    adata.obs_names  = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names  = peak_names
    adata.obs[BATCH_KEY] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"] = sp.csr_matrix(counts)
    adata.layers["tf_idf"] = sp.csr_matrix(tfidf)
    adata.obsm["X_lsi"]       = rng.standard_normal((n_cells, N_PCS)).astype(np.float32)
    adata.obsm["X_umap_atac"] = rng.standard_normal((n_cells, 2)).astype(np.float32)
    adata.obs["atac_leiden"]  = [str(i % 3) for i in range(n_cells)]
    return adata


def _make_mdata(
    n_cells: int = N_CELLS,
    include_rna: bool = True,
    include_atac: bool = True,
    seed: int = 0,
) -> MuData:
    mods = {}
    if include_rna:
        mods["rna"] = _make_rna(n_cells=n_cells, seed=seed)
    if include_atac:
        mods["atac"] = _make_atac(n_cells=n_cells, seed=seed + 1)
    return MuData(mods)


# ---------------------------------------------------------------------------
# MOFA+ — happy path
# ---------------------------------------------------------------------------

class TestRunMofaHappyPath:
    def test_returns_tuple(self):
        out = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert isinstance(out, tuple) and len(out) == 2

    def test_returns_mdata_and_dict(self):
        out, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert isinstance(out, MuData)
        assert isinstance(metrics, dict)

    def test_x_mofa_in_obsm(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "X_mofa" in out.obsm

    def test_x_mofa_shape_cells(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.obsm["X_mofa"].shape[0] == N_CELLS

    def test_x_mofa_shape_factors_positive(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.obsm["X_mofa"].shape[1] >= 1

    def test_x_umap_mofa_in_obsm(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "X_umap_mofa" in out.obsm

    def test_x_umap_mofa_shape(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.obsm["X_umap_mofa"].shape == (N_CELLS, 2)

    def test_no_x_umap_left_on_mdata(self):
        """X_umap must be renamed to X_umap_mofa — not left as X_umap."""
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "X_umap" not in out.obsm

    def test_provenance_written(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "omicsage_mofa" in out.uns

    def test_modalities_still_present(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "rna" in out.mod and "atac" in out.mod


# ---------------------------------------------------------------------------
# MOFA+ — metrics dict
# ---------------------------------------------------------------------------

class TestRunMofaMetrics:
    _REQUIRED = {
        "n_cells", "n_factors", "batch_key", "n_batches",
        "batch_values", "mofa_key", "umap_key",
        "umap_computed", "random_state", "method",
    }

    def test_required_keys_present(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert self._REQUIRED.issubset(set(metrics.keys()))

    def test_n_cells(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["n_cells"] == N_CELLS

    def test_batch_key_recorded(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["batch_key"] == BATCH_KEY

    def test_n_batches(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["n_batches"] == 2

    def test_mofa_key(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["mofa_key"] == "X_mofa"

    def test_umap_key(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["umap_key"] == "X_umap_mofa"

    def test_umap_computed_true(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["umap_computed"] is True

    def test_method_is_mofa(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["method"] == "mofa"

    def test_n_factors_matches_obsm(self):
        out, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["n_factors"] == out.obsm["X_mofa"].shape[1]


# ---------------------------------------------------------------------------
# MOFA+ — provenance
# ---------------------------------------------------------------------------

class TestRunMofaProvenance:
    def test_module_name(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["module"] == "multiome_integration.run_mofa"

    def test_timestamp_non_empty(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["timestamp"]

    def test_params_batch_key(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["batch_key"] == BATCH_KEY

    def test_params_n_factors(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["n_factors"] == 3

    def test_params_random_state(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["random_state"] == 0

    def test_params_use_layer_none(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["use_layer"] is None

    def test_outputs_mofa_key(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["outputs"]["mofa_key"] == "X_mofa"

    def test_outputs_umap_key(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["outputs"]["umap_key"] == "X_umap_mofa"

    def test_outputs_n_factors_actual_matches_obsm(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        actual = out.uns["omicsage_mofa"]["outputs"]["n_factors_actual"]
        assert actual == out.obsm["X_mofa"].shape[1]


# ---------------------------------------------------------------------------
# MOFA+ — inplace
# ---------------------------------------------------------------------------

class TestRunMofaInplace:
    def test_inplace_false_original_untouched(self):
        mdata = _make_mdata()
        run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=False)
        assert "X_mofa" not in mdata.obsm

    def test_inplace_true_modifies_mdata(self):
        mdata = _make_mdata()
        run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=True)
        assert "X_mofa" in mdata.obsm

    def test_inplace_false_returns_different_object(self):
        mdata = _make_mdata()
        out, _ = run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=False)
        assert out is not mdata

    def test_inplace_true_returns_same_object(self):
        mdata = _make_mdata()
        out, _ = run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=True)
        assert out is mdata


# ---------------------------------------------------------------------------
# MOFA+ — validation
# ---------------------------------------------------------------------------

class TestRunMofaValidation:
    def test_missing_rna_raises_keyerror(self):
        with pytest.raises(KeyError, match="rna"):
            run_mofa(_make_mdata(include_rna=False), batch_key=BATCH_KEY)

    def test_missing_atac_raises_keyerror(self):
        with pytest.raises(KeyError, match="atac"):
            run_mofa(_make_mdata(include_atac=False), batch_key=BATCH_KEY)

    def test_missing_batch_key_raises_keyerror(self):
        with pytest.raises(KeyError, match="bad_key"):
            run_mofa(_make_mdata(), batch_key="bad_key")


# ---------------------------------------------------------------------------
# MultiVI — validation only (training tests commented out)
# ---------------------------------------------------------------------------

class TestRunMultiVIValidation:
    """
    Validation tests run without scvi training — they only exercise the
    _validate_multivi_inputs() guard before the scvi import is attempted.
    """

    def test_missing_rna_raises_keyerror(self):
        with pytest.raises(KeyError, match="rna"):
            run_multivi(_make_mdata(include_rna=False), batch_key=BATCH_KEY)

    def test_missing_atac_raises_keyerror(self):
        with pytest.raises(KeyError, match="atac"):
            run_multivi(_make_mdata(include_atac=False), batch_key=BATCH_KEY)

    def test_missing_rna_counts_layer_raises_keyerror(self):
        mdata = _make_mdata()
        del mdata["rna"].layers["counts"]
        with pytest.raises(KeyError, match="counts"):
            run_multivi(mdata, batch_key=BATCH_KEY)

    def test_missing_atac_counts_layer_raises_keyerror(self):
        mdata = _make_mdata()
        del mdata["atac"].layers["counts"]
        with pytest.raises(KeyError, match="counts"):
            run_multivi(mdata, batch_key=BATCH_KEY)

    def test_missing_batch_key_raises_keyerror(self):
        with pytest.raises(KeyError, match="bad_key"):
            run_multivi(_make_mdata(), batch_key="bad_key")


# ---------------------------------------------------------------------------
# MultiVI — happy path (commented out: requires scvi-tools + GPU environment)
# ---------------------------------------------------------------------------

# @skipif_no_scvi
# class TestRunMultiVIHappyPath:
#     def test_returns_tuple(self):
#         out = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert isinstance(out, tuple) and len(out) == 2
#
#     def test_returns_mdata_and_dict(self):
#         out, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert isinstance(out, MuData)
#         assert isinstance(metrics, dict)
#
#     def test_x_multivi_in_obsm(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_multivi" in out.obsm
#
#     def test_x_multivi_shape_cells(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_multivi"].shape[0] == N_CELLS
#
#     def test_x_multivi_shape_latent_positive(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_multivi"].shape[1] > 0
#
#     def test_x_umap_multivi_in_obsm(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_umap_multivi" in out.obsm
#
#     def test_x_umap_multivi_shape(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_umap_multivi"].shape == (N_CELLS, 2)
#
#     def test_no_x_umap_left_on_mdata(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_umap" not in out.obsm
#
#     def test_provenance_written(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "omicsage_multivi" in out.uns
#
#     def test_provenance_params_batch_key(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.uns["omicsage_multivi"]["params"]["batch_key"] == BATCH_KEY
#
#     def test_provenance_params_max_epochs(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.uns["omicsage_multivi"]["params"]["max_epochs"] == 2
#
#     def test_provenance_n_latent_matches_obsm(self):
#         out, _ = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         n_latent = out.uns["omicsage_multivi"]["outputs"]["n_latent_actual"]
#         assert n_latent == out.obsm["X_multivi"].shape[1]
#
#
# @skipif_no_scvi
# class TestRunMultiVIMetrics:
#     def test_required_keys_present(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         required = {"n_cells", "n_latent", "batch_key", "n_batches",
#                     "batch_values", "max_epochs", "multivi_key",
#                     "umap_key", "umap_computed", "random_state", "method"}
#         assert required.issubset(set(metrics.keys()))
#
#     def test_n_cells(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["n_cells"] == N_CELLS
#
#     def test_max_epochs(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["max_epochs"] == 2
#
#     def test_n_batches(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["n_batches"] == 2
#
#     def test_umap_computed_true(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["umap_computed"] is True
#
#     def test_multivi_key(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["multivi_key"] == "X_multivi"
#
#     def test_umap_key(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["umap_key"] == "X_umap_multivi"
#
#     def test_method_is_multivi(self):
#         _, metrics = run_multivi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["method"] == "multivi"
#
#
# @skipif_no_scvi
# class TestRunMultiVIInplace:
#     def test_inplace_false_original_untouched(self):
#         mdata = _make_mdata()
#         run_multivi(mdata, batch_key=BATCH_KEY, max_epochs=2, inplace=False)
#         assert "X_multivi" not in mdata.obsm
#
#     def test_inplace_true_modifies_mdata(self):
#         mdata = _make_mdata()
#         run_multivi(mdata, batch_key=BATCH_KEY, max_epochs=2, inplace=True)
#         assert "X_multivi" in mdata.obsm
