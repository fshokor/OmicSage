from __future__ import annotations

import numpy as np
import pytest
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData

from pipeline.modules.cite.cite_integration import run_mofa, run_totalvi

# Skip all totalVI tests until scvi import issue is resolved
pytest.importorskip("scvi")

# Mark all totalVI tests to skip if scvi-tools is not importable
try:
    import scvi as _scvi_check
    _SCVI_INSTALLED = True
except Exception:
    _SCVI_INSTALLED = False

skipif_no_scvi = pytest.mark.skipif(
    not _SCVI_INSTALLED,
    reason="scvi-tools not installed or not importable in this environment"
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_CELLS   = 80
N_GENES   = 60
N_ADTS    = 20
N_PCS     = 15
BATCH_KEY = "donor"


def _make_rna(n_cells: int = N_CELLS, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    raw  = rng.integers(1, 200, size=(n_cells, N_GENES)).astype(np.float32)
    logX = np.log1p(raw)
    adata = ad.AnnData(X=sp.csr_matrix(logX))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(N_GENES)]
    adata.obs[BATCH_KEY] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"] = sp.csr_matrix(raw)
    adata.obsm["X_pca"] = rng.standard_normal((n_cells, N_PCS)).astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=5, random_state=0)
    return adata


def _make_adt(n_cells: int = N_CELLS, seed: int = 1) -> AnnData:
    rng  = np.random.default_rng(seed)
    raw  = rng.integers(1, 100, size=(n_cells, N_ADTS)).astype(np.float32)
    clr  = rng.normal(0, 2, size=(n_cells, N_ADTS)).astype(np.float32)
    adata = ad.AnnData(X=clr.copy())
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"ADT_{i}" for i in range(N_ADTS)]
    adata.obs[BATCH_KEY] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"]  = sp.csr_matrix(raw)
    adata.layers["adt_clr"] = clr.copy()
    n_pcs_adt = min(N_PCS, N_ADTS - 1)
    adata.obsm["X_pca_adt"]         = rng.standard_normal((n_cells, n_pcs_adt)).astype(np.float32)
    adata.obsm["X_pca_harmony_adt"] = rng.standard_normal((n_cells, n_pcs_adt)).astype(np.float32)
    adata.obsm["X_umap_adt"]        = rng.standard_normal((n_cells, 2)).astype(np.float32)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony_adt", n_neighbors=5, random_state=0)
    return adata


def _make_mdata(n_cells: int = N_CELLS,
                include_rna: bool = True,
                include_adt: bool = True,
                seed: int = 0) -> MuData:
    mods = {}
    if include_rna:
        mods["rna"] = _make_rna(n_cells=n_cells, seed=seed)
    if include_adt:
        mods["adt"] = _make_adt(n_cells=n_cells, seed=seed + 1)
    return MuData(mods)


# ---------------------------------------------------------------------------
# MOFA+ tests
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
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "X_umap" not in out.obsm

    def test_provenance_written(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert "omicsage_mofa" in out.uns

    def test_provenance_params_batch_key(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["batch_key"] == BATCH_KEY

    def test_provenance_params_n_factors(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert out.uns["omicsage_mofa"]["params"]["n_factors"] == 3

    def test_provenance_n_factors_actual_matches_obsm(self):
        out, _ = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        actual = out.uns["omicsage_mofa"]["outputs"]["n_factors_actual"]
        assert actual == out.obsm["X_mofa"].shape[1]


class TestRunMofaMetrics:
    def test_metrics_required_keys(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        required = {"n_cells", "n_factors", "batch_key", "n_batches",
                    "batch_values", "mofa_key", "umap_key",
                    "umap_computed", "random_state"}
        assert required.issubset(set(metrics.keys()))

    def test_metrics_n_cells(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["n_cells"] == N_CELLS

    def test_metrics_batch_key(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["batch_key"] == BATCH_KEY

    def test_metrics_n_batches(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["n_batches"] == 2

    def test_metrics_umap_computed_true(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["umap_computed"] is True

    def test_metrics_mofa_key(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["mofa_key"] == "X_mofa"

    def test_metrics_umap_key(self):
        _, metrics = run_mofa(_make_mdata(), batch_key=BATCH_KEY, n_factors=3)
        assert metrics["umap_key"] == "X_umap_mofa"


class TestRunMofaInplace:
    def test_inplace_false_original_untouched(self):
        mdata = _make_mdata()
        run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=False)
        assert "X_mofa" not in mdata.obsm

    def test_inplace_true_modifies_mdata(self):
        mdata = _make_mdata()
        run_mofa(mdata, batch_key=BATCH_KEY, n_factors=3, inplace=True)
        assert "X_mofa" in mdata.obsm


class TestRunMofaValidation:
    def test_missing_rna_raises_keyerror(self):
        with pytest.raises(KeyError, match="rna"):
            run_mofa(_make_mdata(include_rna=False), batch_key=BATCH_KEY)

    def test_missing_adt_raises_keyerror(self):
        with pytest.raises(KeyError, match="adt"):
            run_mofa(_make_mdata(include_adt=False), batch_key=BATCH_KEY)

    def test_missing_batch_key_raises_keyerror(self):
        with pytest.raises(KeyError, match="bad_key"):
            run_mofa(_make_mdata(), batch_key="bad_key")


# ---------------------------------------------------------------------------
# totalVI tests
# ---------------------------------------------------------------------------

# @skipif_no_scvi
# class TestRunTotalVIHappyPath:
#     def test_returns_tuple(self):
#         out = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert isinstance(out, tuple) and len(out) == 2

#     def test_returns_mdata_and_dict(self):
#         out, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert isinstance(out, MuData)
#         assert isinstance(metrics, dict)

#     def test_x_totalvi_in_obsm(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_totalVI" in out.obsm

#     def test_x_totalvi_shape_cells(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_totalVI"].shape[0] == N_CELLS

#     def test_x_totalvi_shape_latent_positive(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_totalVI"].shape[1] > 0

#     def test_x_umap_totalvi_in_obsm(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_umap_totalVI" in out.obsm

#     def test_x_umap_totalvi_shape(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.obsm["X_umap_totalVI"].shape == (N_CELLS, 2)

#     def test_no_x_umap_left_on_mdata(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "X_umap" not in out.obsm

#     def test_provenance_written(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert "omicsage_totalVI" in out.uns

#     def test_provenance_params_batch_key(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.uns["omicsage_totalVI"]["params"]["batch_key"] == BATCH_KEY

#     def test_provenance_params_max_epochs(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert out.uns["omicsage_totalVI"]["params"]["max_epochs"] == 2

#     def test_provenance_n_latent_matches_obsm(self):
#         out, _ = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         n_latent = out.uns["omicsage_totalVI"]["outputs"]["n_latent"]
#         assert n_latent == out.obsm["X_totalVI"].shape[1]


# @skipif_no_scvi
# class TestRunTotalVIMetrics:
#     def test_metrics_required_keys(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         required = {"n_cells", "n_latent", "batch_key", "n_batches",
#                     "batch_values", "max_epochs", "totalvi_key",
#                     "umap_key", "umap_computed", "random_state"}
#         assert required.issubset(set(metrics.keys()))

#     def test_metrics_n_cells(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["n_cells"] == N_CELLS

#     def test_metrics_max_epochs(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["max_epochs"] == 2

#     def test_metrics_n_batches(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["n_batches"] == 2

#     def test_metrics_umap_computed_true(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["umap_computed"] is True

#     def test_metrics_totalvi_key(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["totalvi_key"] == "X_totalVI"

#     def test_metrics_umap_key(self):
#         _, metrics = run_totalvi(_make_mdata(), batch_key=BATCH_KEY, max_epochs=2)
#         assert metrics["umap_key"] == "X_umap_totalVI"


# @skipif_no_scvi
# class TestRunTotalVIInplace:
#     def test_inplace_false_original_untouched(self):
#         mdata = _make_mdata()
#         run_totalvi(mdata, batch_key=BATCH_KEY, max_epochs=2, inplace=False)
#         assert "X_totalVI" not in mdata.obsm

#     def test_inplace_true_modifies_mdata(self):
#         mdata = _make_mdata()
#         run_totalvi(mdata, batch_key=BATCH_KEY, max_epochs=2, inplace=True)
#         assert "X_totalVI" in mdata.obsm


# @skipif_no_scvi
# class TestRunTotalVIValidation:
#     def test_missing_rna_raises_keyerror(self):
#         with pytest.raises(KeyError, match="rna"):
#             run_totalvi(_make_mdata(include_rna=False), batch_key=BATCH_KEY)

#     def test_missing_adt_raises_keyerror(self):
#         with pytest.raises(KeyError, match="adt"):
#             run_totalvi(_make_mdata(include_adt=False), batch_key=BATCH_KEY)

#     def test_missing_rna_counts_layer_raises_keyerror(self):
#         mdata = _make_mdata()
#         del mdata["rna"].layers["counts"]
#         with pytest.raises(KeyError, match="counts"):
#             run_totalvi(mdata, batch_key=BATCH_KEY)

#     def test_missing_adt_counts_layer_raises_keyerror(self):
#         mdata = _make_mdata()
#         del mdata["adt"].layers["counts"]
#         with pytest.raises(KeyError, match="counts"):
#             run_totalvi(mdata, batch_key=BATCH_KEY)

#     def test_missing_batch_key_raises_keyerror(self):
#         with pytest.raises(KeyError, match="bad_key"):
#             run_totalvi(_make_mdata(), batch_key="bad_key")
