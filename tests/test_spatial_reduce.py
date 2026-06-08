"""
test_spatial_reduce.py — OmicSage Phase 7, Session 2
Tests for spatial_reduce() and generate_spatial_reduce_report().

Conventions:
  - pytest.importorskip() for optional dependencies
  - Fixtures build minimal valid AnnData objects
  - Tests are independent (no shared mutable state)
  - inplace=False default tested first; inplace=True tested separately
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")
scanpy = pytest.importorskip("scanpy")
squidpy = pytest.importorskip("squidpy")

import anndata as ad

from pipeline.modules.spatial.spatial_reduce import spatial_reduce
from reports.templates.spatial.spatial_reduce_report import (
    generate_spatial_reduce_report,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_SPOTS = 80
N_GENES = 200
RNG = np.random.default_rng(42)


def _make_raw_adata(n_obs: int = N_SPOTS, n_vars: int = N_GENES) -> ad.AnnData:
    """Minimal AnnData mimicking raw Visium output from spatial_ingest."""
    import pandas as pd

    # Raw counts (negative binomial-like)
    X = sp.csr_matrix(
        RNG.negative_binomial(5, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    )
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs_names = [f"spot{i}" for i in range(n_obs)]

    adata = ad.AnnData(
        X=X,
        obs={"total_counts": np.array(X.sum(axis=1)).ravel()},
        var={"gene_ids": var_names},
    )
    adata.obs_names = obs_names
    adata.var_names = var_names

    # Spatial coordinates (grid pattern)
    coords = np.column_stack([
        np.tile(np.arange(10), n_obs // 10 + 1)[:n_obs],
        np.repeat(np.arange(n_obs // 10 + 1), 10)[:n_obs],
    ]).astype(float)
    adata.obsm["spatial"] = coords

    # uns["spatial"] required by squidpy
    adata.uns["spatial"] = {"visium": {"images": {}, "scalefactors": {}}}

    # Store raw in layers["counts"] (as spatial_ingest does for visium)
    adata.layers["counts"] = X.copy()

    # Ingest provenance
    adata.uns["omicsage_spatial_ingest"] = {
        "spatial_type": "visium",
        "n_obs": n_obs,
        "n_vars": n_vars,
        "timestamp": "2026-01-01T00:00:00",
    }

    return adata


def _make_benchmark_adata(n_obs: int = N_SPOTS, n_vars: int = N_GENES) -> ad.AnnData:
    """AnnData simulating pre-processed benchmark (already normalized)."""
    import pandas as pd

    X = sp.csr_matrix(
        RNG.random(size=(n_obs, n_vars)).astype(np.float32)
    )
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs_names = [f"spot{i}" for i in range(n_obs)]

    adata = ad.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obsm["spatial"] = np.column_stack([
        np.tile(np.arange(10), n_obs // 10 + 1)[:n_obs],
        np.repeat(np.arange(n_obs // 10 + 1), 10)[:n_obs],
    ]).astype(float)
    adata.uns["spatial"] = {"benchmark": {"images": {}, "scalefactors": {}}}
    # NO layers["counts"] — benchmark path
    adata.uns["omicsage_spatial_ingest"] = {
        "spatial_type": "benchmark",
        "n_obs": n_obs,
        "n_vars": n_vars,
        "timestamp": "2026-01-01T00:00:00",
    }
    return adata


# ---------------------------------------------------------------------------
# Contract tests — return type and provenance
# ---------------------------------------------------------------------------

class TestSpatialReduceContract:

    def test_returns_tuple(self):
        adata = _make_raw_adata()
        result = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_raw_adata()
        out, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert isinstance(out, ad.AnnData)
        assert isinstance(params, dict)

    def test_provenance_written_to_uns(self):
        adata = _make_raw_adata()
        out, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert "omicsage_spatial_reduce" in out.uns
        assert out.uns["omicsage_spatial_reduce"] is params

    def test_provenance_has_required_keys(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        for key in ("module", "timestamp", "params", "outputs"):
            assert key in params, f"Missing provenance key: {key}"

    def test_provenance_params_keys(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        expected = {
            "n_top_genes", "n_comps", "n_neighbors", "coord_type",
            "normalize_total", "target_sum", "log1p", "flavor",
            "skipped_normalization",
        }
        assert expected.issubset(set(params["params"].keys()))

    def test_provenance_outputs_keys(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        expected = {
            "n_hvgs", "n_comps_computed", "pca_variance_ratio_top10",
            "pca_cumulative_variance_top10", "spatial_graph_n_edges",
            "spatial_graph_mean_neighbors",
        }
        assert expected.issubset(set(params["outputs"].keys()))


# ---------------------------------------------------------------------------
# Inplace behaviour
# ---------------------------------------------------------------------------

class TestSpatialReduceInplace:

    def test_inplace_false_returns_new_object(self):
        adata = _make_raw_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10, inplace=False)
        assert out is not adata

    def test_inplace_false_does_not_modify_input(self):
        adata = _make_raw_adata()
        original_uns_keys = set(adata.uns.keys())
        spatial_reduce(adata, n_top_genes=50, n_comps=10, inplace=False)
        assert "omicsage_spatial_reduce" not in adata.uns

    def test_inplace_true_modifies_input(self):
        adata = _make_raw_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10, inplace=True)
        assert out is adata
        assert "omicsage_spatial_reduce" in adata.uns


# ---------------------------------------------------------------------------
# AnnData state after spatial_reduce
# ---------------------------------------------------------------------------

class TestSpatialReduceOutputState:

    @pytest.fixture
    def reduced_adata(self):
        adata = _make_raw_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10, inplace=False)
        return out

    def test_layers_counts_preserved(self, reduced_adata):
        assert "counts" in reduced_adata.layers

    def test_hvg_column_in_var(self, reduced_adata):
        assert "highly_variable" in reduced_adata.var.columns
        assert reduced_adata.var["highly_variable"].dtype == bool

    def test_hvg_count_reasonable(self, reduced_adata):
        n_hvg = int(reduced_adata.var["highly_variable"].sum())
        assert n_hvg > 0
        assert n_hvg <= N_GENES

    def test_pca_in_obsm(self, reduced_adata):
        assert "X_pca" in reduced_adata.obsm

    def test_pca_shape(self, reduced_adata):
        pca = reduced_adata.obsm["X_pca"]
        assert pca.shape[0] == N_SPOTS
        assert pca.shape[1] <= 10  # requested n_comps

    def test_pca_variance_ratio_in_uns(self, reduced_adata):
        assert "pca" in reduced_adata.uns
        assert "variance_ratio" in reduced_adata.uns["pca"]

    def test_spatial_connectivities_in_obsp(self, reduced_adata):
        assert "spatial_connectivities" in reduced_adata.obsp

    def test_spatial_distances_in_obsp(self, reduced_adata):
        assert "spatial_distances" in reduced_adata.obsp

    def test_spatial_neighbors_in_uns(self, reduced_adata):
        assert "spatial_neighbors" in reduced_adata.uns

    def test_spatial_connectivities_shape(self, reduced_adata):
        conn = reduced_adata.obsp["spatial_connectivities"]
        assert conn.shape == (N_SPOTS, N_SPOTS)

    def test_spatial_connectivities_symmetric(self, reduced_adata):
        conn = reduced_adata.obsp["spatial_connectivities"]
        diff = (conn - conn.T).data
        assert len(diff) == 0 or np.allclose(diff, 0, atol=1e-6)

    def test_means_dispersions_in_var(self, reduced_adata):
        assert "means" in reduced_adata.var.columns
        assert "dispersions_norm" in reduced_adata.var.columns

    def test_obs_shape_preserved(self, reduced_adata):
        assert reduced_adata.n_obs == N_SPOTS


# ---------------------------------------------------------------------------
# Normalization behaviour
# ---------------------------------------------------------------------------

class TestNormalizationBehaviour:

    def test_raw_data_normalized_inplace(self):
        adata = _make_raw_adata()
        out, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert not bool(params["params"]["skipped_normalization"])

    def test_benchmark_data_skips_normalization(self):
        adata = _make_benchmark_adata()
        out, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert bool(params["params"]["skipped_normalization"])

    def test_normalize_false_skips_normalization(self):
        adata = _make_raw_adata()
        out, params = spatial_reduce(
            adata, n_top_genes=50, n_comps=10, normalize_total=False, log1p=False
        )
        assert not bool(params["params"]["normalize_total"])

    def test_layers_counts_present_after_reduce_raw(self):
        adata = _make_raw_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert "counts" in out.layers

    def test_layers_counts_present_after_reduce_benchmark(self):
        adata = _make_benchmark_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert "counts" in out.layers


# ---------------------------------------------------------------------------
# Parameter passthrough
# ---------------------------------------------------------------------------

class TestParameterPassthrough:

    def test_n_top_genes_respected(self):
        adata = _make_raw_adata(n_vars=200)
        out, params = spatial_reduce(adata, n_top_genes=30, n_comps=5)
        assert params["outputs"]["n_hvgs"] <= 30

    def test_n_comps_respected(self):
        adata = _make_raw_adata()
        out, params = spatial_reduce(adata, n_top_genes=50, n_comps=8)
        assert out.obsm["X_pca"].shape[1] <= 8

    def test_n_comps_capped_at_min_obs_vars(self):
        # Very small dataset — n_comps cannot exceed min(n_obs, n_vars) - 1
        adata = _make_raw_adata(n_obs=20, n_vars=30)
        out, params = spatial_reduce(adata, n_top_genes=20, n_comps=50)
        # Should not raise — n_comps is capped internally
        assert out.obsm["X_pca"].shape[1] < 50

    def test_coord_type_none_default(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert params["params"]["coord_type"] is None

    def test_flavor_stored_in_provenance(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(
            adata, n_top_genes=50, n_comps=10, flavor="cell_ranger"
        )
        assert params["params"]["flavor"] == "cell_ranger"


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestSpatialReduceValidation:

    def test_rejects_non_anndata(self):
        with pytest.raises(TypeError, match="AnnData"):
            spatial_reduce({"not": "anndata"})

    def test_rejects_missing_spatial_obsm(self):
        adata = _make_raw_adata()
        del adata.obsm["spatial"]
        with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
            spatial_reduce(adata)

    def test_rejects_empty_obs(self):
        adata = _make_raw_adata()
        adata = adata[:0].copy()
        with pytest.raises(ValueError, match="0 observations"):
            spatial_reduce(adata)

    def test_rejects_empty_vars(self):
        adata = _make_raw_adata()
        adata = adata[:, :0].copy()
        with pytest.raises(ValueError, match="0 variables"):
            spatial_reduce(adata)


# ---------------------------------------------------------------------------
# Provenance content correctness
# ---------------------------------------------------------------------------

class TestSpatialReduceProvenanceValues:

    def test_n_hvgs_positive(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert params["outputs"]["n_hvgs"] > 0

    def test_spatial_graph_n_edges_positive(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert params["outputs"]["spatial_graph_n_edges"] > 0

    def test_mean_neighbors_reasonable(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        mean_nb = params["outputs"]["spatial_graph_mean_neighbors"]
        assert 1.0 <= mean_nb <= 6.0

    def test_pca_variance_ratio_sums_to_leq_1(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        total = sum(params["outputs"]["pca_variance_ratio_top10"])
        assert total <= 1.0 + 1e-6

    def test_timestamp_present(self):
        adata = _make_raw_adata()
        _, params = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        assert params["timestamp"]
        assert "T" in params["timestamp"]  # ISO format


# ---------------------------------------------------------------------------
# Report tests
# ---------------------------------------------------------------------------

class TestSpatialReduceReport:

    @pytest.fixture
    def reduced_adata(self):
        adata = _make_raw_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        return out

    def test_report_generates_file(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_reduce_report(
            reduced_adata, out_path, dataset_id="test"
        )
        assert result.endswith(".html")
        import os
        assert os.path.isfile(result)

    def test_report_not_empty(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_reduce_report(reduced_adata, out_path, dataset_id="test")
        import os
        assert os.path.getsize(out_path) > 1000

    def test_report_is_valid_html(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_reduce_report(reduced_adata, out_path, dataset_id="test")
        with open(out_path) as fh:
            content = fh.read()
        assert "<!DOCTYPE html>" in content
        assert "</html>" in content

    def test_report_contains_dataset_id(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_reduce_report(
            reduced_adata, out_path, dataset_id="my_dataset"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "my_dataset" in content

    def test_report_contains_hvg_count(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_reduce_report(reduced_adata, out_path, dataset_id="test")
        with open(out_path) as fh:
            content = fh.read()
        n_hvgs = str(reduced_adata.uns["omicsage_spatial_reduce"]["outputs"]["n_hvgs"])
        assert n_hvgs in content

    def test_report_raises_without_provenance(self, tmp_path):
        adata = _make_raw_adata()
        with pytest.raises(ValueError, match="omicsage_spatial_reduce"):
            generate_spatial_reduce_report(
                adata, str(tmp_path / "x.html"), dataset_id="test"
            )

    def test_report_returns_absolute_path(self, reduced_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_reduce_report(
            reduced_adata, out_path, dataset_id="test"
        )
        import os
        assert os.path.isabs(result)

    def test_report_benchmark_note(self, tmp_path):
        adata = _make_benchmark_adata()
        out, _ = spatial_reduce(adata, n_top_genes=50, n_comps=10)
        out_path = str(tmp_path / "bench_report.html")
        generate_spatial_reduce_report(out, out_path, dataset_id="bench")
        with open(out_path) as fh:
            content = fh.read()
        assert "Normalization skipped" in content
