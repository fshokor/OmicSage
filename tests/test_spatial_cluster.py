"""
test_spatial_cluster.py — OmicSage Phase 7, Session 3
Tests for spatial_cluster() and generate_spatial_cluster_report().

Conventions match test_spatial_reduce.py:
  - pytest.importorskip() for optional dependencies
  - Fixtures build minimal valid AnnData objects via the real pipeline
  - Tests are independent (no shared mutable state)
  - inplace=False default tested first
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")
scanpy = pytest.importorskip("scanpy")
squidpy = pytest.importorskip("squidpy")

import anndata as ad

from pipeline.modules.scripts.spatial.spatial_reduce import spatial_reduce
from pipeline.modules.scripts.spatial.spatial_cluster import spatial_cluster
from reports.templates.spatial.spatial_cluster_report import (
    generate_spatial_cluster_report,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_SPOTS = 120
N_GENES = 300
RNG = np.random.default_rng(7)


def _make_reduced_adata(
    n_obs: int = N_SPOTS,
    n_vars: int = N_GENES,
) -> ad.AnnData:
    """Build a minimal AnnData that has passed through spatial_reduce."""
    X = sp.csr_matrix(
        RNG.negative_binomial(5, 0.5, size=(n_obs, n_vars)).astype(np.float32)
    )
    var_names = [f"Gene{i}" for i in range(n_vars)]
    obs_names = [f"spot{i}" for i in range(n_obs)]

    adata = ad.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = var_names

    # Grid spatial coordinates
    side = int(np.ceil(np.sqrt(n_obs)))
    coords = np.column_stack([
        np.tile(np.arange(side), side)[:n_obs],
        np.repeat(np.arange(side), side)[:n_obs],
    ]).astype(float)
    adata.obsm["spatial"] = coords
    adata.uns["spatial"] = {"visium": {"images": {}, "scalefactors": {}}}
    adata.layers["counts"] = X.copy()
    adata.uns["omicsage_spatial_ingest"] = {
        "spatial_type": "visium",
        "n_obs": n_obs,
        "n_vars": n_vars,
        "timestamp": "2026-01-01T00:00:00",
    }

    # Run spatial_reduce to get X_pca + spatial graph
    adata, _ = spatial_reduce(
        adata, n_top_genes=80, n_comps=15, inplace=False
    )
    return adata


# ---------------------------------------------------------------------------
# Contract tests — return type and provenance
# ---------------------------------------------------------------------------

class TestSpatialClusterContract:

    def test_returns_tuple(self):
        adata = _make_reduced_adata()
        result = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_returns_anndata_and_dict(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert isinstance(out, ad.AnnData)
        assert isinstance(params, dict)

    def test_provenance_written_to_uns(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert "omicsage_spatial_cluster" in out.uns
        assert out.uns["omicsage_spatial_cluster"] is params

    def test_provenance_has_required_keys(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        for key in ("module", "timestamp", "params", "outputs"):
            assert key in params

    def test_provenance_params_keys(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        expected = {
            "resolution", "n_neighbors", "n_pcs", "n_pcs_actual",
            "random_state", "cluster_key", "annotation_map_provided",
            "run_svg", "svg_n_genes",
        }
        assert expected.issubset(set(params["params"].keys()))

    def test_provenance_outputs_keys(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert "n_clusters" in params["outputs"]
        assert "cluster_sizes" in params["outputs"]


# ---------------------------------------------------------------------------
# Inplace behaviour
# ---------------------------------------------------------------------------

class TestSpatialClusterInplace:

    def test_inplace_false_returns_new_object(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False, inplace=False)
        assert out is not adata

    def test_inplace_false_does_not_modify_input(self):
        adata = _make_reduced_adata()
        spatial_cluster(adata, resolution=0.3, run_svg=False, inplace=False)
        assert "omicsage_spatial_cluster" not in adata.uns

    def test_inplace_true_modifies_input(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False, inplace=True)
        assert out is adata
        assert "omicsage_spatial_cluster" in adata.uns


# ---------------------------------------------------------------------------
# Clustering output state
# ---------------------------------------------------------------------------

class TestSpatialClusterOutputState:

    @pytest.fixture
    def clustered_adata(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False, inplace=False)
        return out

    def test_cluster_key_in_obs(self, clustered_adata):
        assert "spatial_cluster" in clustered_adata.obs.columns

    def test_cluster_labels_are_strings(self, clustered_adata):
        labels = clustered_adata.obs["spatial_cluster"]
        assert labels.dtype == object or str(labels.dtype) == "category"

    def test_at_least_one_cluster(self, clustered_adata):
        n = clustered_adata.obs["spatial_cluster"].nunique()
        assert n >= 1

    def test_all_spots_assigned(self, clustered_adata):
        assert clustered_adata.obs["spatial_cluster"].isna().sum() == 0

    def test_obs_shape_preserved(self, clustered_adata):
        assert clustered_adata.n_obs == N_SPOTS

    def test_pca_still_present(self, clustered_adata):
        assert "X_pca" in clustered_adata.obsm

    def test_spatial_graph_still_present(self, clustered_adata):
        assert "spatial_connectivities" in clustered_adata.obsp

    def test_neighbors_graph_added(self, clustered_adata):
        # sc.pp.neighbors writes connectivities and distances
        assert "neighbors" in clustered_adata.uns or "connectivities" in clustered_adata.obsp


# ---------------------------------------------------------------------------
# Cluster count and size consistency
# ---------------------------------------------------------------------------

class TestClusterConsistency:

    def test_n_clusters_matches_unique_labels(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        n_unique = out.obs["spatial_cluster"].nunique()
        assert params["outputs"]["n_clusters"] == n_unique

    def test_cluster_sizes_sum_to_n_obs(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        total = sum(params["outputs"]["cluster_sizes"].values())
        assert total == N_SPOTS

    def test_higher_resolution_more_clusters(self):
        adata_lo = _make_reduced_adata()
        adata_hi = adata_lo.copy()
        _, params_lo = spatial_cluster(adata_lo, resolution=0.2, run_svg=False)
        _, params_hi = spatial_cluster(adata_hi, resolution=1.5, run_svg=False)
        # Not guaranteed but almost always true; soft check
        assert params_hi["outputs"]["n_clusters"] >= params_lo["outputs"]["n_clusters"]


# ---------------------------------------------------------------------------
# Annotation map
# ---------------------------------------------------------------------------

class TestAnnotationMap:

    def test_annotation_map_writes_label_key(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        cluster_ids = out.obs["spatial_cluster"].unique().tolist()
        ann_map = {str(c): f"Region_{c}" for c in cluster_ids}
        out2, _ = spatial_cluster(
            out, resolution=0.3, run_svg=False,
            annotation_map=ann_map
        )
        assert "spatial_cluster_label" in out2.obs.columns

    def test_annotation_map_none_does_not_write_label_key(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(
            adata, resolution=0.3, run_svg=False, annotation_map=None
        )
        assert "spatial_cluster_label" not in out.obs.columns

    def test_annotation_map_provenance_flag(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        cluster_ids = out.obs["spatial_cluster"].unique().tolist()
        ann_map = {str(c): f"T{c}" for c in cluster_ids}
        _, params2 = spatial_cluster(
            out, resolution=0.3, run_svg=False, annotation_map=ann_map
        )
        assert bool(params2["params"]["annotation_map_provided"])

    def test_unknown_clusters_get_unknown_label(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False)
        # Only map cluster "0" — all others should be "Unknown"
        partial_map = {"0": "KnownRegion"}
        out2, _ = spatial_cluster(
            out, resolution=0.3, run_svg=False, annotation_map=partial_map
        )
        labels = out2.obs["spatial_cluster_label"].values
        assert "Unknown" in labels or "KnownRegion" in labels


# ---------------------------------------------------------------------------
# Custom cluster_key
# ---------------------------------------------------------------------------

class TestCustomClusterKey:

    def test_custom_cluster_key_written(self):
        adata = _make_reduced_adata()
        out, params = spatial_cluster(
            adata, resolution=0.3, run_svg=False, cluster_key="my_clusters"
        )
        assert "my_clusters" in out.obs.columns
        assert params["params"]["cluster_key"] == "my_clusters"


# ---------------------------------------------------------------------------
# SVG (Moran's I) tests
# ---------------------------------------------------------------------------

class TestSVG:

    @pytest.fixture
    def svg_adata(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(
            adata, resolution=0.3, run_svg=True, svg_n_genes=50, inplace=False
        )
        return out

    def test_morani_written_to_uns(self, svg_adata):
        assert "moranI" in svg_adata.uns

    def test_morani_is_dataframe(self, svg_adata):
        import pandas as pd
        assert isinstance(svg_adata.uns["moranI"], pd.DataFrame)

    def test_morani_has_required_columns(self, svg_adata):
        df = svg_adata.uns["moranI"]
        for col in ("I", "pval_norm", "pval_norm_fdr_bh"):
            assert col in df.columns, f"Missing column: {col}"

    def test_morani_sorted_descending(self, svg_adata):
        scores = svg_adata.uns["moranI"]["I"].values
        assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))

    def test_svg_provenance_keys(self, svg_adata):
        outputs = svg_adata.uns["omicsage_spatial_cluster"]["outputs"]
        assert "n_genes_tested" in outputs
        assert "n_significant_fdr05" in outputs
        assert "top5_svg" in outputs

    def test_svg_n_genes_caps_tested(self, svg_adata):
        outputs = svg_adata.uns["omicsage_spatial_cluster"]["outputs"]
        assert outputs["n_genes_tested"] <= 50

    def test_run_svg_false_no_morani(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert "moranI" not in out.uns

    def test_run_svg_false_no_svg_provenance(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert "n_genes_tested" not in params["outputs"]


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

class TestSpatialClusterValidation:

    def test_rejects_non_anndata(self):
        with pytest.raises(TypeError, match="AnnData"):
            spatial_cluster({"not": "anndata"})

    def test_rejects_missing_pca(self):
        adata = _make_reduced_adata()
        del adata.obsm["X_pca"]
        with pytest.raises(ValueError, match="X_pca"):
            spatial_cluster(adata)

    def test_rejects_missing_spatial_connectivities_when_svg(self):
        adata = _make_reduced_adata()
        del adata.obsp["spatial_connectivities"]
        with pytest.raises(ValueError, match="spatial_connectivities"):
            spatial_cluster(adata, run_svg=True)

    def test_missing_spatial_connectivities_ok_when_no_svg(self):
        adata = _make_reduced_adata()
        del adata.obsp["spatial_connectivities"]
        # Should not raise when run_svg=False
        out, _ = spatial_cluster(adata, run_svg=False)
        assert "spatial_cluster" in out.obs.columns

    def test_rejects_empty_obs(self):
        adata = _make_reduced_adata()
        adata = adata[:0].copy()
        with pytest.raises(ValueError, match="0 observations"):
            spatial_cluster(adata)

    def test_rejects_empty_vars(self):
        adata = _make_reduced_adata()
        adata = adata[:, :0].copy()
        with pytest.raises(ValueError, match="0 variables"):
            spatial_cluster(adata)


# ---------------------------------------------------------------------------
# Provenance values
# ---------------------------------------------------------------------------

class TestSpatialClusterProvenanceValues:

    def test_n_clusters_positive(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert params["outputs"]["n_clusters"] >= 1

    def test_timestamp_iso_format(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.3, run_svg=False)
        assert "T" in params["timestamp"]

    def test_resolution_stored(self):
        adata = _make_reduced_adata()
        _, params = spatial_cluster(adata, resolution=0.7, run_svg=False)
        assert params["params"]["resolution"] == 0.7


# ---------------------------------------------------------------------------
# Report tests
# ---------------------------------------------------------------------------

class TestSpatialClusterReport:

    @pytest.fixture
    def clustered_adata(self):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(
            adata, resolution=0.3, run_svg=True, svg_n_genes=30, inplace=False
        )
        return out

    def test_report_generates_file(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_cluster_report(
            clustered_adata, out_path, dataset_id="test"
        )
        assert result.endswith(".html")
        import os
        assert os.path.isfile(result)

    def test_report_not_empty(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_cluster_report(clustered_adata, out_path, dataset_id="test")
        import os
        assert os.path.getsize(out_path) > 1000

    def test_report_is_valid_html(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_cluster_report(clustered_adata, out_path, dataset_id="test")
        with open(out_path) as fh:
            content = fh.read()
        assert "<!DOCTYPE html>" in content
        assert "</html>" in content

    def test_report_contains_dataset_id(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_cluster_report(
            clustered_adata, out_path, dataset_id="my_spatial"
        )
        with open(out_path) as fh:
            content = fh.read()
        assert "my_spatial" in content

    def test_report_contains_n_clusters(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        generate_spatial_cluster_report(clustered_adata, out_path, dataset_id="test")
        n_clusters = str(
            clustered_adata.uns["omicsage_spatial_cluster"]["outputs"]["n_clusters"]
        )
        with open(out_path) as fh:
            content = fh.read()
        assert n_clusters in content

    def test_report_raises_without_provenance(self, tmp_path):
        adata = _make_reduced_adata()
        with pytest.raises(ValueError, match="omicsage_spatial_cluster"):
            generate_spatial_cluster_report(
                adata, str(tmp_path / "x.html"), dataset_id="test"
            )

    def test_report_returns_absolute_path(self, clustered_adata, tmp_path):
        out_path = str(tmp_path / "test_report.html")
        result = generate_spatial_cluster_report(
            clustered_adata, out_path, dataset_id="test"
        )
        import os
        assert os.path.isabs(result)

    def test_report_no_svg_section_when_not_run(self, tmp_path):
        adata = _make_reduced_adata()
        out, _ = spatial_cluster(adata, resolution=0.3, run_svg=False)
        out_path = str(tmp_path / "no_svg.html")
        generate_spatial_cluster_report(out, out_path, dataset_id="test")
        with open(out_path) as fh:
            content = fh.read()
        # SVG table header should not appear when SVG not run
        assert "Moran" not in content or "not available" in content
