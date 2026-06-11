"""
tests/test_spatial_ingest_qc.py — OmicSage Phase 7, Session 1
Tests for spatial_ingest.py and spatial_qc.py.

Run with:
    conda activate omicsage
    python -m pytest tests/test_spatial_ingest_qc.py -v
"""

from __future__ import annotations

import os
import tempfile

import numpy as np
import pytest

sq = pytest.importorskip("squidpy", reason="squidpy not installed")

import anndata as ad
import scipy.sparse as sp

from pipeline.modules.spatial.spatial_ingest import (
    _is_benchmark, _is_codex, _is_h5ad,
    _is_merfish, _is_visium, _is_visium_hd, _is_xenium,
    _resolve_spatial_type,
    _validate_spatial_adata,
    list_supported_types,
    spatial_ingest,
)
from pipeline.modules.spatial.spatial_qc import spatial_qc


# ===========================================================================
# Fixtures
# ===========================================================================


@pytest.fixture(scope="module")
def benchmark_adata() -> ad.AnnData:
    adata, _ = spatial_ingest("benchmark")
    return adata


@pytest.fixture
def minimal_spatial_adata() -> ad.AnnData:
    """Minimal synthetic Visium-like AnnData: 50 spots x 100 genes.

    Returns a fresh copy on every call so inplace=True tests cannot
    bleed state into subsequent tests.
    """
    import pandas as pd

    n_spots, n_genes = 50, 100
    rng = np.random.default_rng(42)
    X = sp.random(n_spots, n_genes, density=0.4, format="csr", random_state=42)
    X.data = np.round(X.data * 1000).astype(np.float32) + 1

    gene_names = [f"Gene{i:04d}" for i in range(n_genes - 5)]
    gene_names += [f"MT-Gene{i}" for i in range(5)]

    adata = ad.AnnData(
        X=X.copy(),
        obs=pd.DataFrame(index=[f"spot_{i}" for i in range(n_spots)]),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm["spatial"] = rng.uniform(0, 1000, size=(n_spots, 2))
    adata.uns["spatial"] = {
        "test_library": {
            "images": {},
            "scalefactors": {"spot_diameter_fullres": 89.0,
                             "tissue_hires_scalef": 0.2},
            "metadata": {},
        }
    }
    return adata.copy()


# ===========================================================================
# Auto-detection fingerprint functions
# ===========================================================================


class TestFingerprints:
    def test_is_benchmark(self):
        assert _is_benchmark("benchmark") is True
        assert _is_benchmark("other") is False

    def test_is_h5ad(self):
        assert _is_h5ad("/path/to/data.h5ad") is True
        assert _is_h5ad("/path/to/data.h5") is False

    def test_is_visium(self, tmp_path):
        (tmp_path / "spatial").mkdir()
        assert _is_visium(str(tmp_path)) is True

    def test_is_visium_false_without_spatial(self, tmp_path):
        assert _is_visium(str(tmp_path)) is False

    def test_is_visium_hd(self, tmp_path):
        (tmp_path / "binned_outputs").mkdir()
        assert _is_visium_hd(str(tmp_path)) is True

    def test_is_xenium(self, tmp_path):
        (tmp_path / "transcripts.parquet").touch()
        assert _is_xenium(str(tmp_path)) is True

    def test_is_merfish(self, tmp_path):
        (tmp_path / "cell_by_gene.csv").touch()
        assert _is_merfish(str(tmp_path)) is True

    def test_is_codex(self, tmp_path):
        f = tmp_path / "protein_data.csv"
        f.touch()
        assert _is_codex(str(f)) is True

    def test_is_codex_false_for_directory(self, tmp_path):
        assert _is_codex(str(tmp_path)) is False


# ===========================================================================
# _resolve_spatial_type
# ===========================================================================


class TestResolveSpatialType:
    def test_explicit_benchmark(self):
        assert _resolve_spatial_type("anything", "benchmark") == "benchmark"

    def test_explicit_visium(self):
        assert _resolve_spatial_type("anything", "visium") == "visium"

    def test_explicit_h5ad(self):
        assert _resolve_spatial_type("anything", "h5ad") == "h5ad"

    def test_explicit_xenium(self):
        assert _resolve_spatial_type("anything", "xenium") == "xenium"

    def test_explicit_merfish(self):
        assert _resolve_spatial_type("anything", "merfish") == "merfish"

    def test_explicit_codex(self):
        assert _resolve_spatial_type("anything", "codex") == "codex"

    def test_auto_benchmark(self):
        assert _resolve_spatial_type("benchmark", "auto") == "benchmark"

    def test_auto_h5ad(self):
        assert _resolve_spatial_type("/data/file.h5ad", "auto") == "h5ad"

    def test_auto_visium(self, tmp_path):
        (tmp_path / "spatial").mkdir()
        assert _resolve_spatial_type(str(tmp_path), "auto") == "visium"

    def test_auto_visium_hd(self, tmp_path):
        (tmp_path / "binned_outputs").mkdir()
        assert _resolve_spatial_type(str(tmp_path), "auto") == "visium_hd"

    def test_auto_xenium_beats_visium(self, tmp_path):
        # xenium dir may also have other subdirs — xenium should win
        (tmp_path / "transcripts.parquet").touch()
        (tmp_path / "spatial").mkdir()
        assert _resolve_spatial_type(str(tmp_path), "auto") == "xenium"

    def test_auto_merfish(self, tmp_path):
        (tmp_path / "cell_by_gene.csv").touch()
        assert _resolve_spatial_type(str(tmp_path), "auto") == "merfish"

    def test_auto_unknown_raises(self):
        with pytest.raises(ValueError, match="Cannot auto-detect"):
            _resolve_spatial_type("/nonexistent/unknown_file.xyz", "auto")


# ===========================================================================
# list_supported_types
# ===========================================================================


class TestListSupportedTypes:
    def test_returns_dict(self):
        result = list_supported_types()
        assert isinstance(result, dict)

    def test_implemented_types_present(self):
        result = list_supported_types()
        for t in ["benchmark", "h5ad", "visium"]:
            assert result[t] == "implemented"

    def test_planned_types_present(self):
        result = list_supported_types()
        # visium_hd and xenium are now implemented (Session B)
        for t in ["merfish", "codex"]:
            assert result[t] == "planned"
        # confirm the Session B formats are implemented
        assert result["visium_hd"] == "implemented"
        assert result["xenium"] == "implemented"


# ===========================================================================
# NotImplementedError for planned types
# ===========================================================================


class TestPlannedTypesRaiseNotImplemented:
    @pytest.mark.parametrize("stype", ["merfish", "codex"])  # visium_hd/xenium now implemented
    def test_planned_type_raises(self, stype):
        with pytest.raises(NotImplementedError, match="planned for a future"):
            spatial_ingest("/some/path", spatial_type=stype)

    @pytest.mark.parametrize("stype", ["merfish"])  # xenium/visium_hd now implemented
    def test_planned_type_error_includes_h5ad_workaround(self, stype):
        with pytest.raises(NotImplementedError, match="h5ad"):
            spatial_ingest("/some/path", spatial_type=stype)


# ===========================================================================
# _validate_spatial_adata
# ===========================================================================


class TestValidateSpatialAdata:
    def test_valid_passes(self, minimal_spatial_adata):
        _validate_spatial_adata(minimal_spatial_adata, "test")

    def test_missing_obsm_spatial(self, minimal_spatial_adata):
        del minimal_spatial_adata.obsm["spatial"]
        with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
            _validate_spatial_adata(minimal_spatial_adata, "test")

    def test_missing_uns_spatial(self, minimal_spatial_adata):
        del minimal_spatial_adata.uns["spatial"]
        with pytest.raises(ValueError, match="uns\\['spatial'\\]"):
            _validate_spatial_adata(minimal_spatial_adata, "test")

    def test_wrong_coord_shape(self, minimal_spatial_adata):
        minimal_spatial_adata.obsm["spatial"] = np.zeros((50, 3))
        with pytest.raises(ValueError, match="shape"):
            _validate_spatial_adata(minimal_spatial_adata, "test")


# ===========================================================================
# spatial_ingest — benchmark
# ===========================================================================


class TestSpatialIngestBenchmark:
    def test_returns_anndata(self, benchmark_adata):
        assert isinstance(benchmark_adata, ad.AnnData)

    def test_has_spatial_obsm(self, benchmark_adata):
        assert "spatial" in benchmark_adata.obsm
        assert benchmark_adata.obsm["spatial"].shape[1] == 2

    def test_has_spatial_uns(self, benchmark_adata):
        assert "spatial" in benchmark_adata.uns

    def test_provenance_keys(self, benchmark_adata):
        prov = benchmark_adata.uns["omicsage_spatial_ingest"]
        for key in ["source", "spatial_type", "technology_notes",
                    "n_obs", "n_vars", "timestamp"]:
            assert key in prov

    def test_provenance_spatial_type(self, benchmark_adata):
        prov = benchmark_adata.uns["omicsage_spatial_ingest"]
        assert prov["spatial_type"] == "benchmark"

    def test_technology_notes_non_empty(self, benchmark_adata):
        prov = benchmark_adata.uns["omicsage_spatial_ingest"]
        assert len(prov["technology_notes"]) > 0

    def test_n_obs_positive(self, benchmark_adata):
        assert benchmark_adata.n_obs > 0

    def test_n_vars_positive(self, benchmark_adata):
        assert benchmark_adata.n_vars > 0

    def test_provenance_n_obs_matches(self, benchmark_adata):
        prov = benchmark_adata.uns["omicsage_spatial_ingest"]
        assert prov["n_obs"] == benchmark_adata.n_obs


# ===========================================================================
# spatial_ingest — h5ad
# ===========================================================================


class TestSpatialIngestH5ad:
    def test_load_explicit(self, minimal_spatial_adata):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test.h5ad")
            minimal_spatial_adata.write_h5ad(path)
            adata, params = spatial_ingest(path, spatial_type="h5ad")
            assert isinstance(adata, ad.AnnData)
            assert params["spatial_type"] == "h5ad"

    def test_load_auto_detected(self, minimal_spatial_adata):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "data.h5ad")
            minimal_spatial_adata.write_h5ad(path)
            adata, params = spatial_ingest(path)   # auto
            assert params["spatial_type"] == "h5ad"

    def test_file_not_found_raises(self):
        with pytest.raises(FileNotFoundError):
            spatial_ingest("/nonexistent/file.h5ad", spatial_type="h5ad")

    def test_provenance_source_path(self, minimal_spatial_adata):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "data.h5ad")
            minimal_spatial_adata.write_h5ad(path)
            _, params = spatial_ingest(path, spatial_type="h5ad")
            assert params["source"] == path


# ===========================================================================
# spatial_ingest — visium directory errors
# ===========================================================================


class TestSpatialIngestVisium:
    def test_missing_dir_raises(self):
        with pytest.raises(NotADirectoryError):
            spatial_ingest("/nonexistent/spaceranger/", spatial_type="visium")

    def test_auto_detects_visium(self, tmp_path):
        (tmp_path / "spatial").mkdir()
        # Will fail at sq.read.visium level (no actual data),
        # but should resolve type correctly first
        with pytest.raises(Exception):
            spatial_ingest(str(tmp_path), spatial_type="auto")


# ===========================================================================
# spatial_ingest — unknown type
# ===========================================================================


class TestSpatialIngestUnknownType:
    def test_unknown_type_raises_value_error(self):
        with pytest.raises(ValueError, match="Unknown spatial_type"):
            spatial_ingest("/some/path", spatial_type="foobar")


# ===========================================================================
# spatial_qc — contract
# ===========================================================================


class TestSpatialQcContract:
    def test_returns_tuple(self, minimal_spatial_adata):
        result = spatial_qc(minimal_spatial_adata)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self, minimal_spatial_adata):
        adata_out, params = spatial_qc(minimal_spatial_adata)
        assert isinstance(adata_out, ad.AnnData)
        assert isinstance(params, dict)

    def test_provenance_stored(self, minimal_spatial_adata):
        adata_out, _ = spatial_qc(minimal_spatial_adata)
        assert "omicsage_spatial_qc" in adata_out.uns

    def test_provenance_keys(self, minimal_spatial_adata):
        _, params = spatial_qc(minimal_spatial_adata)
        assert "params" in params
        assert "outputs" in params
        assert "summary_stats" in params
        assert "timestamp" in params

    def test_inplace_false_no_side_effects(self, minimal_spatial_adata):
        n_obs = minimal_spatial_adata.n_obs
        spatial_qc(minimal_spatial_adata, inplace=False)
        assert minimal_spatial_adata.n_obs == n_obs
        assert "omicsage_spatial_qc" not in minimal_spatial_adata.uns

    def test_inplace_true_modifies_input(self, minimal_spatial_adata):
        spatial_qc(minimal_spatial_adata, inplace=True)
        assert "omicsage_spatial_qc" in minimal_spatial_adata.uns


# ===========================================================================
# spatial_qc — metrics
# ===========================================================================


class TestSpatialQcMetrics:
    def test_qc_columns_created(self, minimal_spatial_adata):
        adata_out, _ = spatial_qc(minimal_spatial_adata, filter_spots=False)
        for col in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
            assert col in adata_out.obs.columns

    def test_qc_pass_column_boolean(self, minimal_spatial_adata):
        adata_out, _ = spatial_qc(minimal_spatial_adata, filter_spots=False)
        assert "qc_pass" in adata_out.obs.columns
        assert adata_out.obs["qc_pass"].dtype == bool

    def test_total_counts_non_negative(self, minimal_spatial_adata):
        adata_out, _ = spatial_qc(minimal_spatial_adata, filter_spots=False)
        assert (adata_out.obs["total_counts"] >= 0).all()

    def test_pct_mt_in_range(self, minimal_spatial_adata):
        adata_out, _ = spatial_qc(minimal_spatial_adata, filter_spots=False)
        assert (adata_out.obs["pct_counts_mt"] >= 0).all()
        assert (adata_out.obs["pct_counts_mt"] <= 100).all()

    def test_mt_genes_detected(self, minimal_spatial_adata):
        _, params = spatial_qc(minimal_spatial_adata, mt_prefix="MT-",
                               filter_spots=False)
        assert params["outputs"]["n_mt_genes_detected"] == 5

    def test_mouse_mt_prefix_zero_genes(self, minimal_spatial_adata):
        _, params = spatial_qc(minimal_spatial_adata, mt_prefix="mt-",
                               filter_spots=False)
        assert params["outputs"]["n_mt_genes_detected"] == 0


# ===========================================================================
# spatial_qc — filtering
# ===========================================================================


class TestSpatialQcFiltering:
    def test_filter_reduces_spots(self, minimal_spatial_adata):
        # Get actual total_counts distribution, then cap with max_counts at
        # p25 - 1. This removes the top ~75% of spots and keeps ~25%,
        # using the upper cap rather than a lower bound to avoid the
        # right-skewed fixture distribution problem.
        adata_metrics, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=10_000_000,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=False,
        )
        p25 = float(adata_metrics.obs["total_counts"].quantile(0.25))
        cap = int(p25) - 1

        adata_out, params = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=cap,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=True,
        )
        assert adata_out.n_obs > 0, (
            f"All spots removed with max_counts={cap} "
            f"(p25={p25:.0f}). Check fixture count distribution."
        )
        assert adata_out.n_obs < minimal_spatial_adata.n_obs
        assert adata_out.n_obs == params["outputs"]["n_spots_after"]

    def test_no_filter_keeps_all(self, minimal_spatial_adata):
        n = minimal_spatial_adata.n_obs
        adata_out, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=10_000_000,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=False,
        )
        assert adata_out.n_obs == n

    def test_permissive_keeps_all(self, minimal_spatial_adata):
        n = minimal_spatial_adata.n_obs
        adata_out, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=10_000_000,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=True,
        )
        assert adata_out.n_obs == n

    def test_impossible_threshold_removes_all(self, minimal_spatial_adata):
        adata_out, params = spatial_qc(
            minimal_spatial_adata, min_counts=999_999_999, filter_spots=True
        )
        assert adata_out.n_obs == 0
        assert params["outputs"]["n_spots_after"] == 0

    def test_removed_count_consistent(self, minimal_spatial_adata):
        adata_metrics, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=10_000_000,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=False,
        )
        cap = int(adata_metrics.obs["total_counts"].quantile(0.25)) - 1
        _, params = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=cap,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=True,
        )
        o = params["outputs"]
        assert o["n_spots_before"] == o["n_spots_after"] + o["n_spots_removed"]

    def test_spatial_coords_preserved(self, minimal_spatial_adata):
        adata_metrics, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=10_000_000,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=False,
        )
        cap = int(adata_metrics.obs["total_counts"].quantile(0.25)) - 1
        adata_out, _ = spatial_qc(
            minimal_spatial_adata, min_counts=0, max_counts=cap,
            min_genes=0, max_genes=1_000_000, max_mt_pct=100.0,
            filter_spots=True,
        )
        assert "spatial" in adata_out.obsm
        assert adata_out.obsm["spatial"].shape == (adata_out.n_obs, 2)


# ===========================================================================
# spatial_qc — validation
# ===========================================================================


class TestSpatialQcValidation:
    def test_missing_spatial_raises(self, minimal_spatial_adata):
        del minimal_spatial_adata.obsm["spatial"]
        with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
            spatial_qc(minimal_spatial_adata)

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError):
            spatial_qc("not_an_anndata")

    def test_empty_obs_raises(self):
        import pandas as pd
        empty = ad.AnnData(
            X=np.zeros((0, 10)),
            obs=pd.DataFrame(index=[]),
            var=pd.DataFrame(index=[f"G{i}" for i in range(10)]),
        )
        empty.obsm["spatial"] = np.zeros((0, 2))
        empty.uns["spatial"] = {}
        with pytest.raises(ValueError, match="0 observations"):
            spatial_qc(empty)


# ===========================================================================
# spatial_qc — summary stats
# ===========================================================================


class TestSpatialQcSummaryStats:
    def test_summary_stats_keys(self, minimal_spatial_adata):
        _, params = spatial_qc(minimal_spatial_adata, filter_spots=False)
        for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
            s = params["summary_stats"][metric]
            for stat in ["mean", "median", "std", "min", "max"]:
                assert stat in s

    def test_summary_stats_are_floats(self, minimal_spatial_adata):
        _, params = spatial_qc(minimal_spatial_adata, filter_spots=False)
        for metric_stats in params["summary_stats"].values():
            for val in metric_stats.values():
                assert isinstance(val, float)


# ===========================================================================
# Integration: ingest → qc
# ===========================================================================


class TestIngestQcPipeline:
    def test_benchmark_pipeline(self):
        adata, ingest_params = spatial_ingest("benchmark")
        assert ingest_params["spatial_type"] == "benchmark"

        adata_qc, qc_params = spatial_qc(
            adata, mt_prefix="mt-", min_counts=100, filter_spots=True
        )
        assert isinstance(adata_qc, ad.AnnData)
        assert adata_qc.n_obs > 0
        assert "spatial" in adata_qc.obsm
        assert "omicsage_spatial_ingest" in adata_qc.uns
        assert "omicsage_spatial_qc" in adata_qc.uns

    def test_h5ad_pipeline(self, minimal_spatial_adata):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "spot_data.h5ad")
            minimal_spatial_adata.write_h5ad(path)

            adata, _ = spatial_ingest(path)
            adata_qc, params = spatial_qc(
                adata, mt_prefix="MT-", min_counts=0, filter_spots=False
            )
            assert adata_qc.n_obs == minimal_spatial_adata.n_obs
            assert "total_counts" in adata_qc.obs.columns

    def test_provenance_chain_intact(self, minimal_spatial_adata):
        """Both provenance keys survive the full ingest → qc pipeline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "data.h5ad")
            minimal_spatial_adata.write_h5ad(path)
            adata, _ = spatial_ingest(path)
            adata_qc, _ = spatial_qc(adata, filter_spots=False)
            assert "omicsage_spatial_ingest" in adata_qc.uns
            assert "omicsage_spatial_qc" in adata_qc.uns
