"""
tests/test_spatial_ingest_formats.py — OmicSage Phase 7 extension (Session B)

Tests for the Visium HD and Xenium loader implementations in
pipeline/modules/scripts/spatial/spatial_ingest.py.

Strategy
--------
- spatialdata-io is NOT required in CI — all tests mock the import.
- The SpatialData object returned by spatialdata-io is mocked as a simple
  namespace with a ``tables`` dict containing a minimal AnnData.
- Tests verify:
    - Output contract: obsm["spatial"], uns["spatial"], layers["counts"]
    - Provenance: uns["omicsage_spatial_ingest"] keys
    - Auto-detection fingerprints for visium_hd and xenium
    - bin_size routing: correct table key selected from SpatialData
    - Missing table key raises KeyError with a useful message
    - ImportError raised cleanly when spatialdata-io is absent
    - NotImplementedError still raised for merfish / codex (not implemented)
    - list_supported_types() now shows visium_hd and xenium as "implemented"
"""

from __future__ import annotations

import os
import sys
import types
from pathlib import Path
from unittest.mock import MagicMock, patch

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def rng():
    return np.random.default_rng(42)


def _minimal_adata(n_obs: int = 20, n_genes: int = 50, rng=None) -> ad.AnnData:
    """Minimal AnnData with obsm['spatial'] set — simulates what spatialdata-io returns."""
    if rng is None:
        rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.integers(0, 100, (n_obs, n_genes)).astype(np.float32))
    gene_names = [f"GENE{i:04d}" for i in range(n_genes)]
    obs_names = [str(i + 1) for i in range(n_obs)]  # Xenium uses integer-string cell IDs
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=gene_names),
    )
    adata.obsm["spatial"] = rng.random((n_obs, 2)).astype(np.float32)
    return adata


def _make_sdata_mock(table_key: str, adata: ad.AnnData) -> MagicMock:
    """Minimal SpatialData mock with a tables dict."""
    sdata = MagicMock()
    sdata.tables = {table_key: adata}
    return sdata


def _make_spatialdata_io_mock(func_name: str, table_key: str, adata: ad.AnnData):
    """Return a mock spatialdata_io module with the named function patched."""
    mod = types.ModuleType("spatialdata_io")
    sdata = _make_sdata_mock(table_key, adata)
    setattr(mod, func_name, lambda *a, **kw: sdata)
    return mod


# ---------------------------------------------------------------------------
# Import helper
# ---------------------------------------------------------------------------

def _import_ingest():
    if "pipeline.modules.scripts.spatial.spatial_ingest" in sys.modules:
        import importlib
        return importlib.reload(sys.modules["pipeline.modules.scripts.spatial.spatial_ingest"])
    import pipeline.modules.scripts.spatial.spatial_ingest as m
    return m


# ---------------------------------------------------------------------------
# Auto-detection fingerprint tests
# ---------------------------------------------------------------------------

class TestAutoDetection:
    def test_visium_hd_detected_from_binned_outputs_dir(self, tmp_path):
        m = _import_ingest()
        src = tmp_path / "visium_hd_sample"
        (src / "binned_outputs").mkdir(parents=True)
        assert m._is_visium_hd(str(src)) is True

    def test_visium_hd_not_detected_for_plain_visium(self, tmp_path):
        m = _import_ingest()
        src = tmp_path / "visium_sample"
        (src / "spatial").mkdir(parents=True)
        # has spatial/ but NOT binned_outputs/ → not visium_hd
        assert m._is_visium_hd(str(src)) is False

    def test_xenium_detected_from_transcripts_parquet(self, tmp_path):
        m = _import_ingest()
        src = tmp_path / "xenium_sample"
        src.mkdir()
        (src / "transcripts.parquet").touch()
        assert m._is_xenium(str(src)) is True

    def test_xenium_not_detected_without_transcripts_parquet(self, tmp_path):
        m = _import_ingest()
        src = tmp_path / "xenium_sample"
        src.mkdir()
        assert m._is_xenium(str(src)) is False

    def test_visium_hd_wins_over_visium_in_auto(self, tmp_path):
        """A dir with both spatial/ and binned_outputs/ must be detected as visium_hd."""
        m = _import_ingest()
        src = tmp_path / "hd_sample"
        (src / "spatial").mkdir(parents=True)
        (src / "binned_outputs").mkdir()
        # xenium check (transcripts.parquet absent) → fails
        # merfish check (cell_by_gene.csv absent) → fails
        # visium_hd check → wins (has binned_outputs)
        detected = m._resolve_spatial_type(str(src), "auto")
        assert detected == "visium_hd"

    def test_xenium_wins_over_visium_in_auto(self, tmp_path):
        """A dir with both spatial/ and transcripts.parquet → xenium."""
        m = _import_ingest()
        src = tmp_path / "xenium_sample"
        (src / "spatial").mkdir(parents=True)
        (src / "transcripts.parquet").touch()
        detected = m._resolve_spatial_type(str(src), "auto")
        assert detected == "xenium"


# ---------------------------------------------------------------------------
# Visium HD loader tests
# ---------------------------------------------------------------------------

class TestVisiumHD:
    def _run(self, tmp_path, rng, bin_size=8):
        """Run _load_visium_hd with mocked spatialdata-io."""
        m = _import_ingest()
        n_obs = 25
        adata = _minimal_adata(n_obs=n_obs, rng=rng)
        table_key = f"square_{bin_size:03d}um"
        mock_sio = _make_spatialdata_io_mock("visium_hd", table_key, adata)

        # Create a minimal directory structure
        src = tmp_path / "visium_hd_data"
        (src / "binned_outputs").mkdir(parents=True)

        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            adata_out, lib_id, source_repr = m._load_visium_hd(
                str(src),
                counts_file="filtered_feature_bc_matrix.h5",
                library_id=None,
                load_images=True,
                bin_size=bin_size,
            )
        return adata_out, lib_id, source_repr

    def test_output_obsm_spatial(self, tmp_path, rng):
        adata_out, _, _ = self._run(tmp_path, rng)
        assert "spatial" in adata_out.obsm
        assert adata_out.obsm["spatial"].shape[1] == 2

    def test_output_uns_spatial_contract(self, tmp_path, rng):
        adata_out, lib_id, _ = self._run(tmp_path, rng)
        assert "spatial" in adata_out.uns
        assert lib_id in adata_out.uns["spatial"]
        assert "scalefactors" in adata_out.uns["spatial"][lib_id]

    def test_output_layers_counts(self, tmp_path, rng):
        adata_out, _, _ = self._run(tmp_path, rng)
        assert "counts" in adata_out.layers

    def test_bin_size_8_default(self, tmp_path, rng):
        adata_out, _, _ = self._run(tmp_path, rng, bin_size=8)
        assert adata_out.uns["spatial"][list(adata_out.uns["spatial"].keys())[0]][
            "metadata"
        ]["bin_size_um"] == 8

    def test_bin_size_16_routing(self, tmp_path, rng):
        """bin_size=16 must select table key 'square_016um'."""
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        table_key = "square_016um"
        mock_sio = _make_spatialdata_io_mock("visium_hd", table_key, adata)
        src = tmp_path / "hd16"
        (src / "binned_outputs").mkdir(parents=True)
        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            adata_out, _, _ = m._load_visium_hd(
                str(src), counts_file="", library_id=None, load_images=True, bin_size=16
            )
        assert adata_out is not None

    def test_missing_table_key_raises_key_error(self, tmp_path, rng):
        """If the requested bin_size table is absent, raise KeyError with message."""
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        # Table only has 002um, but we request 008um
        mock_sio = _make_spatialdata_io_mock("visium_hd", "square_002um", adata)
        src = tmp_path / "hd_wrong_bin"
        (src / "binned_outputs").mkdir(parents=True)
        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            with pytest.raises(KeyError, match="bin_size=8"):
                m._load_visium_hd(
                    str(src), counts_file="", library_id=None, load_images=True, bin_size=8
                )

    def test_import_error_without_spatialdata_io(self, tmp_path):
        m = _import_ingest()
        # Remove spatialdata_io from sys.modules to simulate missing package
        with patch.dict(sys.modules, {"spatialdata_io": None}):
            with pytest.raises((ImportError, TypeError)):
                m._load_visium_hd(
                    str(tmp_path), counts_file="", library_id=None,
                    load_images=True, bin_size=8
                )


# ---------------------------------------------------------------------------
# Xenium loader tests
# ---------------------------------------------------------------------------

class TestXenium:
    def _run(self, tmp_path, rng):
        """Run _load_xenium with mocked spatialdata-io."""
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        mock_sio = _make_spatialdata_io_mock("xenium", "table", adata)
        src = tmp_path / "xenium_data"
        src.mkdir()
        (src / "transcripts.parquet").touch()
        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            adata_out, lib_id, source_repr = m._load_xenium(
                str(src),
                counts_file="",
                library_id=None,
                load_images=True,
                bin_size=None,
            )
        return adata_out, lib_id, source_repr

    def test_output_obsm_spatial(self, tmp_path, rng):
        adata_out, _, _ = self._run(tmp_path, rng)
        assert "spatial" in adata_out.obsm
        assert adata_out.obsm["spatial"].shape[1] == 2

    def test_output_uns_spatial_contract(self, tmp_path, rng):
        adata_out, lib_id, _ = self._run(tmp_path, rng)
        assert "spatial" in adata_out.uns
        assert lib_id in adata_out.uns["spatial"]
        assert "scalefactors" in adata_out.uns["spatial"][lib_id]

    def test_output_layers_counts(self, tmp_path, rng):
        adata_out, _, _ = self._run(tmp_path, rng)
        assert "counts" in adata_out.layers

    def test_uns_metadata_platform_xenium(self, tmp_path, rng):
        adata_out, lib_id, _ = self._run(tmp_path, rng)
        meta = adata_out.uns["spatial"][lib_id].get("metadata", {})
        assert meta.get("platform") == "xenium"

    def test_missing_table_raises_key_error(self, tmp_path, rng):
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        # Provide wrong key so lookup fails
        mock_sio = _make_spatialdata_io_mock("xenium", "wrong_key", adata)
        src = tmp_path / "xenium_bad"
        src.mkdir()
        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            with pytest.raises(KeyError, match="table"):
                m._load_xenium(
                    str(src), counts_file="", library_id=None,
                    load_images=True, bin_size=None
                )

    def test_import_error_without_spatialdata_io(self, tmp_path):
        m = _import_ingest()
        with patch.dict(sys.modules, {"spatialdata_io": None}):
            with pytest.raises((ImportError, TypeError)):
                m._load_xenium(
                    str(tmp_path), counts_file="", library_id=None,
                    load_images=True, bin_size=None
                )


# ---------------------------------------------------------------------------
# Registry / list_supported_types tests
# ---------------------------------------------------------------------------

class TestRegistry:
    def test_visium_hd_now_implemented(self):
        m = _import_ingest()
        supported = m.list_supported_types()
        assert supported["visium_hd"] == "implemented"

    def test_xenium_now_implemented(self):
        m = _import_ingest()
        supported = m.list_supported_types()
        assert supported["xenium"] == "implemented"

    def test_merfish_still_planned(self):
        m = _import_ingest()
        assert m.list_supported_types()["merfish"] == "planned"

    def test_codex_still_planned(self):
        m = _import_ingest()
        assert m.list_supported_types()["codex"] == "planned"


# ---------------------------------------------------------------------------
# Still-planned format error tests
# ---------------------------------------------------------------------------

class TestStillPlanned:
    def test_merfish_raises_not_implemented(self, tmp_path):
        m = _import_ingest()
        with pytest.raises(NotImplementedError, match="MERFISH"):
            m._load_merfish(
                str(tmp_path), counts_file="", library_id=None,
                load_images=True, bin_size=None
            )

    def test_codex_raises_not_implemented(self, tmp_path):
        m = _import_ingest()
        with pytest.raises(NotImplementedError, match="CODEX"):
            m._load_codex(
                str(tmp_path), counts_file="", library_id=None,
                load_images=True, bin_size=None
            )


# ---------------------------------------------------------------------------
# Public API integration: spatial_ingest() routes correctly
# ---------------------------------------------------------------------------

class TestPublicAPIRouting:
    def test_spatial_ingest_visium_hd_routed(self, tmp_path, rng):
        """spatial_ingest() with spatial_type='visium_hd' calls _load_visium_hd."""
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        mock_sio = _make_spatialdata_io_mock("visium_hd", "square_008um", adata)

        src = tmp_path / "hd_routed"
        (src / "binned_outputs").mkdir(parents=True)

        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            adata_out, prov = m.spatial_ingest(
                str(src),
                spatial_type="visium_hd",
                bin_size=8,
            )

        assert prov["spatial_type"] == "visium_hd"
        assert prov["bin_size"] == 8
        assert "spatial" in adata_out.obsm
        assert "omicsage_spatial_ingest" in adata_out.uns

    def test_spatial_ingest_xenium_routed(self, tmp_path, rng):
        """spatial_ingest() with spatial_type='xenium' calls _load_xenium."""
        m = _import_ingest()
        adata = _minimal_adata(rng=rng)
        mock_sio = _make_spatialdata_io_mock("xenium", "table", adata)

        src = tmp_path / "xenium_routed"
        src.mkdir()
        (src / "transcripts.parquet").touch()

        with patch.dict(sys.modules, {"spatialdata_io": mock_sio}):
            adata_out, prov = m.spatial_ingest(
                str(src),
                spatial_type="xenium",
            )

        assert prov["spatial_type"] == "xenium"
        assert "spatial" in adata_out.obsm
        assert "omicsage_spatial_ingest" in adata_out.uns
