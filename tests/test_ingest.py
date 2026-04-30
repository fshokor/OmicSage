"""
test_ingest.py
==============
Tests for pipeline/modules/qc/ingest.py

Test coverage
-------------
- Format detection (.h5ad, .h5, MTX directory, unsupported)
- MTX directory loading (GSE166635 HCC1 real data)
- h5ad loading with raw count extraction (GSE194122 CITE)
- AnnData output contract (shape, dtype, uns metadata)
- Edge cases (missing file, bad format)

Run
---
    pytest tests/test_ingest.py -v
    pytest tests/test_ingest.py -v -k "mtx"          # MTX tests only
    pytest tests/test_ingest.py -v -k "h5ad"         # h5ad tests only
    pytest tests/test_ingest.py -v --no-header -rN   # minimal output
"""

import numpy as np
import pytest
import scipy.sparse as sp
from pathlib import Path
from unittest.mock import patch, MagicMock
import anndata as ad

from pipeline.modules.qc.ingest import (
    load_dataset,
    detect_format,
    _is_integer_matrix,
    _find_raw_layer,
)

# ── Paths ──────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parent.parent
TEST_MTX    = REPO_ROOT / "data" / "test" / "GSE166635" / "HCC1"
CITE_H5AD   = REPO_ROOT / "data" / "benchmark" / \
              "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"


# ── Fixtures ───────────────────────────────────────────────────────────────

def make_raw_adata(n_obs=50, n_vars=100) -> ad.AnnData:
    """AnnData with integer counts in X (simulates MTX/H5 input)."""
    X = sp.random(n_obs, n_vars, density=0.1, format="csr")
    X.data = np.floor(X.data * 100).astype(np.float32)
    return ad.AnnData(X=X)


def make_processed_adata(n_obs=50, n_vars=100) -> ad.AnnData:
    """AnnData with normalized X and raw counts in layers['counts']."""
    raw   = sp.random(n_obs, n_vars, density=0.1, format="csr")
    raw.data = np.floor(raw.data * 100).astype(np.float32)
    norm  = raw.copy()
    norm.data = norm.data / norm.data.max()          # fake normalization

    adata = ad.AnnData(X=norm)
    adata.layers["counts"] = raw
    return adata


def make_processed_adata_spliced(n_obs=50, n_vars=100) -> ad.AnnData:
    """AnnData with normalized X and raw counts in layers['spliced']."""
    raw  = sp.random(n_obs, n_vars, density=0.1, format="csr")
    raw.data = np.floor(raw.data * 50).astype(np.float32)
    norm = raw.copy()
    norm.data = norm.data / norm.data.max()

    adata = ad.AnnData(X=norm)
    adata.layers["spliced"] = raw
    return adata


# ── Format detection ────────────────────────────────────────────────────────

class TestDetectFormat:

    def test_h5ad_extension(self, tmp_path):
        f = tmp_path / "sample.h5ad"
        f.touch()
        assert detect_format(f) == "h5ad"

    def test_h5_extension(self, tmp_path):
        f = tmp_path / "sample.h5"
        f.touch()
        assert detect_format(f) == "h5"

    def test_hdf5_extension(self, tmp_path):
        f = tmp_path / "sample.hdf5"
        f.touch()
        assert detect_format(f) == "h5"

    def test_mtx_directory(self, tmp_path):
        (tmp_path / "matrix.mtx.gz").touch()
        assert detect_format(tmp_path) == "mtx_dir"

    def test_mtx_uncompressed(self, tmp_path):
        (tmp_path / "matrix.mtx").touch()
        assert detect_format(tmp_path) == "mtx_dir"

    def test_unsupported_extension_raises(self, tmp_path):
        f = tmp_path / "sample.csv"
        f.touch()
        with pytest.raises(ValueError, match="Unsupported file format"):
            detect_format(f)

    def test_empty_directory_raises(self, tmp_path):
        with pytest.raises(ValueError, match="10x MEX"):
            detect_format(tmp_path)


# ── Integer matrix detection ─────────────────────────────────────────────────

class TestIsIntegerMatrix:

    def test_integer_sparse(self):
        X = sp.csr_matrix(np.array([[1, 0, 3], [0, 5, 0]], dtype=np.float32))
        assert _is_integer_matrix(X) is True

    def test_float_sparse(self):
        X = sp.csr_matrix(np.array([[0.1, 0, 0.3]], dtype=np.float32))
        assert _is_integer_matrix(X) is False

    def test_integer_dense(self):
        X = np.array([[2.0, 0.0], [0.0, 4.0]])
        assert _is_integer_matrix(X) is True

    def test_float_dense(self):
        X = np.array([[0.5, 0.1], [0.0, 0.9]])
        assert _is_integer_matrix(X) is False


# ── Raw layer detection ──────────────────────────────────────────────────────

class TestFindRawLayer:

    def test_finds_counts_layer(self):
        adata = make_processed_adata()
        assert _find_raw_layer(adata) == "counts"

    def test_finds_spliced_layer(self):
        adata = make_processed_adata_spliced()
        assert _find_raw_layer(adata) == "spliced"

    def test_returns_none_when_no_raw(self):
        X = sp.random(10, 20, density=0.3, format="csr")
        X.data = np.random.rand(X.nnz).astype(np.float32)
        adata = ad.AnnData(X=X)
        assert _find_raw_layer(adata) is None


# ── load_dataset output contract ─────────────────────────────────────────────

class TestLoadDatasetContract:

    def test_missing_file_raises(self):
        with pytest.raises(FileNotFoundError):
            load_dataset("/nonexistent/path/file.h5ad", verbose=False)

    def test_x_is_sparse(self, tmp_path):
        """Output X must always be sparse CSR."""
        adata_in = make_raw_adata()
        h5ad_path = tmp_path / "test.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, verbose=False)
        assert sp.issparse(adata_out.X)

    def test_x_contains_integers(self, tmp_path):
        """X must contain integer (raw) counts after ingestion."""
        adata_in = make_processed_adata()
        h5ad_path = tmp_path / "processed.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, verbose=False)
        assert _is_integer_matrix(adata_out.X), \
            "adata.X should contain integer counts after ingestion"

    def test_uns_metadata_populated(self, tmp_path):
        """omicsage_source must be in uns with required keys."""
        adata_in = make_raw_adata()
        h5ad_path = tmp_path / "test.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, verbose=False)
        src = adata_out.uns["omicsage_source"]

        assert "file"      in src
        assert "format"    in src
        assert "raw_layer" in src
        assert "sample"    in src
        assert src["format"] == "h5ad"

    def test_sample_column_added(self, tmp_path):
        """obs must contain a 'sample' column."""
        adata_in = make_raw_adata()
        h5ad_path = tmp_path / "test.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, sample_name="MySample", verbose=False)
        assert "sample" in adata_out.obs.columns
        assert (adata_out.obs["sample"] == "MySample").all()

    def test_normalized_layer_preserved(self, tmp_path):
        """When raw extracted from layers, normalized X saved to layers['normalized']."""
        adata_in = make_processed_adata()
        h5ad_path = tmp_path / "processed.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, verbose=False)
        assert "normalized" in adata_out.layers, \
            "layers['normalized'] should contain the original normalized X"

    def test_shape_preserved(self, tmp_path):
        """Cell and gene counts must be preserved."""
        adata_in = make_raw_adata(n_obs=80, n_vars=200)
        h5ad_path = tmp_path / "test.h5ad"
        adata_in.write_h5ad(h5ad_path)

        adata_out = load_dataset(h5ad_path, verbose=False)
        assert adata_out.n_obs == 80
        assert adata_out.n_vars == 200


# ── Real data tests (skipped if files not present) ──────────────────────────

@pytest.mark.skipif(
    not TEST_MTX.exists(),
    reason="GSE166635 test data not downloaded. Run: python scripts/download_test_data.py"
)
class TestRealMTX:

    def test_mtx_loads_successfully(self):
        adata = load_dataset(TEST_MTX, verbose=False)
        assert adata.n_obs > 0
        assert adata.n_vars > 0

    def test_mtx_x_is_integer(self):
        adata = load_dataset(TEST_MTX, verbose=False)
        assert _is_integer_matrix(adata.X), \
            "MTX data should be raw integer counts"

    def test_mtx_format_in_uns(self):
        adata = load_dataset(TEST_MTX, verbose=False)
        assert adata.uns["omicsage_source"]["format"] == "mtx_dir"

    def test_mtx_reasonable_cell_count(self):
        """HCC1 sample should have thousands of cells."""
        adata = load_dataset(TEST_MTX, verbose=False)
        assert adata.n_obs > 100, "Expected at least 100 cells in HCC1"


@pytest.mark.skipif(
    not CITE_H5AD.exists(),
    reason="GSE194122 CITE file not found in data/benchmark/"
)
class TestRealCiteH5AD:

    def test_cite_loads_successfully(self):
        adata = load_dataset(CITE_H5AD, verbose=False)
        assert adata.n_obs > 0
        assert adata.n_vars > 0

    def test_cite_raw_extracted_from_counts_layer(self):
        """GSE194122 stores raw in layers['counts'] — must be moved to X."""
        adata = load_dataset(CITE_H5AD, verbose=False)
        assert _is_integer_matrix(adata.X), \
            "After ingestion, X should be integer counts from layers['counts']"

    def test_cite_normalized_layer_exists(self):
        """The original normalized X should be preserved."""
        adata = load_dataset(CITE_H5AD, verbose=False)
        assert "normalized" in adata.layers

    def test_cite_raw_layer_recorded_in_uns(self):
        """uns should record which layer was used for raw counts."""
        adata = load_dataset(CITE_H5AD, verbose=False)
        src = adata.uns["omicsage_source"]
        assert src["raw_layer"] == "counts"

    def test_cite_cell_count_plausible(self):
        """CITE file has ~70k cells."""
        adata = load_dataset(CITE_H5AD, verbose=False)
        assert adata.n_obs > 50_000, \
            f"Expected ~70k cells, got {adata.n_obs}"
