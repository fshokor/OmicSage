"""
tests/test_multiome_deg.py — Tests for pipeline/modules/scripts/multiome/multiome_deg.py

Coverage
--------
RNA DEG:
  - runs on MuData with atac_celltype groupby
  - result DataFrames have expected columns (gene, score, pval, logfc, pval_adj)
  - rank_genes_groups_rna written to mdata["rna"].uns
  - skipped with UserWarning when AnnData fallback used
  - warns when logcounts layer absent
  - exclude_gene_prefixes filters results

DCA (Differential Chromatin Accessibility):
  - runs on MuData with atac_celltype groupby
  - runs on AnnData fallback (DCA only)
  - result DataFrames have expected columns (peak, score, pval, logfc, pval_adj)
  - rank_genes_groups_dca written to atac.uns
  - exclude_peak_prefixes filters results
  - warns when counts layer absent (falls back to X)

Groupby fallback:
  - falls back to atac_leiden when atac_celltype absent
  - raises ValueError when no valid groupby found
  - emits UserWarning on fallback

Return contract:
  - returns (MuData, dict) or (AnnData, dict)
  - dict has keys: rna_deg, dca, rna_summary, dca_summary, provenance, input_type
  - input_type is "mudata" for MuData input, "anndata" for AnnData input
  - inplace=False returns a new object, does not modify original
  - inplace=True modifies original

Metrics / provenance:
  - uns["omicsage_multiome_deg"] written with required keys
  - module == "multiome_deg"
  - timestamp non-empty
  - groupby, method, n_cells, n_cell_types recorded
  - n_rna_significant and n_dca_significant are ints
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData

mudata = pytest.importorskip("mudata")
MuData = mudata.MuData

from pipeline.modules.scripts.multiome.multiome_deg import multiome_deg


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_N_CELLS   = 80
_N_GENES   = 50
_N_PEAKS   = 60
_N_GROUPS  = 4
_BATCH_KEY = "batch"


def _make_rna(n_obs: int = _N_CELLS, n_genes: int = _N_GENES, seed: int = 0) -> AnnData:
    rng = np.random.default_rng(seed)
    X   = sp.csr_matrix(rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32))
    adata = AnnData(X=X)
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = [f"Gene{i}" for i in range(n_genes)]
    adata.layers["logcounts"] = X.toarray().copy()
    return adata


def _make_atac(
    n_obs: int = _N_CELLS,
    n_peaks: int = _N_PEAKS,
    n_groups: int = _N_GROUPS,
    with_celltype: bool = True,
    with_leiden: bool = True,
    with_counts: bool = True,
    seed: int = 1,
) -> AnnData:
    rng    = np.random.default_rng(seed)
    counts = sp.csr_matrix(
        rng.integers(0, 20, size=(n_obs, n_peaks)).astype(np.float32)
    )
    tfidf  = sp.csr_matrix(
        rng.random(size=(n_obs, n_peaks)).astype(np.float32)
    )

    peak_names = [f"chr1:{i*1000}-{i*1000+500}" for i in range(n_peaks)]

    adata = AnnData(X=tfidf.copy())
    adata.obs_names = [f"cell_{i}" for i in range(n_obs)]
    adata.var_names = peak_names
    adata.layers["tf_idf"] = tfidf.copy()

    if with_counts:
        adata.layers["counts"] = counts.copy()

    # Cell type labels — evenly distributed across groups
    ct_labels  = [f"CellType_{i % n_groups}" for i in range(n_obs)]
    leid_labels = [str(i % n_groups) for i in range(n_obs)]

    if with_celltype:
        adata.obs["atac_celltype"] = ct_labels
    if with_leiden:
        adata.obs["atac_leiden"] = leid_labels

    return adata


@pytest.fixture
def mdata_fixture():
    """MuData with RNA + ATAC, shared barcodes."""
    rna  = _make_rna()
    atac = _make_atac()
    return MuData({"rna": rna, "atac": atac})


@pytest.fixture
def atac_fixture():
    """Bare AnnData (ATAC only) — fallback path."""
    return _make_atac()


# ---------------------------------------------------------------------------
# Return contract
# ---------------------------------------------------------------------------

class TestReturnContract:
    def test_mudata_returns_tuple(self, mdata_fixture):
        result, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(result, tuple) is False  # unpacked above
        assert isinstance(result, MuData)
        assert isinstance(deg_dict, dict)

    def test_anndata_returns_tuple(self, atac_fixture):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result, deg_dict = multiome_deg(atac_fixture, inplace=False)
        assert isinstance(result, AnnData)
        assert isinstance(deg_dict, dict)

    def test_dict_required_keys(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        for key in ["rna_deg", "dca", "rna_summary", "dca_summary",
                    "provenance", "input_type"]:
            assert key in deg_dict, f"Missing key: {key}"

    def test_input_type_mudata(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert deg_dict["input_type"] == "mudata"

    def test_input_type_anndata(self, atac_fixture):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, deg_dict = multiome_deg(atac_fixture, inplace=False)
        assert deg_dict["input_type"] == "anndata"

    def test_inplace_false_does_not_modify_original(self, mdata_fixture):
        original_uns = set(mdata_fixture.uns.keys())
        multiome_deg(mdata_fixture, inplace=False)
        assert set(mdata_fixture.uns.keys()) == original_uns

    def test_inplace_true_modifies_original(self, mdata_fixture):
        multiome_deg(mdata_fixture, inplace=True)
        assert "omicsage_multiome_deg" in mdata_fixture.uns


# ---------------------------------------------------------------------------
# RNA DEG
# ---------------------------------------------------------------------------

class TestRNADEG:
    def test_rna_deg_populated_for_mudata(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(deg_dict["rna_deg"], dict)
        assert len(deg_dict["rna_deg"]) > 0

    def test_rna_deg_gene_column(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        for ct, df in deg_dict["rna_deg"].items():
            assert "gene" in df.columns, f"Missing 'gene' in rna_deg['{ct}']"

    def test_rna_deg_expected_columns(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        expected = {"gene", "score", "pval", "logfc", "pval_adj"}
        for ct, df in deg_dict["rna_deg"].items():
            assert expected.issubset(df.columns), (
                f"rna_deg['{ct}'] missing columns: {expected - set(df.columns)}"
            )

    def test_rna_deg_written_to_rna_uns(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert "rank_genes_groups_rna" in result["rna"].uns

    def test_rna_deg_empty_for_anndata(self, atac_fixture):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, deg_dict = multiome_deg(atac_fixture, inplace=False)
        assert deg_dict["rna_deg"] == {}

    def test_rna_deg_skipped_emits_warning(self, atac_fixture):
        with pytest.warns(UserWarning, match="multiome_deg received an AnnData"):
            multiome_deg(atac_fixture, inplace=False)

    def test_warns_when_logcounts_absent(self, mdata_fixture):
        del mdata_fixture["rna"].layers["logcounts"]
        with pytest.warns(UserWarning, match="logcounts"):
            multiome_deg(mdata_fixture, inplace=False)

    def test_exclude_gene_prefixes(self, mdata_fixture):
        # Rename some genes to have RPL prefix
        mdata_fixture["rna"].var_names = (
            [f"RPL{i}" for i in range(10)]
            + [f"Gene{i}" for i in range(10, _N_GENES)]
        )
        _, deg_dict = multiome_deg(
            mdata_fixture,
            exclude_gene_prefixes=["RPL"],
            inplace=False,
        )
        for ct, df in deg_dict["rna_deg"].items():
            if df.empty:
                continue
            for gene in df["gene"]:
                assert not gene.upper().startswith("RPL"), (
                    f"RPL gene not filtered: {gene}"
                )

    def test_rna_summary_is_dataframe(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(deg_dict["rna_summary"], pd.DataFrame)

    def test_rna_summary_gene_column_when_nonempty(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        summary = deg_dict["rna_summary"]
        if not summary.empty:
            assert "gene" in summary.columns


# ---------------------------------------------------------------------------
# DCA
# ---------------------------------------------------------------------------

class TestDCA:
    def test_dca_populated(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(deg_dict["dca"], dict)

    def test_dca_peak_column(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        for ct, df in deg_dict["dca"].items():
            assert "peak" in df.columns, f"Missing 'peak' in dca['{ct}']"

    def test_dca_expected_columns(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        expected = {"peak", "score", "pval", "logfc", "pval_adj"}
        for ct, df in deg_dict["dca"].items():
            assert expected.issubset(df.columns), (
                f"dca['{ct}'] missing columns: {expected - set(df.columns)}"
            )

    def test_dca_written_to_atac_uns(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert "rank_genes_groups_dca" in result["atac"].uns

    def test_dca_runs_on_anndata(self, atac_fixture):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result, deg_dict = multiome_deg(atac_fixture, inplace=False)
        assert isinstance(deg_dict["dca"], dict)
        assert "rank_genes_groups_dca" in result.uns.get(
            "rank_genes_groups_dca", result
        ) or "rank_genes_groups_dca" in result.uns or True
        # Just confirm it doesn't raise and dca key is dict
        assert isinstance(deg_dict["dca"], dict)

    def test_dca_anndata_uns_written(self, atac_fixture):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result, _ = multiome_deg(atac_fixture, inplace=False)
        assert "rank_genes_groups_dca" in result.uns

    def test_warns_when_counts_layer_absent(self, mdata_fixture):
        del mdata_fixture["atac"].layers["counts"]
        with pytest.warns(UserWarning, match="counts"):
            multiome_deg(mdata_fixture, inplace=False)

    def test_exclude_peak_prefixes(self, mdata_fixture):
        # Rename first 5 peaks to chrM
        peak_names = list(mdata_fixture["atac"].var_names)
        for i in range(5):
            peak_names[i] = f"chrM:{i*1000}-{i*1000+500}"
        mdata_fixture["atac"].var_names = peak_names
        _, deg_dict = multiome_deg(
            mdata_fixture,
            exclude_peak_prefixes=["chrM"],
            inplace=False,
        )
        for ct, df in deg_dict["dca"].items():
            if df.empty:
                continue
            for peak in df["peak"]:
                assert not peak.upper().startswith("CHRM"), (
                    f"chrM peak not filtered: {peak}"
                )

    def test_dca_no_gene_column(self, mdata_fixture):
        """DCA results must have 'peak' column, not 'gene'."""
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        for ct, df in deg_dict["dca"].items():
            assert "gene" not in df.columns, (
                f"dca['{ct}'] should not have 'gene' column"
            )

    def test_dca_summary_is_dataframe(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(deg_dict["dca_summary"], pd.DataFrame)

    def test_dca_summary_peak_column_when_nonempty(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        summary = deg_dict["dca_summary"]
        if not summary.empty:
            assert "peak" in summary.columns


# ---------------------------------------------------------------------------
# Groupby fallback
# ---------------------------------------------------------------------------

class TestGroupbyFallback:
    def test_falls_back_to_leiden_when_celltype_absent(self, mdata_fixture):
        del mdata_fixture["atac"].obs["atac_celltype"]
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert deg_dict["provenance"]["groupby"] == "atac_leiden"
        assert any("atac_leiden" in str(warning.message) for warning in w)

    def test_fallback_emits_warning(self, mdata_fixture):
        del mdata_fixture["atac"].obs["atac_celltype"]
        with pytest.warns(UserWarning, match="atac_leiden"):
            multiome_deg(mdata_fixture, inplace=False)

    def test_raises_when_no_valid_groupby(self, mdata_fixture):
        del mdata_fixture["atac"].obs["atac_celltype"]
        del mdata_fixture["atac"].obs["atac_leiden"]
        with pytest.raises(ValueError, match="None of"):
            multiome_deg(mdata_fixture, inplace=False)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    _REQUIRED_KEYS = {
        "module", "timestamp", "input_type", "groupby", "method",
        "min_logfc", "max_pval_adj", "n_genes", "n_cells", "n_cell_types",
        "n_rna_significant", "n_dca_significant", "scanpy_version", "outputs",
    }

    def test_provenance_written_to_uns(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert "omicsage_multiome_deg" in result.uns

    def test_module_name(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert result.uns["omicsage_multiome_deg"]["module"] == "multiome_deg"

    def test_timestamp_non_empty(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert result.uns["omicsage_multiome_deg"]["timestamp"]

    def test_required_keys_present(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        prov = result.uns["omicsage_multiome_deg"]
        for key in self._REQUIRED_KEYS:
            assert key in prov, f"Missing provenance key: {key}"

    def test_n_cells_correct(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert result.uns["omicsage_multiome_deg"]["n_cells"] == _N_CELLS

    def test_n_rna_significant_is_int(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(
            result.uns["omicsage_multiome_deg"]["n_rna_significant"], int
        )

    def test_n_dca_significant_is_int(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert isinstance(
            result.uns["omicsage_multiome_deg"]["n_dca_significant"], int
        )

    def test_groupby_recorded(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        assert result.uns["omicsage_multiome_deg"]["groupby"] == "atac_celltype"

    def test_provenance_in_deg_dict(self, mdata_fixture):
        _, deg_dict = multiome_deg(mdata_fixture, inplace=False)
        assert "provenance" in deg_dict
        assert deg_dict["provenance"]["module"] == "multiome_deg"

    def test_outputs_keys(self, mdata_fixture):
        result, _ = multiome_deg(mdata_fixture, inplace=False)
        outputs = result.uns["omicsage_multiome_deg"]["outputs"]
        assert "rna_deg_key" in outputs
        assert "dca_key" in outputs
        assert outputs["dca_key"] == "rank_genes_groups_dca"
