"""
tests/test_multiome_grn.py
Tests for pipeline/modules/multiome/multiome_grn.py

Coverage:
  Happy path
  - Returns (MuData, dict) tuple
  - uns["omicsage_grn"] written
  - uns["grn_network"] written
  - grn_network has required keys
  - n_grn_edges >= 0
  - inplace=False leaves original untouched
  - inplace=True modifies in place

  Metrics dict
  - Required keys present: provenance, n_tfs_rna, n_tfs_atac, n_grn_edges, grn_df
  - n_grn_edges matches grn_df length
  - grn_df is a DataFrame
  - grn_df has expected columns

  Provenance
  - module == "multiome_grn"
  - timestamp present
  - params recorded (motif_db, groupby, n_top_peaks, min_cells)
  - outputs recorded (n_tfs_rna, n_tfs_atac, n_grn_edges)

  Validation
  - Missing "rna" modality raises KeyError
  - Missing "atac" modality raises KeyError
  - Missing groupby column raises KeyError
  - Fewer cells than min_cells raises ValueError

  GRN table
  - grn_df columns: tf, target_gene, rna_score, atac_score, combined_score, cell_type
  - combined_score >= 0 for all rows
  - cell_type values are subset of groupby unique values

  Internal helpers
  - _extract_top_peaks returns list
  - _get_tf_list returns list (handles missing decoupler)
  - _regulons_correlation returns list of dicts with tf/targets/weights keys
  - _build_grn_table returns DataFrame even when no scores available
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from anndata import AnnData
from mudata import MuData

# ---------------------------------------------------------------------------
# Optional dependency availability flags
# ---------------------------------------------------------------------------

try:
    import pyscenic         # noqa: F401
    _PYSCENIC = True
except ModuleNotFoundError:
    _PYSCENIC = False

try:
    import decoupler        # noqa: F401
    _DECOUPLER = True
except ModuleNotFoundError:
    _DECOUPLER = False

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

N_CELLS   = 60
N_GENES   = 40
N_PEAKS   = 80
N_BATCHES = 2
GROUPBY   = "atac_celltype"


def _make_rna(n_cells: int = N_CELLS, seed: int = 0) -> AnnData:
    rng  = np.random.default_rng(seed)
    raw  = rng.integers(1, 100, size=(n_cells, N_GENES)).astype(np.float32)
    logX = np.log1p(raw)
    adata = AnnData(X=sp.csr_matrix(logX))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = [f"gene_{i}" for i in range(N_GENES)]
    adata.obs["batch"] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.obs[GROUPBY] = (["typeA", "typeB"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"] = sp.csr_matrix(raw)
    return adata


def _make_atac(n_cells: int = N_CELLS, seed: int = 1) -> AnnData:
    rng    = np.random.default_rng(seed)
    counts = rng.integers(0, 10, size=(n_cells, N_PEAKS)).astype(np.float32)
    peak_names = [f"chr1:{i*500}-{i*500+499}" for i in range(N_PEAKS)]
    adata = AnnData(X=sp.csr_matrix(counts))
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]
    adata.var_names = peak_names
    adata.obs["batch"] = (["d0", "d1"] * (n_cells // 2 + 1))[:n_cells]
    adata.obs[GROUPBY] = (["typeA", "typeB"] * (n_cells // 2 + 1))[:n_cells]
    adata.layers["counts"] = sp.csr_matrix(counts)

    # Add mock DCA results (rank_genes_groups_dca)
    groups = ["typeA", "typeB"]
    names_arr = np.array(
        [[peak_names[i] for _ in groups] for i in range(N_PEAKS)],
        dtype=object,
    )
    import numpy.lib.recfunctions as rfn
    names_rec = np.core.records.fromarrays(
        [names_arr[:, j] for j in range(len(groups))],
        names=",".join(groups),
    )
    adata.uns["rank_genes_groups_dca"] = {
        "params": {"groupby": GROUPBY, "method": "wilcoxon"},
        "names": names_rec,
    }
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


def _make_deg_dict() -> dict:
    return {
        "provenance": {
            "module": "multiome_deg",
            "timestamp": "2026-01-01T00:00:00",
            "params": {},
            "outputs": {},
        },
        "n_rna_significant": 10,
        "n_dca_significant": 20,
        "n_cell_types": 2,
    }


# ---------------------------------------------------------------------------
# Import the module under test
# ---------------------------------------------------------------------------

from pipeline.modules.multiome.multiome_grn import (   # noqa: E402
    multiome_grn,
    _extract_top_peaks,
    _get_tf_list,
    _regulons_correlation,
    _build_grn_table,
)


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestMultiomeGrnHappyPath:
    def test_returns_tuple(self):
        out = multiome_grn(_make_mdata(), _make_deg_dict())
        assert isinstance(out, tuple) and len(out) == 2

    def test_returns_mdata_and_dict(self):
        mdata_out, result = multiome_grn(_make_mdata(), _make_deg_dict())
        assert isinstance(mdata_out, MuData)
        assert isinstance(result, dict)

    def test_omicsage_grn_written(self):
        mdata_out, _ = multiome_grn(_make_mdata(), _make_deg_dict())
        assert "omicsage_grn" in mdata_out.uns

    def test_grn_network_written(self):
        mdata_out, _ = multiome_grn(_make_mdata(), _make_deg_dict())
        assert "grn_network" in mdata_out.uns

    def test_grn_network_is_dict(self):
        mdata_out, _ = multiome_grn(_make_mdata(), _make_deg_dict())
        assert isinstance(mdata_out.uns["grn_network"], dict)

    def test_grn_network_has_tf_key(self):
        mdata_out, _ = multiome_grn(_make_mdata(), _make_deg_dict())
        assert "tf" in mdata_out.uns["grn_network"]

    def test_n_grn_edges_nonnegative(self):
        _, result = multiome_grn(_make_mdata(), _make_deg_dict())
        assert result["n_grn_edges"] >= 0

    def test_inplace_false_leaves_original(self):
        mdata = _make_mdata()
        original_uns_keys = set(mdata.uns.keys())
        multiome_grn(mdata, _make_deg_dict(), inplace=False)
        assert "omicsage_grn" not in mdata.uns
        assert set(mdata.uns.keys()) == original_uns_keys

    def test_inplace_true_modifies_original(self):
        mdata = _make_mdata()
        multiome_grn(mdata, _make_deg_dict(), inplace=True)
        assert "omicsage_grn" in mdata.uns


# ---------------------------------------------------------------------------
# Metrics dict
# ---------------------------------------------------------------------------

class TestMultiomeGrnMetrics:
    def setup_method(self):
        self.mdata_out, self.result = multiome_grn(
            _make_mdata(), _make_deg_dict()
        )

    def test_provenance_key_present(self):
        assert "provenance" in self.result

    def test_n_tfs_rna_key_present(self):
        assert "n_tfs_rna" in self.result

    def test_n_tfs_atac_key_present(self):
        assert "n_tfs_atac" in self.result

    def test_n_grn_edges_key_present(self):
        assert "n_grn_edges" in self.result

    def test_grn_df_key_present(self):
        assert "grn_df" in self.result

    def test_grn_df_is_dataframe(self):
        assert isinstance(self.result["grn_df"], pd.DataFrame)

    def test_n_grn_edges_matches_grn_df(self):
        assert self.result["n_grn_edges"] == len(self.result["grn_df"])

    def test_n_tfs_rna_nonneg(self):
        assert self.result["n_tfs_rna"] >= 0

    def test_n_tfs_atac_nonneg(self):
        assert self.result["n_tfs_atac"] >= 0


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestMultiomeGrnProvenance:
    def setup_method(self):
        _, self.result = multiome_grn(_make_mdata(), _make_deg_dict())
        self.prov = self.result["provenance"]

    def test_module_name(self):
        assert self.prov["module"] == "multiome_grn"

    def test_timestamp_present(self):
        assert "timestamp" in self.prov
        assert isinstance(self.prov["timestamp"], str)
        assert len(self.prov["timestamp"]) > 0

    def test_params_recorded(self):
        assert "params" in self.prov

    def test_params_motif_db(self):
        assert "motif_db" in self.prov["params"]

    def test_params_groupby(self):
        assert "groupby" in self.prov["params"]

    def test_params_n_top_peaks(self):
        assert "n_top_peaks" in self.prov["params"]

    def test_params_min_cells(self):
        assert "min_cells" in self.prov["params"]

    def test_outputs_recorded(self):
        assert "outputs" in self.prov

    def test_outputs_n_tfs_rna(self):
        assert "n_tfs_rna" in self.prov["outputs"]

    def test_outputs_n_tfs_atac(self):
        assert "n_tfs_atac" in self.prov["outputs"]

    def test_outputs_n_grn_edges(self):
        assert "n_grn_edges" in self.prov["outputs"]


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

class TestMultiomeGrnValidation:
    def test_missing_rna_raises(self):
        mdata = _make_mdata(include_rna=False)
        with pytest.raises(KeyError, match="rna"):
            multiome_grn(mdata, _make_deg_dict())

    def test_missing_atac_raises(self):
        mdata = _make_mdata(include_atac=False)
        with pytest.raises(KeyError, match="atac"):
            multiome_grn(mdata, _make_deg_dict())

    def test_missing_groupby_raises(self):
        mdata = _make_mdata()
        with pytest.raises(KeyError, match="bad_groupby"):
            multiome_grn(mdata, _make_deg_dict(), groupby="bad_groupby")

    def test_too_few_cells_raises(self):
        mdata = _make_mdata(n_cells=5)
        with pytest.raises(ValueError, match="min_cells"):
            multiome_grn(mdata, _make_deg_dict(), min_cells=50)


# ---------------------------------------------------------------------------
# GRN table structure
# ---------------------------------------------------------------------------

class TestGrnTableStructure:
    def setup_method(self):
        _, self.result = multiome_grn(_make_mdata(), _make_deg_dict())
        self.df = self.result["grn_df"]

    def test_has_tf_column(self):
        if not self.df.empty:
            assert "tf" in self.df.columns

    def test_has_target_gene_column(self):
        if not self.df.empty:
            assert "target_gene" in self.df.columns

    def test_has_rna_score_column(self):
        if not self.df.empty:
            assert "rna_score" in self.df.columns

    def test_has_atac_score_column(self):
        if not self.df.empty:
            assert "atac_score" in self.df.columns

    def test_has_combined_score_column(self):
        if not self.df.empty:
            assert "combined_score" in self.df.columns

    def test_has_cell_type_column(self):
        if not self.df.empty:
            assert "cell_type" in self.df.columns

    def test_combined_score_nonneg(self):
        if not self.df.empty:
            assert (self.df["combined_score"] >= 0).all()

    def test_cell_types_are_valid(self):
        if not self.df.empty:
            mdata = _make_mdata()
            valid_cts = set(mdata["rna"].obs[GROUPBY].astype(str).unique())
            assert set(self.df["cell_type"].unique()).issubset(valid_cts)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

class TestInternalHelpers:
    def test_extract_top_peaks_returns_list(self):
        atac = _make_atac()
        peaks = _extract_top_peaks(atac, "rank_genes_groups_dca", n_top_peaks=10)
        assert isinstance(peaks, list)

    def test_extract_top_peaks_length_bounded(self):
        atac = _make_atac()
        peaks = _extract_top_peaks(atac, "rank_genes_groups_dca", n_top_peaks=5)
        # Should return at most n_top_peaks * n_groups unique peaks
        assert len(peaks) <= 5 * 2 + 5  # generous upper bound

    def test_extract_top_peaks_missing_key(self):
        atac = _make_atac()
        atac.uns.pop("rank_genes_groups_dca", None)
        peaks = _extract_top_peaks(atac, "rank_genes_groups_dca", n_top_peaks=10)
        assert peaks == []

    def test_get_tf_list_returns_list(self):
        gene_names = [f"gene_{i}" for i in range(N_GENES)]
        result = _get_tf_list(gene_names)
        assert isinstance(result, list)

    def test_get_tf_list_subset_of_genes(self):
        gene_names = ["TBX21", "GATA3", "RORC", "gene_1", "gene_2"]
        result = _get_tf_list(gene_names)
        assert all(tf in gene_names for tf in result)

    def test_regulons_correlation_returns_list(self):
        rng = np.random.default_rng(0)
        genes = ["GATA3", "gene_1", "gene_2", "gene_3", "gene_4"]
        expr = pd.DataFrame(
            rng.random((30, 5)).astype(np.float32),
            columns=genes,
        )
        result = _regulons_correlation(expr, tf_list=["GATA3"], n_targets=3)
        assert isinstance(result, list)

    def test_regulons_correlation_dict_keys(self):
        rng = np.random.default_rng(0)
        genes = ["GATA3", "gene_1", "gene_2", "gene_3", "gene_4"]
        expr = pd.DataFrame(
            rng.random((30, 5)).astype(np.float32),
            columns=genes,
        )
        result = _regulons_correlation(expr, tf_list=["GATA3"], n_targets=3)
        if result:
            assert "tf" in result[0]
            assert "targets" in result[0]
            assert "weights" in result[0]

    def test_regulons_correlation_no_nonzero_tfs(self):
        rng = np.random.default_rng(0)
        genes = ["gene_1", "gene_2"]
        expr = pd.DataFrame(rng.random((10, 2)).astype(np.float32), columns=genes)
        result = _regulons_correlation(expr, tf_list=["MISSING_TF"], n_targets=3)
        assert result == []

    def test_build_grn_table_returns_dataframe(self):
        rna  = _make_rna()
        atac = _make_atac()
        df = _build_grn_table(
            rna_adata=rna,
            atac_adata=atac,
            groupby=GROUPBY,
            aucell_scores=None,
            rna_regulon_names=[],
            motif_tfs=[],
            n_top_peaks=10,
        )
        assert isinstance(df, pd.DataFrame)

    def test_build_grn_table_empty_when_no_tfs(self):
        rna  = _make_rna()
        atac = _make_atac()
        df = _build_grn_table(
            rna_adata=rna,
            atac_adata=atac,
            groupby=GROUPBY,
            aucell_scores=None,
            rna_regulon_names=[],
            motif_tfs=[],
            n_top_peaks=10,
        )
        assert len(df) == 0
