"""
tests/test_atac_annotate.py — Tests for pipeline/modules/scripts/multiome/atac_annotate.py

Coverage:
  Happy path
  - Returns (AnnData, dict) tuple
  - obsm["gene_activity"] always written
  - uns["gene_activity_var_names"] always written
  - obs["atac_celltype"] always written
  - uns["omicsage_atac_annotate"] always written

  Gene activity matrix (Step A)
  - Shape: (n_cells, n_genes) — n_genes >= 1 when peaks parseable
  - dtype is float32
  - gene_activity_var_names length matches obsm["gene_activity"].shape[1]
  - Values are non-negative (sum of counts)
  - With provided peak_annotation: uses annotation source
  - With uns["atac"]["peak_annotation"]: uses uns source
  - With no annotation: falls back to coordinate grouping + emits UserWarning
  - With un-parseable var_names: emits warning, returns empty matrix (shape n_cells×0)
  - min_peaks_per_gene filters genes correctly
  - peak_annotation_source recorded correctly in metrics and provenance

  Label transfer (Step B)
  - With rna=None: atac_celltype == "Unknown" for all cells + UserWarning
  - With rna provided: labels transferred from obs["cell_type_vote"]
  - Majority vote: all cells in a cluster share the same atac_celltype
  - Cells with no RNA barcode match receive "Unknown"
  - Missing leiden_key raises KeyError
  - Missing rna_label_key raises KeyError
  - n_rna_barcodes_matched > 0 when barcodes overlap
  - n_rna_barcodes_matched == 0 when no barcodes overlap

  Metrics dict
  - Required keys present: n_cells, n_peaks, n_genes_activity,
    atac_celltype_counts, n_rna_barcodes_matched, leiden_key,
    rna_label_key, promoter_upstream_bp, min_peaks_per_gene,
    peak_annotation_source, rna_provided
  - n_cells correct
  - n_peaks correct
  - n_genes_activity matches gene_activity shape
  - atac_celltype_counts is dict[str, int]
  - atac_celltype_counts values sum to n_cells
  - rna_provided True/False correct
  - leiden_key and rna_label_key recorded correctly

  Provenance
  - uns["omicsage_atac_annotate"] present
  - module == "atac_annotate"
  - timestamp non-empty string
  - params: promoter_upstream_bp, min_peaks_per_gene, leiden_key,
            rna_label_key, peak_annotation_source, rna_provided
  - outputs: gene_activity_key, gene_activity_var_names_key,
             n_genes_activity, celltype_key, n_rna_barcodes_matched

  Namespace guards
  - obs["cell_type_vote"] never written (RNA collision guard)
  - obs["cell_type"] never written
  - obs["celltype"] never written

  inplace behaviour
  - inplace=False: copy returned, original obsm untouched
  - inplace=True: same object returned, obsm modified in place

  Edge cases
  - n_cells=200 with structured Leiden clusters
  - Large promoter window (200000 bp) still runs without error
  - min_peaks_per_gene=999 produces empty gene activity matrix
  - Empty peak_annotation DataFrame handled gracefully
"""

from __future__ import annotations

import sys
import os
import warnings

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from anndata import AnnData

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from pipeline.modules.scripts.multiome.atac_annotate import annotate_atac


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_N_CELLS  = 120
_N_PEAKS  = 80
_N_GENES  = 15
_N_LEIDEN = 3   # clusters

# Synthetic peak var_names in 'chr:start-end' format
_PEAK_CHROMS = ["chr1", "chr2", "chrX"]


def _make_peak_names(n: int) -> list[str]:
    """Generate n synthetic peak var_names: chr1:0-500, chr1:600-1100, ..."""
    names = []
    chrom_cycle = _PEAK_CHROMS * (n // len(_PEAK_CHROMS) + 1)
    for i in range(n):
        chrom = chrom_cycle[i]
        start = i * 1000
        end   = start + 500
        names.append(f"{chrom}:{start}-{end}")
    return names


def _make_peak_annotation(peak_names: list[str], gene_names: list[str]) -> pd.DataFrame:
    """
    Create a minimal peak_annotation DataFrame mapping peaks to genes.
    Each gene gets peaks_per_gene = n_peaks // n_genes peaks assigned.
    Uses peak_type="promoter" and distance=0.
    """
    records = []
    n_peaks = len(peak_names)
    n_genes = len(gene_names)
    peaks_per_gene = max(1, n_peaks // n_genes)

    for gi, gene in enumerate(gene_names):
        start_i = gi * peaks_per_gene
        end_i   = min(start_i + peaks_per_gene, n_peaks)
        for pi in range(start_i, end_i):
            records.append({
                "peak":      peak_names[pi],
                "gene_name": gene,
                "distance":  0,
                "peak_type": "promoter",
            })

    return pd.DataFrame(records)


def _make_atac_adata(
    n_cells: int = _N_CELLS,
    n_peaks: int = _N_PEAKS,
    n_leiden: int = _N_LEIDEN,
    with_peak_annotation: bool = False,
    parseable_peaks: bool = True,
    seed: int = 42,
) -> AnnData:
    """
    Create a minimal ATAC AnnData resembling atac_reduce.py output.

    Writes:
      layers["counts"]
      obsm["X_umap_atac"]
      obs["atac_leiden"]            — assigned in round-robin to n_leiden clusters
      uns["atac"]["peak_annotation"] — if with_peak_annotation=True
    """
    rng = np.random.default_rng(seed)

    counts = rng.integers(0, 20, size=(n_cells, n_peaks)).astype(np.float32)
    adata  = ad.AnnData(X=counts.copy())
    adata.obs_names = [f"cell_{i}" for i in range(n_cells)]

    if parseable_peaks:
        adata.var_names = _make_peak_names(n_peaks)
    else:
        adata.var_names = [f"peak_{i}" for i in range(n_peaks)]  # not chr:start-end

    adata.layers["counts"] = counts.copy()
    adata.obsm["X_umap_atac"] = rng.standard_normal((n_cells, 2)).astype(np.float32)

    # Round-robin Leiden assignment
    leiden = np.array([str(i % n_leiden) for i in range(n_cells)])
    adata.obs["atac_leiden"] = leiden

    if with_peak_annotation:
        gene_names = [f"Gene{i}" for i in range(_N_GENES)]
        ann_df = _make_peak_annotation(list(adata.var_names), gene_names)
        adata.uns["atac"] = {"peak_annotation": ann_df}

    return adata


def _make_rna_adata(
    n_cells: int = _N_CELLS,
    n_genes: int = 200,
    n_leiden: int = _N_LEIDEN,
    label_key: str = "cell_type_vote",
    seed: int = 99,
) -> AnnData:
    """
    Minimal RNA AnnData with obs[label_key] and matching barcodes.
    Barcodes match the atac fixture: cell_0 … cell_{n_cells-1}.
    """
    rng  = np.random.default_rng(seed)
    X    = rng.integers(0, 100, size=(n_cells, n_genes)).astype(np.float32)
    rna  = ad.AnnData(X=X)
    rna.obs_names = [f"cell_{i}" for i in range(n_cells)]
    rna.var_names = [f"gene_{i}" for i in range(n_genes)]

    # Assign cell type labels matching Leiden clusters
    ct_map = {str(ci): f"CellType_{ci}" for ci in range(n_leiden)}
    leiden = np.array([str(i % n_leiden) for i in range(n_cells)])
    rna.obs[label_key] = [ct_map[c] for c in leiden]

    return rna


# ---------------------------------------------------------------------------
# Happy path
# ---------------------------------------------------------------------------

class TestAnnotateAtacHappyPath:
    def test_returns_tuple(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        result = annotate_atac(atac)
        assert isinstance(result, tuple) and len(result) == 2

    def test_returns_anndata_and_dict(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac)
        assert isinstance(adata, AnnData)
        assert isinstance(metrics, dict)

    def test_gene_activity_always_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "gene_activity" in adata.obsm

    def test_gene_activity_var_names_always_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "gene_activity_var_names" in adata.uns

    def test_atac_celltype_always_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None)
        assert "atac_celltype" in adata.obs.columns

    def test_provenance_always_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "omicsage_atac_annotate" in adata.uns


# ---------------------------------------------------------------------------
# Gene activity matrix (Step A)
# ---------------------------------------------------------------------------

class TestGeneActivityMatrix:
    def test_shape_n_cells(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.obsm["gene_activity"].shape[0] == atac.n_obs

    def test_shape_n_genes_positive(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.obsm["gene_activity"].shape[1] > 0

    def test_dtype_float32(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.obsm["gene_activity"].dtype == np.float32

    def test_values_non_negative(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert (adata.obsm["gene_activity"] >= 0).all()

    def test_var_names_length_matches_shape(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        mat   = adata.obsm["gene_activity"]
        names = adata.uns["gene_activity_var_names"]
        assert len(names) == mat.shape[1]

    def test_provided_peak_annotation_used(self):
        """Passing peak_annotation= directly should use 'provided' source."""
        atac  = _make_atac_adata(with_peak_annotation=False)
        genes = [f"Gene{i}" for i in range(5)]
        ann   = _make_peak_annotation(list(atac.var_names), genes)
        adata, metrics = annotate_atac(atac, peak_annotation=ann)
        assert metrics["peak_annotation_source"] == "provided"
        assert adata.obsm["gene_activity"].shape[1] > 0

    def test_uns_peak_annotation_used(self):
        """Peak annotation from uns["atac"]["peak_annotation"] should be picked up."""
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac)
        assert metrics["peak_annotation_source"] == "uns"

    def test_no_annotation_fallback_emits_warning(self):
        """When no annotation available, should warn about fallback."""
        atac = _make_atac_adata(with_peak_annotation=False, parseable_peaks=True)
        with pytest.warns(UserWarning, match="peak_annotation not found"):
            annotate_atac(atac)

    def test_no_annotation_fallback_source_recorded(self):
        atac = _make_atac_adata(with_peak_annotation=False, parseable_peaks=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, metrics = annotate_atac(atac)
        assert metrics["peak_annotation_source"] == "coordinate_fallback"

    def test_unparseable_peaks_emits_warning(self):
        """Non-chr:start-end var_names should warn about empty matrix."""
        atac = _make_atac_adata(with_peak_annotation=False, parseable_peaks=False)
        with pytest.warns(UserWarning):
            annotate_atac(atac)

    def test_unparseable_peaks_empty_matrix(self):
        """gene_activity should have 0 genes when var_names are not parseable."""
        atac = _make_atac_adata(with_peak_annotation=False, parseable_peaks=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, metrics = annotate_atac(atac)
        assert adata.obsm["gene_activity"].shape[1] == 0
        assert metrics["n_genes_activity"] == 0

    def test_min_peaks_per_gene_filters(self):
        """min_peaks_per_gene=999 should produce empty matrix."""
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac, min_peaks_per_gene=999)
        assert adata.obsm["gene_activity"].shape[1] == 0
        assert metrics["n_genes_activity"] == 0

    def test_provided_overrides_uns(self):
        """Explicitly provided peak_annotation should override uns."""
        atac  = _make_atac_adata(with_peak_annotation=True)  # uns has annotation
        genes = ["OnlyGene"]
        ann   = _make_peak_annotation(list(atac.var_names)[:5], genes)
        adata, metrics = annotate_atac(atac, peak_annotation=ann)
        assert metrics["peak_annotation_source"] == "provided"
        assert "OnlyGene" in adata.uns["gene_activity_var_names"]

    def test_empty_peak_annotation_df_falls_back(self):
        """Empty DataFrame should fall through to coordinate fallback."""
        atac = _make_atac_adata(with_peak_annotation=False, parseable_peaks=True)
        empty_df = pd.DataFrame(columns=["peak", "gene_name", "distance", "peak_type"])
        with pytest.warns(UserWarning):
            annotate_atac(atac, peak_annotation=empty_df)

    def test_n_genes_activity_in_metrics(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac)
        assert metrics["n_genes_activity"] == adata.obsm["gene_activity"].shape[1]


# ---------------------------------------------------------------------------
# Label transfer (Step B)
# ---------------------------------------------------------------------------

class TestLabelTransfer:
    def test_without_rna_all_unknown(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None)
        assert (adata.obs["atac_celltype"] == "Unknown").all()

    def test_without_rna_emits_warning(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with pytest.warns(UserWarning, match="No RNA AnnData"):
            annotate_atac(atac, rna=None)

    def test_with_rna_labels_transferred(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        adata, _ = annotate_atac(atac, rna=rna)
        unique_labels = set(adata.obs["atac_celltype"].unique())
        # Should have real cell type labels, not all Unknown
        assert unique_labels != {"Unknown"}

    def test_majority_vote_consistent_per_cluster(self):
        """All cells in an ATAC Leiden cluster must share the same atac_celltype."""
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        adata, _ = annotate_atac(atac, rna=rna)
        leiden   = adata.obs["atac_leiden"].astype(str)
        celltype = adata.obs["atac_celltype"].astype(str)
        for cid in leiden.unique():
            mask   = (leiden == cid).values
            labels = celltype[mask].unique()
            assert len(labels) == 1, (
                f"Cluster {cid} has multiple atac_celltype values: {labels}"
            )

    def test_n_rna_barcodes_matched_positive(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        _, metrics = annotate_atac(atac, rna=rna)
        assert metrics["n_rna_barcodes_matched"] == _N_CELLS  # all match

    def test_n_rna_barcodes_matched_zero_no_overlap(self):
        """RNA with completely different barcodes → 0 matched."""
        atac = _make_atac_adata(with_peak_annotation=True)
        # RNA with non-matching barcodes
        rna  = _make_rna_adata()
        rna.obs_names = [f"rna_cell_{i}" for i in range(_N_CELLS)]
        _, metrics = annotate_atac(atac, rna=rna)
        assert metrics["n_rna_barcodes_matched"] == 0

    def test_no_overlap_all_unknown(self):
        """No barcode overlap → all cells get 'Unknown'."""
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        rna.obs_names = [f"rna_cell_{i}" for i in range(_N_CELLS)]
        adata, _ = annotate_atac(atac, rna=rna)
        assert (adata.obs["atac_celltype"] == "Unknown").all()

    def test_missing_leiden_key_raises(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        del atac.obs["atac_leiden"]
        rna  = _make_rna_adata()
        with pytest.raises(KeyError, match="atac_leiden"):
            annotate_atac(atac, rna=rna)

    def test_missing_rna_label_key_raises(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        del rna.obs["cell_type_vote"]
        with pytest.raises(KeyError, match="cell_type_vote"):
            annotate_atac(atac, rna=rna)

    def test_custom_leiden_key(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        atac.obs["my_leiden"] = atac.obs["atac_leiden"].copy()
        del atac.obs["atac_leiden"]
        rna  = _make_rna_adata()
        adata, metrics = annotate_atac(atac, rna=rna, leiden_key="my_leiden")
        assert metrics["leiden_key"] == "my_leiden"
        assert "atac_celltype" in adata.obs.columns

    def test_custom_rna_label_key(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        rna.obs["my_label"] = rna.obs["cell_type_vote"].copy()
        del rna.obs["cell_type_vote"]
        adata, metrics = annotate_atac(atac, rna=rna, rna_label_key="my_label")
        assert metrics["rna_label_key"] == "my_label"

    def test_rna_provided_flag_true(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        _, metrics = annotate_atac(atac, rna=rna)
        assert metrics["rna_provided"] is True

    def test_rna_provided_flag_false(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _, metrics = annotate_atac(atac, rna=None)
        assert metrics["rna_provided"] is False


# ---------------------------------------------------------------------------
# Metrics dict
# ---------------------------------------------------------------------------

class TestMetricsDict:
    _REQUIRED_KEYS = [
        "n_cells", "n_peaks", "n_genes_activity", "atac_celltype_counts",
        "n_rna_barcodes_matched", "leiden_key", "rna_label_key",
        "promoter_upstream_bp", "min_peaks_per_gene",
        "peak_annotation_source", "rna_provided",
    ]

    def test_required_keys_present(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        for key in self._REQUIRED_KEYS:
            assert key in metrics, f"Missing metrics key: {key}"

    def test_n_cells_correct(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert metrics["n_cells"] == _N_CELLS

    def test_n_peaks_correct(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert metrics["n_peaks"] == _N_PEAKS

    def test_n_genes_activity_matches_obsm(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac)
        assert metrics["n_genes_activity"] == adata.obsm["gene_activity"].shape[1]

    def test_celltype_counts_is_dict(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert isinstance(metrics["atac_celltype_counts"], dict)

    def test_celltype_counts_values_are_int(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        for v in metrics["atac_celltype_counts"].values():
            assert isinstance(v, int)

    def test_celltype_counts_sum_to_n_cells(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert sum(metrics["atac_celltype_counts"].values()) == metrics["n_cells"]

    def test_leiden_key_recorded(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert metrics["leiden_key"] == "atac_leiden"

    def test_rna_label_key_recorded(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert metrics["rna_label_key"] == "cell_type_vote"

    def test_promoter_upstream_bp_default(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac)
        assert metrics["promoter_upstream_bp"] == 2000

    def test_promoter_upstream_bp_custom(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac, promoter_upstream_bp=5000)
        assert metrics["promoter_upstream_bp"] == 5000

    def test_min_peaks_per_gene_recorded(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        _, metrics = annotate_atac(atac, min_peaks_per_gene=2)
        assert metrics["min_peaks_per_gene"] == 2


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------

class TestProvenance:
    def test_provenance_present(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "omicsage_atac_annotate" in adata.uns

    def test_module_name(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.uns["omicsage_atac_annotate"]["module"] == "atac_annotate"

    def test_timestamp_non_empty(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.uns["omicsage_atac_annotate"]["timestamp"]

    def test_params_keys_present(self):
        expected = {
            "promoter_upstream_bp", "min_peaks_per_gene", "leiden_key",
            "rna_label_key", "peak_annotation_source", "rna_provided",
        }
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        params = adata.uns["omicsage_atac_annotate"]["params"]
        for key in expected:
            assert key in params, f"Missing provenance params key: {key}"

    def test_outputs_keys_present(self):
        expected = {
            "gene_activity_key", "gene_activity_var_names_key",
            "n_genes_activity", "celltype_key", "n_rna_barcodes_matched",
        }
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        outputs = adata.uns["omicsage_atac_annotate"]["outputs"]
        for key in expected:
            assert key in outputs, f"Missing provenance outputs key: {key}"

    def test_outputs_gene_activity_key(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.uns["omicsage_atac_annotate"]["outputs"]["gene_activity_key"] == "gene_activity"

    def test_outputs_celltype_key(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.uns["omicsage_atac_annotate"]["outputs"]["celltype_key"] == "atac_celltype"

    def test_outputs_n_genes_matches_obsm(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        prov = adata.uns["omicsage_atac_annotate"]["outputs"]["n_genes_activity"]
        assert prov == adata.obsm["gene_activity"].shape[1]

    def test_params_rna_provided_false_without_rna(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None)
        assert adata.uns["omicsage_atac_annotate"]["params"]["rna_provided"] is False

    def test_params_rna_provided_true_with_rna(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        adata, _ = annotate_atac(atac, rna=rna)
        assert adata.uns["omicsage_atac_annotate"]["params"]["rna_provided"] is True

    def test_peak_annotation_source_in_params(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert adata.uns["omicsage_atac_annotate"]["params"]["peak_annotation_source"] == "uns"


# ---------------------------------------------------------------------------
# Namespace guards
# ---------------------------------------------------------------------------

class TestNamespaceGuards:
    def test_no_cell_type_vote_written(self):
        """obs['cell_type_vote'] must never be written — RNA namespace."""
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        adata, _ = annotate_atac(atac, rna=rna)
        assert "cell_type_vote" not in adata.obs.columns

    def test_no_cell_type_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "cell_type" not in adata.obs.columns

    def test_no_celltype_written(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, _ = annotate_atac(atac)
        assert "celltype" not in adata.obs.columns


# ---------------------------------------------------------------------------
# inplace behaviour
# ---------------------------------------------------------------------------

class TestInplace:
    def test_inplace_false_original_untouched(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        original_obsm_keys = set(atac.obsm.keys())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            _ = annotate_atac(atac, rna=None, inplace=False)
        assert set(atac.obsm.keys()) == original_obsm_keys

    def test_inplace_false_returns_copy(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None, inplace=False)
        assert adata is not atac

    def test_inplace_true_modifies_original(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None, inplace=True)
        assert "gene_activity" in atac.obsm

    def test_inplace_true_returns_same_object(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata, _ = annotate_atac(atac, rna=None, inplace=True)
        assert adata is atac


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    def test_n_cells_200(self):
        atac = _make_atac_adata(n_cells=200, with_peak_annotation=True)
        rna  = _make_rna_adata(n_cells=200)
        adata, metrics = annotate_atac(atac, rna=rna)
        assert metrics["n_cells"] == 200

    def test_large_promoter_window(self):
        atac = _make_atac_adata(with_peak_annotation=True)
        adata, metrics = annotate_atac(atac, promoter_upstream_bp=200_000)
        assert metrics["promoter_upstream_bp"] == 200_000

    def test_idempotent(self):
        """Running twice should overwrite cleanly."""
        atac = _make_atac_adata(with_peak_annotation=True)
        rna  = _make_rna_adata()
        adata1, _ = annotate_atac(atac, rna=rna, inplace=False)
        adata2, _ = annotate_atac(adata1, rna=rna, inplace=False)
        assert "gene_activity" in adata2.obsm
        assert "atac_celltype" in adata2.obs.columns

    def test_partial_barcode_overlap(self):
        """Half of RNA barcodes match; other half → 'Unknown' per cluster."""
        atac = _make_atac_adata(n_cells=60, n_leiden=2, with_peak_annotation=True)
        rna  = _make_rna_adata(n_cells=30)  # only cell_0 … cell_29
        adata, metrics = annotate_atac(atac, rna=rna)
        assert 0 < metrics["n_rna_barcodes_matched"] < 60

    def test_single_leiden_cluster(self):
        """All cells in one cluster — should still work."""
        atac = _make_atac_adata(n_leiden=1, with_peak_annotation=True)
        atac.obs["atac_leiden"] = "0"
        rna  = _make_rna_adata()
        adata, metrics = annotate_atac(atac, rna=rna)
        assert set(adata.obs["atac_celltype"].unique()).issubset(
            {"CellType_0", "Unknown"}
        )
