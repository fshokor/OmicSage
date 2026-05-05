"""
tests/test_qc.py
================
Test suite for pipeline/modules/qc/qc.py

All tests use real benchmark data — no synthetic fixtures.

Datasets used:
  CITE-seq  : data/benchmark/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad
  MTX dir   : data/benchmark/GSE166635/  (HCC scRNA-seq)

Both datasets are skipped automatically if the file is not present.

Run:
    python -m pytest tests/test_qc.py -v

API contract (v2):
    run_qc() returns (MuData, dict).
    RNA QC metrics live on mdata["rna"].obs.
    Multi-modal data additionally exposes mdata["adt"] or mdata["atac"].
"""

import numpy as np
import pytest
import scipy.sparse as sp
from pathlib import Path

from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc, _detect_mt_genes, _detect_modality


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CITE_H5AD = Path(
    "data/benchmark/"
    "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
)

MTX_DIR = Path("data/test/GSE166635/HCC1")


# ---------------------------------------------------------------------------
# Shared fixtures — loaded once per session, reused across all tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def cite_raw():
    """Load CITE-seq h5ad, subsample to 5000 cells to stay within RAM limits.

    The full file is 90k cells x 14k genes — too large for Scrublet on 7.6 GB RAM.
    5000 cells is enough for all statistical tests to be meaningful.
    Subsampling is seeded for reproducibility.
    """
    if not CITE_H5AD.exists():
        pytest.skip(f"CITE-seq file not found: {CITE_H5AD}")
    import anndata

    adata_backed = anndata.read_h5ad(CITE_H5AD, backed="r")
    n_cells = adata_backed.n_obs
    rng = np.random.default_rng(42)
    idx = rng.choice(n_cells, size=min(5000, n_cells), replace=False)
    idx_sorted = np.sort(idx)
    adata_sub = adata_backed[idx_sorted].to_memory()
    adata_backed.file.close()

    from pipeline.modules.qc.ingest import _extract_raw_counts
    adata_sub, _ = _extract_raw_counts(adata_sub)

    if not sp.issparse(adata_sub.X):
        adata_sub.X = sp.csr_matrix(adata_sub.X)

    adata_sub.obs["sample"] = "GSE194122_cite_sub5k"
    return adata_sub


@pytest.fixture(scope="session")
def cite_qc(cite_raw):
    """Run QC on CITE-seq data once, reuse result across all tests.

    Returns (cite_raw, mdata, metrics).
    mdata["rna"] is the filtered RNA AnnData.
    mdata["adt"] is the filtered ADT AnnData (if present).
    """
    mdata, metrics = run_qc(
        cite_raw,
        min_genes=200,
        max_genes=6000,
        max_mt_pct=20.0,
        remove_doublets=True,
        generate_report=False,
    )
    return cite_raw, mdata, metrics


@pytest.fixture(scope="session")
def mtx_raw():
    """Load HCC1 MTX sample, subsample to 3000 cells if larger."""
    if not MTX_DIR.exists():
        pytest.skip(f"MTX directory not found: {MTX_DIR}")
    adata = load_dataset(MTX_DIR, verbose=False)
    if adata.n_obs > 3000:
        import scanpy as sc
        sc.pp.subsample(adata, n_obs=3000, random_state=42)
    return adata


# ---------------------------------------------------------------------------
# 1. MT gene detection
# ---------------------------------------------------------------------------

class TestMtGeneDetection:

    def test_mt_genes_detected_in_cite(self, cite_raw):
        """CITE-seq file should have MT- genes in the RNA features."""
        # Subset to RNA features before calling _detect_mt_genes directly
        if "feature_types" in cite_raw.var.columns:
            rna_mask = cite_raw.var["feature_types"].astype(str).isin(
                ["Gene Expression", "GEX"]
            )
            adata_rna = cite_raw[:, rna_mask]
        else:
            adata_rna = cite_raw
        mt_mask = _detect_mt_genes(adata_rna)
        assert mt_mask.sum() > 0, (
            "No MT- genes detected in CITE-seq RNA features — check var_names or gene_ids"
        )

    def test_mt_gene_count_realistic(self, cite_raw):
        """Human cells typically have 13 mt-encoded genes + rRNA/tRNA.
        Expect at least 10 and no more than 100 MT genes."""
        if "feature_types" in cite_raw.var.columns:
            rna_mask = cite_raw.var["feature_types"].astype(str).isin(
                ["Gene Expression", "GEX"]
            )
            adata_rna = cite_raw[:, rna_mask]
        else:
            adata_rna = cite_raw
        mt_mask = _detect_mt_genes(adata_rna)
        n_mt = int(mt_mask.sum())
        assert 10 <= n_mt <= 100, (
            f"Detected {n_mt} MT genes — expected between 10 and 100"
        )

    def test_mt_genes_detected_in_mtx(self, mtx_raw):
        """HCC MTX data (human) should also have MT- genes."""
        mt_mask = _detect_mt_genes(mtx_raw)
        assert mt_mask.sum() > 0, (
            "No MT- genes detected in HCC MTX data"
        )

    def test_mt_mask_length_matches_rna_n_vars(self, cite_qc):
        """MT mask length must match the RNA modality n_vars, not the full mixed object."""
        _, mdata, _ = cite_qc
        adata_rna = mdata["rna"]
        mt_mask = _detect_mt_genes(adata_rna)
        assert len(mt_mask) == adata_rna.n_vars


# ---------------------------------------------------------------------------
# 2. QC metrics computed correctly
# ---------------------------------------------------------------------------

class TestMetricsComputed:

    def test_n_genes_by_counts_in_obs(self, cite_qc):
        _, mdata, _ = cite_qc
        assert "n_genes_by_counts" in mdata["rna"].obs.columns

    def test_total_counts_in_obs(self, cite_qc):
        _, mdata, _ = cite_qc
        assert "total_counts" in mdata["rna"].obs.columns

    def test_pct_counts_mt_in_obs(self, cite_qc):
        _, mdata, _ = cite_qc
        assert "pct_counts_mt" in mdata["rna"].obs.columns

    def test_pct_counts_mt_between_0_and_100(self, cite_qc):
        _, mdata, _ = cite_qc
        mt = mdata["rna"].obs["pct_counts_mt"]
        assert (mt >= 0).all() and (mt <= 100).all(), (
            f"MT% values out of range: min={mt.min():.2f}, max={mt.max():.2f}"
        )

    def test_total_counts_positive(self, cite_qc):
        """All cells that passed QC must have at least 1 UMI."""
        _, mdata, _ = cite_qc
        assert (mdata["rna"].obs["total_counts"] > 0).all()

    def test_n_genes_positive(self, cite_qc):
        """All cells that passed QC must have at least 1 gene detected."""
        _, mdata, _ = cite_qc
        assert (mdata["rna"].obs["n_genes_by_counts"] > 0).all()

    def test_metrics_dict_has_required_keys(self, cite_qc):
        _, _, metrics = cite_qc
        required = [
            "n_cells_input", "n_cells_output", "n_cells_removed",
            "n_removed_low_genes", "n_removed_high_genes",
            "n_removed_high_mt", "n_removed_doublets",
            "median_genes_per_cell", "median_umi_per_cell",
            "median_mt_pct", "median_ribo_pct", "median_hb_pct",
            "n_mt_genes", "n_ribo_genes", "n_hb_genes",
            "thresholds", "modality",
        ]
        for key in required:
            assert key in metrics, f"Missing key in metrics dict: '{key}'"

    def test_n_mt_genes_realistic(self, cite_qc):
        _, _, metrics = cite_qc
        assert 10 <= metrics["n_mt_genes"] <= 100, (
            f"n_mt_genes={metrics['n_mt_genes']} — expected 10-100"
        )

    def test_median_genes_realistic(self, cite_qc):
        """Median genes per cell for BMMC should be biologically plausible."""
        _, _, metrics = cite_qc
        assert 500 <= metrics["median_genes_per_cell"] <= 5000, (
            f"median_genes_per_cell={metrics['median_genes_per_cell']:.0f}"
        )

    def test_median_umi_realistic(self, cite_qc):
        _, _, metrics = cite_qc
        assert 500 <= metrics["median_umi_per_cell"] <= 30000, (
            f"median_umi_per_cell={metrics['median_umi_per_cell']:.0f}"
        )


# ---------------------------------------------------------------------------
# 3. Filters remove low-quality cells
# ---------------------------------------------------------------------------

class TestFilterRemovesLowQualityCells:

    def test_output_fewer_cells_than_input(self, cite_qc):
        cite_raw, mdata, _ = cite_qc
        assert mdata["rna"].n_obs < cite_raw.n_obs, (
            "No cells were removed — QC filters may not be working"
        )

    def test_cell_counts_consistent(self, cite_qc):
        """n_cells_input = n_cells_output + n_cells_removed."""
        _, _, metrics = cite_qc
        assert metrics["n_cells_removed"] == (
            metrics["n_cells_input"] - metrics["n_cells_output"]
        )

    def test_pass_rate_realistic(self, cite_qc):
        """Between 70% and 99% of cells should pass on a clean public dataset."""
        _, _, metrics = cite_qc
        pass_rate = metrics["n_cells_output"] / metrics["n_cells_input"]
        assert 0.70 <= pass_rate <= 1.0, (
            f"Pass rate {pass_rate:.2%} outside expected range [70%, 100%]"
        )

    def test_min_genes_filter_respected(self, cite_qc):
        """No cell in filtered output should have fewer genes than min_genes."""
        _, mdata, metrics = cite_qc
        min_genes = metrics["thresholds"]["min_genes"]
        assert (mdata["rna"].obs["n_genes_by_counts"] >= min_genes).all()

    def test_max_genes_filter_respected(self, cite_qc):
        """No cell in filtered output should have more genes than max_genes."""
        _, mdata, metrics = cite_qc
        max_genes = metrics["thresholds"]["max_genes"]
        assert (mdata["rna"].obs["n_genes_by_counts"] <= max_genes).all()

    def test_max_mt_pct_filter_respected(self, cite_qc):
        """No cell in filtered output should exceed the MT% threshold."""
        _, mdata, metrics = cite_qc
        max_mt = metrics["thresholds"]["max_mt_pct"]
        assert (mdata["rna"].obs["pct_counts_mt"] <= max_mt).all()

    def test_original_adata_not_mutated(self, cite_raw):
        """run_qc must not modify the caller's AnnData."""
        n_obs_before = cite_raw.n_obs
        cols_before  = set(cite_raw.obs.columns)
        run_qc(cite_raw, remove_doublets=False, generate_report=False)
        assert cite_raw.n_obs == n_obs_before, "n_obs changed — adata was mutated"
        assert set(cite_raw.obs.columns) == cols_before, (
            "obs columns changed — adata was mutated"
        )

    def test_remove_doublets_false_keeps_more_cells(self, cite_raw):
        """remove_doublets=False should keep >= cells compared to True."""
        _, metrics_with    = run_qc(cite_raw, remove_doublets=True,  generate_report=False)
        _, metrics_without = run_qc(cite_raw, remove_doublets=False, generate_report=False)
        assert metrics_without["n_cells_output"] >= metrics_with["n_cells_output"]


# ---------------------------------------------------------------------------
# 4. Doublet detection
# ---------------------------------------------------------------------------

class TestDoubletScoresAdded:

    def test_doublet_score_column_present(self, cite_qc):
        _, mdata, _ = cite_qc
        assert "doublet_score" in mdata["rna"].obs.columns

    def test_predicted_doublet_column_present(self, cite_qc):
        _, mdata, _ = cite_qc
        assert "predicted_doublet" in mdata["rna"].obs.columns

    def test_doublet_score_range(self, cite_qc):
        """Doublet scores must be in [0, 1] or NaN if Scrublet failed."""
        _, mdata, _ = cite_qc
        scores = mdata["rna"].obs["doublet_score"].dropna()
        assert (scores >= 0).all() and (scores <= 1).all(), (
            f"Scores out of [0,1]: min={scores.min():.3f}, max={scores.max():.3f}"
        )

    def test_doublets_removed_count_in_metrics(self, cite_qc):
        _, _, metrics = cite_qc
        assert 0 <= metrics["n_removed_doublets"] < metrics["n_cells_input"]


# ---------------------------------------------------------------------------
# 5. Ribosomal and hemoglobin gene metrics
# ---------------------------------------------------------------------------

class TestRiboHbMetrics:

    def test_pct_counts_ribo_in_obs(self, cite_qc):
        """pct_counts_ribo must be present in mdata['rna'].obs after QC."""
        _, mdata, _ = cite_qc
        assert "pct_counts_ribo" in mdata["rna"].obs.columns

    def test_pct_counts_hb_in_obs(self, cite_qc):
        """pct_counts_hb must be present in mdata['rna'].obs after QC."""
        _, mdata, _ = cite_qc
        assert "pct_counts_hb" in mdata["rna"].obs.columns

    def test_pct_counts_ribo_range(self, cite_qc):
        """Ribo% must be in [0, 100] for all cells."""
        _, mdata, _ = cite_qc
        ribo = mdata["rna"].obs["pct_counts_ribo"]
        assert (ribo >= 0).all() and (ribo <= 100).all(), (
            f"Ribo% out of range: min={ribo.min():.2f}, max={ribo.max():.2f}"
        )

    def test_pct_counts_hb_range(self, cite_qc):
        """HB% must be in [0, 100] for all cells."""
        _, mdata, _ = cite_qc
        hb = mdata["rna"].obs["pct_counts_hb"]
        assert (hb >= 0).all() and (hb <= 100).all(), (
            f"HB% out of range: min={hb.min():.2f}, max={hb.max():.2f}"
        )

    def test_ribo_genes_detected(self, cite_qc):
        """At least 1 ribosomal gene (RPS/RPL) must be detected in BMMC data."""
        _, _, metrics = cite_qc
        assert metrics["n_ribo_genes"] > 0, (
            "No ribosomal genes detected — check that var_names contain RPS/RPL symbols"
        )

    def test_ribo_gene_count_realistic(self, cite_qc):
        """Human cells have ~80 cytoplasmic ribosomal protein genes (RPS + RPL).
        Expect at least 50 and no more than 200."""
        _, _, metrics = cite_qc
        assert 50 <= metrics["n_ribo_genes"] <= 200, (
            f"n_ribo_genes={metrics['n_ribo_genes']} — expected 50-200"
        )

    def test_median_ribo_pct_in_metrics(self, cite_qc):
        """metrics dict must contain median_ribo_pct."""
        _, _, metrics = cite_qc
        assert "median_ribo_pct" in metrics

    def test_median_hb_pct_in_metrics(self, cite_qc):
        """metrics dict must contain median_hb_pct."""
        _, _, metrics = cite_qc
        assert "median_hb_pct" in metrics

    def test_n_ribo_genes_in_metrics(self, cite_qc):
        """metrics dict must contain n_ribo_genes."""
        _, _, metrics = cite_qc
        assert "n_ribo_genes" in metrics

    def test_n_hb_genes_in_metrics(self, cite_qc):
        """metrics dict must contain n_hb_genes."""
        _, _, metrics = cite_qc
        assert "n_hb_genes" in metrics

    def test_median_ribo_pct_realistic(self, cite_qc):
        """Ribosomal reads typically account for 10-60% in immune cells.
        Expect median ribo% > 5% for BMMC data."""
        _, _, metrics = cite_qc
        assert metrics["median_ribo_pct"] > 5.0, (
            f"median_ribo_pct={metrics['median_ribo_pct']:.2f}% — expected > 5% for BMMC"
        )

    def test_median_hb_pct_low_in_bmmc(self, cite_qc):
        """BMMC cells are not erythrocytes — median HB% should be near 0."""
        _, _, metrics = cite_qc
        assert metrics["median_hb_pct"] < 5.0, (
            f"median_hb_pct={metrics['median_hb_pct']:.2f}% — expected < 5% for BMMC "
            f"(high HB% would suggest erythrocyte contamination)"
        )

    def test_ribo_var_column_present(self, cite_qc):
        """mdata['rna'].var must contain a 'ribo' boolean column."""
        _, mdata, _ = cite_qc
        assert "ribo" in mdata["rna"].var.columns

    def test_hb_var_column_present(self, cite_qc):
        """mdata['rna'].var must contain a 'hb' boolean column."""
        _, mdata, _ = cite_qc
        assert "hb" in mdata["rna"].var.columns


# ---------------------------------------------------------------------------
# 6. Filtered MuData shape  (was 5 — renumbered after ribo/hb section added)
# ---------------------------------------------------------------------------

class TestFilteredMudataShape:

    def test_filtered_rna_n_vars_is_gex_only(self, cite_qc):
        """mdata['rna'].n_vars must equal the number of GEX features in cite_raw,
        not the total n_vars of the mixed object (RNA + ADT)."""
        cite_raw, mdata, _ = cite_qc
        if "feature_types" in cite_raw.var.columns:
            gex_mask = cite_raw.var["feature_types"].astype(str).isin(
                ["Gene Expression", "GEX"]
            )
            expected_n_vars = int(gex_mask.sum())
        else:
            expected_n_vars = cite_raw.n_vars
        assert mdata["rna"].n_vars == expected_n_vars, (
            f"mdata['rna'].n_vars={mdata['rna'].n_vars} != "
            f"expected GEX n_vars={expected_n_vars}"
        )

    def test_filtered_shape_matches_metrics(self, cite_qc):
        _, mdata, metrics = cite_qc
        assert mdata["rna"].n_obs == metrics["n_cells_output"]

    def test_obs_qc_columns_in_filtered_output(self, cite_qc):
        _, mdata, _ = cite_qc
        for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
            assert col in mdata["rna"].obs.columns

    def test_mtx_qc_runs_without_error(self, mtx_raw):
        """QC should complete without error on MTX-loaded HCC data."""
        mdata, metrics = run_qc(
            mtx_raw, remove_doublets=False, generate_report=False
        )
        assert mdata["rna"].n_obs > 0
        assert metrics["n_cells_input"] > 0

    def test_mtx_returns_mudata_with_rna_key(self, mtx_raw):
        """Plain RNA input must still return a MuData with a 'rna' key."""
        from mudata import MuData
        mdata, _ = run_qc(mtx_raw, remove_doublets=False, generate_report=False)
        assert isinstance(mdata, MuData), "run_qc did not return a MuData object"
        assert "rna" in mdata.mod, "MuData missing 'rna' key for plain RNA input"


# ---------------------------------------------------------------------------
# 6. Ground-truth MT% validation
# ---------------------------------------------------------------------------

class TestQcMetricsMatchGroundTruth:

    def test_ground_truth_column_present(self, cite_raw):
        """GSE194122 CITE-seq file must contain GEX_pct_counts_mt."""
        assert "GEX_pct_counts_mt" in cite_raw.obs.columns

    def _run_permissive_qc(self, cite_raw):
        """Run permissive QC via the new modality-aware run_qc().

        run_qc() now handles GEX subsetting internally — no manual
        subsetting required here. Filters are set wide open so that
        nearly all cells pass, allowing MT% comparison against ground truth.
        """
        mdata, _ = run_qc(
            cite_raw,
            min_genes=1,
            max_genes=999_999,
            max_mt_pct=100.0,
            remove_doublets=False,
            generate_report=False,
        )
        return mdata["rna"]

    def test_mt_pct_correlation_above_threshold(self, cite_raw):
        """OmicSage MT% must correlate with ground truth at r > 0.99.

        Validates that:
        - raw counts are in adata.X (not normalized)
        - MT% is computed on RNA features only (not inflated by ADT counts)
        """
        adata_rna = self._run_permissive_qc(cite_raw)
        our_mt = adata_rna.obs["pct_counts_mt"].values
        gt_mt  = adata_rna.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        assert valid.sum() > 100, f"Too few valid cells: {valid.sum()}"
        corr = float(np.corrcoef(our_mt[valid], gt_mt[valid])[0, 1])
        assert corr > 0.99, (
            f"MT% correlation = {corr:.4f} — expected > 0.99. "
            f"Check that adata.X contains raw counts, not normalized values."
        )

    def test_mt_pct_mean_absolute_error(self, cite_raw):
        """Mean absolute error vs ground truth MT% must be < 0.5%."""
        adata_rna = self._run_permissive_qc(cite_raw)
        our_mt = adata_rna.obs["pct_counts_mt"].values
        gt_mt  = adata_rna.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        mae    = float(np.mean(np.abs(our_mt[valid] - gt_mt[valid])))
        assert mae < 0.5, (
            f"Mean absolute error vs ground truth: {mae:.4f}% — expected < 0.5%"
        )


# ---------------------------------------------------------------------------
# 7. Multi-modal QC — new tests (target: 37/37)
# ---------------------------------------------------------------------------

class TestMultiModalQC:

    # --- _detect_modality() ---

    def test_detect_modality_rna_no_feature_types(self, mtx_raw):
        """Plain RNA AnnData with no feature_types column → 'rna'."""
        adata = mtx_raw.copy()
        if "feature_types" in adata.var.columns:
            del adata.var["feature_types"]
        result = _detect_modality(adata)
        assert result == "rna", f"Expected 'rna', got '{result}'"

    def test_detect_modality_cite(self, cite_raw):
        """CITE-seq AnnData with ADT features → 'cite'."""
        result = _detect_modality(cite_raw)
        assert result == "cite", (
            f"Expected 'cite' for CITE-seq data, got '{result}'. "
            f"feature_types present: {'feature_types' in cite_raw.var.columns}"
        )

    def test_detect_modality_multiome_synthetic(self):
        """Synthetic Multiome AnnData with 'Peaks' feature_types → 'multiome'."""
        import anndata
        import pandas as pd
        n_cells, n_rna, n_atac = 100, 200, 500
        X = sp.random(n_cells, n_rna + n_atac, density=0.1, format="csr")
        var = pd.DataFrame(
            {"feature_types": ["Gene Expression"] * n_rna + ["Peaks"] * n_atac},
            index=[f"gene_{i}" for i in range(n_rna)]
                  + [f"peak_{i}" for i in range(n_atac)],
        )
        adata = anndata.AnnData(X=X, var=var)
        result = _detect_modality(adata)
        assert result == "multiome", f"Expected 'multiome', got '{result}'"

    # --- MuData structure after CITE-seq QC ---

    def test_cite_qc_returns_mudata(self, cite_qc):
        """run_qc on CITE-seq data must return a MuData object."""
        from mudata import MuData
        _, mdata, _ = cite_qc
        assert isinstance(mdata, MuData), (
            f"Expected MuData, got {type(mdata)}"
        )

    def test_cite_qc_mudata_has_rna_and_adt_keys(self, cite_qc):
        """MuData from CITE-seq QC must have both 'rna' and 'adt' keys."""
        _, mdata, _ = cite_qc
        assert "rna" in mdata.mod, "MuData missing 'rna' key"
        assert "adt" in mdata.mod, "MuData missing 'adt' key — ADT features were lost"

    def test_cite_qc_adt_features_preserved(self, cite_raw, cite_qc):
        """ADT feature count in mdata['adt'] must match the raw ADT feature count."""
        _, mdata, _ = cite_qc
        if "feature_types" not in cite_raw.var.columns:
            pytest.skip("feature_types column not present — cannot count ADT features")
        adt_mask = cite_raw.var["feature_types"].astype(str).isin(
            ["ADT", "Antibody Capture"]
        )
        expected_n_adt = int(adt_mask.sum())
        actual_n_adt   = mdata["adt"].n_vars
        assert actual_n_adt == expected_n_adt, (
            f"ADT features after QC: {actual_n_adt} != expected {expected_n_adt}. "
            f"ADT features were dropped or incorrectly counted."
        )

    def test_cite_qc_total_counts_based_on_rna_only(self, cite_raw, cite_qc):
        """total_counts in mdata['rna'].obs must reflect RNA UMIs only.

        Computes the expected RNA-only total_counts directly from the raw
        GEX submatrix and checks that the QC-computed values match.
        """
        _, mdata, _ = cite_qc

        if "feature_types" not in cite_raw.var.columns:
            pytest.skip("feature_types not present — cannot isolate RNA counts")

        # Build expected RNA-only counts from raw data
        gex_mask = cite_raw.var["feature_types"].astype(str).isin(
            ["Gene Expression", "GEX"]
        )
        adata_rna_raw = cite_raw[:, gex_mask]
        rna_totals = np.asarray(adata_rna_raw.X.sum(axis=1)).flatten()

        # Look up the same cells in mdata["rna"]
        filtered_barcodes = mdata["rna"].obs_names
        raw_barcodes      = cite_raw.obs_names
        shared_idx        = raw_barcodes.get_indexer(filtered_barcodes)
        expected_totals   = rna_totals[shared_idx]

        actual_totals = mdata["rna"].obs["total_counts"].values

        # Allow floating-point rounding tolerance
        np.testing.assert_allclose(
            actual_totals, expected_totals, rtol=1e-5,
            err_msg=(
                "total_counts does not match RNA-only sum — "
                "ADT counts may have been included in the metric."
            ),
        )

    def test_cite_qc_all_modalities_have_same_cells(self, cite_qc):
        """After QC, mdata['rna'] and mdata['adt'] must contain the same cell barcodes."""
        _, mdata, _ = cite_qc
        if "adt" not in mdata.mod:
            pytest.skip("No ADT modality — CITE-seq detection may have failed")
        rna_barcodes = set(mdata["rna"].obs_names)
        adt_barcodes = set(mdata["adt"].obs_names)
        assert rna_barcodes == adt_barcodes, (
            f"Barcode mismatch: "
            f"{len(rna_barcodes - adt_barcodes)} barcodes in RNA but not ADT, "
            f"{len(adt_barcodes - rna_barcodes)} in ADT but not RNA."
        )
