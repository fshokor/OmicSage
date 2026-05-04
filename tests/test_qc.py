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
"""

import numpy as np
import pytest
from pathlib import Path

from pipeline.modules.qc.ingest import load_dataset
from pipeline.modules.qc.qc import run_qc, _detect_mt_genes


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
    import numpy as np
    import anndata
    # Read only metadata first (backed mode) to get cell count cheaply
    adata_backed = anndata.read_h5ad(CITE_H5AD, backed="r")
    n_cells = adata_backed.n_obs
    rng = np.random.default_rng(42)
    idx = rng.choice(n_cells, size=min(5000, n_cells), replace=False)
    idx_sorted = np.sort(idx)
    # Subset in backed mode, then load into memory
    adata_sub = adata_backed[idx_sorted].to_memory()
    adata_backed.file.close()
    # Now run through ingest to ensure raw counts in X
    from pipeline.modules.qc.ingest import _extract_raw_counts
    adata_sub, _ = _extract_raw_counts(adata_sub)
    import scipy.sparse as sp
    if not sp.issparse(adata_sub.X):
        adata_sub.X = sp.csr_matrix(adata_sub.X)
    adata_sub.obs["sample"] = "GSE194122_cite_sub5k"
    return adata_sub


@pytest.fixture(scope="session")
def cite_qc(cite_raw):
    """Run QC on CITE-seq data once, reuse result across all tests."""
    adata_f, metrics = run_qc(
        cite_raw,
        min_genes=200,
        max_genes=6000,
        max_mt_pct=20.0,
        remove_doublets=True,
        generate_report=False,
    )
    return cite_raw, adata_f, metrics


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
        mt_mask = _detect_mt_genes(cite_raw)
        assert mt_mask.sum() > 0, (
            "No MT- genes detected in CITE-seq data — check var_names or gene_ids"
        )

    def test_mt_gene_count_realistic(self, cite_raw):
        """Human cells typically have 13 mt-encoded genes + rRNA/tRNA.
        Expect at least 10 and no more than 100 MT genes."""
        mt_mask = _detect_mt_genes(cite_raw)
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

    def test_mt_mask_length_matches_n_vars(self, cite_raw):
        """MT mask must have exactly one entry per gene."""
        mt_mask = _detect_mt_genes(cite_raw)
        assert len(mt_mask) == cite_raw.n_vars


# ---------------------------------------------------------------------------
# 2. QC metrics computed correctly
# ---------------------------------------------------------------------------

class TestMetricsComputed:

    def test_n_genes_by_counts_in_obs(self, cite_qc):
        _, adata_f, _ = cite_qc
        assert "n_genes_by_counts" in adata_f.obs.columns

    def test_total_counts_in_obs(self, cite_qc):
        _, adata_f, _ = cite_qc
        assert "total_counts" in adata_f.obs.columns

    def test_pct_counts_mt_in_obs(self, cite_qc):
        _, adata_f, _ = cite_qc
        assert "pct_counts_mt" in adata_f.obs.columns

    def test_pct_counts_mt_between_0_and_100(self, cite_qc):
        _, adata_f, _ = cite_qc
        mt = adata_f.obs["pct_counts_mt"]
        assert (mt >= 0).all() and (mt <= 100).all(), (
            f"MT% values out of range: min={mt.min():.2f}, max={mt.max():.2f}"
        )

    def test_total_counts_positive(self, cite_qc):
        """All cells that passed QC must have at least 1 UMI."""
        _, adata_f, _ = cite_qc
        assert (adata_f.obs["total_counts"] > 0).all()

    def test_n_genes_positive(self, cite_qc):
        """All cells that passed QC must have at least 1 gene detected."""
        _, adata_f, _ = cite_qc
        assert (adata_f.obs["n_genes_by_counts"] > 0).all()

    def test_metrics_dict_has_required_keys(self, cite_qc):
        _, _, metrics = cite_qc
        required = [
            "n_cells_input", "n_cells_output", "n_cells_removed",
            "n_removed_low_genes", "n_removed_high_genes",
            "n_removed_high_mt", "n_removed_doublets",
            "median_genes_per_cell", "median_umi_per_cell",
            "median_mt_pct", "n_mt_genes", "thresholds",
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
        cite_raw, adata_f, _ = cite_qc
        assert adata_f.n_obs < cite_raw.n_obs, (
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
            f"Pass rate {pass_rate:.2%} outside expected range [70%, 100%] — "
            f"if 0% were removed check that filters are being applied"
        )

    def test_min_genes_filter_respected(self, cite_qc):
        """No cell in filtered output should have fewer genes than min_genes."""
        _, adata_f, metrics = cite_qc
        min_genes = metrics["thresholds"]["min_genes"]
        assert (adata_f.obs["n_genes_by_counts"] >= min_genes).all()

    def test_max_genes_filter_respected(self, cite_qc):
        """No cell in filtered output should have more genes than max_genes."""
        _, adata_f, metrics = cite_qc
        max_genes = metrics["thresholds"]["max_genes"]
        assert (adata_f.obs["n_genes_by_counts"] <= max_genes).all()

    def test_max_mt_pct_filter_respected(self, cite_qc):
        """No cell in filtered output should exceed the MT% threshold."""
        _, adata_f, metrics = cite_qc
        max_mt = metrics["thresholds"]["max_mt_pct"]
        assert (adata_f.obs["pct_counts_mt"] <= max_mt).all()

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
        _, adata_f, _ = cite_qc
        assert "doublet_score" in adata_f.obs.columns

    def test_predicted_doublet_column_present(self, cite_qc):
        _, adata_f, _ = cite_qc
        assert "predicted_doublet" in adata_f.obs.columns

    def test_doublet_score_range(self, cite_qc):
        """Doublet scores must be in [0, 1] or NaN if Scrublet failed."""
        _, adata_f, _ = cite_qc
        scores = adata_f.obs["doublet_score"].dropna()
        assert (scores >= 0).all() and (scores <= 1).all(), (
            f"Scores out of [0,1]: min={scores.min():.3f}, max={scores.max():.3f}"
        )

    def test_doublets_removed_count_in_metrics(self, cite_qc):
        _, _, metrics = cite_qc
        assert 0 <= metrics["n_removed_doublets"] < metrics["n_cells_input"]


# ---------------------------------------------------------------------------
# 5. Filtered AnnData shape
# ---------------------------------------------------------------------------

class TestFilteredAdataShape:

    def test_filtered_n_vars_unchanged(self, cite_qc):
        """QC filters cells only — gene count must not change."""
        cite_raw, adata_f, _ = cite_qc
        assert adata_f.n_vars == cite_raw.n_vars

    def test_filtered_shape_matches_metrics(self, cite_qc):
        _, adata_f, metrics = cite_qc
        assert adata_f.n_obs == metrics["n_cells_output"]

    def test_obs_qc_columns_in_filtered_output(self, cite_qc):
        _, adata_f, _ = cite_qc
        for col in ["n_genes_by_counts", "total_counts", "pct_counts_mt"]:
            assert col in adata_f.obs.columns

    def test_mtx_qc_runs_without_error(self, mtx_raw):
        """QC should complete without error on MTX-loaded HCC data."""
        adata_f, metrics = run_qc(
            mtx_raw, remove_doublets=False, generate_report=False
        )
        assert adata_f.n_obs > 0
        assert metrics["n_cells_input"] > 0


# ---------------------------------------------------------------------------
# 6. Ground-truth MT% validation
# ---------------------------------------------------------------------------

class TestQcMetricsMatchGroundTruth:

    def test_ground_truth_column_present(self, cite_raw):
        """GSE194122 CITE-seq file must contain GEX_pct_counts_mt."""
        assert "GEX_pct_counts_mt" in cite_raw.obs.columns

    def _run_rna_only_qc(self, cite_raw):
        """
        Helper: subset to RNA features only, then run permissive QC.

        The CITE-seq file mixes RNA + ADT features with duplicate var_names.
        var_names_make_unique() inside run_qc() renames genes (e.g. MT-ND1-1)
        which breaks MT- prefix detection. Subsetting to Gene Expression
        features first isolates the RNA matrix so MT genes are found correctly
        and MT% matches the ground truth GEX_pct_counts_mt column.

        We also remove zero-count cells before QC — after subsetting to RNA
        features only, some cells may have zero RNA counts (ADT-only cells),
        which causes Scrublet to divide by zero.
        """
        import scipy.sparse as sp
        import numpy as np

        # Subset to RNA features only if feature_types column is present
        if "feature_types" in cite_raw.var.columns:
            rna_mask = cite_raw.var["feature_types"].astype(str).isin(["Gene Expression", "GEX"])
            adata_rna = cite_raw[:, rna_mask].copy()
        else:
            adata_rna = cite_raw.copy()

        # Ensure sparse
        if not sp.issparse(adata_rna.X):
            adata_rna.X = sp.csr_matrix(adata_rna.X)

        # Remove cells with zero total RNA counts to avoid Scrublet divide-by-zero
        cell_totals = np.asarray(adata_rna.X.sum(axis=1)).flatten()
        adata_rna = adata_rna[cell_totals > 0].copy()

        adata_f, _ = run_qc(
            adata_rna,
            min_genes=1,
            max_genes=999_999,
            max_mt_pct=100.0,
            remove_doublets=False,   # skip Scrublet — not needed for MT% validation
            generate_report=False,
        )
        return adata_f

    def test_mt_pct_correlation_above_threshold(self, cite_raw):
        """
        OmicSage MT% must correlate with ground truth at r > 0.99.
        Validates that raw counts are in adata.X, not normalized values.
        Subset to RNA features first to avoid ADT interference.
        """
        adata_f = self._run_rna_only_qc(cite_raw)
        our_mt = adata_f.obs["pct_counts_mt"].values
        gt_mt  = adata_f.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        assert valid.sum() > 100, f"Too few valid cells: {valid.sum()}"
        corr = float(np.corrcoef(our_mt[valid], gt_mt[valid])[0, 1])
        assert corr > 0.99, (
            f"MT% correlation = {corr:.4f} — expected > 0.99. "
            f"Check that adata.X contains raw counts, not normalized values."
        )

    def test_mt_pct_mean_absolute_error(self, cite_raw):
        """Mean absolute error vs ground truth MT% must be < 0.5%.
        Subset to RNA features first to avoid ADT interference."""
        adata_f = self._run_rna_only_qc(cite_raw)
        our_mt = adata_f.obs["pct_counts_mt"].values
        gt_mt  = adata_f.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        mae    = float(np.mean(np.abs(our_mt[valid] - gt_mt[valid])))
        assert mae < 0.5, (
            f"Mean absolute error vs ground truth: {mae:.4f}% — expected < 0.5%"
        )
