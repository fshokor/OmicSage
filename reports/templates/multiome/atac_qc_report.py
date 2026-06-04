"""
OmicSage — ATAC QC Report
reports/templates/multiome/atac_qc_report.py

Generated after atac_qc step.
Output: multiome_01_atac_qc_report.html

Usage
-----
    from reports.templates.multiome.atac_qc_report import run_atac_qc_report
    run_atac_qc_report(
        atac=atac_qcd,
        metrics=metrics,
        report_path="reports/GSE194122_multiome/multiome_01_atac_qc_report.html",
        dataset_name="BMMC Multiome (NeurIPS 2021)",
    )
"""

from __future__ import annotations

import base64
import io
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData

_DPI = 130


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _get_dense(x) -> np.ndarray:
    if hasattr(x, "toarray"):
        return x.toarray()
    return np.asarray(x, dtype=float)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _plot_n_peaks_histogram(atac: AnnData, threshold_low: int, threshold_high: int) -> str:
    """Histogram of n_peaks_by_counts per cell with QC thresholds."""
    vals = atac.obs["n_peaks_by_counts"].values
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=60, color="#4C78A8", alpha=0.8, edgecolor="none")
    ax.axvline(threshold_low,  color="#e74c3c", linewidth=1.5,
               linestyle="--", label=f"min = {threshold_low:,}")
    ax.axvline(threshold_high, color="#e74c3c", linewidth=1.5,
               linestyle=":",  label=f"max = {threshold_high:,}")
    ax.set_xlabel("Peaks per cell (n_peaks_by_counts)", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Distribution of Peaks per Cell", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_total_counts_histogram(atac: AnnData, threshold_low: int, threshold_high: int) -> str:
    """Histogram of total_peak_counts per cell with QC thresholds."""
    col = "total_peak_counts"
    if col not in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, f"'{col}' not found in obs", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    vals = atac.obs[col].values
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(np.log10(vals + 1), bins=60, color="#e07b3a", alpha=0.8, edgecolor="none")
    ax.axvline(np.log10(threshold_low  + 1), color="#e74c3c", linewidth=1.5,
               linestyle="--", label=f"min = {threshold_low:,}")
    ax.axvline(np.log10(threshold_high + 1), color="#e74c3c", linewidth=1.5,
               linestyle=":",  label=f"max = {threshold_high:,}")
    ax.set_xlabel("log10(total peak counts + 1)", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Distribution of Total Peak Counts per Cell", fontsize=12, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_nucleosome_signal(atac: AnnData, threshold: float) -> str:
    """Histogram of nucleosome signal per cell."""
    col = "nucleosome_signal"
    if col not in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "nucleosome_signal not available\n(no fragment file)",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    vals = atac.obs[col].values
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=60, color="#9b59b6", alpha=0.8, edgecolor="none")
    ax.axvline(threshold, color="#e74c3c", linewidth=1.5,
               linestyle="--", label=f"threshold = {threshold}")
    ax.set_xlabel("Nucleosome signal", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title(
        "Nucleosome Signal Distribution\n"
        "(< 2 = good signal-to-noise; CellRanger-ARC pre-computed)",
        fontsize=11, fontweight="bold",
    )
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_reads_in_peaks(atac: AnnData) -> str:
    """Histogram of fraction of reads in peaks per cell."""
    col = "reads_in_peaks_frac"
    if col not in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, f"'{col}' not found in obs", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    vals = atac.obs[col].values
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(vals, bins=60, color="#27ae60", alpha=0.8, edgecolor="none")
    ax.axvline(0.15, color="#e74c3c", linewidth=1.2, linestyle="--",
               label="typical low-quality threshold (0.15)")
    ax.set_xlabel("Fraction of reads in peaks", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Reads in Peaks Fraction\n(low values = poor signal-to-noise)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_counts_vs_peaks_scatter(atac: AnnData, metrics: dict) -> str:
    """Scatter: total_peak_counts vs n_peaks_by_counts, coloured by QC pass/fail."""
    if "atac_qc_pass" not in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "atac_qc_pass not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    qc_pass = atac.obs["atac_qc_pass"].values.astype(bool)
    counts  = atac.obs["total_peak_counts"].values
    peaks   = atac.obs["n_peaks_by_counts"].values

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(counts[~qc_pass], peaks[~qc_pass], s=2, alpha=0.4,
               color="#e74c3c", label=f"QC fail (n={int((~qc_pass).sum()):,})",
               rasterized=True)
    ax.scatter(counts[qc_pass],  peaks[qc_pass],  s=2, alpha=0.3,
               color="#2980b9", label=f"QC pass (n={int(qc_pass.sum()):,})",
               rasterized=True)
    ax.set_xscale("log")
    ax.set_xlabel("Total peak counts (log scale)", fontsize=10)
    ax.set_ylabel("n_peaks_by_counts", fontsize=10)
    ax.set_title("Total Counts vs Peaks per Cell\n(red = QC fail)",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, markerscale=4)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_doublet_scores(atac: AnnData) -> str:
    """Histogram of Scrublet doublet scores."""
    col = "atac_doublet_score"
    if col not in atac.obs.columns or atac.obs[col].isna().all():
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "Doublet scores not available", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    vals = atac.obs[col].dropna().values
    n_flagged = int(atac.obs.get("atac_predicted_doublet", 0).sum()) \
        if "atac_predicted_doublet" in atac.obs.columns else 0

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.hist(vals, bins=50, color="#f39c12", alpha=0.8, edgecolor="none")
    ax.set_xlabel("Scrublet doublet score", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title(
        f"ATAC Doublet Score Distribution\n"
        f"({n_flagged:,} cells flagged — flag only, not removed)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML sections
# ---------------------------------------------------------------------------

def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — ATAC QC Report — {dataset_name}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }}
    header {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}
    section {{ background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }}
    section h2 {{ font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }}
    section h3 {{ font-size: 0.98rem; font-weight: 600; color: #16213e; margin: 18px 0 8px; }}
    section p  {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    code {{ font-family: "SFMono-Regular", Consolas, monospace;
            background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em; }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{ background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }}
    .stat-card.warn {{ background: #fff3cd; }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-card.warn .stat-value {{ color: #856404; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }}
    th {{ background: #f0f2ff; color: #0f3460; font-weight: 600;
          padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0; }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}
    .fig-grid {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }}
    .fig-wrap {{ flex: 1 1 300px; max-width: 520px; }}
    .fig-wrap.wide {{ flex: 1 1 100%; max-width: 100%; }}
    .fig-wrap h3 {{ font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }}
    .fig-wrap img {{ width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }}
    .callout {{ background: #e8f4fd; border-left: 4px solid #3498db;
                padding: 12px 16px; border-radius: 0 6px 6px 0;
                margin: 12px 0; font-size: 0.88rem; color: #1a3a5c; }}
    .callout.warn {{ background: #fff8e1; border-left-color: #f39c12; color: #5d3a00; }}
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage — ATAC QC Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    thr = metrics.get("thresholds", {})
    n_before  = metrics.get("n_cells_before", 0)
    n_after   = metrics.get("n_cells_after",  0)
    n_pass    = metrics.get("n_qc_pass",      0)
    n_fail    = metrics.get("n_qc_fail",      0)
    pass_rate = 100 * n_pass / n_before if n_before > 0 else 0.0

    cards = [
        ("Cells (input)",      f"{n_before:,}",                       False),
        ("Cells (QC pass)",    f"{n_pass:,}",                         False),
        ("Cells (QC fail)",    f"{n_fail:,}",                         n_fail > 0),
        ("Pass rate",          f"{pass_rate:.1f}%",                   False),
        ("Peaks (before)",     f"{metrics.get('n_peaks_before',0):,}",False),
        ("Peaks (after)",      f"{metrics.get('n_peaks_after',0):,}", False),
    ]
    stat_html = "".join(
        f'<div class="stat-card{"  warn" if is_warn else ""}">'
        f'<div class="stat-value">{val}</div>'
        f'<div class="stat-label">{lbl}</div></div>'
        for lbl, val, is_warn in cards
    )

    cr_note = (
        "CellRanger-ARC pre-computed metrics used for nucleosome signal and fragment counts."
        if metrics.get("cellranger_metrics_available")
        else "No CellRanger-ARC obs columns found — only peak-count-based metrics computed."
    )

    rows = [
        ("Cells before QC",           f"{n_before:,}"),
        ("Cells passing QC",          f"{n_pass:,}"),
        ("Cells failing QC",          f"{n_fail:,}"),
        ("filter_cells applied",      str(metrics.get("filter_cells_applied", False))),
        ("Peaks before filter",       f"{metrics.get('n_peaks_before',0):,}"),
        ("Peaks after filter",        f"{metrics.get('n_peaks_after',0):,}"),
        ("Peaks removed (< min_cells)",f"{metrics.get('n_peaks_removed_feature_filter',0):,}"),
        ("Median n_peaks_by_counts",  f"{metrics.get('median_n_peaks_by_counts',0):.0f}"),
        ("Median total_peak_counts",  f"{metrics.get('median_total_peak_counts',0):.0f}"),
        ("Nucleosome filter applied", str(metrics.get("nucleosome_filter_applied", False))),
        ("Fragment file available",   str(metrics.get("fragment_file_available", False))),
        ("min_peaks threshold",       f"{thr.get('min_peaks','?'):,}"),
        ("max_peaks threshold",       f"{thr.get('max_peaks','?'):,}"),
        ("min_peak_counts threshold", f"{thr.get('min_peak_counts','?'):,}"),
        ("max_peak_counts threshold", f"{thr.get('max_peak_counts','?'):,}"),
        ("max_nucleosome_signal",     str(thr.get("max_nucleosome_signal", "?"))),
        ("min_cells per peak",        str(thr.get("min_cells_per_peak", "?"))),
        ("Fail — low peaks",          str(metrics.get("n_fail_low_peaks",  0))),
        ("Fail — high peaks",         str(metrics.get("n_fail_high_peaks", 0))),
        ("Fail — low counts",         str(metrics.get("n_fail_low_counts", 0))),
        ("Fail — high counts",        str(metrics.get("n_fail_high_counts",0))),
        ("Fail — high nucleosome",    str(metrics.get("n_fail_high_nucleosome", 0))),
    ]
    rows_html = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="callout warn">
        <strong>filter_cells=False (default).</strong> {cr_note}
        Cell filtering was not applied because RNA QC already removed low-quality barcodes.
        The <code>atac_qc_pass</code> column in <code>.obs</code> records which cells
        would be removed if stricter ATAC-only filtering were desired.
      </div>
      <div class="stat-grid">{stat_html}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows_html}</tbody></table>
    </section>"""


def _section_figures(
    fig_n_peaks: str,
    fig_counts: str,
    fig_nucleosome: str,
    fig_reads_in_peaks: str,
    fig_scatter: str,
    fig_doublets: str,
) -> str:
    return f"""
    <section>
      <h2>QC Metric Distributions</h2>
      <div class="callout">
        Thresholds shown as dashed red lines. These are calibrated for the
        NeurIPS 2021 BMMC multiome dataset following the sc-best-practices
        ATAC QC chapter (Lance &amp; Martens 2022).
      </div>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Peaks per Cell</h3>
          <img src="data:image/png;base64,{fig_n_peaks}" alt="n peaks histogram">
        </div>
        <div class="fig-wrap">
          <h3>Total Peak Counts per Cell (log10)</h3>
          <img src="data:image/png;base64,{fig_counts}" alt="total counts histogram">
        </div>
        <div class="fig-wrap">
          <h3>Nucleosome Signal</h3>
          <img src="data:image/png;base64,{fig_nucleosome}" alt="nucleosome signal">
        </div>
        <div class="fig-wrap">
          <h3>Reads in Peaks Fraction</h3>
          <img src="data:image/png;base64,{fig_reads_in_peaks}" alt="reads in peaks">
        </div>
        <div class="fig-wrap wide">
          <h3>Total Counts vs Peaks per Cell (QC pass/fail)</h3>
          <img src="data:image/png;base64,{fig_scatter}" alt="counts vs peaks scatter">
        </div>
        <div class="fig-wrap">
          <h3>Scrublet Doublet Scores (ATAC)</h3>
          <img src="data:image/png;base64,{fig_doublets}" alt="doublet scores">
        </div>
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_atac_qc_report(
    atac: AnnData,
    metrics: dict,
    report_path: str = "reports/multiome_01_atac_qc_report.html",
    dataset_name: str = "dataset",
) -> str:
    """Generate a self-contained HTML ATAC QC report.

    Parameters
    ----------
    atac : AnnData
        Output of ``atac_qc()``, with QC metrics in ``.obs``.
    metrics : dict
        Metrics dict returned by ``atac_qc()``.
    report_path : str
        Output path for the HTML file.
    dataset_name : str
        Dataset label shown in the report header.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    thr = metrics.get("thresholds", {})

    print(f"Building ATAC QC report for '{dataset_name}' ...", flush=True)

    print("  Rendering n_peaks histogram ...", flush=True)
    fig_n_peaks = _plot_n_peaks_histogram(
        atac,
        threshold_low=thr.get("min_peaks", 750),
        threshold_high=thr.get("max_peaks", 500_000),
    )
    print("  Rendering total counts histogram ...", flush=True)
    fig_counts = _plot_total_counts_histogram(
        atac,
        threshold_low=thr.get("min_peak_counts", 1_500),
        threshold_high=thr.get("max_peak_counts", 100_000),
    )
    print("  Rendering nucleosome signal histogram ...", flush=True)
    fig_nucleosome = _plot_nucleosome_signal(
        atac,
        threshold=thr.get("max_nucleosome_signal", 2.0),
    )
    print("  Rendering reads-in-peaks histogram ...", flush=True)
    fig_reads = _plot_reads_in_peaks(atac)
    print("  Rendering counts vs peaks scatter ...", flush=True)
    fig_scatter = _plot_counts_vs_peaks_scatter(atac, metrics)
    print("  Rendering doublet score histogram ...", flush=True)
    fig_doublets = _plot_doublet_scores(atac)

    sections = [
        _section_summary(metrics, dataset_name, timestamp),
        _section_figures(
            fig_n_peaks, fig_counts, fig_nucleosome,
            fig_reads, fig_scatter, fig_doublets,
        ),
    ]

    html = _render_page(sections=sections, timestamp=timestamp, dataset_name=dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
