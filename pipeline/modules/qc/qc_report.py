"""
OmicSage — QC Report Generator
pipeline/modules/qc/qc_report.py

Generates a self-contained HTML report after QC is complete.

Usage (called automatically by run_qc if generate_report=True):
    from pipeline.modules.qc.qc_report import generate_qc_report

    generate_qc_report(
        adata_raw=adata_before_qc,
        adata_filtered=adata_after_qc,
        metrics=metrics,
        output_path="reports/qc_report.html",
        sample_name="GSE194122_BMMC",
    )
"""

from __future__ import annotations

import base64
import logging
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _violin_plot(values_before, values_after, label,
                 threshold_low=None, threshold_high=None, colour="#4C78A8"):
    fig, ax = plt.subplots(figsize=(6, 4))
    parts = ax.violinplot([values_before, values_after], positions=[1, 2],
                          showmedians=True, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(colour)
        pc.set_alpha(0.6)
    parts["cmedians"].set_colors("#333333")
    parts["cmedians"].set_linewidth(2)
    if threshold_low is not None:
        ax.axhline(threshold_low, color="#e07b3a", linewidth=1.2,
                   linestyle="--", label=f"min {threshold_low}")
    if threshold_high is not None:
        ax.axhline(threshold_high, color="#e05252", linewidth=1.2,
                   linestyle="--", label=f"max {threshold_high}")
    if threshold_low or threshold_high:
        ax.legend(fontsize=8, frameon=False)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Before QC", "After QC"])
    ax.set_ylabel(label, fontsize=10)
    ax.set_title(label, fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _scatter_umi_genes(adata_raw: AnnData, mt_threshold: float) -> str:
    fig, ax = plt.subplots(figsize=(6, 5))
    x = adata_raw.obs["total_counts"].values
    y = adata_raw.obs["n_genes_by_counts"].values
    c = adata_raw.obs["pct_counts_mt"].values
    sc = ax.scatter(x, y, c=c, cmap="RdYlGn_r", s=1.5, alpha=0.5,
                    vmin=0, vmax=mt_threshold * 1.5, rasterized=True)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("MT%", fontsize=9)
    ax.set_xlabel("Total UMI counts", fontsize=10)
    ax.set_ylabel("Genes detected", fontsize=10)
    ax.set_title("UMI vs Genes (colour = MT%)", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _doublet_histogram(scores: np.ndarray, threshold: float = 0.25) -> str:
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scores, bins=60, color="#d2a8ff", alpha=0.8, edgecolor="white", linewidth=0.4)
    ax.axvline(threshold, color="#e05252", linewidth=1.5, linestyle="--",
               label=f"threshold ≈ {threshold:.2f}")
    ax.set_xlabel("Doublet score", fontsize=10)
    ax.set_ylabel("Cell count", fontsize=10)
    ax.set_title("Scrublet Doublet Score Distribution", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _mt_correlation_plot(our_mt: np.ndarray, gt_mt: np.ndarray, corr: float) -> str:
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(gt_mt, our_mt, s=1.5, alpha=0.4, color="#54a868", rasterized=True)
    lims = [min(gt_mt.min(), our_mt.min()), max(gt_mt.max(), our_mt.max())]
    ax.plot(lims, lims, "--", color="#e07b3a", linewidth=1, label="y = x")
    ax.set_xlabel("Ground truth MT% (GEX_pct_counts_mt)", fontsize=10)
    ax.set_ylabel("OmicSage MT%", fontsize=10)
    ax.set_title(f"MT% Validation  (r = {corr:.4f})", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML renderer  — matches deg_report._render_page exactly
# ---------------------------------------------------------------------------

def _render_page(title: str, sections: list[str], timestamp: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      font-size: 14px;
      line-height: 1.6;
      color: #1a1a2e;
      background: #f7f8fc;
    }}

    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white;
      padding: 32px 40px 24px;
    }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}

    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}

    section {{
      background: white;
      border-radius: 10px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.07);
      padding: 28px 32px;
      margin-bottom: 24px;
    }}
    section h2 {{
      font-size: 1.15rem;
      font-weight: 700;
      color: #0f3460;
      border-bottom: 2px solid #e8eaf6;
      padding-bottom: 10px;
      margin-bottom: 18px;
    }}
    section h3 {{
      font-size: 1rem;
      font-weight: 600;
      color: #16213e;
      margin: 18px 0 10px;
    }}
    section p {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}

    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}

    .note {{
      font-size: 0.82rem;
      color: #7a5c00;
      background: #fffbe6;
      border-left: 3px solid #f0c040;
      padding: 8px 12px;
      border-radius: 4px;
      margin-bottom: 14px;
    }}

    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 0.85em;
    }}

    .stat-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      margin-bottom: 24px;
    }}
    .stat-card {{
      background: #f0f2ff;
      border-radius: 8px;
      padding: 14px 20px;
      min-width: 130px;
      text-align: center;
      flex: 1 1 130px;
    }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}

    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.88rem;
      margin-top: 8px;
    }}
    th {{
      background: #f0f2ff;
      color: #0f3460;
      font-weight: 600;
      padding: 9px 12px;
      text-align: left;
      border-bottom: 2px solid #d0d4f0;
    }}
    td {{
      padding: 8px 12px;
      border-bottom: 1px solid #eee;
      vertical-align: middle;
    }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}

    .pass {{ color: #155724; font-weight: 600; }}
    .warn {{ color: #856404; font-weight: 600; }}

    .fig-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 18px;
      margin-top: 12px;
    }}
    .fig-wrap {{
      flex: 1 1 300px;
      max-width: 520px;
    }}
    .fig-wrap.wide {{
      flex: 1 1 100%;
      max-width: 100%;
    }}
    .fig-wrap h3 {{
      font-size: 0.9rem;
      margin-bottom: 6px;
      color: #16213e;
    }}
    .fig-wrap img {{
      width: 100%;
      border-radius: 6px;
      border: 1px solid #e8eaf6;
    }}

    footer {{
      text-align: center;
      font-size: 0.78rem;
      color: #aaa;
      padding: 24px 0 32px;
    }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage — Quality Control Report</h1>
    <p>Generated {timestamp}</p>
  </header>
  <main>
    {body}
  </main>
  <footer>
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    &middot; MIT License
  </footer>
</body>
</html>
"""


def _section_summary(metrics: dict, sample_name: str, timestamp: str) -> str:
    t = metrics["thresholds"]
    pct_kept = 100 * metrics["n_cells_output"] / metrics["n_cells_input"]

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells input",      f"{metrics['n_cells_input']:,}"),
            ("Cells kept",       f"{metrics['n_cells_output']:,}"),
            ("Cells removed",    f"{metrics['n_cells_removed']:,}"),
            ("Pass rate",        f"{pct_kept:.1f}%"),
            ("MT genes found",   f"{metrics['n_mt_genes']:,}"),
            ("Doublets removed", f"{metrics['n_removed_doublets']:,}"),
        ]
    )

    filter_rows = "".join(
        f"<tr><td>{param}</td><td>{thresh}</td><td>{removed:,}</td></tr>"
        for param, thresh, removed in [
            ("min_genes",       t["min_genes"],         metrics["n_removed_low_genes"]),
            ("max_genes",       t["max_genes"],         metrics["n_removed_high_genes"]),
            ("max_mt_pct",      f"{t['max_mt_pct']}%",  metrics["n_removed_high_mt"]),
            ("remove_doublets", str(t["remove_doublets"]), metrics["n_removed_doublets"]),
        ]
    )

    median_rows = (
        f"<tr><td>Median genes / cell (pre-QC)</td><td>{metrics['median_genes_per_cell']:.0f}</td></tr>"
        f"<tr><td>Median UMI / cell (pre-QC)</td><td>{metrics['median_umi_per_cell']:.0f}</td></tr>"
        f"<tr><td>Median MT% (pre-QC)</td><td>{metrics['median_mt_pct']:.2f}%</td></tr>"
    )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Sample: <strong>{sample_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>Pre-QC Distribution Medians</h3>
      <table>
        <thead><tr><th>Metric</th><th>Value</th></tr></thead>
        <tbody>{median_rows}</tbody>
      </table>
      <h3>Filter Thresholds Applied</h3>
      <table>
        <thead><tr><th>Parameter</th><th>Threshold</th><th>Cells removed</th></tr></thead>
        <tbody>{filter_rows}</tbody>
      </table>
    </section>
    """


def _section_distributions(plots: dict) -> str:
    figs_html = "".join(
        f'<div class="fig-wrap"><h3>{title}</h3>'
        f'<img src="data:image/png;base64,{plots[key]}" alt="{title}"></div>'
        for key, title in [
            ("violin_genes", "Genes per Cell"),
            ("violin_umi",   "Total UMI Counts"),
            ("violin_mt",    "Mitochondrial %"),
        ]
    )
    return f"""
    <section>
      <h2>Distributions — Before vs After QC</h2>
      <div class="fig-grid">{figs_html}</div>
    </section>
    """


def _section_scatter_doublet(plots: dict) -> str:
    return f"""
    <section>
      <h2>Cell-Level Diagnostics</h2>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>UMI vs Genes (colour = MT%)</h3>
          <img src="data:image/png;base64,{plots['scatter']}" alt="UMI vs genes scatter">
        </div>
        <div class="fig-wrap">
          <h3>Scrublet Doublet Score Distribution</h3>
          <img src="data:image/png;base64,{plots['doublet_hist']}" alt="Doublet histogram">
        </div>
      </div>
    </section>
    """


def _section_validation(plots: dict, corr_value: float) -> str:
    corr_class  = "pass" if corr_value >= 0.99 else "warn"
    corr_symbol = "PASS" if corr_value >= 0.99 else "CHECK"
    target_note = "(target >= 0.99)" if corr_value >= 0.99 else "(target >= 0.99 -- investigate)"
    return f"""
    <section>
      <h2>Ground-Truth MT% Validation</h2>
      <p>Pearson r between OmicSage <code>pct_counts_mt</code> and
         <code>obs['GEX_pct_counts_mt']</code>:
         <span class="{corr_class}">{corr_symbol} r = {corr_value:.4f}</span>
         {target_note}
      </p>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>MT% Correlation Plot</h3>
          <img src="data:image/png;base64,{plots['mt_corr']}" alt="MT% correlation">
        </div>
      </div>
    </section>
    """


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_qc_report(
    adata_raw: AnnData,
    adata_filtered: AnnData,
    metrics: dict[str, Any],
    output_path: str | Path = "reports/qc_report.html",
    sample_name: str = "sample",
) -> Path:
    """
    Generate a self-contained HTML QC report.

    Parameters
    ----------
    adata_raw:
        AnnData before QC filtering, after QC metrics were computed.
        Must have obs: n_genes_by_counts, total_counts, pct_counts_mt, doublet_score.
    adata_filtered:
        AnnData after all QC filters were applied.
    metrics:
        Dict returned by run_qc().
    output_path:
        Where to write the HTML file.
    sample_name:
        Label shown in the report header.

    Returns
    -------
    Path
        Resolved path to the written HTML file.
    """
    timestamp   = datetime.now().strftime("%Y-%m-%d %H:%M")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    t = metrics["thresholds"]

    # Build plots
    logger.info("Generating QC report plots ...")
    plots = {}

    plots["violin_genes"] = _violin_plot(
        adata_raw.obs["n_genes_by_counts"].values,
        adata_filtered.obs["n_genes_by_counts"].values,
        "Genes per cell",
        threshold_low=t["min_genes"], threshold_high=t["max_genes"],
        colour="#4C78A8",
    )
    plots["violin_umi"] = _violin_plot(
        adata_raw.obs["total_counts"].values,
        adata_filtered.obs["total_counts"].values,
        "Total UMI counts",
        colour="#54a868",
    )
    plots["violin_mt"] = _violin_plot(
        adata_raw.obs["pct_counts_mt"].values,
        adata_filtered.obs["pct_counts_mt"].values,
        "Mitochondrial %",
        threshold_high=t["max_mt_pct"],
        colour="#e05252",
    )
    plots["scatter"] = _scatter_umi_genes(adata_raw, t["max_mt_pct"])

    if "doublet_score" in adata_raw.obs.columns:
        scores = adata_raw.obs["doublet_score"].values
        threshold = float(np.nanpercentile(scores, 97)) if not np.all(np.isnan(scores)) else 0.25
        plots["doublet_hist"] = _doublet_histogram(scores, threshold)
    else:
        plots["doublet_hist"] = _doublet_histogram(np.zeros(10), 0.25)

    # Ground-truth validation
    corr_value = None
    if "GEX_pct_counts_mt" in adata_raw.obs.columns:
        our_mt = adata_raw.obs["pct_counts_mt"].values
        gt_mt  = adata_raw.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        if valid.sum() > 10:
            corr_value = float(np.corrcoef(our_mt[valid], gt_mt[valid])[0, 1])
            plots["mt_corr"] = _mt_correlation_plot(our_mt[valid], gt_mt[valid], corr_value)

    # Assemble sections
    sections = [
        _section_summary(metrics, sample_name, timestamp),
        _section_distributions(plots),
        _section_scatter_doublet(plots),
    ]
    if corr_value is not None:
        sections.append(_section_validation(plots, corr_value))

    html = _render_page(
        title=f"OmicSage -- QC Report -- {sample_name}",
        sections=sections,
        timestamp=timestamp,
    )
    output_path.write_text(html, encoding="utf-8")
    size_kb = output_path.stat().st_size / 1024
    logger.info("QC report written -> %s (%.1f KB)", output_path, size_kb)
    print(f"  QC report -> {output_path}  ({size_kb:.0f} KB)")
    return output_path
