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

The report is a single .html file — no external dependencies.
It contains:
  - Summary cards (cells in/out, genes, MT genes, doublets)
  - Violin plots: genes per cell, UMI counts, MT%
  - Scatter plot: UMI vs genes (coloured by MT%)
  - Doublet score histogram
  - Filter thresholds table
  - Ground-truth MT% correlation (if GEX_pct_counts_mt is present)
"""

from __future__ import annotations

import base64
import json
import logging
from io import BytesIO
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend — must be set before pyplot import
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from anndata import AnnData

logger = logging.getLogger(__name__)

# ── Colour palette (dark scientific theme) ──────────────────────────────────
_PALETTE = {
    "bg":        "#0d1117",
    "surface":   "#161b22",
    "border":    "#30363d",
    "accent":    "#58a6ff",
    "accent2":   "#3fb950",
    "accent3":   "#f78166",
    "accent4":   "#d2a8ff",
    "text":      "#e6edf3",
    "subtext":   "#8b949e",
    "warn":      "#e3b341",
}


# ── Figure helpers ───────────────────────────────────────────────────────────

def _fig_to_b64(fig: plt.Figure) -> str:
    """Encode a matplotlib figure as a base64 PNG string."""
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _base_fig(figsize=(7, 4)):
    fig, ax = plt.subplots(figsize=figsize, facecolor=_PALETTE["bg"])
    ax.set_facecolor(_PALETTE["surface"])
    for spine in ax.spines.values():
        spine.set_edgecolor(_PALETTE["border"])
    ax.tick_params(colors=_PALETTE["subtext"], labelsize=9)
    ax.xaxis.label.set_color(_PALETTE["subtext"])
    ax.yaxis.label.set_color(_PALETTE["subtext"])
    ax.title.set_color(_PALETTE["text"])
    return fig, ax


# ── Individual plots ─────────────────────────────────────────────────────────

def _violin_plot(values_before: np.ndarray, values_after: np.ndarray,
                 label: str, threshold_low=None, threshold_high=None,
                 colour: str = "#58a6ff") -> str:
    fig, ax = _base_fig(figsize=(6, 4))
    data = [values_before, values_after]
    parts = ax.violinplot(data, positions=[1, 2], showmedians=True,
                          showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(colour)
        pc.set_alpha(0.6)
    parts["cmedians"].set_colors(_PALETTE["text"])
    parts["cmedians"].set_linewidth(2)

    if threshold_low is not None:
        ax.axhline(threshold_low, color=_PALETTE["warn"], linewidth=1,
                   linestyle="--", label=f"min {threshold_low}")
    if threshold_high is not None:
        ax.axhline(threshold_high, color=_PALETTE["accent3"], linewidth=1,
                   linestyle="--", label=f"max {threshold_high}")
    if threshold_low or threshold_high:
        ax.legend(fontsize=8, facecolor=_PALETTE["surface"],
                  labelcolor=_PALETTE["subtext"], framealpha=0.8)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Before QC", "After QC"], color=_PALETTE["text"])
    ax.set_ylabel(label)
    ax.set_title(label)
    ax.grid(axis="y", color=_PALETTE["border"], linewidth=0.5, alpha=0.5)
    return _fig_to_b64(fig)


def _scatter_umi_genes(adata_raw: AnnData, mt_threshold: float) -> str:
    fig, ax = _base_fig(figsize=(6, 5))
    x = adata_raw.obs["total_counts"].values
    y = adata_raw.obs["n_genes_by_counts"].values
    c = adata_raw.obs["pct_counts_mt"].values

    sc = ax.scatter(x, y, c=c, cmap="RdYlGn_r", s=1.5, alpha=0.5,
                    vmin=0, vmax=mt_threshold * 1.5)
    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.ax.tick_params(colors=_PALETTE["subtext"], labelsize=8)
    cbar.set_label("MT%", color=_PALETTE["subtext"], fontsize=9)
    cbar.outline.set_edgecolor(_PALETTE["border"])

    ax.set_xlabel("Total UMI counts")
    ax.set_ylabel("Genes detected")
    ax.set_title("UMI vs Genes (colour = MT%)")
    ax.grid(color=_PALETTE["border"], linewidth=0.4, alpha=0.4)
    return _fig_to_b64(fig)


def _doublet_histogram(scores: np.ndarray, threshold: float = 0.25) -> str:
    fig, ax = _base_fig(figsize=(6, 4))
    ax.hist(scores, bins=60, color=_PALETTE["accent4"], alpha=0.7,
            edgecolor=_PALETTE["border"], linewidth=0.4)
    ax.axvline(threshold, color=_PALETTE["accent3"], linewidth=1.5,
               linestyle="--", label=f"threshold ≈ {threshold:.2f}")
    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Cell count")
    ax.set_title("Scrublet Doublet Score Distribution")
    ax.legend(fontsize=8, facecolor=_PALETTE["surface"],
              labelcolor=_PALETTE["subtext"])
    ax.grid(axis="y", color=_PALETTE["border"], linewidth=0.5, alpha=0.5)
    return _fig_to_b64(fig)


def _mt_correlation_plot(our_mt: np.ndarray, ground_truth_mt: np.ndarray,
                         corr: float) -> str:
    fig, ax = _base_fig(figsize=(5, 5))
    ax.scatter(ground_truth_mt, our_mt, s=1.5, alpha=0.4,
               color=_PALETTE["accent2"])
    lims = [min(ground_truth_mt.min(), our_mt.min()),
            max(ground_truth_mt.max(), our_mt.max())]
    ax.plot(lims, lims, "--", color=_PALETTE["warn"], linewidth=1,
            label="y = x")
    ax.set_xlabel("Ground truth MT% (GEX_pct_counts_mt)")
    ax.set_ylabel("OmicSage MT%")
    ax.set_title(f"MT% Validation  (r = {corr:.4f})")
    ax.legend(fontsize=8, facecolor=_PALETTE["surface"],
              labelcolor=_PALETTE["subtext"])
    ax.grid(color=_PALETTE["border"], linewidth=0.4, alpha=0.4)
    return _fig_to_b64(fig)


# ── HTML assembly ────────────────────────────────────────────────────────────

def _card(value: str, label: str, colour: str = "#58a6ff") -> str:
    return f"""
    <div class="card">
      <div class="card-value" style="color:{colour}">{value}</div>
      <div class="card-label">{label}</div>
    </div>"""


def _img_block(b64: str, caption: str) -> str:
    return f"""
    <figure class="plot-block">
      <img src="data:image/png;base64,{b64}" alt="{caption}">
      <figcaption>{caption}</figcaption>
    </figure>"""


CSS = """
:root {{
  --bg:      {bg};
  --surface: {surface};
  --border:  {border};
  --accent:  {accent};
  --text:    {text};
  --sub:     {subtext};
}}
* {{ box-sizing: border-box; margin: 0; padding: 0; }}
body {{
  font-family: 'JetBrains Mono', 'Fira Mono', monospace;
  background: var(--bg);
  color: var(--text);
  padding: 2rem;
  line-height: 1.6;
}}
header {{
  border-bottom: 1px solid var(--border);
  padding-bottom: 1.2rem;
  margin-bottom: 2rem;
}}
header h1 {{
  font-size: 1.6rem;
  letter-spacing: 0.05em;
  color: var(--accent);
}}
header p {{ color: var(--sub); font-size: 0.85rem; margin-top: 0.3rem; }}
h2 {{
  font-size: 1rem;
  letter-spacing: 0.08em;
  text-transform: uppercase;
  color: var(--sub);
  border-left: 3px solid var(--accent);
  padding-left: 0.75rem;
  margin: 2rem 0 1rem;
}}
.cards {{
  display: flex;
  flex-wrap: wrap;
  gap: 1rem;
  margin-bottom: 1.5rem;
}}
.card {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 1rem 1.5rem;
  min-width: 160px;
  flex: 1;
}}
.card-value {{ font-size: 1.8rem; font-weight: 700; }}
.card-label  {{ font-size: 0.75rem; color: var(--sub); margin-top: 0.2rem; }}
.plots {{
  display: flex;
  flex-wrap: wrap;
  gap: 1.5rem;
  margin-bottom: 1.5rem;
}}
.plot-block {{
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: 6px;
  padding: 1rem;
  flex: 1;
  min-width: 300px;
}}
.plot-block img {{ width: 100%; border-radius: 4px; }}
.plot-block figcaption {{
  font-size: 0.72rem;
  color: var(--sub);
  text-align: center;
  margin-top: 0.5rem;
}}
table {{
  width: 100%;
  border-collapse: collapse;
  font-size: 0.82rem;
  background: var(--surface);
  border-radius: 6px;
  overflow: hidden;
  border: 1px solid var(--border);
}}
th {{
  background: var(--border);
  color: var(--sub);
  text-align: left;
  padding: 0.6rem 1rem;
  font-weight: 600;
  letter-spacing: 0.05em;
}}
td {{ padding: 0.5rem 1rem; border-top: 1px solid var(--border); }}
tr:hover td {{ background: rgba(88,166,255,0.04); }}
.pass {{ color: {accent2}; }}
.warn {{ color: {warn}; }}
footer {{
  margin-top: 3rem;
  padding-top: 1rem;
  border-top: 1px solid var(--border);
  font-size: 0.72rem;
  color: var(--sub);
}}
""".format(**_PALETTE)


def _build_html(
    sample_name: str,
    metrics: dict[str, Any],
    plots: dict[str, str],
    corr_value: float | None,
    timestamp: str,
) -> str:
    t = metrics["thresholds"]

    # Summary cards
    pct_kept = 100 * metrics["n_cells_output"] / metrics["n_cells_input"]
    cards_html = (
        _card(f"{metrics['n_cells_input']:,}",  "Cells input",    _PALETTE["subtext"]) +
        _card(f"{metrics['n_cells_output']:,}",  "Cells kept",     _PALETTE["accent2"]) +
        _card(f"{metrics['n_cells_removed']:,}", "Cells removed",  _PALETTE["accent3"]) +
        _card(f"{pct_kept:.1f}%",                "Pass rate",      _PALETTE["accent"]) +
        _card(f"{metrics['n_mt_genes']:,}",       "MT genes found", _PALETTE["accent4"]) +
        _card(f"{metrics['n_removed_doublets']:,}","Doublets removed",_PALETTE["warn"])
    )

    # Distribution medians
    median_rows = f"""
    <tr><td>Median genes / cell (pre-QC)</td><td>{metrics['median_genes_per_cell']:.0f}</td></tr>
    <tr><td>Median UMI / cell (pre-QC)</td><td>{metrics['median_umi_per_cell']:.0f}</td></tr>
    <tr><td>Median MT% (pre-QC)</td><td>{metrics['median_mt_pct']:.2f}%</td></tr>
    """

    # Filter table
    def _row(param, value, removed):
        return f"<tr><td>{param}</td><td>{value}</td><td>{removed:,} cells</td></tr>"

    filter_rows = (
        _row("min_genes",       t["min_genes"],     metrics["n_removed_low_genes"]) +
        _row("max_genes",       t["max_genes"],     metrics["n_removed_high_genes"]) +
        _row("max_mt_pct",      f"{t['max_mt_pct']}%", metrics["n_removed_high_mt"]) +
        _row("remove_doublets", str(t["remove_doublets"]), metrics["n_removed_doublets"])
    )

    # Validation section
    if corr_value is not None:
        corr_class = "pass" if corr_value >= 0.99 else "warn"
        corr_symbol = "✓" if corr_value >= 0.99 else "⚠"
        validation_section = f"""
        <h2>Ground-Truth Validation</h2>
        <p style="font-size:0.85rem; color:var(--sub); margin-bottom:1rem;">
          Pearson r between OmicSage MT% and
          <code>obs['GEX_pct_counts_mt']</code>:
          <span class="{corr_class}">{corr_symbol} r = {corr_value:.4f}</span>
          {'(target ≥ 0.99 ✓)' if corr_value >= 0.99 else '(target ≥ 0.99 — investigate)'}
        </p>
        <div class="plots">
          {_img_block(plots["mt_corr"], "MT% — OmicSage vs Ground Truth")}
        </div>
        """
    else:
        validation_section = ""

    violin_plots = (
        _img_block(plots["violin_genes"], "Genes per cell") +
        _img_block(plots["violin_umi"],   "Total UMI counts") +
        _img_block(plots["violin_mt"],    "Mitochondrial %")
    )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage QC Report — {sample_name}</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link href="https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;600;700&display=swap" rel="stylesheet">
  <style>{CSS}</style>
</head>
<body>
  <header>
    <h1>⬡ OmicSage · QC Report</h1>
    <p>Sample: <strong>{sample_name}</strong> &nbsp;|&nbsp; Generated: {timestamp}</p>
  </header>

  <h2>Summary</h2>
  <div class="cards">{cards_html}</div>

  <h2>Distributions — Before vs After QC</h2>
  <div class="plots">{violin_plots}</div>

  <h2>UMI vs Genes Scatter</h2>
  <div class="plots">
    {_img_block(plots["scatter"], "Each point = one cell; colour = MT%")}
  </div>

  <h2>Doublet Detection</h2>
  <div class="plots">
    {_img_block(plots["doublet_hist"], "Scrublet doublet score distribution")}
  </div>

  <h2>Pre-QC Distribution Medians</h2>
  <table>
    <tr><th>Metric</th><th>Value</th></tr>
    {median_rows}
  </table>

  <h2>Filter Thresholds Applied</h2>
  <table>
    <tr><th>Parameter</th><th>Threshold</th><th>Cells removed</th></tr>
    {filter_rows}
  </table>

  {validation_section}

  <footer>
    OmicSage · Phase 1 · pipeline/modules/qc/qc_report.py ·
    Report auto-generated — do not edit manually.
  </footer>
</body>
</html>"""


# ── Public API ───────────────────────────────────────────────────────────────

def generate_qc_report(
    adata_raw: AnnData,
    adata_filtered: AnnData,
    metrics: dict[str, Any],
    output_path: str | Path = "reports/qc_report.html",
    sample_name: str = "sample",
) -> Path:
    """Generate a self-contained HTML QC report.

    Parameters
    ----------
    adata_raw:
        AnnData *before* QC filtering, but *after* QC metrics were computed.
        Must have obs columns: n_genes_by_counts, total_counts, pct_counts_mt,
        doublet_score, predicted_doublet.
    adata_filtered:
        AnnData after all QC filters were applied.
    metrics:
        Dict returned by ``run_qc()``.
    output_path:
        Where to write the HTML file.
    sample_name:
        Label shown in the report header.

    Returns
    -------
    Path
        Resolved path to the written HTML file.
    """
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    t = metrics["thresholds"]

    # ── Build all plots ──────────────────────────────────────────────────
    logger.info("Generating QC report plots...")

    plots = {}

    before_genes = adata_raw.obs["n_genes_by_counts"].values
    after_genes  = adata_filtered.obs["n_genes_by_counts"].values
    plots["violin_genes"] = _violin_plot(
        before_genes, after_genes, "Genes per cell",
        threshold_low=t["min_genes"], threshold_high=t["max_genes"],
        colour=_PALETTE["accent"],
    )

    before_umi = adata_raw.obs["total_counts"].values
    after_umi  = adata_filtered.obs["total_counts"].values
    plots["violin_umi"] = _violin_plot(
        before_umi, after_umi, "Total UMI counts",
        colour=_PALETTE["accent2"],
    )

    before_mt = adata_raw.obs["pct_counts_mt"].values
    after_mt  = adata_filtered.obs["pct_counts_mt"].values
    plots["violin_mt"] = _violin_plot(
        before_mt, after_mt, "Mitochondrial %",
        threshold_high=t["max_mt_pct"],
        colour=_PALETTE["accent3"],
    )

    plots["scatter"] = _scatter_umi_genes(adata_raw, t["max_mt_pct"])

    # Doublet histogram — use scores if available
    if "doublet_score" in adata_raw.obs.columns:
        scores = adata_raw.obs["doublet_score"].values
        if not np.all(np.isnan(scores)):
            # Estimate threshold as midpoint between score peaks
            threshold = float(np.nanpercentile(scores, 97))
            plots["doublet_hist"] = _doublet_histogram(scores, threshold)
        else:
            plots["doublet_hist"] = _doublet_histogram(np.zeros(10), 0.25)
    else:
        plots["doublet_hist"] = _doublet_histogram(np.zeros(10), 0.25)

    # ── Ground-truth MT% correlation ─────────────────────────────────────
    corr_value = None
    if "GEX_pct_counts_mt" in adata_raw.obs.columns:
        our_mt = adata_raw.obs["pct_counts_mt"].values
        gt_mt  = adata_raw.obs["GEX_pct_counts_mt"].values
        valid  = np.isfinite(our_mt) & np.isfinite(gt_mt)
        if valid.sum() > 10:
            corr_value = float(np.corrcoef(our_mt[valid], gt_mt[valid])[0, 1])
            logger.info("MT%% ground-truth correlation: %.4f", corr_value)
            plots["mt_corr"] = _mt_correlation_plot(our_mt[valid], gt_mt[valid], corr_value)

    # ── Write HTML ───────────────────────────────────────────────────────
    html = _build_html(sample_name, metrics, plots, corr_value, timestamp)
    output_path.write_text(html, encoding="utf-8")

    size_kb = output_path.stat().st_size / 1024
    logger.info("QC report written → %s (%.1f KB)", output_path, size_kb)
    print(f"\n  ✓ QC report → {output_path}  ({size_kb:.0f} KB)\n")

    return output_path
