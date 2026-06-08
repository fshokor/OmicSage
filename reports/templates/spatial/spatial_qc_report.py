"""
spatial_qc_report.py — OmicSage Phase 7, Session 1
HTML report for spatial QC results.

Generates a self-contained HTML report with:
  - QC summary table
  - Violin plots: total_counts, n_genes_by_counts, pct_counts_mt
  - Spatial scatter plots: total_counts and pct_counts_mt overlaid on tissue
  - Before/after spot count comparison
"""

from __future__ import annotations

import base64
import io
import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")

try:
    import squidpy as sq

    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_spatial_qc_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
) -> str:
    """Generate a self-contained HTML QC report for Visium data.

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_qc` (contains
        ``uns["omicsage_spatial_qc"]``).
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label used in the report title.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    if "omicsage_spatial_qc" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_qc'] not found. "
            "Run spatial_qc() before generating the report."
        )

    qc_info = adata.uns["omicsage_spatial_qc"]
    figures = _build_figures(adata, qc_info)
    html = _render_html(adata, qc_info, figures, dataset_id)

    output_path = str(output_path)
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    return os.path.abspath(output_path)


# ---------------------------------------------------------------------------
# Figure generation
# ---------------------------------------------------------------------------


def _fig_to_base64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _build_figures(adata: ad.AnnData, qc_info: dict) -> dict[str, str]:
    """Build all figures and return them as base64-encoded PNG strings."""
    figures = {}

    # ---- 1. Violin plots -------------------------------------------------
    figures["violin"] = _violin_plots(adata, qc_info)

    # ---- 2. Spatial scatter plots (if squidpy available + spatial ok) ----
    if _SQUIDPY_AVAILABLE and "spatial" in adata.obsm and "spatial" in adata.uns:
        figures["spatial_counts"] = _spatial_scatter(adata, "total_counts")
        figures["spatial_mt"] = _spatial_scatter(adata, "pct_counts_mt")
    else:
        figures["spatial_counts"] = None
        figures["spatial_mt"] = None

    # ---- 3. Threshold bar chart -----------------------------------------
    figures["threshold_bar"] = _threshold_bar(qc_info)

    return figures


def _violin_plots(adata: ad.AnnData, qc_info: dict) -> str:
    params = qc_info.get("params", {})
    metrics = [
        ("total_counts", "Total UMI counts / spot",
         params.get("min_counts"), params.get("max_counts")),
        ("n_genes_by_counts", "Genes detected / spot",
         params.get("min_genes"), params.get("max_genes")),
        ("pct_counts_mt", "MT gene % / spot",
         None, params.get("max_mt_pct")),
    ]
    available = [m for m in metrics if m[0] in adata.obs.columns]

    fig, axes = plt.subplots(1, len(available), figsize=(5 * len(available), 4))
    if len(available) == 1:
        axes = [axes]

    for ax, (col, label, lo, hi) in zip(axes, available):
        vals = adata.obs[col].values
        parts = ax.violinplot(vals, positions=[0], showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor("#4C72B0")
            pc.set_alpha(0.7)
        # threshold lines
        if lo is not None:
            ax.axhline(lo, color="red", linestyle="--", linewidth=1,
                       label=f"min={lo}")
        if hi is not None:
            ax.axhline(hi, color="orange", linestyle="--", linewidth=1,
                       label=f"max={hi}")
        ax.set_xticks([])
        ax.set_ylabel(label)
        ax.set_title(label)
        if lo is not None or hi is not None:
            ax.legend(fontsize=8)

    fig.suptitle("QC metrics per spot", fontsize=12, fontweight="bold")
    fig.tight_layout()
    return _fig_to_base64(fig)


def _spatial_scatter(adata: ad.AnnData, color_key: str) -> Optional[str]:
    """Spatial scatter coloured by a QC metric using squidpy."""
    try:
        fig, ax = plt.subplots(figsize=(5, 5))
        sq.pl.spatial_scatter(
            adata,
            color=color_key,
            ax=ax,
            show=False,
            frameon=False,
        )
        ax.set_title(color_key.replace("_", " ").title())
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _threshold_bar(qc_info: dict) -> str:
    outputs = qc_info.get("outputs", {})
    categories = [
        ("Low counts", outputs.get("removed_low_counts", 0)),
        ("High counts", outputs.get("removed_high_counts", 0)),
        ("Low genes", outputs.get("removed_low_genes", 0)),
        ("High genes", outputs.get("removed_high_genes", 0)),
        ("High MT%", outputs.get("removed_high_mt", 0)),
    ]
    labels, values = zip(*categories)

    fig, ax = plt.subplots(figsize=(6, 3))
    colors = ["#e74c3c" if v > 0 else "#95a5a6" for v in values]
    bars = ax.bar(labels, values, color=colors, edgecolor="white")
    ax.set_ylabel("Spots removed")
    ax.set_title("Spots removed per filter criterion\n(spots may overlap)")
    for bar, val in zip(bars, values):
        if val > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5,
                str(val),
                ha="center",
                va="bottom",
                fontsize=9,
            )
    ax.tick_params(axis="x", labelsize=8)
    fig.tight_layout()
    return _fig_to_base64(fig)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


def _render_html(
    adata: ad.AnnData,
    qc_info: dict,
    figures: dict,
    dataset_id: str,
) -> str:
    outputs = qc_info.get("outputs", {})
    params = qc_info.get("params", {})
    stats = qc_info.get("summary_stats", {})
    timestamp = qc_info.get("timestamp", datetime.now().isoformat())

    n_before = outputs.get("n_spots_before", "?")
    n_after = outputs.get("n_spots_after", "?")
    n_removed = outputs.get("n_spots_removed", "?")
    pct_kept = (
        f"{100 * n_after / n_before:.1f}%"
        if isinstance(n_before, int) and n_before > 0
        else "?"
    )

    def stat_row(metric_key, label):
        s = stats.get(metric_key, {})
        if not s:
            return ""
        return (
            f"<tr><td>{label}</td>"
            f"<td>{s['mean']:.1f}</td>"
            f"<td>{s['median']:.1f}</td>"
            f"<td>{s['std']:.1f}</td>"
            f"<td>{s['min']:.1f}</td>"
            f"<td>{s['max']:.1f}</td></tr>"
        )

    def img_section(title, b64, alt="figure"):
        if b64 is None:
            return f"<p><em>{title} — not available (squidpy or spatial data missing)</em></p>"
        return (
            f"<h3>{title}</h3>"
            f'<img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="max-width:100%;border:1px solid #e0e0e0;border-radius:4px;">'
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage — Spatial QC Report: {dataset_id}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          margin: 0; padding: 0; background: #f8f9fa; color: #212529; }}
  .header {{ background: linear-gradient(135deg, #2c3e50, #3498db);
             color: white; padding: 2rem 2.5rem; }}
  .header h1 {{ margin: 0 0 0.3rem; font-size: 1.6rem; }}
  .header p  {{ margin: 0; opacity: 0.85; font-size: 0.9rem; }}
  .container {{ max-width: 1100px; margin: 0 auto; padding: 1.5rem 2rem; }}
  .card {{ background: white; border-radius: 8px; padding: 1.5rem 2rem;
           margin-bottom: 1.5rem; box-shadow: 0 1px 4px rgba(0,0,0,0.08); }}
  .card h2 {{ margin-top: 0; font-size: 1.15rem; color: #2c3e50;
              border-bottom: 2px solid #3498db; padding-bottom: 0.4rem; }}
  .card h3 {{ font-size: 1rem; color: #34495e; margin: 1rem 0 0.4rem; }}
  .kpi-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
               gap: 1rem; margin-bottom: 0.5rem; }}
  .kpi {{ background: #f0f4f8; border-radius: 6px; padding: 1rem;
          text-align: center; border-left: 4px solid #3498db; }}
  .kpi .value {{ font-size: 1.8rem; font-weight: 700; color: #2c3e50; }}
  .kpi .label {{ font-size: 0.8rem; color: #7f8c8d; margin-top: 0.2rem; }}
  .kpi.warn {{ border-left-color: #e74c3c; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
  th {{ background: #f0f4f8; padding: 0.5rem 0.8rem; text-align: left;
        font-weight: 600; color: #2c3e50; }}
  td {{ padding: 0.4rem 0.8rem; border-bottom: 1px solid #f0f0f0; }}
  tr:last-child td {{ border-bottom: none; }}
  .param-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                 gap: 0.4rem; font-size: 0.88rem; }}
  .param-item {{ display: flex; justify-content: space-between;
                 background: #f8f9fa; padding: 0.3rem 0.6rem; border-radius: 4px; }}
  .param-item .key {{ color: #7f8c8d; }}
  .param-item .val {{ font-weight: 600; color: #2c3e50; }}
  .footer {{ text-align: center; padding: 1.5rem; color: #95a5a6; font-size: 0.8rem; }}
</style>
</head>
<body>
<div class="header">
  <h1>🔬 OmicSage — Spatial QC Report</h1>
  <p>Dataset: <strong>{dataset_id}</strong> &nbsp;|&nbsp; Generated: {timestamp[:19].replace("T", " ")}</p>
</div>

<div class="container">

  <!-- KPI cards -->
  <div class="card">
    <h2>Spot Summary</h2>
    <div class="kpi-grid">
      <div class="kpi">
        <div class="value">{n_before:,}</div>
        <div class="label">Spots (input)</div>
      </div>
      <div class="kpi {'warn' if isinstance(n_removed, int) and n_removed > 0 else ''}">
        <div class="value">{n_removed:,}</div>
        <div class="label">Spots removed</div>
      </div>
      <div class="kpi">
        <div class="value">{n_after:,}</div>
        <div class="label">Spots retained</div>
      </div>
      <div class="kpi">
        <div class="value">{pct_kept}</div>
        <div class="label">% kept</div>
      </div>
      <div class="kpi">
        <div class="value">{adata.n_vars:,}</div>
        <div class="label">Genes</div>
      </div>
    </div>
  </div>

  <!-- QC metrics summary table -->
  <div class="card">
    <h2>QC Metric Summary (retained spots)</h2>
    <table>
      <thead>
        <tr><th>Metric</th><th>Mean</th><th>Median</th><th>Std</th><th>Min</th><th>Max</th></tr>
      </thead>
      <tbody>
        {stat_row("total_counts", "Total UMI counts")}
        {stat_row("n_genes_by_counts", "Genes detected")}
        {stat_row("pct_counts_mt", "MT gene %")}
      </tbody>
    </table>
  </div>

  <!-- Threshold parameters -->
  <div class="card">
    <h2>Filter Thresholds Applied</h2>
    <div class="param-grid">
      <div class="param-item"><span class="key">min_counts</span><span class="val">{params.get('min_counts','?'):,}</span></div>
      <div class="param-item"><span class="key">max_counts</span><span class="val">{params.get('max_counts','?'):,}</span></div>
      <div class="param-item"><span class="key">min_genes</span><span class="val">{params.get('min_genes','?')}</span></div>
      <div class="param-item"><span class="key">max_genes</span><span class="val">{params.get('max_genes','?'):,}</span></div>
      <div class="param-item"><span class="key">max_mt_pct</span><span class="val">{params.get('max_mt_pct','?')}%</span></div>
      <div class="param-item"><span class="key">mt_prefix</span><span class="val">{params.get('mt_prefix','?')}</span></div>
      <div class="param-item"><span class="key">filter_spots</span><span class="val">{params.get('filter_spots','?')}</span></div>
    </div>
  </div>

  <!-- Figures -->
  <div class="card">
    <h2>QC Figures</h2>
    {img_section("Violin plots — QC metrics per spot", figures.get("violin"), "violin plots")}
    {img_section("Spots removed per filter criterion", figures.get("threshold_bar"), "threshold bar chart")}
    {img_section("Spatial distribution — Total UMI counts", figures.get("spatial_counts"), "spatial counts")}
    {img_section("Spatial distribution — MT gene %", figures.get("spatial_mt"), "spatial MT%")}
  </div>

</div>
<div class="footer">
  Generated by OmicSage &nbsp;|&nbsp; Phase 7 — Spatial Transcriptomics
</div>
</body>
</html>"""

    return html
