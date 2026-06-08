"""
spatial_qc_report.py — OmicSage Phase 7
HTML report for spatial QC results.
Matches the _render_page / section structure of the RNA pipeline reports.
"""

from __future__ import annotations

import base64
import logging
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


def _fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _violin_plots(adata, params):
    metrics = [
        ("total_counts",      "Total UMI counts / spot",
         params.get("min_counts"), params.get("max_counts")),
        ("n_genes_by_counts", "Genes detected / spot",
         params.get("min_genes"),  params.get("max_genes")),
        ("pct_counts_mt",     "MT gene % / spot",
         None,                     params.get("max_mt_pct")),
    ]
    available = [m for m in metrics if m[0] in adata.obs.columns]
    if not available:
        return None
    fig, axes = plt.subplots(1, len(available), figsize=(5 * len(available), 4))
    if len(available) == 1:
        axes = [axes]
    for ax, (col, label, lo, hi) in zip(axes, available):
        vals = adata.obs[col].values
        parts = ax.violinplot(vals, positions=[0], showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor("#4C78A8"); pc.set_alpha(0.7)
        if lo is not None:
            ax.axhline(lo, color="#e07b3a", linestyle="--", linewidth=1.2,
                       label=f"min={lo}")
        if hi is not None:
            ax.axhline(hi, color="#e05252", linestyle="--", linewidth=1.2,
                       label=f"max={hi}")
        ax.set_xticks([])
        ax.set_ylabel(label)
        ax.set_title(label, fontsize=10, fontweight="bold")
        if lo is not None or hi is not None:
            ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("QC metrics per spot", fontsize=11, fontweight="bold")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _spatial_scatter(adata, color_key):
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    try:
        fig, ax = plt.subplots(figsize=(5, 5))
        sq.pl.spatial_scatter(adata, color=color_key, ax=ax, show=False, frameon=False)
        ax.set_title(color_key.replace("_", " ").title(), fontsize=10, fontweight="bold")
        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception:
        return None


def _threshold_bar(outputs):
    categories = [
        ("Low counts",  outputs.get("removed_low_counts",  0)),
        ("High counts", outputs.get("removed_high_counts", 0)),
        ("Low genes",   outputs.get("removed_low_genes",   0)),
        ("High genes",  outputs.get("removed_high_genes",  0)),
        ("High MT%",    outputs.get("removed_high_mt",     0)),
    ]
    labels, values = zip(*categories)
    fig, ax = plt.subplots(figsize=(6, 3))
    colors = ["#e05252" if v > 0 else "#95a5a6" for v in values]
    bars = ax.bar(labels, values, color=colors, edgecolor="white")
    ax.set_ylabel("Spots removed")
    ax.set_title("Spots removed per filter criterion", fontsize=10, fontweight="bold")
    for bar, val in zip(bars, values):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5, str(val),
                    ha="center", va="bottom", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _section_summary(adata, qc_info, dataset_id, timestamp):
    outputs = qc_info.get("outputs", {})
    params  = qc_info.get("params",  {})
    stats   = qc_info.get("summary_stats", {})
    n_before  = outputs.get("n_spots_before", 0)
    n_after   = outputs.get("n_spots_after",  0)
    n_removed = outputs.get("n_spots_removed", 0)
    pct_kept  = f"{100 * n_after / n_before:.1f}%" if n_before > 0 else "?"

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Spots input",   f"{n_before:,}"),
            ("Spots kept",    f"{n_after:,}"),
            ("Spots removed", f"{n_removed:,}"),
            ("Pass rate",     pct_kept),
            ("Genes",         f"{adata.n_vars:,}"),
        ]
    )
    threshold_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in params.items()
    )
    stat_rows = ""
    for col, label in [
        ("total_counts", "Total UMI counts"),
        ("n_genes_by_counts", "Genes detected"),
        ("pct_counts_mt", "MT gene %"),
    ]:
        s = stats.get(col, {})
        if s:
            stat_rows += (
                f"<tr><td>{label}</td>"
                f"<td>{s['mean']:.1f}</td><td>{s['median']:.1f}</td>"
                f"<td>{s['std']:.1f}</td><td>{s['min']:.1f}</td>"
                f"<td>{s['max']:.1f}</td></tr>"
            )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>QC Metric Summary (retained spots)</h3>
      <table>
        <thead>
          <tr><th>Metric</th><th>Mean</th><th>Median</th><th>Std</th><th>Min</th><th>Max</th></tr>
        </thead>
        <tbody>{stat_rows}</tbody>
      </table>
      <h3>Filter Thresholds Applied</h3>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{threshold_rows}</tbody>
      </table>
    </section>
    """


def _section_distributions(adata, params):
    b64 = _violin_plots(adata, params)
    if not b64:
        return ""
    return f"""
    <section>
      <h2>QC Metric Distributions</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Violin plots per spot</h3>
          <img src="data:image/png;base64,{b64}" alt="violin plots">
        </div>
      </div>
    </section>
    """


def _section_spatial(adata):
    b64_counts = _spatial_scatter(adata, "total_counts")
    b64_mt     = _spatial_scatter(adata, "pct_counts_mt")
    if not b64_counts and not b64_mt:
        return ""
    figs = ""
    if b64_counts:
        figs += (
            '<div class="fig-wrap"><h3>Total UMI counts</h3>'
            f'<img src="data:image/png;base64,{b64_counts}" alt="spatial counts"></div>'
        )
    if b64_mt:
        figs += (
            '<div class="fig-wrap"><h3>MT gene %</h3>'
            f'<img src="data:image/png;base64,{b64_mt}" alt="spatial MT%"></div>'
        )
    return f"""
    <section>
      <h2>Spatial Distribution of QC Metrics</h2>
      <div class="fig-grid">{figs}</div>
    </section>
    """


def _section_filter_breakdown(outputs):
    b64 = _threshold_bar(outputs)
    return f"""
    <section>
      <h2>Filter Breakdown</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Spots removed per filter criterion</h3>
          <img src="data:image/png;base64,{b64}" alt="threshold bar">
        </div>
      </div>
      <p class="note">Spots may fail multiple filters simultaneously &mdash; counts may overlap.</p>
    </section>
    """


_PAGE_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
             color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
              box-shadow: 0 1px 4px rgba(0,0,0,0.07);
              padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                 border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
    section h3 { font-size: 1rem; font-weight: 600; color: #16213e; margin: 18px 0 10px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
            border-left: 3px solid #f0c040; padding: 8px 12px;
            border-radius: 4px; margin-bottom: 14px; }
    code { font-family: "SFMono-Regular", Consolas, monospace;
           background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em; }
    .stat-grid { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }
    .stat-card { background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
                 min-width: 130px; text-align: center; flex: 1 1 130px; }
    .stat-value { font-size: 1.4rem; font-weight: 700; color: #0f3460; }
    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }
    table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }
    th { background: #f0f2ff; color: #0f3460; font-weight: 600;
         padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9ff; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 560px; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
"""


def _render_page(title, sections, timestamp):
    body = "\n".join(sections)
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "  <meta charset=\"UTF-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
        f"  <title>{title}</title>\n"
        f"  <style>{_PAGE_CSS}</style>\n"
        "</head>\n<body>\n"
        "  <header>\n"
        "    <h1>OmicSage &#8212; Spatial QC Report</h1>\n"
        f"    <p>Generated {timestamp}</p>\n"
        "  </header>\n"
        "  <main>\n"
        f"    {body}\n"
        "  </main>\n"
        "  <footer>\n"
        "    Generated by <a href=\"https://github.com/fshokor/OmicSage\">OmicSage</a>\n"
        "    &middot; MIT License\n"
        "  </footer>\n"
        "</body>\n</html>"
    )


def generate_spatial_qc_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
) -> str:
    if "omicsage_spatial_qc" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_qc'] not found. "
            "Run spatial_qc() before generating the report."
        )
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    qc_info   = adata.uns["omicsage_spatial_qc"]
    params    = qc_info.get("params",  {})
    outputs   = qc_info.get("outputs", {})

    sections = [
        _section_summary(adata, qc_info, dataset_id, timestamp),
        _section_distributions(adata, params),
        _section_spatial(adata),
        _section_filter_breakdown(outputs),
    ]
    html = _render_page(
        title=f"OmicSage -- Spatial QC -- {dataset_id}",
        sections=sections,
        timestamp=timestamp,
    )
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    size_kb = Path(output_path).stat().st_size / 1024
    logger.info("Spatial QC report -> %s (%.1f KB)", output_path, size_kb)
    return str(Path(output_path).resolve())
