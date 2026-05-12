"""
OmicSage — Normalization Report
pipeline/modules/qc/normalization_report.py

Usage
-----
    from pipeline.modules.qc.normalization_report import run_normalization_report
    run_normalization_report(
        adata_norm=adata_norm,
        metrics=metrics,
        report_path="reports/GSE194122/02_normalization_report.html",
        dataset_name="GSE194122_CITE",
    )
"""

from __future__ import annotations

import base64
import io
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_hvg_scatter(adata_norm: AnnData) -> str:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    var = adata_norm.var.copy()
    if "variances_norm" in var.columns:
        y_col, y_label = "variances_norm", "Normalised variance"
        x_col = "means" if "means" in var.columns else var.columns[0]
    elif "dispersions_norm" in var.columns:
        y_col, y_label = "dispersions_norm", "Normalised dispersion"
        x_col = "means" if "means" in var.columns else var.columns[0]
    else:
        n_hvg = int(var["highly_variable"].sum())
        n_non = len(var) - n_hvg
        fig2, ax2 = plt.subplots(figsize=(4, 4))
        ax2.bar(["HVG", "Non-HVG"], [n_hvg, n_non], color=["#e07b3a", "#aaaaaa"])
        ax2.set_ylabel("Gene count"); ax2.set_title("Highly Variable Genes")
        ax2.spines[["top", "right"]].set_visible(False)
        return _fig_to_b64(fig2)

    hvg_mask = var["highly_variable"]
    ax.scatter(var.loc[~hvg_mask, x_col], var.loc[~hvg_mask, y_col],
               s=4, alpha=0.4, color="#aaaaaa", label="Non-HVG", rasterized=True)
    ax.scatter(var.loc[hvg_mask, x_col], var.loc[hvg_mask, y_col],
               s=6, alpha=0.7, color="#e07b3a",
               label=f"HVG (n={hvg_mask.sum()})", rasterized=True)
    ax.set_xlabel("Mean expression", fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title("Highly Variable Genes", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_library_size(adata_norm: AnnData) -> str:
    import scipy.sparse as sp
    raw  = adata_norm.layers["counts"]
    norm = adata_norm.layers["logcounts"]
    raw_arr  = raw.toarray()  if sp.issparse(raw)  else np.asarray(raw)
    norm_arr = norm.toarray() if sp.issparse(norm) else np.asarray(norm)
    raw_sums  = raw_arr.sum(axis=1)
    norm_sums = norm_arr.sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    for ax, data, title, color in zip(
        axes,
        [raw_sums, norm_sums],
        ["Raw counts per cell", "log1p-normalised sum per cell"],
        ["#4C78A8", "#54a868"],
    ):
        ax.violinplot(data, positions=[0], showmedians=True, widths=0.6)
        for pc in ax.collections:
            pc.set_facecolor(color); pc.set_alpha(0.6)
        ax.set_xticks([])
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel("Sum per cell", fontsize=10)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))
        ax.spines[["top", "right"]].set_visible(False)
        median = np.median(data)
        ax.axhline(median, color="#333", linewidth=1, linestyle="--", alpha=0.5)
        ax.text(0.55, median, f"median={median:,.1f}", va="center",
                fontsize=8, transform=ax.get_yaxis_transform())
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_top_hvgs(adata_norm: AnnData, n: int = 20) -> str:
    var = adata_norm.var[adata_norm.var["highly_variable"]].copy()
    rank_col = None
    for col in ["variances_norm", "dispersions_norm", "highly_variable_rank"]:
        if col in var.columns:
            rank_col = col; break
    if rank_col is None or rank_col == "highly_variable_rank":
        top = var.head(n); values = np.ones(len(top))
        xlabel = "HVG (no variance metric available)"
    else:
        top = var.nlargest(n, rank_col); values = top[rank_col].values
        xlabel = rank_col.replace("_", " ").title()
    fig, ax = plt.subplots(figsize=(7, max(3, n * 0.32)))
    ax.barh(range(len(top)), values[::-1], color="#e07b3a", alpha=0.8)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index.tolist()[::-1], fontsize=9)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_title(f"Top {n} Highly Variable Genes", fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_gene_detection(adata_norm: AnnData) -> str:
    import scipy.sparse as sp
    raw = adata_norm.layers["counts"]
    raw_arr = raw.toarray() if sp.issparse(raw) else np.asarray(raw)
    detection_rate = (raw_arr > 0).mean(axis=0)
    sorted_rates = np.sort(detection_rate)[::-1]
    x = np.arange(1, len(sorted_rates) + 1)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, sorted_rates, color="#4C78A8", linewidth=1.5)
    ax.axhline(0.01, color="#aaa", linestyle="--", linewidth=0.8, label="1% cells")
    ax.axhline(0.10, color="#aaa", linestyle=":",  linewidth=0.8, label="10% cells")
    ax.set_xlabel("Genes ranked by detection rate", fontsize=11)
    ax.set_ylabel("Fraction of cells detected", fontsize=11)
    ax.set_title("Gene Detection Rate", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML renderer — matches deg_report._render_page exactly
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
      font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc;
    }}
    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white; padding: 32px 40px 24px;
    }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}
    section {{
      background: white; border-radius: 10px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.07);
      padding: 28px 32px; margin-bottom: 24px;
    }}
    section h2 {{
      font-size: 1.15rem; font-weight: 700; color: #0f3460;
      border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px;
    }}
    section h3 {{ font-size: 1rem; font-weight: 600; color: #16213e; margin: 18px 0 10px; }}
    section p {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em;
    }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{
      background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
      min-width: 130px; text-align: center; flex: 1 1 130px;
    }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }}
    th {{
      background: #f0f2ff; color: #0f3460; font-weight: 600;
      padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0;
    }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}
    .fig-grid {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }}
    .fig-wrap {{ flex: 1 1 300px; max-width: 520px; }}
    .fig-wrap.wide {{ flex: 1 1 100%; max-width: 100%; }}
    .fig-wrap h3 {{ font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }}
    .fig-wrap img {{ width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }}
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage -- Normalization Report</h1>
    <p>Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    &middot; MIT License
  </footer>
</body>
</html>
"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",             f"{metrics.get('n_cells', '?'):,}" if isinstance(metrics.get('n_cells'), int) else metrics.get('n_cells', '?')),
            ("Genes",             f"{metrics.get('n_genes', '?'):,}" if isinstance(metrics.get('n_genes'), int) else metrics.get('n_genes', '?')),
            ("HVGs selected",     str(metrics.get("n_hvg_selected", "?"))),
            ("HVG flavor",        metrics.get("hvg_flavor", "?")),
            ("Target sum",        f"{metrics.get('target_sum', 1e4):,.0f}"),
            ("Batch key",         metrics.get("batch_key") or "None"),
        ]
    )
    param_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("log1p applied",          "Yes" if metrics.get("log1p_applied") else "No"),
            ("Raw counts layer",       metrics.get("raw_counts_in_layer", "counts")),
            ("Normalized layer",       metrics.get("normalized_in_layer", "logcounts")),
            ("Mean norm sum per cell", f"{metrics.get('mean_counts_per_cell_after_norm', 0):.3f}"),
        ]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{param_rows}</tbody>
      </table>
    </section>
    """


def _section_figures(plots: dict) -> str:
    return f"""
    <section>
      <h2>Figures</h2>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>Highly Variable Genes — Mean vs Variance</h3>
          <img src="data:image/png;base64,{plots['hvg']}" alt="HVG scatter">
        </div>
        <div class="fig-wrap">
          <h3>Library Size Distribution</h3>
          <img src="data:image/png;base64,{plots['libsize']}" alt="Library size">
        </div>
        <div class="fig-wrap">
          <h3>Gene Detection Rate</h3>
          <img src="data:image/png;base64,{plots['detection']}" alt="Gene detection">
        </div>
        <div class="fig-wrap wide">
          <h3>Top 20 Highly Variable Genes</h3>
          <img src="data:image/png;base64,{plots['top_hvg']}" alt="Top HVGs">
        </div>
      </div>
    </section>
    """


def _section_provenance(adata_norm: AnnData) -> str:
    prov = adata_norm.uns.get("omicsage_normalization", {})
    rows = "".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov.items()
    )
    return f"""
    <section>
      <h2>Provenance</h2>
      <table>
        <thead><tr><th>Key</th><th>Value</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_normalization_report(
    adata_norm: AnnData,
    metrics: dict,
    report_path: str = "reports/normalization_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the normalization HTML report and write it to disk.

    Parameters
    ----------
    adata_norm : AnnData
        Normalized AnnData returned by normalize(). Must have layers['counts'] and ['logcounts'].
    metrics : dict
        Metrics dict returned by normalize().
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Label shown in the report header.

    Returns
    -------
    str
        Absolute path to the written report file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building normalization report for '{dataset_name}' ...", flush=True)
    print("  Rendering HVG scatter ...", flush=True)
    plots = {
        "hvg":       _plot_hvg_scatter(adata_norm),
    }
    print("  Rendering library size ...", flush=True)
    plots["libsize"]   = _plot_library_size(adata_norm)
    print("  Rendering top HVGs ...", flush=True)
    plots["top_hvg"]   = _plot_top_hvgs(adata_norm)
    print("  Rendering gene detection ...", flush=True)
    plots["detection"] = _plot_gene_detection(adata_norm)

    sections = [
        _section_summary(metrics, dataset_name, timestamp),
        _section_figures(plots),
        _section_provenance(adata_norm),
    ]

    html = _render_page(
        title=f"OmicSage -- Normalization Report -- {dataset_name}",
        sections=sections,
        timestamp=timestamp,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
