"""
OmicSage — CITE-seq ADT Normalization Report
reports/templates/cite/cite_normalize_report.py

Generated after normalize_adt step.
Output: cite_01_normalize_report.html

Usage
-----
    from reports.templates.cite.cite_normalize_report import run_cite_normalize_report
    run_cite_normalize_report(
        adt=adt_norm,
        metrics=metrics,
        report_path="reports/GSE194122/cite_01_normalize_report.html",
        dataset_name="BMMC CITE-seq (NeurIPS 2021)",
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


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _plot_clr_violin(adt: AnnData) -> str:
    """Violin plot of CLR values across proteins (up to 40)."""
    import scipy.sparse as sp
    clr = adt.layers.get("adt_clr", adt.X)
    if sp.issparse(clr):
        clr = clr.toarray()
    clr = np.asarray(clr, dtype=float)

    n_proteins = clr.shape[1]
    max_show = 40
    idx = np.linspace(0, n_proteins - 1, min(max_show, n_proteins), dtype=int)
    clr_sub = clr[:, idx]
    labels = [adt.var_names[i] for i in idx]

    fig, ax = plt.subplots(figsize=(max(10, len(idx) * 0.35), 4.5))
    parts = ax.violinplot(
        [clr_sub[:, j] for j in range(clr_sub.shape[1])],
        positions=range(clr_sub.shape[1]),
        showmedians=True, widths=0.7,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#4C78A8"); pc.set_alpha(0.6)
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color("#16213e"); parts[key].set_linewidth(1.0)
    ax.set_xticks(range(clr_sub.shape[1]))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("CLR value", fontsize=10)
    ax.set_title(
        f"CLR-Normalised ADT Distribution"
        + (f" (showing {len(idx)} of {n_proteins})" if n_proteins > max_show else ""),
        fontsize=12, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_clr_heatmap(adt: AnnData) -> str:
    """Per-protein CLR mean bar chart."""
    import scipy.sparse as sp
    clr = adt.layers.get("adt_clr", adt.X)
    if sp.issparse(clr):
        clr = clr.toarray()
    clr = np.asarray(clr, dtype=float)
    means = clr.mean(axis=0)
    order = np.argsort(means)[::-1]
    labels = [adt.var_names[i] for i in order]
    values = means[order]

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.28), 4))
    ax.bar(range(len(labels)), values, color="#e07b3a", alpha=0.8, width=0.7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Mean CLR value", fontsize=10)
    ax.set_title("Mean CLR Expression per Protein (ranked)", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_raw_vs_clr(adt: AnnData) -> str:
    """Scatter: raw total counts vs CLR sum per cell."""
    import scipy.sparse as sp
    raw = adt.layers.get("counts", None)
    clr = adt.layers.get("adt_clr", adt.X)
    if raw is None:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "Raw counts layer not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)
    if sp.issparse(raw): raw = raw.toarray()
    if sp.issparse(clr): clr = clr.toarray()
    raw_sum = np.asarray(raw, dtype=float).sum(axis=1)
    clr_sum = np.asarray(clr, dtype=float).sum(axis=1)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(raw_sum, clr_sum, s=2, alpha=0.4, color="#4C78A8", rasterized=True)
    ax.set_xlabel("Raw total ADT counts per cell", fontsize=10)
    ax.set_ylabel("CLR sum per cell", fontsize=10)
    ax.set_title("Raw Counts vs CLR Sum", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML
# ---------------------------------------------------------------------------

def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq ADT Normalization Report -- {dataset_name}</title>
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
    section p  {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    code {{ font-family: "SFMono-Regular", Consolas, monospace;
            background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em; }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{ background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
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
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq ADT Normalization Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",      f"{metrics.get('n_cells', '?'):,}"),
            ("Proteins",   str(metrics.get("n_proteins", "?"))),
            ("CLR axis",   str(metrics.get("clr_axis", "?"))),
            ("DSB applied", "Yes" if metrics.get("dsb_applied") else "No"),
            ("CLR max",    f"{metrics.get('clr_max', 0):.2f}"),
            ("CLR min",    f"{metrics.get('clr_min', 0):.2f}"),
        ]
    )
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("CLR axis description", "per-protein across cells (axis=0)" if metrics.get("clr_axis") == 0
                                     else "per-cell across proteins (axis=1)"),
            ("Raw counts layer",     metrics.get("raw_counts_in_layer", "counts")),
            ("CLR layer",            metrics.get("clr_in_layer", "adt_clr")),
            ("CLR mean per cell",    f"{metrics.get('clr_mean_per_cell', 0):.3f}"),
            ("Raw median total",     f"{metrics.get('raw_median_total_counts_per_cell', 0):.1f}"),
        ]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{cards}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows}</tbody></table>
    </section>"""


def _section_figures(fig_violin: str, fig_bar: str, fig_scatter: str) -> str:
    return f"""
    <section>
      <h2>Figures</h2>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>CLR Distribution per Protein (violin)</h3>
          <img src="data:image/png;base64,{fig_violin}" alt="CLR violin">
        </div>
        <div class="fig-wrap wide">
          <h3>Mean CLR Expression per Protein (ranked)</h3>
          <img src="data:image/png;base64,{fig_bar}" alt="CLR mean bar">
        </div>
        <div class="fig-wrap">
          <h3>Raw Counts vs CLR Sum per Cell</h3>
          <img src="data:image/png;base64,{fig_scatter}" alt="Raw vs CLR scatter">
        </div>
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_cite_normalize_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_01_normalize_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building CITE normalize report for '{dataset_name}' ...", flush=True)
    print("  Rendering CLR violin ...", flush=True)
    fig_violin  = _plot_clr_violin(adt)
    print("  Rendering CLR mean bar ...", flush=True)
    fig_bar     = _plot_clr_heatmap(adt)
    print("  Rendering raw vs CLR scatter ...", flush=True)
    fig_scatter = _plot_raw_vs_clr(adt)

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_figures(fig_violin, fig_bar, fig_scatter),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
