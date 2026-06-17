"""
OmicSage -- Leiden Clustering Report
pipeline/modules/clustering/cluster_report.py

Usage
-----
    from pipeline.modules.clustering.cluster_report import run_cluster_report
    run_cluster_report(
        adata_clustered=adata_clustered,
        metrics=metrics,
        report_path="reports/GSE194122/04_cluster_report.html",
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
# Figure helpers — plotting logic unchanged from original
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_umap_resolutions(adata_clustered: AnnData, metrics: dict) -> str:
    if "X_umap" not in adata_clustered.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "UMAP not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    resolutions = metrics["resolutions"]
    best_res    = metrics["best_resolution"]
    umap        = adata_clustered.obsm["X_umap"]
    n = len(resolutions)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows), squeeze=False)
    axes_flat  = axes.flatten()

    for idx, res in enumerate(resolutions):
        ax      = axes_flat[idx]
        obs_key = f"leiden_{res:.4g}"
        if obs_key not in adata_clustered.obs.columns:
            ax.text(0.5, 0.5, f"leiden_{res:.4g}\nnot found",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off"); continue
        labels   = adata_clustered.obs[obs_key].astype("category")
        cats     = labels.cat.categories.tolist()
        cmap_    = plt.get_cmap("tab20", max(len(cats), 1))
        col_map  = {l: cmap_(i) for i, l in enumerate(cats)}
        ax.scatter(umap[:, 0], umap[:, 1], c=[col_map[l] for l in labels],
                   s=1.5, alpha=0.6, rasterized=True)
        ax.set_title(f"resolution = {res:.2g}  ({len(cats)} clusters)",
                     fontsize=10, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=8); ax.set_ylabel("UMAP 2", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        if res == best_res:
            for spine in ax.spines.values():
                spine.set_visible(True); spine.set_edgecolor("#e0a800"); spine.set_linewidth(2.5)
            ax.set_title(f"* resolution = {res:.2g}  ({len(cats)} clusters)  [selected]",
                         fontsize=10, fontweight="bold", color="#b8860b")

    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)
    fig.suptitle("UMAP -- Leiden Clustering across Resolutions",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_silhouette_bar(metrics: dict) -> str:
    resolutions = metrics["resolutions"]
    scores      = [metrics["silhouette_scores"][r] for r in resolutions]
    n_clusters  = [metrics["n_clusters"][r] for r in resolutions]
    best_res    = metrics["best_resolution"]
    colors      = ["#e07b3a" if r == best_res else "#4C78A8" for r in resolutions]

    fig, ax1 = plt.subplots(figsize=(max(6, len(resolutions) * 1.4), 4))
    x    = np.arange(len(resolutions))
    bars = ax1.bar(x, scores, color=colors, width=0.6, alpha=0.80, edgecolor="white")
    for bar, score in zip(bars, scores):
        ypos = bar.get_height() + 0.005 if score >= 0 else bar.get_height() - 0.02
        ax1.text(bar.get_x() + bar.get_width() / 2, ypos,
                 f"{score:.3f}", ha="center", va="bottom", fontsize=8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{r:.2g}" for r in resolutions], fontsize=10)
    ax1.set_xlabel("Leiden Resolution", fontsize=11)
    ax1.set_ylabel("Silhouette Score", fontsize=11, color="#4C78A8")
    ax1.axhline(0, color="#bbb", linewidth=0.8, linestyle="--")
    ax1.spines[["top"]].set_visible(False)

    ax2 = ax1.twinx()
    ax2.plot(x, n_clusters, color="#2ca02c", marker="o", linewidth=2,
             markersize=5, label="n clusters", zorder=5)
    for xi, nc in zip(x, n_clusters):
        ax2.text(xi, nc + max(n_clusters) * 0.03, str(nc),
                 ha="center", va="bottom", fontsize=8, color="#2ca02c")
    ax2.set_ylabel("Number of Clusters", fontsize=11, color="#2ca02c")
    ax2.spines[["top"]].set_visible(False)
    ax2.tick_params(axis="y", labelcolor="#2ca02c")

    n_expected = metrics.get("n_clusters_expected")
    if n_expected is not None:
        ax2.axhline(n_expected, color="#d62728", linewidth=1.2, linestyle=":",
                    label=f"Expected: {n_expected}")
        ax2.text(len(resolutions) - 0.5, n_expected,
                 f" expected={n_expected}", va="center", fontsize=8, color="#d62728")

    selection_reason = metrics.get("selection_reason", "")
    ax1.set_title(
        f"Silhouette Score & Cluster Count vs Resolution\n"
        f"Selected: res={best_res:.2g}  |  reason: {selection_reason}",
        fontsize=12, fontweight="bold",
    )
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_elements = [
        Patch(facecolor="#e07b3a", label=f"Selected (res={best_res:.2g})"),
        Patch(facecolor="#4C78A8", label="Other resolutions"),
        Line2D([0], [0], color="#2ca02c", marker="o", label="n clusters"),
    ]
    if n_expected is not None:
        legend_elements.append(
            Line2D([0], [0], color="#d62728", linestyle=":", label=f"Expected ({n_expected})")
        )
    ax1.legend(handles=legend_elements, frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_resolution_selection(metrics: dict) -> str:
    resolutions = metrics["resolutions"]
    n_clust     = [metrics["n_clusters"][r] for r in resolutions]
    deltas      = [metrics["n_clusters_delta"][r] for r in resolutions]
    stability   = [metrics["stability_scores"][r] for r in resolutions]
    silhouette  = [metrics["silhouette_scores"][r] for r in resolutions]
    best_res    = metrics["best_resolution"]
    best_idx    = resolutions.index(best_res)
    n_expected  = metrics.get("n_clusters_expected")
    reason      = metrics.get("selection_reason", "")
    x           = np.arange(len(resolutions))
    xlbls       = [f"{r:.2g}" for r in resolutions]
    vline_kw    = dict(color="#e07b3a", linewidth=1.8, linestyle="--", alpha=0.8)

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    fig.suptitle(
        f"Resolution Selection Diagnostics  --  selected={best_res:.2g}  "
        f"({metrics['best_n_clusters']} clusters, reason='{reason}')",
        fontsize=13, fontweight="bold",
    )

    ax = axes[0, 0]
    ax.plot(x, n_clust, color="#1f77b4", marker="o", linewidth=2, markersize=6)
    ax.axvline(best_idx, **vline_kw, label=f"Selected (res={best_res:.2g})")
    if n_expected is not None:
        ax.axhline(n_expected, color="#d62728", linewidth=1.2, linestyle=":",
                   label=f"Expected ({n_expected})")
    ax.legend(frameon=False, fontsize=8)
    for xi, nc in zip(x, n_clust):
        ax.text(xi, nc + max(n_clust) * 0.03, str(nc), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlbls)
    ax.set_title("Number of Clusters", fontsize=11, fontweight="bold")
    ax.set_xlabel("Resolution"); ax.set_ylabel("n clusters")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[0, 1]
    delta_colors = ["#aec7e8" if r != best_res else "#e07b3a" for r in resolutions]
    ax.bar(x, deltas, color=delta_colors, width=0.6, alpha=0.85, edgecolor="white")
    ax.axvline(best_idx, **vline_kw)
    for xi, d in zip(x, deltas):
        ax.text(xi, d + max(deltas) * 0.03 if d >= 0 else d - max(deltas) * 0.05,
                f"+{d}" if d > 0 else str(d), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlbls)
    ax.set_title("Cluster Count Delta (new clusters per step)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Resolution"); ax.set_ylabel("delta n clusters")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 0]
    stab_colors = ["#98df8a" if r != best_res else "#e07b3a" for r in resolutions]
    ax.bar(x, stability, color=stab_colors, width=0.6, alpha=0.85, edgecolor="white")
    ax.axvline(best_idx, **vline_kw)
    for xi, s in zip(x, stability):
        ax.text(xi, s + 0.01, f"{s:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlbls)
    ax.set_ylim(0, 1.15)
    ax.set_title("Stability Score  (1 = plateau)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Resolution"); ax.set_ylabel("Stability")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1, 1]
    sil_colors = ["#ffbb78" if r != best_res else "#e07b3a" for r in resolutions]
    ax.bar(x, silhouette, color=sil_colors, width=0.6, alpha=0.85, edgecolor="white")
    ax.axvline(best_idx, **vline_kw)
    ax.axhline(0, color="#bbb", linewidth=0.8, linestyle="--")
    for xi, s in zip(x, silhouette):
        yp = s + 0.005 if s >= 0 else s - 0.015
        ax.text(xi, yp, f"{s:.3f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(xlbls)
    ax.set_title("Silhouette Score  (geometric cluster quality)", fontsize=11, fontweight="bold")
    ax.set_xlabel("Resolution"); ax.set_ylabel("Silhouette")
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_cluster_sizes(adata_clustered: AnnData, metrics: dict) -> str:
    best_res = metrics["best_resolution"]
    obs_key  = f"leiden_{best_res:.4g}"
    if obs_key not in adata_clustered.obs.columns:
        obs_key = "leiden"
    if obs_key not in adata_clustered.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "Cluster labels not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    counts     = adata_clustered.obs[obs_key].value_counts().sort_values(ascending=False)
    n_clusters = len(counts)
    cmap_      = plt.get_cmap("tab20", n_clusters)
    fig, ax    = plt.subplots(figsize=(max(8, n_clusters * 0.5), 4))
    ax.bar(range(n_clusters), counts.values,
           color=[cmap_(i) for i in range(n_clusters)],
           width=0.7, alpha=0.85, edgecolor="white")
    ax.set_xticks(range(n_clusters))
    ax.set_xticklabels(counts.index.tolist(), fontsize=8, rotation=45, ha="right")
    ax.set_xlabel("Cluster", fontsize=11)
    ax.set_ylabel("Number of Cells", fontsize=11)
    ax.set_title(
        f"Cluster Size Distribution  (resolution={best_res:.2g}, "
        f"{n_clusters} clusters, {adata_clustered.n_obs:,} cells)",
        fontsize=12, fontweight="bold",
    )
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.spines[["top", "right"]].set_visible(False)
    median_size = int(np.median(counts.values))
    ax.axhline(median_size, color="#c0392b", linewidth=1.2, linestyle="--",
               label=f"Median: {median_size:,} cells")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_cell_type(adata_clustered: AnnData) -> Optional[str]:
    if "X_umap" not in adata_clustered.obsm:
        return None
    if "cell_type" not in adata_clustered.obs.columns:
        return None
    umap        = adata_clustered.obsm["X_umap"]
    cell_types  = adata_clustered.obs["cell_type"].astype(str)
    unique_types = sorted(cell_types.unique())
    cmap_       = plt.get_cmap("tab20", len(unique_types))
    col_map     = {ct: cmap_(i) for i, ct in enumerate(unique_types)}
    n_types     = len(unique_types)
    fig, ax     = plt.subplots(figsize=(8 + max(0, (n_types - 20) // 10) * 1.5,
                                        max(6, n_types * 0.22)))
    ax.scatter(umap[:, 0], umap[:, 1], c=[col_map[ct] for ct in cell_types],
               s=1.5, alpha=0.6, rasterized=True)
    from matplotlib.patches import Patch
    legend_el = [Patch(facecolor=col_map[ct], label=ct) for ct in unique_types]
    ax.legend(handles=legend_el, bbox_to_anchor=(1.02, 1), loc="upper left",
              borderaxespad=0, frameon=False, fontsize=7.5,
              ncol=max(1, n_types // 25), handlelength=1.2, labelspacing=0.4)
    ax.set_title("UMAP -- Ground-Truth Cell Type", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=10); ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
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
    <h1>OmicSage -- Leiden Clustering Report</h1>
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


def _section_summary(adata_clustered: AnnData, metrics: dict,
                     dataset_name: str, timestamp: str) -> str:
    best_res  = metrics["best_resolution"]
    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",           f"{adata_clustered.n_obs:,}"),
            ("Selected res",    f"{best_res:.2g}"),
            ("Clusters",        str(metrics["best_n_clusters"])),
            ("Silhouette",      f"{metrics['best_silhouette']:.4f}"),
            ("Stability",       f"{metrics['best_stability']:.4f}"),
            ("Selection",       metrics.get("selection_reason", "?")),
        ]
    )
    per_res_rows = "".join(
        f"<tr><td>res={r:.2g}{'  *' if r == best_res else ''}</td>"
        f"<td>{metrics['n_clusters'][r]}</td>"
        f"<td>{metrics['n_clusters_delta'][r]:+d}</td>"
        f"<td>{metrics['stability_scores'][r]:.3f}</td>"
        f"<td>{metrics['silhouette_scores'][r]:.4f}</td></tr>"
        for r in metrics["resolutions"]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>Per-Resolution Metrics</h3>
      <table>
        <thead><tr><th>Resolution</th><th>Clusters</th><th>Delta</th>
               <th>Stability</th><th>Silhouette</th></tr></thead>
        <tbody>{per_res_rows}</tbody>
      </table>
    </section>
    """


def _section_figures(fig_umap: str, fig_sil: str, fig_diag: str,
                     fig_sizes: str, fig_ct: Optional[str], best_res: float) -> str:
    cell_type_html = ""
    if fig_ct is not None:
        cell_type_html = f"""
        <div class="fig-wrap wide">
          <h3>UMAP -- Ground-Truth Cell Type</h3>
          <img src="data:image/png;base64,{fig_ct}" alt="UMAP cell type">
        </div>"""
    return f"""
    <section>
      <h2>Figures</h2>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>UMAP -- Leiden Clustering across Resolutions</h3>
          <img src="data:image/png;base64,{fig_umap}" alt="UMAP resolutions">
        </div>
        <div class="fig-wrap wide">
          <h3>Silhouette Score &amp; Cluster Count vs Resolution</h3>
          <img src="data:image/png;base64,{fig_sil}" alt="Silhouette bar chart">
        </div>
        <div class="fig-wrap wide">
          <h3>Resolution Selection Diagnostics (4-panel)</h3>
          <img src="data:image/png;base64,{fig_diag}" alt="Resolution diagnostics">
        </div>
        <div class="fig-wrap wide">
          <h3>Cluster Size Distribution -- Best Resolution ({best_res:.2g})</h3>
          <img src="data:image/png;base64,{fig_sizes}" alt="Cluster sizes">
        </div>
        {cell_type_html}
      </div>
    </section>
    """


def _section_provenance(adata_clustered: AnnData) -> str:
    prov = adata_clustered.uns.get("omicsage_cluster", {})
    prov_display = {}
    for k, v in prov.items():
        if isinstance(v, dict):
            for subk, subv in v.items():
                display_subv = f"{subv:.4f}" if isinstance(subv, float) else subv
                prov_display[f"{k}[{subk}]"] = display_subv
        else:
            prov_display[k] = v
    rows = "".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov_display.items()
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

def run_cluster_report(
    adata_clustered: AnnData,
    metrics: dict,
    report_path: str = "reports/cluster_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the Leiden clustering HTML report and write it to disk.

    Parameters
    ----------
    adata_clustered : AnnData
        Clustered AnnData returned by cluster().
    metrics : dict
        Metrics dict returned by cluster().
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

    print(f"Building clustering report for '{dataset_name}' ...", flush=True)
    print("  Rendering UMAP resolution panels ...", flush=True)
    fig_umap  = _plot_umap_resolutions(adata_clustered, metrics)
    print("  Rendering silhouette bar chart ...", flush=True)
    fig_sil   = _plot_silhouette_bar(metrics)
    print("  Rendering resolution selection diagnostics ...", flush=True)
    fig_diag  = _plot_resolution_selection(metrics)
    print("  Rendering cluster size distribution ...", flush=True)
    fig_sizes = _plot_cluster_sizes(adata_clustered, metrics)
    fig_ct    = _plot_umap_cell_type(adata_clustered)
    if fig_ct is not None:
        print("  Rendering ground-truth cell type UMAP ...", flush=True)

    sections = [
        _section_summary(adata_clustered, metrics, dataset_name, timestamp),
        _section_figures(fig_umap, fig_sil, fig_diag, fig_sizes, fig_ct,
                         metrics["best_resolution"]),
        _section_provenance(adata_clustered),
    ]

    html = _render_page(
        title=f"OmicSage -- Cluster Report -- {dataset_name}",
        sections=sections,
        timestamp=timestamp,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
