"""
OmicSage — ATAC Dimensionality Reduction Report
reports/templates/multiome/atac_reduce_report.py

Generated after atac_reduce step.
Output: multiome_02_atac_reduce_report.html

Usage
-----
    from reports.templates.multiome.atac_reduce_report import run_atac_reduce_report
    run_atac_reduce_report(
        atac=atac_reduced,
        metrics=metrics,
        report_path="reports/GSE194122_multiome/multiome_02_atac_reduce_report.html",
        dataset_name="BMMC Multiome (NeurIPS 2021)",
        batch_key="batch",
    )
"""

from __future__ import annotations

import base64
import io
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def _plot_variance_explained(metrics: dict) -> str:
    """
    Scree-style bar chart of SVD component variance ratios.
    Component 1 (dropped) shown in red; components 2–N (used) in blue.
    Mirrors the Signac DepthCor / variance plot pattern.
    """
    var_all = np.array(metrics.get("variance_ratio_all", []))
    n_used  = metrics.get("n_lsi_components_used", len(var_all) - 1)

    if len(var_all) == 0:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "Variance data not available", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    n_total = len(var_all)
    x = np.arange(1, n_total + 1)
    # Component 1 dropped (red), components 2–N used (blue), rest grey
    colors = ["#e74c3c"] + ["#4C78A8"] * n_used + ["#aaaaaa"] * max(0, n_total - n_used - 1)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # --- Bar chart ---
    ax = axes[0]
    ax.bar(x, var_all * 100, color=colors[:n_total], width=0.8, alpha=0.85)
    ax.axvline(1.5, color="#e74c3c", linewidth=1.2, linestyle="--",
               label="Component 1 dropped (depth correlation)")
    ax.set_xlabel("SVD Component", fontsize=11)
    ax.set_ylabel("Variance explained (%)", fontsize=11)
    ax.set_title("LSI Variance per Component", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))

    # --- Cumulative (components 2–N only) ---
    ax = axes[1]
    var_used = var_all[1:n_used + 1]
    cumvar = np.cumsum(var_used)
    x_used = np.arange(2, len(cumvar) + 2)
    ax.plot(x_used, cumvar * 100, color="#4C78A8", linewidth=2)
    ax.set_xlabel("LSI Component (2–N)", fontsize=11)
    ax.set_ylabel("Cumulative variance (%)", fontsize=11)
    ax.set_title("Cumulative Variance (components 2–N)", fontsize=13, fontweight="bold")
    cum_total = float(var_used.sum()) * 100
    ax.annotate(
        f"{cum_total:.1f}% total",
        xy=(x_used[-1], cumvar[-1] * 100),
        xytext=(max(2, x_used[-1] - len(x_used) // 4), cumvar[-1] * 100 - 5),
        fontsize=8, color="#4C78A8",
        arrowprops=dict(arrowstyle="->", color="#4C78A8", lw=0.8),
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_by_key(atac: AnnData, key: Optional[str], title: Optional[str] = None) -> Optional[str]:
    """UMAP coloured by any categorical or continuous obs column."""
    umap_key = "X_umap_atac"
    if umap_key not in atac.obsm:
        return None
    if key is None or key not in atac.obs.columns:
        return None

    umap   = atac.obsm[umap_key]
    labels = atac.obs[key]
    fig, ax = plt.subplots(figsize=(6, 5))

    if labels.dtype.kind in ("O", "S", "U") or str(labels.dtype) == "category":
        unique = sorted(labels.astype(str).unique())
        cmap   = plt.get_cmap("tab20", max(len(unique), 2))
        for i, lbl in enumerate(unique):
            mask = labels.astype(str) == lbl
            ax.scatter(umap[mask, 0], umap[mask, 1], s=2, alpha=0.6,
                       color=cmap(i), label=lbl, rasterized=True)
        ax.legend(markerscale=4, frameon=False, fontsize=7,
                  loc="upper right", ncol=max(1, len(unique) // 12))
    else:
        sc = ax.scatter(umap[:, 0], umap[:, 1],
                        c=labels.values.astype(float),
                        cmap="viridis", s=2, alpha=0.6, rasterized=True)
        plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)

    ax.set_title(title or f"UMAP — {key}", fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_lsi_scatter(atac: AnnData) -> str:
    """
    Scatter of LSI component 2 vs component 3, coloured by total_peak_counts.
    Shows that components 2+ capture biology, not depth.
    """
    if "X_lsi" not in atac.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "X_lsi not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    lsi = atac.obsm["X_lsi"]
    if lsi.shape[1] < 2:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "< 2 LSI components", ha="center", va="center",
                transform=ax.transAxes)
        ax.axis("off")
        return _fig_to_b64(fig)

    color_col = "total_peak_counts"
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    if color_col in atac.obs.columns:
        sc_plot = ax.scatter(lsi[:, 0], lsi[:, 1],
                             c=atac.obs[color_col].values.astype(float),
                             cmap="YlOrRd", s=2, alpha=0.6, rasterized=True)
        plt.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02, label="Total peak counts")
    else:
        ax.scatter(lsi[:, 0], lsi[:, 1], s=2, alpha=0.4, color="#4C78A8", rasterized=True)

    ax.set_xlabel("LSI 2", fontsize=10)
    ax.set_ylabel("LSI 3", fontsize=10)
    ax.set_title(
        "LSI Components 2 vs 3\n(coloured by total peak counts — should show no depth gradient)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_cluster_sizes(atac: AnnData) -> str:
    """Bar chart of Leiden cluster sizes."""
    col = "atac_leiden"
    if col not in atac.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "atac_leiden not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    counts = atac.obs[col].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(max(5, len(counts) * 0.4), 4))
    ax.bar(counts.index.astype(str), counts.values, color="#4C78A8", alpha=0.8, width=0.7)
    ax.set_xlabel("Leiden cluster", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title(
        f"ATAC Leiden Cluster Sizes ({len(counts)} clusters)",
        fontsize=12, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — ATAC Reduction Report — {dataset_name}</title>
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
    .callout {{ background: #e8f4fd; border-left: 4px solid #3498db;
                padding: 12px 16px; border-radius: 0 6px 6px 0;
                margin: 12px 0; font-size: 0.88rem; color: #1a3a5c; }}
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage — ATAC Dimensionality Reduction Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    cum_var_pct = metrics.get("cumulative_variance_used", 0) * 100
    cards = [
        ("Cells",              f"{metrics.get('n_cells', '?'):,}"),
        ("Peaks",              f"{metrics.get('n_peaks', '?'):,}"),
        ("SVD components",     str(metrics.get("n_components_computed", "?"))),
        ("LSI components used",str(metrics.get("n_lsi_components_used", "?"))),
        ("Cumul. variance",    f"{cum_var_pct:.1f}%"),
        ("Neighbors (k)",      str(metrics.get("n_neighbors", "?"))),
        ("Leiden clusters",    str(metrics.get("n_leiden_clusters", "?"))),
    ]
    stat_html = "".join(
        f'<div class="stat-card"><div class="stat-value">{val}</div>'
        f'<div class="stat-label">{lbl}</div></div>'
        for lbl, val in cards
    )
    rows = [
        ("TF-IDF method",        "log(TF) × log(1 + N/df)"),
        ("LSI method",           "sklearn TruncatedSVD"),
        ("Component 1 dropped",  "Yes — correlates with sequencing depth"),
        ("Post-SVD normalisation","L2 row normalisation"),
        ("Embeddings computed",  ", ".join(metrics.get("embeddings_computed", []))),
        ("Cluster key",          f"<code>{metrics.get('cluster_key', 'atac_leiden')}</code>"),
        ("Leiden resolution",    str(metrics.get("leiden_resolution", "?"))),
    ]
    rows_html = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="callout">
        LSI component 1 was dropped unconditionally — it correlates with sequencing
        depth rather than biology (Signac, sc-best-practices ATAC chapter, Seurat v5
        ATAC vignette all confirm this). Components 2–{metrics.get('n_components_computed','N')}
        are used for the neighbor graph, UMAP, and clustering.
      </div>
      <div class="stat-grid">{stat_html}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows_html}</tbody></table>
    </section>"""


def _section_variance(fig_var: str) -> str:
    return f"""
    <section>
      <h2>LSI Variance Explained</h2>
      <p>Red bar = component 1 (dropped). Blue bars = components used.
         Right panel shows cumulative variance for components 2–N.</p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_var}" alt="LSI variance">
        </div>
      </div>
    </section>"""


def _section_embeddings(
    fig_lsi: str,
    fig_umap_leiden: Optional[str],
    fig_umap_batch: Optional[str],
    fig_umap_gt: Optional[str],
    fig_clusters: str,
    batch_key: Optional[str],
) -> str:
    umap_leiden_html = ""
    if fig_umap_leiden:
        umap_leiden_html = f"""
        <div class="fig-wrap">
          <h3>UMAP — Leiden clusters (atac_leiden)</h3>
          <img src="data:image/png;base64,{fig_umap_leiden}" alt="UMAP Leiden">
        </div>"""
    umap_batch_html = ""
    if fig_umap_batch:
        umap_batch_html = f"""
        <div class="fig-wrap">
          <h3>UMAP — coloured by {batch_key}</h3>
          <img src="data:image/png;base64,{fig_umap_batch}" alt="UMAP batch">
        </div>"""
    umap_gt_html = ""
    if fig_umap_gt:
        umap_gt_html = f"""
        <div class="fig-wrap">
          <h3>UMAP — ground truth cell types</h3>
          <img src="data:image/png;base64,{fig_umap_gt}" alt="UMAP ground truth">
        </div>"""
    return f"""
    <section>
      <h2>Embeddings</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>LSI 2 vs 3 (coloured by total peak counts)</h3>
          <img src="data:image/png;base64,{fig_lsi}" alt="LSI scatter">
        </div>
        {umap_leiden_html}
        {umap_batch_html}
        {umap_gt_html}
        <div class="fig-wrap">
          <h3>Leiden Cluster Sizes</h3>
          <img src="data:image/png;base64,{fig_clusters}" alt="cluster sizes">
        </div>
      </div>
    </section>"""


def _section_provenance(atac: AnnData) -> str:
    prov = atac.uns.get("omicsage_atac_reduce", {})
    rows = "".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov.items()
        if k != "variance_ratio_used"   # skip long list
    )
    return f"""
    <section>
      <h2>Provenance</h2>
      <table><thead><tr><th>Key</th><th>Value</th></tr></thead>
      <tbody>{rows}</tbody></table>
    </section>"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_atac_reduce_report(
    atac: AnnData,
    metrics: dict,
    report_path: str = "reports/multiome_02_atac_reduce_report.html",
    dataset_name: str = "dataset",
    batch_key: Optional[str] = "batch",
) -> str:
    """Generate a self-contained HTML ATAC dimensionality reduction report.

    Parameters
    ----------
    atac : AnnData
        Output of atac_reduce(), with X_lsi, X_umap_atac, atac_leiden in obs.
    metrics : dict
        Metrics dict returned by atac_reduce().
    report_path : str
        Output path for the HTML file.
    dataset_name : str
        Dataset label shown in the report header.
    batch_key : str, optional
        obs column for the batch UMAP panel. Default "batch".

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building ATAC reduce report for '{dataset_name}' ...", flush=True)

    print("  Rendering variance explained ...", flush=True)
    fig_var = _plot_variance_explained(metrics)

    print("  Rendering LSI scatter ...", flush=True)
    fig_lsi = _plot_lsi_scatter(atac)

    print("  Rendering UMAP (Leiden) ...", flush=True)
    fig_umap_leiden = _plot_umap_by_key(atac, "atac_leiden",
                                         "UMAP — ATAC Leiden clusters")

    print(f"  Rendering UMAP ({batch_key}) ...", flush=True)
    fig_umap_batch = _plot_umap_by_key(atac, batch_key,
                                        f"UMAP — {batch_key}") if batch_key else None

    print("  Rendering UMAP (ground truth) ...", flush=True)
    fig_umap_gt = _plot_umap_by_key(atac, "cell_type_groundtruth",
                                     "UMAP — ground truth cell types")

    print("  Rendering cluster sizes ...", flush=True)
    fig_clusters = _plot_cluster_sizes(atac)

    sections = [
        _section_summary(metrics, dataset_name, timestamp),
        _section_variance(fig_var),
        _section_embeddings(
            fig_lsi, fig_umap_leiden, fig_umap_batch, fig_umap_gt,
            fig_clusters, batch_key,
        ),
        _section_provenance(atac),
    ]

    html = _render_page(sections=sections, timestamp=timestamp, dataset_name=dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
