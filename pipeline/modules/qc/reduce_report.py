"""
OmicSage — Dimensionality Reduction Report
pipeline/modules/qc/reduce_report.py

Usage
-----
    from pipeline.modules.qc.reduce_report import run_reduce_report
    run_reduce_report(
        adata_reduced=adata_reduced,
        metrics=metrics,
        report_path="reports/GSE194122/03_reduce_report.html",
        dataset_name="GSE194122_CITE",
        batch_key="batch",
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


def _plot_scree(adata_reduced: AnnData, metrics: dict) -> str:
    variance_ratio    = np.array(metrics["variance_explained_per_pc"])
    cumvar            = np.cumsum(variance_ratio)
    n_pcs_used        = metrics["n_pcs_used"]
    pc_selection_method = metrics.get("pc_selection_method", "")
    n_comps           = len(variance_ratio)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    x = np.arange(1, n_comps + 1)

    ax = axes[0]
    colors = ["#e07b3a" if i < n_pcs_used else "#aaaaaa" for i in range(n_comps)]
    ax.bar(x, variance_ratio * 100, color=colors, width=0.8, alpha=0.85)
    ax.axvline(n_pcs_used + 0.5, color="#c0392b", linewidth=1.5, linestyle="--",
               label=f"Selected: {n_pcs_used} PCs ({pc_selection_method})")
    ax.set_xlabel("Principal Component", fontsize=11)
    ax.set_ylabel("Variance explained (%)", fontsize=11)
    ax.set_title("Scree Plot", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))

    ax = axes[1]
    ax.plot(x, cumvar * 100, color="#4C78A8", linewidth=2)
    ax.axvline(n_pcs_used, color="#c0392b", linewidth=1.5, linestyle="--",
               label=f"Selected: {n_pcs_used} PCs")
    ax.axhline(float(cumvar[n_pcs_used - 1]) * 100, color="#c0392b",
               linewidth=0.8, linestyle=":")
    ax.set_xlabel("Number of PCs", fontsize=11)
    ax.set_ylabel("Cumulative variance explained (%)", fontsize=11)
    ax.set_title("Cumulative Variance", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))
    cum_at_selected = float(cumvar[n_pcs_used - 1]) * 100
    ax.annotate(f"{cum_at_selected:.1f}%",
                xy=(n_pcs_used, cum_at_selected),
                xytext=(n_pcs_used + max(1, n_comps // 10), cum_at_selected - 5),
                fontsize=8, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.8))

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_qc(adata_reduced: AnnData) -> str:
    if "X_umap" not in adata_reduced.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "UMAP not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    umap = adata_reduced.obsm["X_umap"]
    obs  = adata_reduced.obs
    qc_cols = [
        ("n_genes_by_counts", "Genes per cell",   "YlOrRd"),
        ("total_counts",      "Total UMI counts", "Blues"),
        ("pct_counts_mt",     "MT% per cell",     "Reds"),
        ("doublet_score",     "Doublet score",    "Purples"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for ax, (col, title, cmap) in zip(axes.flatten(), qc_cols):
        if col in obs.columns:
            sc_plot = ax.scatter(umap[:, 0], umap[:, 1],
                                 c=obs[col].values.astype(float),
                                 cmap=cmap, s=2, alpha=0.6, rasterized=True)
            plt.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02)
        else:
            ax.scatter(umap[:, 0], umap[:, 1], s=2, alpha=0.4,
                       color="#aaaaaa", rasterized=True)
            ax.text(0.5, 0.02, f"'{col}' not in obs", ha="center",
                    transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.suptitle("UMAP -- QC Metric Overlays", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_pca_qc(adata_reduced: AnnData) -> str:
    if "X_pca" not in adata_reduced.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "PCA not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    pca = adata_reduced.obsm["X_pca"]
    obs = adata_reduced.obs
    qc_cols = [
        ("n_genes_by_counts", "Genes per cell",   "YlOrRd"),
        ("total_counts",      "Total UMI counts", "Blues"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, (col, title, cmap) in zip(axes, qc_cols):
        if col in obs.columns:
            sc_plot = ax.scatter(pca[:, 0], pca[:, 1],
                                 c=obs[col].values.astype(float),
                                 cmap=cmap, s=3, alpha=0.6, rasterized=True)
            plt.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02)
        else:
            ax.scatter(pca[:, 0], pca[:, 1], s=3, alpha=0.4,
                       color="#aaaaaa", rasterized=True)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("PC 1", fontsize=10); ax.set_ylabel("PC 2", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("PCA -- QC Metric Overlays (PC1 vs PC2)", fontsize=13, fontweight="bold")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_batch(adata_reduced: AnnData, batch_key: Optional[str]) -> Optional[str]:
    if batch_key is None or "X_umap" not in adata_reduced.obsm:
        return None
    if batch_key not in adata_reduced.obs.columns:
        logger.warning("batch_key='%s' not in obs -- skipping batch UMAP", batch_key)
        return None
    umap    = adata_reduced.obsm["X_umap"]
    batches = adata_reduced.obs[batch_key].astype(str)
    unique  = sorted(batches.unique())
    cmap    = plt.get_cmap("tab20", len(unique))
    cmap_   = {b: cmap(i) for i, b in enumerate(unique)}
    fig, ax = plt.subplots(figsize=(7, 5))
    for b in unique:
        mask = batches == b
        ax.scatter(umap[mask, 0], umap[mask, 1], s=2, alpha=0.6,
                   color=cmap_[b], label=b, rasterized=True)
    ax.legend(markerscale=4, frameon=False, fontsize=8,
              loc="upper right", ncol=max(1, len(unique) // 10))
    ax.set_title(f"UMAP -- coloured by {batch_key}", fontsize=13, fontweight="bold")
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
    <h1>OmicSage -- Dimensionality Reduction Report</h1>
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
    cum_var_pct = metrics.get("cumulative_variance_explained", 0) * 100
    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",          f"{metrics.get('n_cells', '?'):,}" if isinstance(metrics.get('n_cells'), int) else "?"),
            ("HVGs used",      str(metrics.get("n_hvg_used", "?"))),
            ("PCs computed",   str(metrics.get("n_comps_computed", "?"))),
            ("PCs used",       str(metrics.get("n_pcs_used", "?"))),
            ("Cumul. var",     f"{cum_var_pct:.1f}%"),
            ("Neighbors (k)",  str(metrics.get("n_neighbors", "?"))),
        ]
    )
    param_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("PC selection method", metrics.get("pc_selection_method", "?")),
            ("Embeddings computed", ", ".join(metrics.get("embeddings_computed", []))),
            ("t-SNE computed",      "Yes" if metrics.get("run_tsne") else "No"),
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


def _section_scree(fig_scree: str) -> str:
    return f"""
    <section>
      <h2>Scree Plot</h2>
      <p>Orange bars = PCs selected for the neighbor graph. Grey = computed but not used.
         Red dashed line = selected cut-off.</p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_scree}" alt="Scree plot">
        </div>
      </div>
    </section>
    """


def _section_umap_qc(fig_umap: str, fig_pca: str, fig_batch: Optional[str],
                     batch_key: Optional[str]) -> str:
    batch_html = ""
    if fig_batch is not None:
        batch_html = f"""
        <div class="fig-wrap wide">
          <h3>UMAP -- coloured by {batch_key}</h3>
          <img src="data:image/png;base64,{fig_batch}" alt="UMAP batch">
        </div>"""
    return f"""
    <section>
      <h2>Embeddings</h2>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>UMAP -- QC Metric Overlays</h3>
          <img src="data:image/png;base64,{fig_umap}" alt="UMAP QC">
        </div>
        <div class="fig-wrap wide">
          <h3>PCA -- QC Metric Overlays (PC1 vs PC2)</h3>
          <img src="data:image/png;base64,{fig_pca}" alt="PCA QC">
        </div>
        {batch_html}
      </div>
    </section>
    """


def _section_provenance(adata_reduced: AnnData) -> str:
    prov = adata_reduced.uns.get("omicsage_reduce", {})
    prov_display = {
        k: (f"[{v[0]:.4f}, {v[1]:.4f}, ... ] ({len(v)} values)"
            if k == "variance_explained_per_pc" and isinstance(v, list) else v)
        for k, v in prov.items()
    }
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

def run_reduce_report(
    adata_reduced: AnnData,
    metrics: dict,
    report_path: str = "reports/reduce_report.html",
    dataset_name: str = "dataset",
    batch_key: Optional[str] = None,
) -> str:
    """
    Generate the dimensionality reduction HTML report and write it to disk.

    Parameters
    ----------
    adata_reduced : AnnData
        Reduced AnnData returned by reduce(). Must have obsm['X_pca'] and obsm['X_umap'].
    metrics : dict
        Metrics dict returned by reduce().
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Label shown in the report header.
    batch_key : str, optional
        obs column for the batch UMAP panel. Omitted if None.

    Returns
    -------
    str
        Absolute path to the written report file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building dimensionality reduction report for '{dataset_name}' ...", flush=True)
    print("  Rendering scree plot ...", flush=True)
    fig_scree = _plot_scree(adata_reduced, metrics)
    print("  Rendering UMAP QC overlays ...", flush=True)
    fig_umap  = _plot_umap_qc(adata_reduced)
    print("  Rendering PCA QC overlays ...", flush=True)
    fig_pca   = _plot_pca_qc(adata_reduced)
    fig_batch = _plot_umap_batch(adata_reduced, batch_key)
    if fig_batch is not None:
        print("  Rendering batch UMAP ...", flush=True)

    sections = [
        _section_summary(metrics, dataset_name, timestamp),
        _section_scree(fig_scree),
        _section_umap_qc(fig_umap, fig_pca, fig_batch, batch_key),
        _section_provenance(adata_reduced),
    ]

    html = _render_page(
        title=f"OmicSage -- Reduce Report -- {dataset_name}",
        sections=sections,
        timestamp=timestamp,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
