"""
OmicSage — CITE-seq ADT Dimensionality Reduction Report
reports/templates/cite/cite_reduce_report.py

Generated after reduce_adt step.
Output: cite_03_reduce_report.html
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
import numpy as np
from anndata import AnnData

_DPI = 130
_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px; box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
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
    td { padding: 8px 12px; border-bottom: 1px solid #eee; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9ff; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 520px; }
    .fig-wrap.wide { flex: 1 1 100%; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
"""


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_pca_scatter(adt: AnnData) -> str:
    """ADT PCA — PC1 vs PC2, coloured by doublet score if available."""
    pca_key = "X_pca_adt"
    if pca_key not in adt.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "X_pca_adt not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    pca = adt.obsm[pca_key]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Left: plain scatter
    axes[0].scatter(pca[:, 0], pca[:, 1], s=2, alpha=0.5,
                    color="#4C78A8", rasterized=True)
    axes[0].set_title("ADT PCA — PC1 vs PC2", fontsize=11, fontweight="bold")
    axes[0].set_xlabel("PC 1", fontsize=9); axes[0].set_ylabel("PC 2", fontsize=9)
    axes[0].spines[["top", "right"]].set_visible(False)

    # Right: coloured by batch if available
    batch_col = next((c for c in ["batch", "donor", "sample"] if c in adt.obs.columns), None)
    if batch_col:
        labels = adt.obs[batch_col].astype(str)
        unique = sorted(labels.unique())
        cmap   = plt.get_cmap("tab20", max(len(unique), 2))
        for i, b in enumerate(unique):
            mask = labels == b
            axes[1].scatter(pca[mask, 0], pca[mask, 1], s=2, alpha=0.6,
                            color=cmap(i), label=b, rasterized=True)
        axes[1].legend(markerscale=4, frameon=False, fontsize=8,
                       ncol=max(1, len(unique) // 10))
        axes[1].set_title(f"ADT PCA — coloured by {batch_col}", fontsize=11, fontweight="bold")
    else:
        axes[1].scatter(pca[:, 0], pca[:, 1], s=2, alpha=0.5,
                        color="#e07b3a", rasterized=True)
        axes[1].set_title("ADT PCA — PC1 vs PC2", fontsize=11, fontweight="bold")
    axes[1].set_xlabel("PC 1", fontsize=9); axes[1].set_ylabel("PC 2", fontsize=9)
    axes[1].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap(adt: AnnData) -> str:
    """ADT UMAP coloured by batch."""
    umap_key = "X_umap_adt"
    if umap_key not in adt.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "X_umap_adt not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    umap = adt.obsm[umap_key]
    batch_col = next((c for c in ["batch", "donor", "sample"] if c in adt.obs.columns), None)

    fig, ax = plt.subplots(figsize=(7, 5))
    if batch_col:
        labels = adt.obs[batch_col].astype(str)
        unique = sorted(labels.unique())
        cmap   = plt.get_cmap("tab20", max(len(unique), 2))
        for i, b in enumerate(unique):
            mask = labels == b
            ax.scatter(umap[mask, 0], umap[mask, 1], s=2, alpha=0.6,
                       color=cmap(i), label=b, rasterized=True)
        ax.legend(markerscale=4, frameon=False, fontsize=8,
                  ncol=max(1, len(unique) // 10))
        ax.set_title(f"ADT UMAP (pre-Harmony) — coloured by {batch_col}",
                     fontsize=12, fontweight="bold")
    else:
        ax.scatter(umap[:, 0], umap[:, 1], s=2, alpha=0.5,
                   color="#4C78A8", rasterized=True)
        ax.set_title("ADT UMAP (pre-Harmony)", fontsize=12, fontweight="bold")

    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq ADT Reduction Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq ADT Dimensionality Reduction Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    var_total = metrics.get("variance_explained_total")
    var_str   = f"{var_total * 100:.1f}%" if var_total else "?"
    cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",          f"{metrics.get('n_cells', '?'):,}"),
            ("Proteins",       str(metrics.get("n_vars", "?"))),
            ("PCs computed",   str(metrics.get("n_comps_actual", "?"))),
            ("PCs for graph",  str(metrics.get("n_pcs_used", "?"))),
            ("Variance expl.", var_str),
            ("Neighbors (k)",  str(metrics.get("n_neighbors", "?"))),
        ]
    )
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("PCA key",  metrics.get("pca_key", "X_pca_adt")),
            ("UMAP key", metrics.get("umap_key", "X_umap_adt")),
            ("UMAP computed", "Yes" if metrics.get("umap_computed") else "No"),
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


def _section_figures(fig_pca: str, fig_umap: str) -> str:
    return f"""
    <section>
      <h2>Embeddings</h2>
      <p>PCA and UMAP computed from CLR-normalised ADT values.
         UMAP shown here is pre-Harmony — batch effects may be visible.
         The post-Harmony UMAP is in the next report (cite_04).</p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>ADT PCA — PC1 vs PC2</h3>
          <img src="data:image/png;base64,{fig_pca}" alt="ADT PCA">
        </div>
        <div class="fig-wrap">
          <h3>ADT UMAP (pre-Harmony)</h3>
          <img src="data:image/png;base64,{fig_umap}" alt="ADT UMAP">
        </div>
      </div>
    </section>"""


def run_cite_reduce_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_03_reduce_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building CITE reduce report for '{dataset_name}' ...", flush=True)
    fig_pca  = _plot_pca_scatter(adt)
    fig_umap = _plot_umap(adt)

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_figures(fig_pca, fig_umap),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
