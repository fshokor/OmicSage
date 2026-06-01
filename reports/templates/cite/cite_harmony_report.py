"""
OmicSage — CITE-seq ADT Harmony Batch Correction Report
reports/templates/cite/cite_harmony_report.py

Generated after harmony_adt step.
Output: cite_04_harmony_report.html
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


def _plot_umap_by_batch(adt: AnnData, batch_key: str, title: str,
                         pca_key: str = "X_umap_adt") -> str:
    if pca_key not in adt.obsm or batch_key not in adt.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, f"{pca_key} or {batch_key} not found",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    coords = adt.obsm[pca_key]
    labels = adt.obs[batch_key].astype(str)
    unique = sorted(labels.unique())
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, b in enumerate(unique):
        mask = labels == b
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6,
                   color=cmap(i), label=b, rasterized=True)
    ax.legend(markerscale=4, frameon=False, fontsize=8,
              ncol=max(1, len(unique) // 10))
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1" if "umap" in pca_key.lower() else "PC 1", fontsize=9)
    ax.set_ylabel("UMAP 2" if "umap" in pca_key.lower() else "PC 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_pca_before_after(adt: AnnData, batch_key: str) -> str:
    """Side-by-side PCA: uncorrected (X_pca_adt) vs harmony (X_pca_harmony_adt)."""
    has_pre  = "X_pca_adt" in adt.obsm
    has_post = "X_pca_harmony_adt" in adt.obsm
    has_batch = batch_key in adt.obs.columns

    if not (has_pre or has_post):
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "PCA embeddings not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    titles = ["ADT PCA (uncorrected)", "ADT PCA (Harmony-corrected)"]
    keys   = ["X_pca_adt", "X_pca_harmony_adt"]

    for ax, key, title in zip(axes, keys, titles):
        if key not in adt.obsm:
            ax.text(0.5, 0.5, f"{key} not found", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="#888")
            ax.axis("off"); continue
        coords = adt.obsm[key]
        if has_batch:
            labels = adt.obs[batch_key].astype(str)
            unique = sorted(labels.unique())
            cmap   = plt.get_cmap("tab20", max(len(unique), 2))
            for i, b in enumerate(unique):
                mask = labels == b
                ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6,
                           color=cmap(i), label=b, rasterized=True)
            ax.legend(markerscale=4, frameon=False, fontsize=7,
                      ncol=max(1, len(unique) // 10))
        else:
            ax.scatter(coords[:, 0], coords[:, 1], s=2, alpha=0.5,
                       color="#4C78A8", rasterized=True)
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("PC 1", fontsize=9); ax.set_ylabel("PC 2", fontsize=9)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(f"PCA Before / After Harmony  (coloured by {batch_key})",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq ADT Harmony Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq ADT Harmony Batch Correction Report</h1>
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
            ("Cells",        f"{metrics.get('n_cells', '?'):,}"),
            ("Batches",      str(metrics.get("n_batches", "?"))),
            ("Batch key",    metrics.get("batch_key", "?")),
            ("PCs used",     str(metrics.get("n_pcs_used", "?"))),
            ("Neighbors",    str(metrics.get("n_neighbors", "?"))),
            ("UMAP computed","Yes" if metrics.get("umap_computed") else "No"),
        ]
    )
    batch_vals = metrics.get("batch_values", [])
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Harmony embedding key", metrics.get("harmony_key", "X_pca_harmony_adt")),
            ("UMAP key",              metrics.get("umap_key", "X_umap_adt")),
            ("Batch values",          ", ".join(str(b) for b in batch_vals)),
            ("Random state",          str(metrics.get("random_state", 0))),
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


def _section_figures(fig_pca_compare: str, fig_umap: str, batch_key: str) -> str:
    return f"""
    <section>
      <h2>Batch Correction Results</h2>
      <p>Harmony batch correction applied to ADT PCA embedding (<code>X_pca_adt</code>).
         Corrected embedding stored in <code>X_pca_harmony_adt</code>.
         UMAP recomputed from the corrected embedding (<code>X_umap_adt</code>).</p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>PCA Before / After Harmony (coloured by {batch_key})</h3>
          <img src="data:image/png;base64,{fig_pca_compare}" alt="PCA before after">
        </div>
        <div class="fig-wrap">
          <h3>ADT UMAP (post-Harmony, coloured by {batch_key})</h3>
          <img src="data:image/png;base64,{fig_umap}" alt="ADT UMAP post-Harmony">
        </div>
      </div>
    </section>"""


def run_cite_harmony_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_04_harmony_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    batch_key = metrics.get("batch_key", "batch")

    print(f"Building CITE harmony report for '{dataset_name}' ...", flush=True)
    fig_pca_compare = _plot_pca_before_after(adt, batch_key)
    fig_umap        = _plot_umap_by_batch(adt, batch_key,
                                          f"ADT UMAP (post-Harmony) — {batch_key}",
                                          pca_key="X_umap_adt")

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_figures(fig_pca_compare, fig_umap, batch_key),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
