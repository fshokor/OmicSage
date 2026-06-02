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


def _plot_embedding(
    adt: AnnData,
    embed_key: str,
    color_keys: list[str],
    title_prefix: str,
) -> str:
    """
    Generic multi-panel embedding plot.
    One panel per color_key — mimics sc.pl.umap(color=[...]).
    Categorical obs columns get tab20 colours + legend.
    Numeric obs columns get a viridis scatter + colorbar.
    Falls back gracefully if a key is missing from obs.
    """
    if embed_key not in adt.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, f"{embed_key} not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    coords = adt.obsm[embed_key]

    # Resolve valid keys — warn about missing ones
    valid_keys = []
    for k in color_keys:
        if k in adt.obs.columns:
            valid_keys.append(k)
        else:
            print(f"    [warn] color key '{k}' not in adt.obs — skipped", flush=True)

    # Always add a plain (no colour) panel as first panel
    all_keys = [None] + valid_keys   # None = plain scatter

    n_panels = len(all_keys)
    ncols    = min(3, n_panels)
    nrows    = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 5.5, nrows * 4.5),
                             squeeze=False)

    for idx, key in enumerate(all_keys):
        row, col = divmod(idx, ncols)
        ax = axes[row][col]

        if key is None:
            # Plain scatter
            ax.scatter(coords[:, 0], coords[:, 1],
                       s=2, alpha=0.5, color="#4C78A8", rasterized=True)
            ax.set_title(f"{title_prefix}", fontsize=10, fontweight="bold")

        else:
            vals = adt.obs[key]
            import pandas as pd
            is_numeric = pd.api.types.is_numeric_dtype(vals)

            if is_numeric:
                v = vals.values.astype(float)
                vmin, vmax = np.nanpercentile(v, [2, 98])
                sc = ax.scatter(coords[:, 0], coords[:, 1],
                                c=v, cmap="viridis",
                                vmin=vmin, vmax=vmax,
                                s=2, alpha=0.6, rasterized=True)
                plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
                ax.set_title(f"{title_prefix}\n{key}", fontsize=10, fontweight="bold")
            else:
                labels = vals.astype(str)
                unique = sorted(labels.unique())
                cmap   = plt.get_cmap("tab20", max(len(unique), 2))
                for i, lbl in enumerate(unique):
                    mask = (labels == lbl).values
                    ax.scatter(coords[mask, 0], coords[mask, 1],
                               s=2, alpha=0.7, color=cmap(i),
                               label=lbl, rasterized=True)
                ncol = max(1, len(unique) // 12)
                ax.legend(markerscale=4, frameon=False, fontsize=7,
                          loc="upper right", ncol=ncol,
                          title=key, title_fontsize=7)
                ax.set_title(f"{title_prefix}\n{key}", fontsize=10, fontweight="bold")

        ax.set_xlabel(f"{embed_key.split('_')[-1].upper()} 1", fontsize=8)
        ax.set_ylabel(f"{embed_key.split('_')[-1].upper()} 2", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

    # Hide unused axes
    for idx in range(n_panels, nrows * ncols):
        row, col = divmod(idx, ncols)
        axes[row][col].set_visible(False)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap(adt: AnnData, color_keys: Optional[list[str]] = None) -> str:
    """ADT UMAP with one panel per color key."""
    keys = color_keys or ["batch"]
    return _plot_embedding(adt, "X_umap_adt", keys, "ADT UMAP (pre-Harmony)")


def _plot_pca_multi(adt: AnnData, color_keys: Optional[list[str]] = None) -> str:
    """ADT PCA with one panel per color key."""
    keys = color_keys or ["batch"]
    return _plot_embedding(adt, "X_pca_adt", keys, "ADT PCA")


def _section_obs_columns(adt: AnnData) -> str:
    """
    Table of all obs columns with dtype and n_unique values.
    Use this to decide which columns to add to umap_color_keys in the config.
    """
    rows = ""
    for col in adt.obs.columns:
        dtype  = str(adt.obs[col].dtype)
        n_uniq = adt.obs[col].nunique()
        sample_vals = ", ".join(
            str(v) for v in adt.obs[col].astype(str).unique()[:5]
        )
        if n_uniq > 5:
            sample_vals += " ..."
        rows += (
            f"<tr><td><code>{col}</code></td>"
            f"<td>{dtype}</td>"
            f"<td>{n_uniq}</td>"
            f"<td style='font-size:0.82rem;color:#555;'>{sample_vals}</td></tr>"
        )

    return f"""
    <section>
      <h2>Available obs Columns</h2>
      <p>These are all columns available in <code>adt.obs</code>.
         Add any column name to <code>reduce_adt.params.umap_color_keys</code>
         in the config to colour the PCA and UMAP panels by that variable.</p>
      <table>
        <thead>
          <tr><th>Column</th><th>dtype</th><th>N unique</th><th>Sample values</th></tr>
        </thead>
        <tbody>{rows}</tbody>
      </table>
    </section>"""


def _section_protein_panel(adt: AnnData, metrics: dict) -> str:
    """
    Two-column layout: Kept | Removed.
    Run 1: all proteins in Kept, Removed column empty.
    Run 2: proteins split by isotype filtering.
    """
    all_proteins_before = metrics.get("n_vars_before_filter", adt.n_vars)
    removed = metrics.get("isotype_controls_removed", [])
    kept    = list(adt.var_names)

    def _li_list(items, color):
        if not items:
            return '<span style="color:#aaa;font-size:0.85rem;">—</span>'
        # 3 proteins per row using CSS grid
        cells = "".join(
            f'<div style="padding:3px 6px;font-size:0.83rem;'
            f'font-family:SFMono-Regular,Consolas,monospace;color:{color};">{p}</div>'
            for p in items
        )
        return (
            '<div style="display:grid;grid-template-columns:repeat(3,1fr);gap:2px;">'
            + cells +
            '</div>'
        )

    if not removed:
        callout = (
            '<div style="background:#e8f4fd;border-left:4px solid #3498db;'
            'padding:12px 16px;border-radius:0 6px 6px 0;'
            'margin-bottom:16px;font-size:0.88rem;color:#1a3a5c;">'
            '<strong>First run — no isotype controls removed.</strong> '
            'Identify isotype controls (typically IgG / Ctrl) in the Kept column, '
            'add them to <code>reduce_adt.params.isotype_controls</code> in the config, '
            'then re-run with <code>--step reduce_adt --force</code>.'
            '</div>'
        )
    else:
        removed_str = ", ".join(f"<code>{p}</code>" for p in removed)
        callout = (
            '<div style="background:#d4edda;border-left:4px solid #28a745;'
            'padding:12px 16px;border-radius:0 6px 6px 0;'
            'margin-bottom:16px;font-size:0.88rem;color:#155724;">'
            f'<strong>Isotype controls removed before PCA:</strong> {removed_str}<br>'
            f'{all_proteins_before} proteins → {adt.n_vars} proteins used for PCA.'
            '</div>'
        )

    kept_col    = _li_list(kept,    "#1a1a2e")
    removed_col = _li_list(removed, "#cc3333")

    return (
        '<section>'
        f'<h2>Protein Panel ({all_proteins_before} proteins, {adt.n_vars} used for PCA)</h2>'
        + callout +
        '<div style="display:grid;grid-template-columns:1fr 1fr;gap:24px;">'
        '<div>'
        f'<div style="font-weight:600;color:#0f3460;margin-bottom:8px;'
        f'border-bottom:2px solid #e8eaf6;padding-bottom:6px;">Kept ({len(kept)})</div>'
        + kept_col +
        '</div>'
        '<div>'
        f'<div style="font-weight:600;color:#cc3333;margin-bottom:8px;'
        f'border-bottom:2px solid #fde;padding-bottom:6px;">Removed ({len(removed)})</div>'
        + removed_col +
        '</div>'
        '</div>'
        '</section>'
    )


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


def _section_figures(
    fig_pca: str,
    fig_umap: str,
    active_layer: str = "adt_clr",
) -> str:
    pca_b64  = (
        '<div class="fig-wrap wide"><h3>ADT PCA</h3>'
        f'<img src="data:image/png;base64,{fig_pca}" alt="ADT PCA"></div>'
    )
    umap_b64 = (
        '<div class="fig-wrap wide"><h3>ADT UMAP (pre-Harmony)</h3>'
        f'<img src="data:image/png;base64,{fig_umap}" alt="ADT UMAP"></div>'
    )
    return f"""
    <section>
      <h2>Embeddings</h2>
      <p>PCA and UMAP computed from <code>{active_layer}</code>-normalised ADT values.
         Each panel shows one colour key. UMAP is pre-Harmony — batch effects
         may be visible. The post-Harmony UMAP is in the next report (cite_04).</p>
      <div class="fig-grid">
        {pca_b64}
        {umap_b64}
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
    active_layer = metrics.get("active_layer", "adt_clr")
    color_keys   = metrics.get("umap_color_keys") or ["batch"]

    print("  Rendering PCA multi-panel ...", flush=True)
    fig_pca  = _plot_pca_multi(adt, color_keys)
    print("  Rendering UMAP multi-panel ...", flush=True)
    fig_umap = _plot_umap(adt, color_keys)

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_obs_columns(adt),
            _section_protein_panel(adt, metrics),
            _section_figures(fig_pca, fig_umap, active_layer),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
