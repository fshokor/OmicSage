"""
spatial_cluster_report.py — OmicSage Phase 7, Session 3
HTML report for spatial clustering and spatially variable gene results.

Generates a self-contained HTML report with:
  - Cluster summary: n_clusters, sizes, resolution used
  - Spatial scatter: spots coloured by cluster label
  - UMAP scatter: spots coloured by cluster label (if UMAP available)
  - Top spatially variable genes: bar chart of Moran's I scores
  - Spatial scatter of top 3 SVGs
  - Summary table of top 20 SVGs
"""

from __future__ import annotations

import base64
import io
import os
from datetime import datetime
from typing import Optional

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_spatial_cluster_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
) -> str:
    """Generate a self-contained HTML clustering report for Visium data.

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_cluster` (contains
        ``uns["omicsage_spatial_cluster"]``).
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label used in the report title.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    if "omicsage_spatial_cluster" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_cluster'] not found. "
            "Run spatial_cluster() before generating the report."
        )

    cluster_info = adata.uns["omicsage_spatial_cluster"]
    figures = _build_figures(adata, cluster_info)
    html = _render_html(adata, cluster_info, figures, dataset_id)

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


def _build_figures(adata: ad.AnnData, cluster_info: dict) -> dict[str, Optional[str]]:
    figures: dict[str, Optional[str]] = {}
    cluster_key = cluster_info["params"]["cluster_key"]

    figures["spatial_clusters"] = _spatial_cluster_scatter(adata, cluster_key)
    figures["cluster_sizes"] = _cluster_size_bar(cluster_info)
    figures["svg_bar"] = _svg_bar(adata)
    figures["svg_spatial"] = _svg_spatial_scatter(adata)

    return figures


def _spatial_cluster_scatter(
    adata: ad.AnnData, cluster_key: str
) -> Optional[str]:
    """Spatial scatter coloured by Leiden cluster."""
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    if cluster_key not in adata.obs.columns:
        return None
    try:
        fig, ax = plt.subplots(figsize=(6, 6))
        sq.pl.spatial_scatter(
            adata,
            color=cluster_key,
            ax=ax,
            show=False,
            frameon=False,
        )
        ax.set_title("Leiden clusters (spatial)")
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _cluster_size_bar(cluster_info: dict) -> Optional[str]:
    """Horizontal bar chart of spots per cluster."""
    sizes = cluster_info["outputs"].get("cluster_sizes", {})
    if not sizes:
        return None
    try:
        # Sort by cluster id numerically where possible
        try:
            labels = sorted(sizes.keys(), key=lambda x: int(x))
        except ValueError:
            labels = sorted(sizes.keys())
        values = [sizes[k] for k in labels]

        fig, ax = plt.subplots(figsize=(5, max(3, len(labels) * 0.35)))
        colors = plt.cm.tab20(np.linspace(0, 1, len(labels)))
        ax.barh(labels, values, color=colors, edgecolor="white")
        ax.set_xlabel("Number of spots")
        ax.set_ylabel("Cluster")
        ax.set_title("Spots per cluster")
        ax.invert_yaxis()
        for i, (label, val) in enumerate(zip(labels, values)):
            ax.text(val + max(values) * 0.01, i, str(val), va="center", fontsize=8)
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _svg_bar(adata: ad.AnnData) -> Optional[str]:
    """Bar chart of Moran's I score for top 20 SVGs."""
    if "moranI" not in adata.uns:
        return None
    try:
        moran_df: pd.DataFrame = adata.uns["moranI"]
        top = moran_df.head(20).copy()
        genes = top.index.tolist()
        scores = top["I"].values
        sig = top["pval_norm_fdr_bh"].values < 0.05

        colors = ["#e74c3c" if s else "#95a5a6" for s in sig]

        fig, ax = plt.subplots(figsize=(7, max(3, len(genes) * 0.35)))
        ax.barh(genes, scores, color=colors, edgecolor="white")
        ax.set_xlabel("Moran's I")
        ax.set_title("Top spatially variable genes\n(red = FDR < 0.05)")
        ax.invert_yaxis()
        ax.set_xlim(0, max(scores) * 1.15)
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _svg_spatial_scatter(adata: ad.AnnData) -> Optional[str]:
    """Spatial scatter of top 3 SVGs side by side."""
    if not _SQUIDPY_AVAILABLE:
        return None
    if "moranI" not in adata.uns:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    try:
        moran_df: pd.DataFrame = adata.uns["moranI"]
        top_genes = moran_df.head(3).index.tolist()
        # Only keep genes that are actually in adata
        top_genes = [g for g in top_genes if g in adata.var_names]
        if not top_genes:
            return None

        n = len(top_genes)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
        if n == 1:
            axes = [axes]

        for ax, gene in zip(axes, top_genes):
            sq.pl.spatial_scatter(
                adata,
                color=gene,
                ax=ax,
                show=False,
                frameon=False,
            )
            ax.set_title(gene)

        fig.suptitle("Top spatially variable genes", fontsize=12, fontweight="bold")
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


def _render_html(
    adata: ad.AnnData,
    cluster_info: dict,
    figures: dict,
    dataset_id: str,
) -> str:
    params = cluster_info.get("params", {})
    outputs = cluster_info.get("outputs", {})
    timestamp = cluster_info.get("timestamp", datetime.now().isoformat())

    n_clusters = outputs.get("n_clusters", "?")
    n_annotated = outputs.get("n_annotated_spots", 0)
    n_svg_tested = outputs.get("n_genes_tested", None)
    n_svg_sig = outputs.get("n_significant_fdr05", None)
    top5_svg = outputs.get("top5_svg", [])
    cluster_key = params.get("cluster_key", "spatial_cluster")
    resolution = params.get("resolution", "?")
    run_svg = params.get("run_svg", False)

    def img_section(title: str, b64: Optional[str], alt: str = "figure") -> str:
        if b64 is None:
            return f"<p><em>{title} — not available</em></p>"
        return (
            f"<h3>{title}</h3>"
            f'<img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="max-width:100%;border:1px solid #e0e0e0;border-radius:4px;">'
        )

    # Pre-compute values that would require backslashes or nested quotes
    # inside the main f-string (not allowed in Python < 3.12)
    annotation_map_val = "yes" if params.get("annotation_map_provided") else "no"
    resolution_val = params.get("resolution", "?")
    n_neighbors_val = params.get("n_neighbors", "?")
    n_pcs_actual_val = params.get("n_pcs_actual", "?")

    svg_section_html = ""  # built separately to avoid nested f-string with backslashes

    # SVG table
    svg_table_html = ""
    if "moranI" in adata.uns:
        try:
            moran_df: pd.DataFrame = adata.uns["moranI"]
            top20 = moran_df.head(20).reset_index()
            rows = ""
            for _, row in top20.iterrows():
                sig_badge = (
                    '<span style="color:#e74c3c;font-weight:600;">★</span>'
                    if row.get("pval_norm_fdr_bh", 1.0) < 0.05
                    else ""
                )
                rows += (
                    f"<tr>"
                    f"<td>{row.iloc[0]}</td>"
                    f"<td>{row['I']:.4f}</td>"
                    f"<td>{row.get('pval_norm', float('nan')):.2e}</td>"
                    f"<td>{row.get('pval_norm_fdr_bh', float('nan')):.2e} {sig_badge}</td>"
                    f"</tr>"
                )
            svg_table_html = f"""
            <table>
              <thead>
                <tr>
                  <th>Gene</th><th>Moran's I</th>
                  <th>p-value</th><th>FDR-adjusted p</th>
                </tr>
              </thead>
              <tbody>{rows}</tbody>
            </table>
            <p style="font-size:0.8rem;color:#7f8c8d;">
              ★ = FDR &lt; 0.05 &nbsp;|&nbsp; Showing top 20 of {len(moran_df)} tested genes
            </p>
            """
        except Exception:
            svg_table_html = "<p><em>SVG table could not be rendered.</em></p>"

    # Top 5 SVG badges
    top5_html = " ".join(
        f'<span style="display:inline-block;background:#fdebd0;color:#a04000;'
        f'padding:0.2rem 0.5rem;border-radius:4px;font-size:0.85rem;'
        f'font-weight:600;margin:2px;">{g}</span>'
        for g in top5_svg
    ) or "<em>none</em>"

    # Build SVG section separately — avoids nested f-string with backslash escapes
    if run_svg:
        svg_bar_fig = img_section("Top 20 SVGs \u2014 Moran's I score", figures.get("svg_bar"), "SVG bar chart")
        svg_spatial_fig = img_section("Spatial expression of top 3 SVGs", figures.get("svg_spatial"), "SVG spatial scatter")
        svg_section_html = (
            '<div class="card">'
            "<h2>Spatially Variable Genes (Moran\u2019s I)</h2>"
            + svg_bar_fig
            + svg_spatial_fig
            + "<h3>Top 20 SVGs</h3>"
            + svg_table_html
            + "</div>"
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage — Spatial Cluster Report: {dataset_id}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          margin: 0; padding: 0; background: #f8f9fa; color: #212529; }}
  .header {{ background: linear-gradient(135deg, #2c3e50, #8e44ad);
             color: white; padding: 2rem 2.5rem; }}
  .header h1 {{ margin: 0 0 0.3rem; font-size: 1.6rem; }}
  .header p  {{ margin: 0; opacity: 0.85; font-size: 0.9rem; }}
  .container {{ max-width: 1100px; margin: 0 auto; padding: 1.5rem 2rem; }}
  .card {{ background: white; border-radius: 8px; padding: 1.5rem 2rem;
           margin-bottom: 1.5rem; box-shadow: 0 1px 4px rgba(0,0,0,0.08); }}
  .card h2 {{ margin-top: 0; font-size: 1.15rem; color: #2c3e50;
              border-bottom: 2px solid #8e44ad; padding-bottom: 0.4rem; }}
  .card h3 {{ font-size: 1rem; color: #34495e; margin: 1rem 0 0.4rem; }}
  .kpi-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
               gap: 1rem; margin-bottom: 0.5rem; }}
  .kpi {{ background: #f0f4f8; border-radius: 6px; padding: 1rem;
          text-align: center; border-left: 4px solid #8e44ad; }}
  .kpi .value {{ font-size: 1.8rem; font-weight: 700; color: #2c3e50; }}
  .kpi .label {{ font-size: 0.8rem; color: #7f8c8d; margin-top: 0.2rem; }}
  .kpi.green {{ border-left-color: #27ae60; }}
  .kpi.blue  {{ border-left-color: #3498db; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
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
  <h1>🔬 OmicSage — Spatial Cluster Report</h1>
  <p>Dataset: <strong>{dataset_id}</strong> &nbsp;|&nbsp; Generated: {timestamp[:19].replace("T", " ")}</p>
</div>

<div class="container">

  <!-- KPI summary -->
  <div class="card">
    <h2>Clustering Summary</h2>
    <div class="kpi-grid">
      <div class="kpi">
        <div class="value">{adata.n_obs:,}</div>
        <div class="label">Spots clustered</div>
      </div>
      <div class="kpi">
        <div class="value">{n_clusters}</div>
        <div class="label">Leiden clusters</div>
      </div>
      <div class="kpi blue">
        <div class="value">{resolution}</div>
        <div class="label">Resolution</div>
      </div>
      {f'<div class="kpi green"><div class="value">{n_svg_sig}</div><div class="label">SVGs (FDR&lt;0.05)</div></div>' if n_svg_sig is not None else ""}
      {f'<div class="kpi"><div class="value">{n_svg_tested:,}</div><div class="label">Genes tested</div></div>' if n_svg_tested is not None else ""}
    </div>
    {f"<p><strong>Top SVGs:</strong> {top5_html}</p>" if top5_svg else ""}
  </div>

  <!-- Parameters -->
  <div class="card">
    <h2>Parameters Used</h2>
    <div class="param-grid">
      <div class="param-item"><span class="key">resolution</span><span class="val">{resolution_val}</span></div>
      <div class="param-item"><span class="key">n_neighbors (KNN)</span><span class="val">{n_neighbors_val}</span></div>
      <div class="param-item"><span class="key">n_pcs</span><span class="val">{n_pcs_actual_val}</span></div>
      <div class="param-item"><span class="key">cluster_key</span><span class="val">{cluster_key}</span></div>
      <div class="param-item"><span class="key">annotation_map</span><span class="val">{annotation_map_val}</span></div>
      <div class="param-item"><span class="key">run_svg</span><span class="val">{run_svg}</span></div>
    </div>
  </div>

  <!-- Clustering figures -->
  <div class="card">
    <h2>Clustering Figures</h2>
    {img_section("Clusters in spatial context", figures.get("spatial_clusters"), "spatial clusters")}
    {img_section("Spots per cluster", figures.get("cluster_sizes"), "cluster size bar chart")}
  </div>

  {svg_section_html}

</div>
<div class="footer">
  Generated by OmicSage &nbsp;|&nbsp; Phase 7 — Spatial Transcriptomics
</div>
</body>
</html>"""

    return html
