"""
spatial_reduce_report.py — OmicSage Phase 7, Session 2
HTML report for spatial reduction results.

Generates a self-contained HTML report with:
  - Summary table: HVGs selected, PCA components, spatial graph stats
  - HVG scatter: mean vs dispersion with HVGs highlighted
  - PCA variance explained: elbow plot (cumulative variance per PC)
  - PCA scatter: spots coloured by total_counts (quality check)
  - Spatial connectivity: histogram of neighbours per spot (~6 for Visium)
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

matplotlib.use("Agg")

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_spatial_reduce_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
) -> str:
    """Generate a self-contained HTML reduction report for Visium data.

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_reduce` (contains
        ``uns["omicsage_spatial_reduce"]``).
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label used in the report title.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    if "omicsage_spatial_reduce" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_reduce'] not found. "
            "Run spatial_reduce() before generating the report."
        )

    reduce_info = adata.uns["omicsage_spatial_reduce"]
    figures = _build_figures(adata, reduce_info)
    html = _render_html(adata, reduce_info, figures, dataset_id)

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


def _build_figures(adata: ad.AnnData, reduce_info: dict) -> dict[str, Optional[str]]:
    """Build all figures and return as base64-encoded PNG strings."""
    figures: dict[str, Optional[str]] = {}
    figures["hvg_scatter"] = _hvg_scatter(adata)
    figures["pca_variance"] = _pca_variance(adata, reduce_info)
    figures["pca_scatter"] = _pca_scatter(adata)
    figures["spatial_neighbors"] = _spatial_neighbors_hist(adata)
    return figures


def _hvg_scatter(adata: ad.AnnData) -> Optional[str]:
    """Mean vs dispersions scatter, HVGs highlighted in orange."""
    required = {"means", "dispersions_norm", "highly_variable"}
    if not required.issubset(adata.var.columns):
        return None
    try:
        means = adata.var["means"].values
        disp = adata.var["dispersions_norm"].values
        hvg = adata.var["highly_variable"].values

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(
            means[~hvg], disp[~hvg],
            s=3, alpha=0.3, color="#95a5a6", label="non-HVG", rasterized=True,
        )
        ax.scatter(
            means[hvg], disp[hvg],
            s=5, alpha=0.6, color="#e67e22", label=f"HVG (n={hvg.sum():,})",
            rasterized=True,
        )
        ax.set_xlabel("Mean expression")
        ax.set_ylabel("Normalized dispersion")
        ax.set_title("Highly Variable Genes")
        ax.legend(fontsize=9)
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _pca_variance(adata: ad.AnnData, reduce_info: dict) -> Optional[str]:
    """Elbow plot: individual + cumulative explained variance per PC."""
    try:
        vr = np.array(adata.uns["pca"]["variance_ratio"])
        cumvr = np.cumsum(vr)
        pcs = np.arange(1, len(vr) + 1)

        fig, ax1 = plt.subplots(figsize=(7, 4))
        ax2 = ax1.twinx()

        ax1.bar(pcs, vr * 100, color="#3498db", alpha=0.7, label="Per-PC variance %")
        ax2.plot(pcs, cumvr * 100, color="#e74c3c", linewidth=2,
                 marker="o", markersize=3, label="Cumulative variance %")
        ax2.axhline(90, color="#e74c3c", linestyle="--", linewidth=0.8, alpha=0.5)

        ax1.set_xlabel("Principal component")
        ax1.set_ylabel("Variance explained (%)", color="#3498db")
        ax2.set_ylabel("Cumulative variance (%)", color="#e74c3c")
        ax1.set_title("PCA — Variance Explained")
        ax1.tick_params(axis="y", colors="#3498db")
        ax2.tick_params(axis="y", colors="#e74c3c")

        n_comps = reduce_info["outputs"].get("n_comps_computed", len(vr))
        ax1.set_xlim(0.5, min(n_comps + 0.5, len(vr) + 0.5))
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _pca_scatter(adata: ad.AnnData) -> Optional[str]:
    """PCA embedding scatter coloured by total_counts (spot quality check)."""
    if "X_pca" not in adata.obsm:
        return None
    try:
        pca = adata.obsm["X_pca"]
        color_vals = (
            adata.obs["total_counts"].values
            if "total_counts" in adata.obs.columns
            else np.ones(adata.n_obs)
        )
        label = "total_counts" if "total_counts" in adata.obs.columns else ""

        fig, ax = plt.subplots(figsize=(5, 4))
        sc_plot = ax.scatter(
            pca[:, 0], pca[:, 1],
            c=color_vals, cmap="viridis", s=4, alpha=0.7, rasterized=True,
        )
        plt.colorbar(sc_plot, ax=ax, label=label)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title(f"PCA — coloured by {label}" if label else "PCA")
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


def _spatial_neighbors_hist(adata: ad.AnnData) -> Optional[str]:
    """Histogram of number of spatial neighbours per spot."""
    if "spatial_connectivities" not in adata.obsp:
        return None
    try:
        conn = adata.obsp["spatial_connectivities"]
        n_neighbors_per_spot = np.array(conn.sum(axis=1)).ravel().astype(int)

        fig, ax = plt.subplots(figsize=(5, 3))
        unique, counts = np.unique(n_neighbors_per_spot, return_counts=True)
        ax.bar(unique, counts, color="#2ecc71", edgecolor="white")
        ax.axvline(6, color="#e74c3c", linestyle="--", linewidth=1.2,
                   label="expected 6 (Visium grid)")
        ax.set_xlabel("Neighbours per spot")
        ax.set_ylabel("Number of spots")
        ax.set_title("Spatial neighbours distribution")
        ax.legend(fontsize=8)
        ax.set_xticks(sorted(unique))
        fig.tight_layout()
        return _fig_to_base64(fig)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------


def _render_html(
    adata: ad.AnnData,
    reduce_info: dict,
    figures: dict,
    dataset_id: str,
) -> str:
    params = reduce_info.get("params", {})
    outputs = reduce_info.get("outputs", {})
    timestamp = reduce_info.get("timestamp", datetime.now().isoformat())

    n_hvgs = outputs.get("n_hvgs", "?")
    n_comps = outputs.get("n_comps_computed", "?")
    n_edges = outputs.get("spatial_graph_n_edges", "?")
    mean_nb = outputs.get("spatial_graph_mean_neighbors", "?")
    cum_var = outputs.get("pca_cumulative_variance_top10", None)
    skipped_norm = params.get("skipped_normalization", False)

    cum_var_str = f"{cum_var * 100:.1f}%" if cum_var is not None else "?"

    def img_section(title: str, b64: Optional[str], alt: str = "figure") -> str:
        if b64 is None:
            return (
                f"<p><em>{title} — not available</em></p>"
            )
        return (
            f"<h3>{title}</h3>"
            f'<img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="max-width:100%;border:1px solid #e0e0e0;border-radius:4px;">'
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage — Spatial Reduce Report: {dataset_id}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
          margin: 0; padding: 0; background: #f8f9fa; color: #212529; }}
  .header {{ background: linear-gradient(135deg, #2c3e50, #27ae60);
             color: white; padding: 2rem 2.5rem; }}
  .header h1 {{ margin: 0 0 0.3rem; font-size: 1.6rem; }}
  .header p  {{ margin: 0; opacity: 0.85; font-size: 0.9rem; }}
  .container {{ max-width: 1100px; margin: 0 auto; padding: 1.5rem 2rem; }}
  .card {{ background: white; border-radius: 8px; padding: 1.5rem 2rem;
           margin-bottom: 1.5rem; box-shadow: 0 1px 4px rgba(0,0,0,0.08); }}
  .card h2 {{ margin-top: 0; font-size: 1.15rem; color: #2c3e50;
              border-bottom: 2px solid #27ae60; padding-bottom: 0.4rem; }}
  .card h3 {{ font-size: 1rem; color: #34495e; margin: 1rem 0 0.4rem; }}
  .kpi-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
               gap: 1rem; margin-bottom: 0.5rem; }}
  .kpi {{ background: #f0f4f8; border-radius: 6px; padding: 1rem;
          text-align: center; border-left: 4px solid #27ae60; }}
  .kpi .value {{ font-size: 1.8rem; font-weight: 700; color: #2c3e50; }}
  .kpi .label {{ font-size: 0.8rem; color: #7f8c8d; margin-top: 0.2rem; }}
  .kpi.info {{ border-left-color: #3498db; }}
  .kpi.warn {{ border-left-color: #e67e22; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
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
  .badge {{ display: inline-block; padding: 0.15rem 0.5rem; border-radius: 4px;
            font-size: 0.78rem; font-weight: 600; }}
  .badge-green {{ background: #d5f5e3; color: #1e8449; }}
  .badge-orange {{ background: #fdebd0; color: #a04000; }}
  .footer {{ text-align: center; padding: 1.5rem; color: #95a5a6; font-size: 0.8rem; }}
</style>
</head>
<body>
<div class="header">
  <h1>🔬 OmicSage — Spatial Reduce Report</h1>
  <p>Dataset: <strong>{dataset_id}</strong> &nbsp;|&nbsp; Generated: {timestamp[:19].replace("T", " ")}</p>
</div>

<div class="container">

  <!-- KPI summary -->
  <div class="card">
    <h2>Reduction Summary</h2>
    <div class="kpi-grid">
      <div class="kpi">
        <div class="value">{adata.n_obs:,}</div>
        <div class="label">Spots (input)</div>
      </div>
      <div class="kpi">
        <div class="value">{adata.n_vars:,}</div>
        <div class="label">Genes (total)</div>
      </div>
      <div class="kpi">
        <div class="value">{n_hvgs:,}</div>
        <div class="label">HVGs selected</div>
      </div>
      <div class="kpi info">
        <div class="value">{n_comps}</div>
        <div class="label">PCA components</div>
      </div>
      <div class="kpi info">
        <div class="value">{cum_var_str}</div>
        <div class="label">Variance (top 10 PCs)</div>
      </div>
      <div class="kpi">
        <div class="value">{n_edges:,}</div>
        <div class="label">Spatial graph edges</div>
      </div>
      <div class="kpi {'warn' if isinstance(mean_nb, float) and abs(mean_nb - 6) > 0.5 else ''}">
        <div class="value">{mean_nb}</div>
        <div class="label">Mean neighbours/spot</div>
      </div>
    </div>
    {"<p><span class='badge badge-orange'>⚠ Normalization skipped</span> — benchmark dataset detected (pre-processed).</p>" if skipped_norm else ""}
  </div>

  <!-- Parameters -->
  <div class="card">
    <h2>Parameters Used</h2>
    <div class="param-grid">
      <div class="param-item"><span class="key">n_top_genes</span><span class="val">{params.get('n_top_genes','?'):,}</span></div>
      <div class="param-item"><span class="key">n_comps</span><span class="val">{params.get('n_comps','?')}</span></div>
      <div class="param-item"><span class="key">n_neighbors</span><span class="val">{params.get('n_neighbors','?')}</span></div>
      <div class="param-item"><span class="key">coord_type</span><span class="val">{params.get('coord_type') or 'None (auto)'}</span></div>
      <div class="param-item"><span class="key">normalize_total</span><span class="val">{params.get('normalize_total','?')}</span></div>
      <div class="param-item"><span class="key">target_sum</span><span class="val">{params.get('target_sum','?'):.0f}</span></div>
      <div class="param-item"><span class="key">log1p</span><span class="val">{params.get('log1p','?')}</span></div>
      <div class="param-item"><span class="key">hvg_flavor</span><span class="val">{params.get('flavor','?')}</span></div>
    </div>
  </div>

  <!-- Figures -->
  <div class="card">
    <h2>Figures</h2>
    {img_section("Highly Variable Genes — mean vs dispersion", figures.get("hvg_scatter"), "HVG scatter")}
    {img_section("PCA — Variance Explained (elbow plot)", figures.get("pca_variance"), "PCA variance")}
    {img_section("PCA — Spots coloured by total_counts", figures.get("pca_scatter"), "PCA scatter")}
    {img_section("Spatial neighbours per spot", figures.get("spatial_neighbors"), "neighbours histogram")}
  </div>

</div>
<div class="footer">
  Generated by OmicSage &nbsp;|&nbsp; Phase 7 — Spatial Transcriptomics
</div>
</body>
</html>"""

    return html
