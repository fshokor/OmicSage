"""
spatial_deconvolve_report.py — OmicSage Phase 7
HTML report for cell2location deconvolution results.
Matches the _render_page / section structure of the RNA pipeline reports.
"""

from __future__ import annotations

import base64
import logging
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Optional

import anndata as ad
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


def _fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _mean_abundance_bar(adata, cell_type_names):
    try:
        available = [ct for ct in cell_type_names if ct in adata.obs.columns]
        if not available:
            return None
        means = {ct: float(adata.obs[ct].mean()) for ct in available}
        sorted_cts = sorted(means, key=means.get, reverse=False)
        values = [means[ct] for ct in sorted_cts]
        fig, ax = plt.subplots(figsize=(6, max(3, len(sorted_cts) * 0.4)))
        colors = plt.cm.tab20(np.linspace(0, 1, len(sorted_cts)))
        ax.barh(sorted_cts, values, color=colors, edgecolor="white")
        ax.set_xlabel("Mean cell abundance per spot (5% quantile)")
        ax.set_title("Mean cell type abundance across all spots",
                     fontsize=10, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception:
        return None


def _dominant_celltype_spatial(adata, cell_type_names):
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    try:
        available = [ct for ct in cell_type_names if ct in adata.obs.columns]
        if not available:
            return None
        adata.obs["dominant_cell_type"] = (
            adata.obs[available].idxmax(axis=1).astype("category")
        )
        fig, ax = plt.subplots(figsize=(6, 6))
        sq.pl.spatial_scatter(adata, color="dominant_cell_type",
                              ax=ax, show=False, frameon=False)
        ax.set_title("Dominant cell type per spot", fontsize=10, fontweight="bold")
        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception:
        return None


def _spatial_top6(adata, cell_type_names):
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    try:
        available = [ct for ct in cell_type_names if ct in adata.obs.columns]
        if not available:
            return None
        means = {ct: float(adata.obs[ct].mean()) for ct in available}
        top6  = sorted(means, key=means.get, reverse=True)[:6]
        n     = len(top6)
        ncols = min(3, n)
        nrows = int(np.ceil(n / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 5 * nrows))
        axes = np.array(axes).ravel() if n > 1 else [axes]
        for ax, ct in zip(axes, top6):
            sq.pl.spatial_scatter(adata, color=ct, ax=ax, show=False,
                                  frameon=False, cmap="magma")
            ax.set_title(ct, fontsize=9)
        for ax in axes[n:]:
            ax.set_visible(False)
        fig.suptitle("Top 6 cell types", fontsize=11, fontweight="bold")
        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception:
        return None


def _section_summary(adata, deconv_info, dataset_id, timestamp):
    outputs     = deconv_info.get("outputs", {})
    skipped     = deconv_info.get("skipped", True)
    skip_reason = deconv_info.get("skip_reason", "")
    n_spots     = outputs.get("n_spots", adata.n_obs)
    n_ct        = outputs.get("n_cell_types", 0)
    n_shared    = outputs.get("n_shared_genes", 0)

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Spots", f"{n_spots:,}"),
            ("Cell types", "Skipped" if skipped else n_ct),
        ] + ([] if skipped else [("Shared genes", f"{n_shared:,}")])
    )
    skip_note = (
        f'<p class="note">&#x26A0; Deconvolution skipped &mdash; {skip_reason}</p>'
        if skipped else ""
    )
    ct_names  = outputs.get("cell_type_names", [])
    ct_badges = " ".join(f"<code>{ct}</code>" for ct in ct_names)

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
      {skip_note}
      <div class="stat-grid">{stat_cards}</div>
      {"<h3>Cell Types</h3><p>" + ct_badges + "</p>" if ct_badges else ""}
    </section>
    """


def _section_params(deconv_info):
    params = deconv_info.get("params", {})
    if not params:
        return ""
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in params.items()
    )
    return f"""
    <section>
      <h2>Parameters Used</h2>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_abundance(adata, deconv_info):
    ct_names = deconv_info.get("outputs", {}).get("cell_type_names", [])
    if not ct_names:
        return ""
    b64_bar  = _mean_abundance_bar(adata, ct_names)
    b64_dom  = _dominant_celltype_spatial(adata, ct_names)
    b64_top6 = _spatial_top6(adata, ct_names)
    figs = ""
    if b64_bar:
        figs += (
            '<div class="fig-wrap"><h3>Mean abundance per cell type</h3>'
            f'<img src="data:image/png;base64,{b64_bar}" alt="mean abundance"></div>'
        )
    if b64_dom:
        figs += (
            '<div class="fig-wrap"><h3>Dominant cell type per spot</h3>'
            f'<img src="data:image/png;base64,{b64_dom}" alt="dominant cell type"></div>'
        )
    if b64_top6:
        figs += (
            '<div class="fig-wrap"><h3>Top 6 cell types &mdash; spatial</h3>'
            f'<img src="data:image/png;base64,{b64_top6}" alt="spatial top 6"></div>'
        )
    if not figs:
        return ""
    return f"""
    <section>
      <h2>Cell Type Abundances</h2>
      <div class="fig-grid">{figs}</div>
    </section>
    """


_PAGE_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
             color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
              box-shadow: 0 1px 4px rgba(0,0,0,0.07);
              padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                 border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
    section h3 { font-size: 1rem; font-weight: 600; color: #16213e; margin: 18px 0 10px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
            border-left: 3px solid #f0c040; padding: 8px 12px;
            border-radius: 4px; margin-bottom: 14px; }
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
    td { padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9ff; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 560px; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
"""


def _render_page(title, sections, timestamp):
    body = "\n".join(sections)
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "  <meta charset=\"UTF-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
        f"  <title>{title}</title>\n"
        f"  <style>{_PAGE_CSS}</style>\n"
        "</head>\n<body>\n"
        "  <header>\n"
        "    <h1>OmicSage &#8212; Spatial Deconvolution Report</h1>\n"
        f"    <p>Generated {timestamp}</p>\n"
        "  </header>\n"
        "  <main>\n"
        f"    {body}\n"
        "  </main>\n"
        "  <footer>\n"
        "    Generated by <a href=\"https://github.com/fshokor/OmicSage\">OmicSage</a>\n"
        "    &middot; MIT License\n"
        "  </footer>\n"
        "</body>\n</html>"
    )


def generate_spatial_deconvolve_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
) -> str:
    if "omicsage_spatial_deconvolve" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_deconvolve'] not found. "
            "Run spatial_deconvolve() before generating the report."
        )
    timestamp   = datetime.now().strftime("%Y-%m-%d %H:%M")
    deconv_info = adata.uns["omicsage_spatial_deconvolve"]
    skipped     = deconv_info.get("skipped", True)

    sections = [_section_summary(adata, deconv_info, dataset_id, timestamp)]
    if not skipped:
        sections += [
            _section_abundance(adata, deconv_info),
            _section_params(deconv_info),
        ]
    html = _render_page(
        title=f"OmicSage -- Spatial Deconvolution -- {dataset_id}",
        sections=sections,
        timestamp=timestamp,
    )
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    size_kb = Path(output_path).stat().st_size / 1024
    logger.info("Spatial deconvolve report -> %s (%.1f KB)", output_path, size_kb)
    return str(Path(output_path).resolve())
