"""
spatial_downstream_report.py — OmicSage Phase 7, Session 5
HTML report for spatial downstream analyses.

Style: identical _PAGE_CSS / _render_page pattern as all other spatial reports.

Sections:
  1. Run Summary          — stat cards + analysis status table
  2. Region Clustering    — spatial scatter coloured by region_cluster
  3. Cell-type Markers    — table: top 10 correlated genes per cell type
  4. Cell-type SVGs       — table: top 5 SVGs per cell type
  5. Co-occurrence        — heatmap / note
  6. Neighbourhood Enrichment — heatmap / note
  7. Ligand-Receptor      — top interactions table
  8. SVG Pathway Enrichment — bar chart of top 10 GO terms
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
import pandas as pd

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Shared CSS + page renderer  (identical to all other spatial reports)
# ---------------------------------------------------------------------------

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
    .skip-note { font-size: 0.82rem; color: #555; background: #f4f4f4;
                 border-left: 3px solid #ccc; padding: 8px 12px;
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


def _render_page(title: str, header_subtitle: str, sections: list[str], timestamp: str) -> str:
    body = "\n".join(s for s in sections if s)
    return (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "  <meta charset=\"UTF-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
        f"  <title>{title}</title>\n"
        f"  <style>{_PAGE_CSS}</style>\n"
        "</head>\n<body>\n"
        "  <header>\n"
        "    <h1>OmicSage &#8212; Spatial Downstream Report</h1>\n"
        f"    <p>{header_subtitle} &middot; Generated {timestamp}</p>\n"
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


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _img_tag(b64: str, alt: str = "figure") -> str:
    return f'<img src="data:image/png;base64,{b64}" alt="{alt}">'


def _stat_card(value: str, label: str) -> str:
    return (
        f'<div class="stat-card">'
        f'<div class="stat-value">{value}</div>'
        f'<div class="stat-label">{label}</div>'
        f'</div>'
    )


def _skip_section(title: str, reason: str) -> str:
    return (
        f'<section><h2>{title}</h2>'
        f'<p class="skip-note">Not run / data not available: {reason}</p>'
        f'</section>'
    )


def _squidpy_fig_b64(
    adata: ad.AnnData,
    plot_fn: str,
    figsize: tuple = (6, 5),
    **kwargs,
) -> Optional[str]:
    """Call a sq.pl.* function by name and return base64 PNG or None.

    squidpy plot functions may return an Axes, a list of Axes, or None
    depending on version and whether show=True/False is honoured.
    We capture the figure via the returned object when possible, and fall
    back to plt.gcf() otherwise.  We always close all figures afterwards to
    avoid memory leaks.
    """
    if not _SQUIDPY_AVAILABLE:
        return None
    try:
        plt.close("all")
        fn = getattr(sq.pl, plot_fn)
        # Pass figsize as a keyword; some sq.pl functions ignore it, which is fine.
        result = fn(adata, figsize=figsize, **kwargs)
        # Resolve the figure object from whatever sq.pl.* returned
        if result is None:
            fig = plt.gcf()
        elif hasattr(result, "get_figure"):
            # Single Axes object
            fig = result.get_figure()
        elif hasattr(result, "__len__") and len(result) > 0:
            # List/array of Axes
            first = result[0]
            if hasattr(first, "get_figure"):
                fig = first.get_figure()
            else:
                fig = plt.gcf()
        else:
            fig = plt.gcf()
        # Sanity-check: if gcf() returns an empty figure squidpy may have
        # plotted into a different figure — use the most recently created one.
        if not fig.axes:
            figs = list(map(plt.figure, plt.get_fignums()))
            if figs:
                fig = figs[-1]
        try:
            fig.tight_layout()
        except Exception:
            pass
        return _fig_to_b64(fig)
    except Exception as e:
        logger.warning("[downstream_report] sq.pl.%s failed: %s", plot_fn, e)
        plt.close("all")
        return None


def _spatial_scatter_b64(
    adata: ad.AnnData,
    color: str,
    figsize: tuple = (6, 6),
    library_key: Optional[str] = None,
) -> Optional[str]:
    """Spatial scatter coloured by an obs column.

    Multi-sample strategy: render each sample separately and stitch into a
    2-column grid to avoid extreme aspect ratios with ≥2 samples.
    """
    import numpy as _np
    import io as _io
    import base64 as _b64_mod

    if not _SQUIDPY_AVAILABLE or "spatial" not in adata.obsm:
        return None

    try:
        samples = []
        if library_key and library_key in adata.obs.columns:
            samples = sorted(adata.obs[library_key].astype(str).unique().tolist())

        # Single-sample fast path
        if len(samples) <= 1:
            plt.close("all")
            kwargs: dict = {"color": color, "frameon": False, "figsize": figsize}
            if library_key:
                kwargs["library_key"] = library_key
            result = sq.pl.spatial_scatter(adata, **kwargs)
            if result is None:
                fig = plt.gcf()
            elif hasattr(result, "get_figure"):
                fig = result.get_figure()
            elif hasattr(result, "__len__") and len(result) > 0:
                ax0 = result[0]
                fig = ax0.get_figure() if hasattr(ax0, "get_figure") else plt.gcf()
            else:
                fig = plt.gcf()
            if not fig.axes:
                all_figs = [plt.figure(n) for n in plt.get_fignums()]
                if all_figs:
                    fig = all_figs[-1]
            try:
                fig.tight_layout()
            except Exception:
                pass
            return _fig_to_b64(fig)

        # Multi-sample: per-sample render + 2-column grid stitch
        per_b64s: list[tuple[str, str]] = []
        for sample_id in samples:
            try:
                mask = adata.obs[library_key].astype(str) == sample_id
                ad_sub = adata[mask].copy()
                # uns["spatial"] is NOT filtered by obs subsetting — prune to
                # the current sample so squidpy doesn't demand a library_key.
                if "spatial" in ad_sub.uns and sample_id in ad_sub.uns["spatial"]:
                    ad_sub.uns["spatial"] = {sample_id: ad_sub.uns["spatial"][sample_id]}
                plt.close("all")
                kwargs = {"color": color, "frameon": False, "figsize": figsize}
                result = sq.pl.spatial_scatter(ad_sub, **kwargs)
                if result is None:
                    fig = plt.gcf()
                elif hasattr(result, "get_figure"):
                    fig = result.get_figure()
                elif hasattr(result, "__len__") and len(result) > 0:
                    ax0 = result[0]
                    fig = ax0.get_figure() if hasattr(ax0, "get_figure") else plt.gcf()
                else:
                    fig = plt.gcf()
                if not fig.axes:
                    all_figs = [plt.figure(n) for n in plt.get_fignums()]
                    if all_figs:
                        fig = all_figs[-1]
                try:
                    fig.tight_layout()
                except Exception:
                    pass
                buf = _io.BytesIO()
                fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
                buf.seek(0)
                per_b64s.append((sample_id, _b64_mod.b64encode(buf.read()).decode()))
                plt.close(fig)
            except Exception as e:
                logger.warning("[downstream_report] scatter failed for sample %s: %s", sample_id, e)
                plt.close("all")

        if not per_b64s:
            return None

        imgs = []
        for _, b in per_b64s:
            arr = plt.imread(_io.BytesIO(_b64_mod.b64decode(b)))
            imgs.append(arr)

        ncols = 2
        nrows = (len(imgs) + ncols - 1) // ncols
        max_h = max(a.shape[0] for a in imgs)
        max_w = max(a.shape[1] for a in imgs)
        n_ch  = imgs[0].shape[2] if imgs[0].ndim == 3 else 1

        def _pad(a):
            ph = max_h - a.shape[0]
            pw = max_w - a.shape[1]
            return _np.pad(a, ((0, ph), (0, pw), (0, 0)),
                           mode="constant", constant_values=1.0)

        padded = [_pad(a) for a in imgs]
        while len(padded) < nrows * ncols:
            padded.append(_np.ones((max_h, max_w, n_ch), dtype=padded[0].dtype))

        rows = [_np.hstack(padded[r * ncols:(r + 1) * ncols]) for r in range(nrows)]
        grid = _np.vstack(rows)

        dpi = 100
        fig2, ax2 = plt.subplots(figsize=(grid.shape[1] / dpi, grid.shape[0] / dpi), dpi=dpi)
        ax2.imshow(grid)
        ax2.axis("off")
        for idx, (sample_id, _) in enumerate(per_b64s):
            row_i = idx // ncols
            col_i = idx % ncols
            ax2.text(col_i * max_w + max_w * 0.02, row_i * max_h + max_h * 0.04,
                     sample_id, fontsize=8, fontweight="bold", color="white", va="top",
                     bbox=dict(facecolor="black", alpha=0.45, pad=2, linewidth=0))
        fig2.subplots_adjust(left=0, right=1, top=1, bottom=0)

        buf2 = _io.BytesIO()
        fig2.savefig(buf2, format="png", dpi=dpi, bbox_inches="tight")
        buf2.seek(0)
        result_b64 = _b64_mod.b64encode(buf2.read()).decode()
        plt.close(fig2)
        return result_b64

    except Exception as e:
        logger.warning("[downstream_report] spatial_scatter(%s) failed: %s", color, e)
        plt.close("all")
        return None


def _get_library_key(adata: ad.AnnData) -> Optional[str]:
    ingest_prov = adata.uns.get("omicsage_spatial_ingest", {})
    if "library_key" in ingest_prov and ingest_prov["library_key"] is not None:
        return ingest_prov["library_key"]
    library_ids = list(adata.uns.get("spatial", {}).keys())
    if len(library_ids) <= 1:
        return None
    id_set = set(library_ids)
    for candidate in ("library_id", "sample", "patient", "donor_id", "batch"):
        if candidate in adata.obs.columns:
            if id_set.issubset(set(adata.obs[candidate].astype(str).unique())):
                return candidate
    return None


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _section_summary(
    adata: ad.AnnData,
    prov: dict,
    dataset_id: str,
    timestamp: str,
) -> str:
    analyses = prov.get("analyses", {})

    def _status(key: str) -> str:
        info = analyses.get(key, {})
        if not info:
            return "not requested"
        if info.get("skipped"):
            return f"skipped ({info.get('reason', '')})"
        return "✓"

    n_regions   = analyses.get("region_clustering", {}).get("n_regions", "—")
    n_ct_expr   = analyses.get("celltype_expression", {}).get("n_cell_types", "—")
    n_ct_svg    = analyses.get("celltype_svg", {}).get("n_cell_types", "—")
    n_pathways  = analyses.get("svg_gsea", {}).get("n_pathways", "—")

    stat_cards = (
        _stat_card(str(n_regions),  "region clusters")
        + _stat_card(str(n_ct_expr), "cell types (expr)")
        + _stat_card(str(n_ct_svg),  "cell types (SVG)")
        + _stat_card(str(n_pathways), "SVG pathways")
    )

    status_rows = "".join(
        f"<tr><td>{label}</td><td>{_status(key)}</td></tr>"
        for label, key in [
            ("Region clustering",      "region_clustering"),
            ("Cell-type expression",   "celltype_expression"),
            ("Cell-type SVGs",         "celltype_svg"),
            ("Co-occurrence",          "co_occurrence"),
            ("Neighbourhood enrichment","nhood_enrichment"),
            ("Ligand-receptor",        "ligrec"),
            ("SVG pathway enrichment", "svg_gsea"),
        ]
    )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>Analysis Status</h3>
      <table>
        <thead><tr><th>Analysis</th><th>Status</th></tr></thead>
        <tbody>{status_rows}</tbody>
      </table>
    </section>
    """


def _section_region_clustering(adata: ad.AnnData, prov: dict) -> str:
    info = prov.get("analyses", {}).get("region_clustering", {})
    if info.get("skipped"):
        return _skip_section("Region Clustering", info.get("reason", "not run"))

    library_key = _get_library_key(adata)

    b64_spatial = None
    if "region_cluster" in adata.obs.columns:
        b64_spatial = _spatial_scatter_b64(
            adata, "region_cluster", figsize=(6, 5), library_key=library_key
        )

    b64_umap = None
    if "X_umap_celltype" in adata.obsm and "region_cluster" in adata.obs.columns:
        try:
            prev = adata.obsm.get("X_umap")
            adata.obsm["X_umap"] = adata.obsm["X_umap_celltype"]
            fig, ax = plt.subplots(figsize=(5, 4))
            import scanpy as sc
            sc.pl.umap(adata, color="region_cluster", ax=ax, show=False, frameon=False)
            try:
                fig.tight_layout()
            except Exception:
                pass
            b64_umap = _fig_to_b64(fig)
            if prev is not None:
                adata.obsm["X_umap"] = prev
            else:
                adata.obsm.pop("X_umap", None)
        except Exception as e:
            logger.warning("[downstream_report] region UMAP failed: %s", e)
            plt.close("all")

    n_regions = info.get("n_regions", "?")
    resolution = info.get("resolution", "?")

    figs = ""
    if b64_umap:
        figs += (
            '<div class="fig-wrap"><h3>UMAP of cell type composition</h3>'
            + _img_tag(b64_umap, "region cluster UMAP")
            + "</div>"
        )
    if b64_spatial:
        figs += (
            '<div class="fig-wrap"><h3>Region clusters on tissue</h3>'
            + _img_tag(b64_spatial, "region clusters spatial")
            + "</div>"
        )

    return f"""
    <section>
      <h2>Region Clustering (Cell Type Composition)</h2>
      <p>Spots are clustered by their deconvolved cell type composition using Leiden
         clustering on the cell type abundance KNN graph (resolution={resolution}).
         Each region represents a tissue niche with a characteristic cell type mixture.</p>
      <p class="note"><strong>{n_regions}</strong> tissue regions identified.</p>
      {('<div class="fig-grid">' + figs + "</div>") if figs else ""}
    </section>
    """


def _section_celltype_markers(adata: ad.AnnData, prov: dict) -> str:
    info = prov.get("analyses", {}).get("celltype_expression", {})
    if info.get("skipped"):
        return _skip_section("Cell-type Marker Genes", info.get("reason", "not run"))

    marker_dict: dict = adata.uns.get("celltype_marker_genes", {})
    if not marker_dict:
        return _skip_section("Cell-type Marker Genes", "no results in uns['celltype_marker_genes']")

    rows = ""
    for ct, genes in sorted(marker_dict.items()):
        gene_str = ", ".join(f"<code>{g}</code>" for g in genes[:10])
        rows += f"<tr><td>{ct}</td><td>{gene_str}</td></tr>"

    return f"""
    <section>
      <h2>Cell-type Marker Genes (Spatial Correlation)</h2>
      <p>Genes whose expression most strongly correlates (Spearman) with each cell type's
         abundance across spots. Top 10 of {info.get('n_marker_genes', '?')} stored genes shown.</p>
      <table>
        <thead><tr><th>Cell Type</th><th>Top Correlated Genes</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_celltype_svg(adata: ad.AnnData, prov: dict) -> str:
    info = prov.get("analyses", {}).get("celltype_svg", {})
    if info.get("skipped"):
        return _skip_section("Cell-type Specific SVGs", info.get("reason", "not run"))

    celltype_svg: dict = adata.uns.get("celltype_svg", {})
    if not celltype_svg:
        return _skip_section("Cell-type Specific SVGs", "no results in uns['celltype_svg']")

    rows = ""
    for ct, df in sorted(celltype_svg.items()):
        try:
            top5 = df.head(5).index.tolist()
            gene_str = ", ".join(f"<code>{g}</code>" for g in top5)
            rows += f"<tr><td>{ct}</td><td>{gene_str}</td></tr>"
        except Exception:
            rows += f"<tr><td>{ct}</td><td><em>error</em></td></tr>"

    return f"""
    <section>
      <h2>Cell-type Specific Spatially Variable Genes</h2>
      <p>Moran's I computed on spots where each cell type is above median abundance.
         Top 5 SVGs per cell type shown (analytical p-values).</p>
      <table>
        <thead><tr><th>Cell Type</th><th>Top SVGs</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_co_occurrence(adata: ad.AnnData, prov: dict, dominant_celltype_key: str) -> str:
    info = prov.get("analyses", {}).get("co_occurrence", {})
    if info.get("skipped"):
        return _skip_section("Spatial Co-occurrence", info.get("reason", "not run"))

    co_key = f"{dominant_celltype_key}_co_occurrence"
    if co_key not in adata.uns:
        return _skip_section("Spatial Co-occurrence", f"uns['{co_key}'] not found")

    # Pick the most abundant cluster to illustrate
    b64 = None
    if _SQUIDPY_AVAILABLE and dominant_celltype_key in adata.obs.columns:
        try:
            top_cluster = (
                adata.obs[dominant_celltype_key]
                .value_counts()
                .index[0]
            )
            # squidpy reads uns[f"{cluster_key}_colors"] for the palette.
            # After deconvolution the color list may be sized for all cell types
            # (11) while dominant_cell_type only has 6 — this causes a palette
            # length mismatch.  Strip the stale key so squidpy auto-generates
            # a fresh palette from the actual categories present.
            import copy as _copy
            adata_plot = _copy.copy(adata)
            adata_plot.uns = dict(adata.uns)
            colors_key = f"{dominant_celltype_key}_colors"
            adata_plot.uns.pop(colors_key, None)
            # Also ensure the column is Categorical with only observed categories
            col = adata_plot.obs[dominant_celltype_key]
            if hasattr(col, "cat"):
                adata_plot.obs = adata_plot.obs.copy()
                adata_plot.obs[dominant_celltype_key] = (
                    col.astype(str).astype("category")
                )
            b64 = _squidpy_fig_b64(
                adata_plot,
                "co_occurrence",
                figsize=(8, 5),
                cluster_key=dominant_celltype_key,
                clusters=top_cluster,
            )
        except Exception as e:
            logger.warning("[downstream_report] co_occurrence plot failed: %s", e)

    fig_html = ""
    if b64:
        fig_html = (
            '<div class="fig-grid">'
            '<div class="fig-wrap"><h3>Co-occurrence scores vs distance</h3>'
            + _img_tag(b64, "co-occurrence plot")
            + "</div></div>"
        )

    return f"""
    <section>
      <h2>Spatial Co-occurrence</h2>
      <p>Co-occurrence score (conditional probability ratio) of each cell type
         co-occurring with the most abundant type across increasing radii.
         Scores &gt; 1 indicate spatial enrichment at that distance.</p>
      {fig_html if fig_html else '<p class="note">Co-occurrence data computed but plot could not be rendered (check logs for details).</p>'}
    </section>
    """


def _section_nhood_enrichment(adata: ad.AnnData, prov: dict, dominant_celltype_key: str) -> str:
    info = prov.get("analyses", {}).get("nhood_enrichment", {})
    if info.get("skipped"):
        return _skip_section("Neighbourhood Enrichment", info.get("reason", "not run"))

    nhood_key = f"{dominant_celltype_key}_nhood_enrichment"
    if nhood_key not in adata.uns:
        return _skip_section("Neighbourhood Enrichment", f"uns['{nhood_key}'] not found")

    b64 = None
    if _SQUIDPY_AVAILABLE:
        try:
            # The nhood enrichment z-score matrix is sized for all Categorical
            # levels, but dominant_cell_type may have fewer observed categories
            # than the original deconvolution cell type list.  Casting to string
            # drops unused levels and ensures the matrix size matches the obs column.
            # Also strip any stale _colors key to avoid palette length mismatches.
            import copy as _copy
            adata_plot = _copy.copy(adata)
            adata_plot.uns = dict(adata.uns)
            colors_key = f"{dominant_celltype_key}_colors"
            adata_plot.uns.pop(colors_key, None)
            if dominant_celltype_key in adata_plot.obs.columns:
                col = adata_plot.obs[dominant_celltype_key]
                if hasattr(col, "cat"):
                    adata_plot.obs = adata_plot.obs.copy()
                    adata_plot.obs[dominant_celltype_key] = (
                        col.astype(str).astype("category")
                    )
            b64 = _squidpy_fig_b64(
                adata_plot,
                "nhood_enrichment",
                figsize=(7, 6),
                cluster_key=dominant_celltype_key,
                method="average",
            )
        except Exception as e:
            logger.warning("[downstream_report] nhood_enrichment plot failed: %s", e)

    fig_html = ""
    if b64:
        fig_html = (
            '<div class="fig-grid">'
            '<div class="fig-wrap"><h3>Neighbourhood enrichment z-scores</h3>'
            + _img_tag(b64, "neighbourhood enrichment heatmap")
            + "</div></div>"
        )

    n_perms = info.get("n_perms", "?")
    return f"""
    <section>
      <h2>Neighbourhood Enrichment</h2>
      <p>Permutation-based test ({n_perms} permutations) on the spatial adjacency graph.
         High z-scores indicate cell type pairs that are more often neighbours than
         expected by chance. Negative values indicate spatial exclusion.</p>
      {fig_html if fig_html else '<p class="note">Enrichment data computed but plot could not be rendered (check logs for details).</p>'}
    </section>
    """


def _deserialize_ligrec(raw: dict) -> dict:
    """Reconstruct ligrec dict-of-DataFrames from the serialized checkpoint format.

    _serialize_ligrec_uns stores each DataFrame as three keys:
        <name>_data     JSON records string
        <name>_columns  JSON array of stringified column names (tuple strings)
        <name>_index    JSON array of stringified index names (tuple strings)

    This function reassembles each DataFrame and restores the MultiIndex on
    both axes so sq.pl.ligrec receives the format it expects.
    """
    import json as _json
    from ast import literal_eval

    def _parse_label(s: str):
        try:
            return literal_eval(s)
        except Exception:
            return s

    # Collect data/columns/index triplets
    keys = {k.rsplit("_", 1)[0] for k in raw if k.endswith(("_data", "_columns", "_index"))}
    result: dict = {}
    for k in keys:
        data_json    = raw.get(f"{k}_data")
        columns_json = raw.get(f"{k}_columns")
        index_json   = raw.get(f"{k}_index")
        if data_json and columns_json and index_json:
            cols = [_parse_label(c) for c in _json.loads(columns_json)]
            idx  = [_parse_label(i) for i in _json.loads(index_json)]
            df   = pd.read_json(data_json, orient="records")
            df.columns = pd.MultiIndex.from_tuples(cols) if isinstance(cols[0], tuple) else cols
            df.index   = pd.MultiIndex.from_tuples(idx)  if isinstance(idx[0],  tuple) else idx
            result[k] = df
        else:
            result[k] = raw.get(f"{k}_data", raw.get(k))
    # pass through any non-DataFrame keys (metadata scalars etc.)
    for k, v in raw.items():
        base = k.rsplit("_", 1)[0]
        if not k.endswith(("_data", "_columns", "_index")) and base not in result:
            result[k] = v
    return result


def _lr_bar_chart(means_df: "pd.DataFrame", pvals_df: "pd.DataFrame",
                  n_top: int = 20, alpha: float = 0.001) -> Optional[str]:
    """Option 2 — ranked horizontal bar chart of top LR pairs by -log10(p).

    Rows = LR pairs (ligand_receptor label).
    Colour = mean expression (viridis).
    X-axis = -log10(min p-value across all cell type pairs).
    Only pairs significant in at least one cell type combination are shown.
    """
    try:
        import numpy as _np

        # Flatten MultiIndex columns → per LR pair: min p-value and max mean
        min_p  = pvals_df.min(axis=1)          # Series indexed by (ligand, receptor)
        max_mu = means_df.max(axis=1)

        sig = min_p[min_p < alpha].nsmallest(n_top)
        if sig.empty:
            return None

        labels   = [f"{a} → {b}" for a, b in sig.index]
        neg_logp = -_np.log10(sig.values.clip(1e-300))
        colours  = max_mu.loc[sig.index].values

        fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.32)))
        sc = ax.barh(
            range(len(labels)), neg_logp,
            color=plt.cm.viridis(colours / (colours.max() or 1)),
            edgecolor="none", height=0.7,
        )
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel(r"$-\log_{10}(p)$", fontsize=9)
        ax.set_title(
            f"Top {len(labels)} LR pairs  (p < {alpha})", fontsize=10, fontweight="bold"
        )
        ax.axvline(-_np.log10(alpha), color="#c0392b", lw=1, ls="--", alpha=0.6)
        ax.spines[["top", "right"]].set_visible(False)

        # Colorbar for mean expression
        sm = plt.cm.ScalarMappable(
            cmap="viridis",
            norm=plt.Normalize(vmin=0, vmax=float(colours.max() or 1)),
        )
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
        cbar.set_label("max mean expr", fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception as e:
        logger.warning("[downstream_report] LR bar chart failed: %s", e)
        return None


def _lr_focused_dotplot(adata: "ad.AnnData", ligrec_key: str,
                        dominant_celltype_key: str,
                        pvals_df: "pd.DataFrame",
                        n_cell_types: int = 6,
                        n_top_pairs: int = 40,
                        alpha: float = 0.001) -> Optional[str]:
    """Focused squidpy dotplot restricted to top N cell type pairs × top N LR pairs.

    squidpy has no built-in row cap — the only way to limit the y-axis is to
    physically subset uns[ligrec_key]["means"] and ["pvalues"] to the rows we
    want *before* passing to sq.pl.ligrec.  We write the subset into a shallow
    copy of adata.uns so the original data is never mutated.

    Row selection: top `n_top_pairs` LR pairs ranked by minimum p-value across
    the selected cell-type column pairs.
    """
    if not _SQUIDPY_AVAILABLE:
        return None
    try:
        import copy as _copy
        import numpy as _np

        # ── Select top source / target cell types ────────────────────────────
        sig_mask = pvals_df < alpha
        src_counts = sig_mask.sum(axis=0).groupby(level=0).sum()
        tgt_counts = sig_mask.sum(axis=0).groupby(level=1).sum()
        top_src = src_counts.nlargest(n_cell_types).index.tolist()
        top_tgt = tgt_counts.nlargest(n_cell_types).index.tolist()
        if not top_src or not top_tgt:
            return None

        # ── Retrieve full means/pvals from uns ───────────────────────────────
        ligrec_data = adata.uns[ligrec_key]
        means_full = ligrec_data.get("means")
        pvals_full = ligrec_data.get("pvalues")
        if means_full is None or pvals_full is None:
            return None

        # ── Column subset: keep only top_src × top_tgt pairs ────────────────
        keep_cols = [
            c for c in pvals_full.columns
            if c[0] in top_src and c[1] in top_tgt
        ]
        if not keep_cols:
            return None
        pvals_sub = pvals_full[keep_cols]
        means_sub = means_full[keep_cols] if all(c in means_full.columns for c in keep_cols) else means_full

        # ── Row subset: top n_top_pairs by min p-value across kept columns ───
        row_min_p = pvals_sub.min(axis=1)
        # Keep only rows that are significant in at least one kept column
        sig_rows = row_min_p < alpha
        if sig_rows.sum() == 0:
            # Relax to show something
            sig_rows = row_min_p < 0.05
        if sig_rows.sum() == 0:
            return None

        top_rows = (
            row_min_p[sig_rows]
            .nsmallest(n_top_pairs)
            .index
        )
        pvals_trimmed = pvals_sub.loc[top_rows]
        means_trimmed = means_sub.loc[top_rows]

        # ── Build a shallow adata copy with trimmed ligrec in uns ────────────
        adata_plot = _copy.copy(adata)
        adata_plot.uns = dict(adata.uns)
        adata_plot.uns[ligrec_key] = dict(ligrec_data)
        adata_plot.uns[ligrec_key]["means"]   = means_trimmed
        adata_plot.uns[ligrec_key]["pvalues"] = pvals_trimmed

        # ── Figsize from actual trimmed row/col counts ───────────────────────
        n_rows = len(pvals_trimmed)
        n_cols = len(keep_cols)
        fig_w = min(4 + n_cols * 0.9, 22)
        fig_h = min(3 + n_rows * 0.35, 20)

        return _squidpy_fig_b64(
            adata_plot,
            "ligrec",
            figsize=(fig_w, fig_h),
            cluster_key=dominant_celltype_key,
            source_groups=top_src,
            target_groups=top_tgt,
            pvalue_threshold=1.0,   # don't let squidpy filter further
            alpha=alpha,
        )
    except Exception as e:
        logger.warning("[downstream_report] LR focused dotplot failed: %s", e)
        return None


def _section_ligrec(adata: ad.AnnData, prov: dict, dominant_celltype_key: str) -> str:
    info = prov.get("analyses", {}).get("ligrec", {})
    if info.get("skipped"):
        return _skip_section("Ligand-Receptor Communication", info.get("reason", "not run"))

    ligrec_key = f"{dominant_celltype_key}_ligrec"
    if ligrec_key not in adata.uns:
        return _skip_section("Ligand-Receptor Communication", f"uns['{ligrec_key}'] not found")

    # Deserialize into a temporary local copy — never mutate adata.uns.
    raw = adata.uns[ligrec_key]
    if any(k.endswith("_data") for k in raw):
        import copy as _copy
        adata = _copy.copy(adata)
        adata.uns = dict(adata.uns)
        adata.uns[ligrec_key] = _deserialize_ligrec(raw)

    ligrec_data = adata.uns[ligrec_key]
    means_df = ligrec_data.get("means")
    pvals_df = ligrec_data.get("pvalues")

    # ── Summary note ──────────────────────────────────────────────────────────
    summary_html = ""
    if means_df is not None and pvals_df is not None:
        n_sig_001 = int((pvals_df < 0.001).values.sum())
        n_sig_005 = int((pvals_df < 0.05).values.sum())
        n_pairs   = pvals_df.shape[1]
        summary_html = (
            f'<p class="note">'
            f'{n_sig_001:,} interactions at p&nbsp;&lt;&nbsp;0.001 '
            f'({n_sig_005:,} at p&nbsp;&lt;&nbsp;0.05) '
            f'across {n_pairs} cell-type pairs.</p>'
        )

    # ── Plot 1: ranked bar chart (primary, always shown) ──────────────────────
    bar_b64 = None
    if means_df is not None and pvals_df is not None:
        bar_b64 = _lr_bar_chart(means_df, pvals_df, n_top=20, alpha=0.001)

    # ── Plot 2: focused squidpy dotplot (detail, lightbox-friendly) ───────────
    dot_b64 = None
    if means_df is not None and pvals_df is not None and _SQUIDPY_AVAILABLE:
        dot_b64 = _lr_focused_dotplot(
            adata, ligrec_key, dominant_celltype_key,
            pvals_df, n_cell_types=6, n_top_pairs=40, alpha=0.001,
        )

    # ── Assemble HTML ─────────────────────────────────────────────────────────
    figs_html = ""
    if bar_b64 or dot_b64:
        bar_wrap = dot_wrap = ""
        if bar_b64:
            bar_wrap = (
                '<div class="fig-wrap">'
                '<h3>Top 20 LR pairs by significance</h3>'
                + _img_tag(bar_b64, "LR ranked bar chart")
                + "</div>"
            )
        if dot_b64:
            dot_wrap = (
                '<div class="fig-wrap">'
                '<h3>Focused dotplot — top 6 cell types (click to expand)</h3>'
                + _img_tag(dot_b64, "LR focused dotplot")
                + "</div>"
            )
        figs_html = f'<div class="fig-grid">{bar_wrap}{dot_wrap}</div>'

    return f"""
    <section>
      <h2>Ligand-Receptor Communication</h2>
      <p>Permutation test (CellPhoneDB-like) for ligand-receptor interactions between
         spatially co-localised cell types, using the OmniPath database.
         Bar chart shows the top 20 pairs ranked by &minus;log&#8321;&#8320;(p).
         Dotplot is restricted to the 6 cell types with the most significant interactions.</p>
      {summary_html}
      {figs_html}
    </section>
    """


def _section_svg_gsea(adata: ad.AnnData, prov: dict) -> str:
    info = prov.get("analyses", {}).get("svg_gsea", {})
    if info.get("skipped"):
        return _skip_section("SVG Pathway Enrichment", info.get("reason", "not run"))

    gsea_df: pd.DataFrame = adata.uns.get("svg_gsea")
    if gsea_df is None or len(gsea_df) == 0:
        return _skip_section("SVG Pathway Enrichment", "no GSEA results")

    # Defensive coercion: gseapy may return numeric columns as object dtype.
    # nlargest/nsmallest require a real numeric dtype — coerce before sorting.
    _FLOAT_COLS = {"ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val"}
    gsea_df = gsea_df.copy()
    for _col in _FLOAT_COLS:
        if _col in gsea_df.columns:
            gsea_df[_col] = pd.to_numeric(gsea_df[_col], errors="coerce")

    # Bar chart of top 10 pathways by NES
    b64 = None
    try:
        nes_col = next(
            (c for c in ["NES", "nes"] if c in gsea_df.columns), None
        )
        term_col = next(
            (c for c in ["Term", "term"] if c in gsea_df.columns), None
        )
        if nes_col and term_col:
            top10 = gsea_df.nlargest(10, nes_col).copy()
            labels = [t[:45] + "…" if len(t) > 45 else t for t in top10[term_col].tolist()]
            nes_vals = top10[nes_col].values

            fig, ax = plt.subplots(figsize=(7, 4))
            colors = ["#e74c3c" if v > 0 else "#3498db" for v in nes_vals]
            ax.barh(range(len(labels)), nes_vals, color=colors, height=0.6)
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=8)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_xlabel("Normalised Enrichment Score (NES)")
            ax.set_title("Top 10 enriched pathways (SVG-ranked GSEA)")
            ax.invert_yaxis()
            try:
                fig.tight_layout()
            except Exception:
                pass
            b64 = _fig_to_b64(fig)
    except Exception as e:
        logger.warning("[downstream_report] GSEA bar chart failed: %s", e)

    # Table of top 20 results
    table_html = ""
    try:
        fdr_col = next(
            (c for c in ["FDR q-val", "fdr"] if c in gsea_df.columns), None
        )
        nes_col = next(
            (c for c in ["NES", "nes"] if c in gsea_df.columns), None
        )
        term_col = next(
            (c for c in ["Term", "term"] if c in gsea_df.columns), None
        )
        if all(c is not None for c in [fdr_col, nes_col, term_col]):
            top20 = gsea_df.nsmallest(20, fdr_col)
            rows = ""
            for _, row in top20.iterrows():
                fdr_val = float(row[fdr_col])
                sig = (
                    '<span style="color:#e74c3c;font-weight:700;">★</span>'
                    if fdr_val < 0.05 else ""
                )
                rows += (
                    f"<tr><td>{row[term_col]}</td>"
                    f"<td>{float(row[nes_col]):.3f}</td>"
                    f"<td>{fdr_val:.3e} {sig}</td></tr>"
                )
            table_html = f"""
            <h3>Top 20 enriched pathways</h3>
            <table>
              <thead><tr><th>Pathway</th><th>NES</th><th>FDR q-val</th></tr></thead>
              <tbody>{rows}</tbody>
            </table>
            <p class="note">&#9733; = FDR &lt; 0.05 &nbsp;&middot;&nbsp;
              Gene sets: <code>{info.get('gene_sets', '?')}</code> &nbsp;&middot;&nbsp;
              {info.get('n_significant', 0)} significant of {info.get('n_pathways', '?')} tested</p>
            """
    except Exception as e:
        logger.warning("[downstream_report] GSEA table failed: %s", e)

    fig_html = ""
    if b64:
        fig_html = (
            '<div class="fig-grid">'
            '<div class="fig-wrap">'
            + _img_tag(b64, "SVG GSEA bar chart")
            + "</div></div>"
        )

    return f"""
    <section>
      <h2>SVG Pathway Enrichment (GSEA)</h2>
      <p>Pre-ranked GSEA using Moran's I scores as the gene ranking.
         Genes with high spatial autocorrelation are ranked first.
         Red bars indicate positively enriched pathways; blue = depleted.</p>
      {fig_html}
      {table_html}
    </section>
    """


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_spatial_downstream_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
    dominant_celltype_key: str = "dominant_cell_type",
) -> str:
    """Generate a self-contained HTML downstream analysis report.

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_downstream`.
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label shown in the report header.
    dominant_celltype_key
        ``obs`` column with dominant cell type per spot.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    if "omicsage_spatial_downstream" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_downstream'] not found. "
            "Run spatial_downstream() before generating the report."
        )

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    prov = adata.uns["omicsage_spatial_downstream"]

    sections = [
        _section_summary(adata, prov, dataset_id, timestamp),
        _section_region_clustering(adata, prov),
        _section_celltype_markers(adata, prov),
        _section_celltype_svg(adata, prov),
        _section_co_occurrence(adata, prov, dominant_celltype_key),
        _section_nhood_enrichment(adata, prov, dominant_celltype_key),
        _section_ligrec(adata, prov, dominant_celltype_key),
        _section_svg_gsea(adata, prov),
    ]

    html = _render_page(
        title=f"OmicSage -- Spatial Downstream -- {dataset_id}",
        header_subtitle=f"Dataset: {dataset_id}",
        sections=sections,
        timestamp=timestamp,
    )

    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    size_kb = Path(output_path).stat().st_size / 1024
    logger.info("Spatial downstream report -> %s (%.1f KB)", output_path, size_kb)
    return str(Path(output_path).resolve())
