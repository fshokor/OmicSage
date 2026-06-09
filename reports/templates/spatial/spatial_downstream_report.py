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
    """Call a sq.pl.* function by name and return base64 PNG or None."""
    if not _SQUIDPY_AVAILABLE:
        return None
    try:
        fn = getattr(sq.pl, plot_fn)
        result = fn(adata, figsize=figsize, **kwargs)
        # sq.pl.* may return axes, list of axes, or None (writes directly to gcf)
        if result is None:
            fig = plt.gcf()
        elif hasattr(result, "__len__") and not hasattr(result, "get_figure"):
            first = result[0] if len(result) > 0 else None
            fig = first.get_figure() if first is not None else plt.gcf()
        elif hasattr(result, "get_figure"):
            fig = result.get_figure()
        else:
            fig = plt.gcf()
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
    figsize: tuple = (6, 5),
    library_key: Optional[str] = None,
) -> Optional[str]:
    """Spatial scatter coloured by an obs column."""
    if not _SQUIDPY_AVAILABLE or "spatial" not in adata.obsm:
        return None
    try:
        kwargs: dict = {
            "color": color,
            "frameon": False,
            "return_ax": True,
            "figsize": figsize,
        }
        if library_key:
            kwargs["library_key"] = library_key
        axes = sq.pl.spatial_scatter(adata, **kwargs)
        if axes is None:
            fig = plt.gcf()
        elif hasattr(axes, "__len__") and not hasattr(axes, "get_figure"):
            first = axes[0] if len(axes) > 0 else None
            fig = first.get_figure() if first is not None else plt.gcf()
        else:
            fig = axes.get_figure()
        try:
            fig.tight_layout()
        except Exception:
            pass
        return _fig_to_b64(fig)
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
            b64 = _squidpy_fig_b64(
                adata,
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
      {fig_html if fig_html else '<p class="note">Co-occurrence data computed. Install squidpy for visualisation.</p>'}
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
        b64 = _squidpy_fig_b64(
            adata,
            "nhood_enrichment",
            figsize=(7, 6),
            cluster_key=dominant_celltype_key,
            method="average",
        )

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
      {fig_html if fig_html else '<p class="note">Enrichment data computed. Install squidpy for visualisation.</p>'}
    </section>
    """


def _section_ligrec(adata: ad.AnnData, prov: dict, dominant_celltype_key: str) -> str:
    info = prov.get("analyses", {}).get("ligrec", {})
    if info.get("skipped"):
        return _skip_section("Ligand-Receptor Communication", info.get("reason", "not run"))

    ligrec_key = f"{dominant_celltype_key}_ligrec"
    if ligrec_key not in adata.uns:
        return _skip_section("Ligand-Receptor Communication", f"uns['{ligrec_key}'] not found")

    # Attempt to render the squidpy dotplot
    b64 = None
    if _SQUIDPY_AVAILABLE:
        try:
            b64 = _squidpy_fig_b64(
                adata,
                "ligrec",
                figsize=(10, 6),
                cluster_key=dominant_celltype_key,
                pvalue_threshold=0.05,
                alpha=0.001,
            )
        except Exception as e:
            logger.warning("[downstream_report] ligrec dotplot failed: %s", e)

    # Fallback: summary table from raw means/pvalues
    table_html = ""
    try:
        ligrec_data = adata.uns[ligrec_key]
        means_df = ligrec_data.get("means")
        pvals_df = ligrec_data.get("pvalues")

        if means_df is not None and pvals_df is not None:
            # Multi-index columns: (source_cluster, target_cluster)
            # Flatten to find top significant interactions
            n_sig = int((pvals_df < 0.05).values.sum())
            table_html = (
                f'<p class="note">'
                f'{n_sig:,} significant ligand-receptor interactions at p &lt; 0.05 '
                f'across {pvals_df.shape[1]} cell type pairs.</p>'
            )
    except Exception as e:
        logger.warning("[downstream_report] ligrec table failed: %s", e)

    fig_html = ""
    if b64:
        fig_html = (
            '<div class="fig-grid">'
            '<div class="fig-wrap" style="max-width:900px;"><h3>Significant LR pairs (p &lt; 0.001)</h3>'
            + _img_tag(b64, "ligand-receptor dotplot")
            + "</div></div>"
        )

    return f"""
    <section>
      <h2>Ligand-Receptor Communication</h2>
      <p>Permutation test (CellPhoneDB-like) for ligand-receptor interactions between
         spatially co-localised cell types, using the OmniPath database.</p>
      {table_html}
      {fig_html if fig_html else ""}
    </section>
    """


def _section_svg_gsea(adata: ad.AnnData, prov: dict) -> str:
    info = prov.get("analyses", {}).get("svg_gsea", {})
    if info.get("skipped"):
        return _skip_section("SVG Pathway Enrichment", info.get("reason", "not run"))

    gsea_df: pd.DataFrame = adata.uns.get("svg_gsea")
    if gsea_df is None or len(gsea_df) == 0:
        return _skip_section("SVG Pathway Enrichment", "no GSEA results")

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
