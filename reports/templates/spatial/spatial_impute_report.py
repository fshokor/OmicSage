"""
spatial_impute_report.py — OmicSage Phase 7 extension

HTML report for the spatial imputation step.
Style: identical _PAGE_CSS / _render_page pattern as all other spatial reports.

Sections
--------
1. Run Summary         — stat cards: n_genes_imputed, method, mean mapping
                         score, n_spots, sc_reference filename
2. Mapping Score Dist  — histogram of per-spot Tangram mapping scores
                         (spots < 0.1 highlighted; not shown for gimVI)
3. Top Imputed Genes   — spatial scatter of top 5 imputed genes by variance
                         (the visual payoff — genes not in original HVG set)
4. Imputation Validation — scatter: measured vs imputed for genes present in
                           both datasets, Spearman r in title
"""

from __future__ import annotations

import base64
import io
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
import scipy.sparse as sp
from scipy.stats import spearmanr

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
    .fig-wrap { flex: 1 1 300px; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6;
                    cursor: zoom-in; transition: box-shadow 0.15s; display: block; }
    .fig-wrap img:hover { box-shadow: 0 4px 16px rgba(15,52,96,0.18); }
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
        "    <h1>OmicSage &#8212; Spatial Imputation Report</h1>\n"
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


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------

def _section_run_summary(
    prov: dict,
    sc_ref_label: str,
    timestamp: str,
) -> str:
    out = prov.get("outputs", {})
    method          = prov.get("method", "—")
    n_genes         = out.get("n_genes_imputed", 0)
    n_spots         = out.get("n_spots", 0)
    mean_score      = out.get("mean_mapping_score", float("nan"))
    n_poor          = out.get("n_poor_spots", 0)

    score_str = f"{mean_score:.3f}" if not np.isnan(mean_score) else "N/A"
    poor_str  = str(n_poor) if method == "tangram" else "N/A"

    cards = (
        _stat_card(str(n_genes), "Genes imputed")
        + _stat_card(method.capitalize(), "Method")
        + _stat_card(score_str, "Mean mapping score")
        + _stat_card(f"{n_spots:,}", "Spots")
        + _stat_card(poor_str, "Poor-score spots (<0.1)")
    )

    note = ""
    if method == "tangram":
        note = (
            '<p class="note">Mapping score interpretation: values ≥ 0.1 indicate '
            'reliable sc→spot assignment. Spots below 0.1 may be in tissue regions '
            'not well-represented by the scRNA-seq reference.</p>'
        )

    return (
        f'<section>'
        f'<h2>1. Run Summary</h2>'
        f'<p class="timestamp">Generated {timestamp}</p>'
        f'<div class="stat-grid">{cards}</div>'
        f'<p><strong>Method:</strong> <code>{method}</code> &nbsp; '
        f'<strong>SC reference:</strong> <code>{sc_ref_label}</code></p>'
        f'{note}'
        f'</section>'
    )


def _section_mapping_score_hist(adata: ad.AnnData, prov: dict) -> str:
    method       = prov.get("method", "tangram")
    tangram_mode = prov.get("outputs", {}).get("tangram_mode", "clusters")

    if method != "tangram":
        return _skip_section(
            "2. Mapping Score Distribution",
            "Not available for gimVI (no per-spot score).",
        )

    if tangram_mode == "clusters":
        return (
            "<section><h2>2. Mapping Score Distribution</h2>"
            "<p>Per-spot mapping scores are not available in <code>clusters</code> mode. "
            "In clusters mode Tangram maps per-cell-type signatures onto spots rather than "
            "individual cells, so there is no per-spot quality score. "
            "To obtain per-spot mapping scores, set <code>tangram_mode: cells</code> in your config "
            "(note: cells mode requires significantly more memory).</p></section>"
        )

    if "tangram_mapping_score" not in adata.obs.columns:
        return _skip_section(
            "2. Mapping Score Distribution",
            "tangram_mapping_score not found in adata.obs.",
        )

    scores = adata.obs["tangram_mapping_score"].values.astype(float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(scores, bins=40, color="#4C78A8", edgecolor="white", linewidth=0.4)
    ax.axvline(0.1, color="#e05c5c", linestyle="--", linewidth=1.4,
               label="Poor threshold (0.1)")
    ax.set_xlabel("Tangram mapping score", fontsize=11)
    ax.set_ylabel("Number of spots", fontsize=11)
    ax.set_title("Per-spot mapping score distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    b64 = _fig_to_b64(fig)

    n_poor = int((scores < 0.1).sum())
    pct_poor = 100 * n_poor / len(scores) if len(scores) > 0 else 0

    return (
        f'<section>'
        f'<h2>2. Mapping Score Distribution</h2>'
        f'<p>Per-spot Tangram mapping score: how well each Visium spot is '
        f'explained by the scRNA-seq reference cells. '
        f'<strong>{n_poor} spots ({pct_poor:.1f}%)</strong> score below 0.1 '
        f'(red dashed line) and may have unreliable imputed values.</p>'
        f'<div class="fig-grid">'
        f'<div class="fig-wrap"><h3>Mapping score histogram</h3>'
        f'{_img_tag(b64, "Mapping score distribution")}</div>'
        f'</div>'
        f'</section>'
    )


def _section_top_imputed_genes(adata: ad.AnnData, prov: dict, library_key: str = None) -> str:
    """Spatial scatter of top 5 imputed genes by variance."""
    if "imputed_expression" not in adata.obsm:
        return _skip_section(
            "3. Top Imputed Genes on H&E",
            "imputed_expression not found in adata.obsm.",
        )

    raw = adata.obsm["imputed_expression"]
    gene_names = adata.uns.get("omicsage_spatial_impute", {}).get(
        "outputs", {}
    ).get("genes_imputed", None)

    if isinstance(raw, np.ndarray):
        if gene_names is None or len(gene_names) != raw.shape[1]:
            return _skip_section(
                "3. Top Imputed Genes on H&E",
                "Gene names missing from uns — cannot reconstruct imputed DataFrame.",
            )
        imputed = pd.DataFrame(raw, index=adata.obs_names, columns=gene_names)
    elif isinstance(raw, pd.DataFrame):
        imputed = raw
    else:
        return _skip_section(
            "3. Top Imputed Genes on H&E",
            f"Unexpected type for imputed_expression: {type(raw).__name__}.",
        )

    if imputed.shape[1] == 0:
        return _skip_section("3. Top Imputed Genes on H&E", "Imputed expression is empty.")

    # Top 5 genes by variance in imputed values
    variances = imputed.var(axis=0)
    top5 = variances.nlargest(5).index.tolist()

    if not _SQUIDPY_AVAILABLE:
        return _skip_section(
            "3. Top Imputed Genes on H&E",
            "squidpy not available — cannot render spatial scatter.",
        )

    # Add imputed expression to adata.obs temporarily for sq.pl.spatial_scatter
    figs_html = ""
    for gene in top5:
        key = f"_imputed_{gene}"
        adata.obs[key] = imputed[gene].values if len(imputed) == adata.n_obs else np.nan
        try:
            fig, ax = plt.subplots(figsize=(5, 5))
            import squidpy as sq
            scatter_kwargs = dict(
                color=key,
                ax=ax,
                title=f"Imputed: {gene}",
                colormap="magma",
                size=1.2,
                show=False,
            )
            if library_key:
                scatter_kwargs["library_key"] = library_key
            sq.pl.spatial_scatter(adata, **scatter_kwargs)
            ax.set_title(f"Imputed: {gene}", fontsize=10, fontweight="bold")
            b64 = _fig_to_b64(fig)
            figs_html += (
                f'<div class="fig-wrap"><h3>{gene}</h3>'
                f'{_img_tag(b64, gene)}</div>'
            )
        except Exception as exc:
            logger.warning(f"Could not plot imputed gene {gene}: {exc}", exc_info=True)
            plt.close("all")
        finally:
            if key in adata.obs.columns:
                del adata.obs[key]

    if not figs_html:
        return _skip_section("3. Top Imputed Genes on H&E", "Could not render spatial plots.")

    return (
        f'<section>'
        f'<h2>3. Top Imputed Genes on Tissue</h2>'
        f'<p>Top 5 genes by variance in imputed expression — these are the '
        f'genes not in the original Visium panel whose spatial expression '
        f'patterns are predicted from the scRNA-seq reference. Brighter spots '
        f'indicate higher predicted expression.</p>'
        f'<div class="fig-grid">{figs_html}</div>'
        f'</section>'
    )


def _section_validation(adata: ad.AnnData, prov: dict) -> str:
    """Scatter of measured vs imputed expression for overlapping genes."""
    if "imputed_expression" not in adata.obsm:
        return _skip_section(
            "4. Imputation Validation",
            "imputed_expression not found in adata.obsm.",
        )

    raw = adata.obsm["imputed_expression"]
    gene_names = adata.uns.get("omicsage_spatial_impute", {}).get(
        "outputs", {}
    ).get("genes_imputed", None)

    if isinstance(raw, np.ndarray):
        if gene_names is None or len(gene_names) != raw.shape[1]:
            return _skip_section(
                "4. Imputation Validation",
                "Gene names missing from uns — cannot reconstruct imputed DataFrame.",
            )
        imputed = pd.DataFrame(raw, index=adata.obs_names, columns=gene_names)
    elif isinstance(raw, pd.DataFrame):
        imputed = raw
    else:
        return _skip_section(
            "4. Imputation Validation",
            f"Unexpected type for imputed_expression: {type(raw).__name__}.",
        )

    # Find genes present in both spatial measured data and imputed DataFrame
    measured_genes = set(adata.var_names)
    imputed_genes  = set(imputed.columns)
    overlap = list(measured_genes & imputed_genes)

    if len(overlap) < 5:
        return _skip_section(
            "4. Imputation Validation",
            f"Only {len(overlap)} overlapping genes between measured and "
            "imputed data — insufficient for validation scatter.",
        )

    # Use up to 50 genes for the scatter to keep it readable
    np.random.seed(42)
    sample_genes = overlap if len(overlap) <= 50 else list(
        np.random.choice(overlap, 50, replace=False)
    )

    # Log-normalise measured counts before comparing to imputed values.
    # Tangram imputed expression is already on a normalised scale, so
    # comparing raw counts vs normalised imputed causes artificially low
    # Spearman r. We use log1p(counts / sum * 10000) to match the typical
    # normalisation applied before Tangram training.
    x_arr = adata[:, sample_genes].X
    if sp.issparse(x_arr):
        x_arr = x_arr.toarray()
    x_arr = np.asarray(x_arr, dtype=np.float64)
    lib_sizes = x_arr.sum(axis=1, keepdims=True)
    lib_sizes[lib_sizes == 0] = 1  # avoid divide-by-zero
    x_norm = np.log1p(x_arr / lib_sizes * 10000)
    measured_mean = x_norm.mean(axis=0)
    imputed_mean  = imputed[sample_genes].values.mean(axis=0)

    rho, pval = spearmanr(measured_mean, imputed_mean)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.scatter(measured_mean, imputed_mean, alpha=0.65, s=22,
               color="#4C78A8", edgecolors="none")
    ax.set_xlabel("Mean measured expression (log-normalised)", fontsize=11)
    ax.set_ylabel("Mean imputed expression", fontsize=11)
    ax.set_title(
        f"Measured vs imputed (n={len(sample_genes)} genes)\n"
        f"Spearman r = {rho:.3f}  (p = {pval:.2e})",
        fontsize=10, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    b64 = _fig_to_b64(fig)

    quality_note = ""
    if rho >= 0.7:
        quality_note = (
            '<p class="note" style="border-left-color:#54a868; background:#f0fff4; color:#1a5c2a;">'
            f'Good imputation quality: Spearman r = {rho:.3f}. '
            'Imputed values track measured expression well.</p>'
        )
    elif rho >= 0.4:
        quality_note = (
            '<p class="note">'
            f'Moderate imputation quality: Spearman r = {rho:.3f}. '
            'Imputation may be less reliable for lowly-expressed genes.</p>'
        )
    else:
        quality_note = (
            '<p class="note" style="border-left-color:#e05c5c; background:#fff5f5; color:#7a1a1a;">'
            f'Low imputation correlation: Spearman r = {rho:.3f}. '
            'Check that the scRNA-seq reference is from a matched cell type '
            'composition and that gene IDs are consistent.</p>'
        )

    return (
        f'<section>'
        f'<h2>4. Imputation Validation</h2>'
        f'<p>Scatter of mean measured vs mean imputed expression across '
        f'{len(sample_genes)} genes shared between the spatial panel and the '
        f'imputed gene set. High Spearman correlation indicates that the '
        f'imputation model has learned biologically realistic gene relationships.</p>'
        f'{quality_note}'
        f'<div class="fig-grid">'
        f'<div class="fig-wrap"><h3>Measured vs imputed mean expression</h3>'
        f'{_img_tag(b64, "Validation scatter")}</div>'
        f'</div>'
        f'</section>'
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def generate_spatial_impute_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
    sc_ref_label: str = "",
) -> str:
    """Generate the spatial imputation HTML report.

    Parameters
    ----------
    adata
        AnnData post-imputation (contains ``uns["omicsage_spatial_impute"]``).
    output_path
        Path to write the HTML report.
    dataset_id
        Human-readable dataset name for the report header.
    sc_ref_label
        Filename or description of the scRNA-seq reference, shown in Run Summary.

    Returns
    -------
    str
        Absolute path to the written report.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    prov      = adata.uns.get("omicsage_spatial_impute", {})

    if not sc_ref_label:
        # Try to derive from provenance or checkpoint path
        sc_ref_label = prov.get("sc_reference", "unknown")

    if prov.get("skipped"):
        reason = prov.get("skip_reason", "unknown")
        html = _render_page(
            title=f"OmicSage — Spatial Imputation — {dataset_id}",
            header_subtitle=dataset_id,
            sections=[
                f'<section><h2>Imputation Skipped</h2>'
                f'<p class="skip-note">{reason}</p></section>'
            ],
            timestamp=timestamp,
        )
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(html, encoding="utf-8")
        return str(Path(output_path).resolve())

    # Extract library_key from ingest provenance — required by squidpy
    # sq.pl.spatial_scatter when the AnnData contains multiple library IDs.
    library_key = (
        adata.uns
        .get("omicsage_spatial_ingest", {})
        .get("library_key", None)
    )

    sections = [
        _section_run_summary(prov, sc_ref_label, timestamp),
        _section_mapping_score_hist(adata, prov),
        _section_top_imputed_genes(adata, prov, library_key=library_key),
        _section_validation(adata, prov),
    ]

    html = _render_page(
        title=f"OmicSage — Spatial Imputation — {dataset_id}",
        header_subtitle=dataset_id,
        sections=sections,
        timestamp=timestamp,
    )

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    logger.info(f"[spatial_impute_report] -> {output_path}")
    return str(Path(output_path).resolve())


if __name__ == "__main__":
    import argparse
    import scanpy as sc

    parser = argparse.ArgumentParser()
    parser.add_argument("--adata",   required=True, help="Path to imputed .h5ad")
    parser.add_argument("--output",  required=True, help="Output HTML path")
    parser.add_argument("--dataset", default="spatial")
    parser.add_argument("--sc-ref",  default="", help="SC reference label")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata)
    generate_spatial_impute_report(
        adata, args.output, dataset_id=args.dataset, sc_ref_label=args.sc_ref
    )
