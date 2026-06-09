"""
spatial_deconvolve_report.py — OmicSage Phase 7, Session 5 (multi-method)
HTML report for spatial deconvolution results.

Works with both deconvolution methods:
  - method="nnls"          : proportions, sum to 1 per spot
  - method="cell2location" : absolute abundances

The report reads provenance from uns["omicsage_spatial_deconvolve"]:
  - outputs["method"], outputs["per_sample"], outputs["library_key"]
  - outputs["cell_type_names"], outputs["n_spots"], outputs["n_shared_genes"]
  - skipped + skip_reason for graceful-skip mode

Style matches spatial_qc_report.py (_PAGE_CSS / _render_page pattern).
Spatial figures use the canonical squidpy API via _squidpy_scatter_b64.
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

matplotlib.use("Agg")
logger = logging.getLogger(__name__)

try:
    import squidpy as sq
    _SQUIDPY_AVAILABLE = True
except ImportError:
    _SQUIDPY_AVAILABLE = False


# ---------------------------------------------------------------------------
# Squidpy + figure helpers (matches the other spatial reports)
# ---------------------------------------------------------------------------


def _fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _resolve_img_key(adata: ad.AnnData, img_key: Optional[str]) -> Optional[str]:
    """Return img_key if the image exists in uns['spatial'], else try aliases, else None."""
    if img_key is None:
        return None
    spatial_uns = adata.uns.get("spatial", {})
    for _sample, sample_data in spatial_uns.items():
        images = sample_data.get("images", {})
        if img_key in images:
            return img_key
        for alias in ("hires", "lowres"):
            if alias in images:
                return alias
    return None


def _get_library_key(adata: ad.AnnData) -> Optional[str]:
    """Return the obs column that maps spots to their library/sample ID.

    Reading order (first match wins):
    1. uns["omicsage_spatial_ingest"]["library_key"] -- set at ingest time
    2. Auto-detection fallback by matching uns['spatial'] keys vs obs columns
    Returns None for single-sample data (squidpy does not need library_key).
    """
    ingest_prov = adata.uns.get("omicsage_spatial_ingest", {})
    if "library_key" in ingest_prov and ingest_prov["library_key"] is not None:
        return ingest_prov["library_key"]

    library_ids = list(adata.uns.get("spatial", {}).keys())
    if len(library_ids) <= 1:
        return None
    id_set = set(library_ids)
    for candidate in ("library_id", "sample", "patient", "donor_id", "batch", "slide"):
        if candidate in adata.obs.columns:
            if id_set.issubset(set(adata.obs[candidate].astype(str).unique())):
                return candidate
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name in ("object", "category"):
            if id_set.issubset(set(adata.obs[col].astype(str).unique())):
                return col
    return None


def _squidpy_scatter_b64(
    adata,
    color,
    img_key=None,
    library_key=None,
    cmap=None,
    title=None,
    figsize=(6, 5),
):
    """Render sq.pl.spatial_scatter and return a base64 PNG string.

    Uses the canonical squidpy API:
      - ``img_res_key`` (NOT ``img_key`` -- the docstring shorthand is misleading)
      - ``figsize`` is top-level
      - ``return_ax=True`` returns axes so we can grab the figure cleanly
      - NO ``show=`` parameter exists; passing it crashes inside matplotlib
      - When ``library_key`` is set, squidpy makes one panel per library and
        owns the figure entirely.
    """
    kwargs = dict(color=color, frameon=False, return_ax=True, figsize=figsize)
    if img_key:
        kwargs["img_res_key"] = img_key
    if library_key:
        kwargs["library_key"] = library_key
    if cmap:
        kwargs["cmap"] = cmap

    axes = sq.pl.spatial_scatter(adata, **kwargs)

    if axes is None:
        fig = plt.gcf()
    elif hasattr(axes, "__len__") and not hasattr(axes, "get_figure"):
        first = axes[0] if len(axes) > 0 else None
        fig = first.get_figure() if first is not None else plt.gcf()
    else:
        fig = axes.get_figure()

    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)
    try:
        fig.tight_layout()
    except Exception:
        pass
    return _fig_to_b64(fig)


def _stitch_panels(b64_panels: list, ncols: int = 3) -> Optional[str]:
    """Combine N base64 PNGs into one image grid. Returns the stitched base64."""
    if not b64_panels:
        return None
    if len(b64_panels) == 1:
        return b64_panels[0]
    imgs = [plt.imread(io.BytesIO(base64.b64decode(b))) for b in b64_panels]
    ncols = min(len(imgs), ncols)
    nrows = (len(imgs) + ncols - 1) // ncols
    max_h = max(a.shape[0] for a in imgs)
    max_w = max(a.shape[1] for a in imgs)
    ch    = imgs[0].shape[2] if imgs[0].ndim == 3 else 1

    def _pad(a):
        ph = max_h - a.shape[0]
        pw = max_w - a.shape[1]
        return np.pad(a, ((0, ph), (0, pw), (0, 0)),
                      mode="constant", constant_values=1)

    padded = [_pad(a) for a in imgs]
    while len(padded) < nrows * ncols:
        padded.append(np.ones((max_h, max_w, ch), dtype=padded[0].dtype))

    rows = [np.hstack(padded[r * ncols:(r + 1) * ncols]) for r in range(nrows)]
    grid = np.vstack(rows)
    fig, ax = plt.subplots(figsize=(grid.shape[1] / 100, grid.shape[0] / 100))
    ax.imshow(grid); ax.axis("off")
    fig.tight_layout(pad=0)
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Individual figure functions
# ---------------------------------------------------------------------------


def _fig_mean_abundance_bar(adata: ad.AnnData, cell_type_names: list) -> Optional[str]:
    """Horizontal bar of mean abundance per cell type, sorted ascending."""
    try:
        available = [ct for ct in cell_type_names if ct in adata.obs.columns]
        if not available:
            return None
        means = {ct: float(adata.obs[ct].mean()) for ct in available}
        sorted_cts = sorted(means, key=means.get)
        values = [means[ct] for ct in sorted_cts]

        fig, ax = plt.subplots(figsize=(6, max(3, len(sorted_cts) * 0.4)))
        colors = plt.cm.tab20(np.linspace(0, 1, len(sorted_cts)))
        ax.barh(sorted_cts, values, color=colors, edgecolor="white")
        ax.set_xlabel("Mean per-spot abundance / proportion")
        ax.set_title("Mean cell-type abundance across all spots",
                     fontsize=10, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        return _fig_to_b64(fig)
    except Exception as e:
        logger.warning("mean abundance bar failed: %s", e)
        return None


def _fig_dominant_celltype(
    adata: ad.AnnData,
    img_key: Optional[str],
) -> Optional[str]:
    """Spatial scatter coloured by dominant cell type per spot."""
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm:
        return None
    if "dominant_cell_type" not in adata.obs.columns:
        return None
    resolved    = _resolve_img_key(adata, img_key)
    library_key = _get_library_key(adata)
    suffix = " (on H&E)" if resolved else ""
    try:
        return _squidpy_scatter_b64(
            adata,
            color="dominant_cell_type",
            img_key=resolved,
            library_key=library_key,
            title=f"Dominant cell type per spot{suffix}",
            figsize=(6, 6),
        )
    except Exception as e:
        logger.warning("dominant cell type figure failed: %s", e)
        return None


def _fig_top_celltypes_spatial(
    adata: ad.AnnData,
    cell_type_names: list,
    img_key: Optional[str],
    n_top: int = 6,
) -> Optional[str]:
    """Spatial scatter for top-N cell types by mean abundance, in a grid."""
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm:
        return None

    available = [ct for ct in cell_type_names if ct in adata.obs.columns]
    if not available:
        return None

    means = {ct: float(adata.obs[ct].mean()) for ct in available}
    top   = sorted(means, key=means.get, reverse=True)[:n_top]

    resolved    = _resolve_img_key(adata, img_key)
    library_key = _get_library_key(adata)
    try:
        b64_panels = []
        for ct in top:
            b = _squidpy_scatter_b64(
                adata,
                color=ct,
                img_key=resolved,
                library_key=library_key,
                cmap="magma",
                title=ct,
                figsize=(5, 5),
            )
            if b:
                b64_panels.append(b)
        return _stitch_panels(b64_panels, ncols=3)
    except Exception as e:
        logger.warning("top cell types figure failed: %s", e)
        return None


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------


def _section_summary(
    adata: ad.AnnData,
    deconv_info: dict,
    dataset_id: str,
    timestamp: str,
) -> str:
    outputs     = deconv_info.get("outputs", {})
    skipped     = deconv_info.get("skipped", True)
    skip_reason = deconv_info.get("skip_reason", "")
    method      = outputs.get("method", "?")
    per_sample  = outputs.get("per_sample", False)
    library_key = outputs.get("library_key")
    n_spots     = outputs.get("n_spots", adata.n_obs)
    n_ct        = outputs.get("n_cell_types", 0)
    n_shared    = outputs.get("n_shared_genes", 0)

    if skipped:
        stat_cards = (
            f'<div class="stat-card"><div class="stat-value">{n_spots:,}</div>'
            f'<div class="stat-label">Spots</div></div>'
            f'<div class="stat-card"><div class="stat-value">Skipped</div>'
            f'<div class="stat-label">Deconvolution</div></div>'
        )
        skip_note = (
            f'<p class="note">&#x26A0; Deconvolution skipped &mdash; '
            f'{skip_reason}</p>'
        )
        return f"""
        <section>
          <h2>Run Summary</h2>
          <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
          {skip_note}
          <div class="stat-grid">{stat_cards}</div>
        </section>
        """

    method_label = {
        "nnls":          "NNLS (proportions, memory-safe)",
        "cell2location": "cell2location (Bayesian, gold standard)",
    }.get(method, method)

    cards = [
        ("Method",         method_label),
        ("Spots",          f"{n_spots:,}"),
        ("Cell types",     str(n_ct)),
        ("Shared genes",   f"{n_shared:,}"),
    ]
    if per_sample and library_key:
        cards.append(("Per-sample loop", f"by <code>{library_key}</code>"))
    elif per_sample:
        cards.append(("Per-sample loop", "yes"))

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in cards
    )

    ct_names  = outputs.get("cell_type_names", [])
    ct_badges = " ".join(f"<code>{ct}</code>" for ct in ct_names)

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      {"<h3>Cell Types</h3><p>" + ct_badges + "</p>" if ct_badges else ""}
    </section>
    """


def _section_params(deconv_info: dict) -> str:
    params = deconv_info.get("params", {})
    if not params:
        return ""

    # Render scalars first, then lists/dicts as their str representation
    def _fmt(v):
        if v is None:
            return "<em>None</em>"
        if isinstance(v, bool):
            return "yes" if v else "no"
        if isinstance(v, (list, tuple)):
            return ", ".join(str(x) for x in v) if v else "<em>None</em>"
        return str(v)

    rows = "".join(
        f"<tr><td><code>{k}</code></td><td>{_fmt(v)}</td></tr>"
        for k, v in params.items()
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


def _section_abundance(
    adata: ad.AnnData,
    deconv_info: dict,
    img_key: Optional[str],
) -> str:
    ct_names = deconv_info.get("outputs", {}).get("cell_type_names", [])
    if not ct_names:
        return ""

    b64_bar  = _fig_mean_abundance_bar(adata, ct_names)
    b64_dom  = _fig_dominant_celltype(adata, img_key)
    b64_top  = _fig_top_celltypes_spatial(adata, ct_names, img_key, n_top=6)

    if not (b64_bar or b64_dom or b64_top):
        return ""

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
    if b64_top:
        figs += (
            '<div class="fig-wrap" style="max-width:900px;">'
            '<h3>Top 6 cell types &mdash; spatial abundance</h3>'
            f'<img src="data:image/png;base64,{b64_top}" alt="top cell types"></div>'
        )

    return f"""
    <section>
      <h2>Cell Type Abundances</h2>
      <p>NNLS outputs proportions (sum to 1 per spot). cell2location outputs
         absolute abundances. The dominant-cell-type panel uses an argmax
         across cell types per spot &mdash; written to <code>obs["dominant_cell_type"]</code>.</p>
      <div class="fig-grid">{figs}</div>
    </section>
    """


# ---------------------------------------------------------------------------
# Shared CSS + page renderer
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


def _render_page(title: str, sections: list, timestamp: str) -> str:
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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_spatial_deconvolve_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
    img_key: Optional[str] = "hires",
) -> str:
    """Generate a self-contained HTML deconvolution report.

    Works for both NNLS and cell2location outputs (same provenance contract).

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_deconvolve`.
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label shown in the report header.
    img_key
        H&E image key in ``uns["spatial"][sample]["images"]``.
        ``"hires"`` (default) or ``"lowres"``. Figures degrade gracefully
        to spots-only when the image is absent.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
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
            _section_abundance(adata, deconv_info, img_key),
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
