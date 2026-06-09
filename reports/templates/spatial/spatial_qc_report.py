"""
spatial_qc_report.py — OmicSage Phase 7, Session 5 (image update)
HTML report for spatial QC results.
Matches the _render_page / section structure of the RNA pipeline reports.

Changes vs Session 4:
  - _spatial_scatter() now accepts img_key to overlay H&E image (tutorial standard)
  - New figure: tissue overview coloured by in_tissue / array position
  - generate_spatial_qc_report() gains img_key parameter (default "hires")
  - _resolve_img_key() helper for graceful fallback when image absent
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


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def _fig_to_b64(fig):
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return b64


def _squidpy_scatter_b64(
    adata,
    color,
    img_key=None,
    library_key=None,
    title=None,
    figsize=(6, 5),
):
    """Render ``sq.pl.spatial_scatter`` and return a base64 PNG string.

    Uses the canonical squidpy API (per scverse/squidpy source):
      - ``img_res_key`` (NOT ``img_key`` -- the docstring shorthand is misleading)
      - ``figsize`` is a top-level parameter
      - ``return_ax=True`` returns axes so we can grab the figure cleanly
      - NO ``show=`` parameter exists -- passing it crashes deep in matplotlib
      - When ``library_key`` is set, squidpy creates one panel per library
        automatically; we let it own the figure entirely.

    Returns
    -------
    str or None
        Base64-encoded PNG of the rendered figure, or None if rendering fails.
    """
    kwargs = dict(color=color, frameon=False, return_ax=True, figsize=figsize)
    if img_key:
        kwargs["img_res_key"] = img_key
    if library_key:
        kwargs["library_key"] = library_key

    axes = sq.pl.spatial_scatter(adata, **kwargs)

    # axes can be a single Axes or a Sequence[Axes] (multi-sample)
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


def _resolve_img_key(adata: ad.AnnData, img_key: Optional[str]) -> Optional[str]:
    """
    Return img_key if the image exists in adata.uns['spatial'], otherwise try
    common aliases ('hires', 'lowres'), otherwise return None so figures
    degrade gracefully to spots-only.
    """
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


# ---------------------------------------------------------------------------
# Figure functions
# ---------------------------------------------------------------------------

def _get_library_key(adata: ad.AnnData) -> Optional[str]:
    """Return the obs column that maps spots to their library/sample ID.

    Reading order (first match wins):
    1. ``uns["omicsage_spatial_ingest"]["library_key"]`` -- set by spatial_ingest()
       using either the user-supplied value or auto-detection at load time.
    2. Auto-detection fallback -- for h5ad files loaded outside OmicSage that
       lack the provenance key.

    Returns ``None`` for single-sample data (squidpy does not need library_key).
    """
    # 1. Read from ingest provenance (authoritative)
    ingest_prov = adata.uns.get("omicsage_spatial_ingest", {})
    if "library_key" in ingest_prov and ingest_prov["library_key"] is not None:
        return ingest_prov["library_key"]

    # 2. Auto-detection fallback
    library_ids = list(adata.uns.get("spatial", {}).keys())
    if len(library_ids) <= 1:
        return None
    id_set = set(library_ids)
    for candidate in ("library_id", "sample", "patient", "donor_id", "batch", "slide"):
        if candidate in adata.obs.columns:
            if id_set.issubset(set(adata.obs[candidate].astype(str).unique())):
                return candidate
    # Fallback: scan all string/category columns
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name in ("object", "category"):
            if id_set.issubset(set(adata.obs[col].astype(str).unique())):
                return col
    return None


def _violin_plots(adata, params):
    metrics = [
        ("total_counts",      "Total UMI counts / spot",
         params.get("min_counts"), params.get("max_counts")),
        ("n_genes_by_counts", "Genes detected / spot",
         params.get("min_genes"),  params.get("max_genes")),
        ("pct_counts_mt",     "MT gene % / spot",
         None,                     params.get("max_mt_pct")),
    ]
    available = [m for m in metrics if m[0] in adata.obs.columns]
    if not available:
        return None
    fig, axes = plt.subplots(1, len(available), figsize=(5 * len(available), 4))
    if len(available) == 1:
        axes = [axes]
    for ax, (col, label, lo, hi) in zip(axes, available):
        vals = adata.obs[col].values
        parts = ax.violinplot(vals, positions=[0], showmedians=True)
        for pc in parts["bodies"]:
            pc.set_facecolor("#4C78A8"); pc.set_alpha(0.7)
        if lo is not None:
            ax.axhline(lo, color="#e07b3a", linestyle="--", linewidth=1.2,
                       label=f"min={lo}")
        if hi is not None:
            ax.axhline(hi, color="#e05252", linestyle="--", linewidth=1.2,
                       label=f"max={hi}")
        ax.set_xticks([])
        ax.set_ylabel(label)
        ax.set_title(label, fontsize=10, fontweight="bold")
        if lo is not None or hi is not None:
            ax.legend(fontsize=8, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("QC metrics per spot", fontsize=11, fontweight="bold")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _spatial_scatter(adata, color_key, img_key: Optional[str] = None):
    """Spatial scatter coloured by color_key, delegating figure creation to squidpy."""
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None
    if color_key not in adata.obs.columns:
        return None
    resolved = _resolve_img_key(adata, img_key)
    library_key = _get_library_key(adata)
    title = color_key.replace("_", " ").title()
    if resolved:
        title += " (on H&E)"
    try:
        return _squidpy_scatter_b64(
            adata, color=color_key,
            img_key=resolved, library_key=library_key,
            title=title, figsize=(5, 5),
        )
    except Exception as e:
        logger.warning("figure failed (%s): %s", __name__, e)
        return None


def _tissue_overview(adata, img_key: Optional[str]) -> Optional[str]:
    """
    Two-panel overview of the tissue array:
      Left : spots coloured by in_tissue flag (shows capture area layout)
      Right: spots coloured by array_row (shows spatial gradient — QC check)
    This is the first plot in most Visium tutorials.
    Falls back gracefully when in_tissue / array_row columns are absent.
    """
    if not _SQUIDPY_AVAILABLE:
        return None
    if "spatial" not in adata.obsm or "spatial" not in adata.uns:
        return None

    panels = []
    for col in ("in_tissue", "array_row"):
        if col in adata.obs.columns:
            panels.append(col)

    if not panels:
        return None

    resolved = _resolve_img_key(adata, img_key)
    library_key = _get_library_key(adata)
    # Render one panel at a time so each gets a clean title, then stitch into
    # one row.  When library_key is set squidpy makes N subplots per call (one
    # per sample) — we collect those figures and concatenate them horizontally.
    try:
        b64s = []
        for col in panels:
            title = col.replace("_", " ").title()
            if resolved:
                title += " (on H&E)"
            b64 = _squidpy_scatter_b64(
                adata, color=col,
                img_key=resolved, library_key=library_key,
                title=title, figsize=(5, 5),
            )
            if b64:
                b64s.append(b64)
        if not b64s:
            return None
        if len(b64s) == 1:
            return b64s[0]
        # Stitch multiple panels side-by-side using numpy/PIL via matplotlib
        import io as _io
        import numpy as _np
        imgs = []
        for b in b64s:
            import base64 as _b64
            arr = plt.imread(_io.BytesIO(_b64.b64decode(b)))
            imgs.append(arr)
        max_h = max(a.shape[0] for a in imgs)
        # Pad shorter images to same height
        padded = []
        for a in imgs:
            if a.shape[0] < max_h:
                pad = _np.ones((max_h - a.shape[0], a.shape[1], a.shape[2]),
                               dtype=a.dtype)
                a = _np.vstack([a, pad])
            padded.append(a)
        combined = _np.hstack(padded)
        fig2, ax2 = plt.subplots(figsize=(combined.shape[1] / 100,
                                          combined.shape[0] / 100))
        ax2.imshow(combined)
        ax2.axis("off")
        fig2.tight_layout(pad=0)
        return _fig_to_b64(fig2)
    except Exception as e:
        logger.warning("figure failed (%s): %s", __name__, e)
        return None


def _threshold_bar(outputs):
    categories = [
        ("Low counts",  outputs.get("removed_low_counts",  0)),
        ("High counts", outputs.get("removed_high_counts", 0)),
        ("Low genes",   outputs.get("removed_low_genes",   0)),
        ("High genes",  outputs.get("removed_high_genes",  0)),
        ("High MT%",    outputs.get("removed_high_mt",     0)),
    ]
    labels, values = zip(*categories)
    fig, ax = plt.subplots(figsize=(6, 3))
    colors = ["#e05252" if v > 0 else "#95a5a6" for v in values]
    bars = ax.bar(labels, values, color=colors, edgecolor="white")
    ax.set_ylabel("Spots removed")
    ax.set_title("Spots removed per filter criterion", fontsize=10, fontweight="bold")
    for bar, val in zip(bars, values):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5, str(val),
                    ha="center", va="bottom", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

def _section_summary(adata, qc_info, dataset_id, timestamp):
    outputs = qc_info.get("outputs", {})
    params  = qc_info.get("params",  {})
    stats   = qc_info.get("summary_stats", {})
    n_before  = outputs.get("n_spots_before", 0)
    n_after   = outputs.get("n_spots_after",  0)
    n_removed = outputs.get("n_spots_removed", 0)
    pct_kept  = f"{100 * n_after / n_before:.1f}%" if n_before > 0 else "?"

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Spots input",   f"{n_before:,}"),
            ("Spots kept",    f"{n_after:,}"),
            ("Spots removed", f"{n_removed:,}"),
            ("Pass rate",     pct_kept),
            ("Genes",         f"{adata.n_vars:,}"),
        ]
    )
    threshold_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in params.items()
    )
    stat_rows = ""
    for col, label in [
        ("total_counts", "Total UMI counts"),
        ("n_genes_by_counts", "Genes detected"),
        ("pct_counts_mt", "MT gene %"),
    ]:
        s = stats.get(col, {})
        if s:
            stat_rows += (
                f"<tr><td>{label}</td>"
                f"<td>{s['mean']:.1f}</td><td>{s['median']:.1f}</td>"
                f"<td>{s['std']:.1f}</td><td>{s['min']:.1f}</td>"
                f"<td>{s['max']:.1f}</td></tr>"
            )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_id}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>QC Metric Summary (retained spots)</h3>
      <table>
        <thead>
          <tr><th>Metric</th><th>Mean</th><th>Median</th><th>Std</th><th>Min</th><th>Max</th></tr>
        </thead>
        <tbody>{stat_rows}</tbody>
      </table>
      <h3>Filter Thresholds Applied</h3>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{threshold_rows}</tbody>
      </table>
    </section>
    """


def _section_distributions(adata, params):
    b64 = _violin_plots(adata, params)
    if not b64:
        return ""
    return f"""
    <section>
      <h2>QC Metric Distributions</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Violin plots per spot</h3>
          <img src="data:image/png;base64,{b64}" alt="violin plots">
        </div>
      </div>
    </section>
    """


def _section_tissue_overview(adata, img_key: Optional[str]) -> str:
    """
    Tissue array overview — first plot in all Visium tutorials.
    Shows in_tissue capture mask and array_row gradient on the H&E image.
    """
    b64 = _tissue_overview(adata, img_key)
    if not b64:
        return ""
    resolved = _resolve_img_key(adata, img_key)
    caption = (
        "Spots overlaid on the H&amp;E image. "
        "Left: capture-area mask (<code>in_tissue</code>). "
        "Right: array row gradient confirming correct orientation."
        if resolved else
        "Spot array layout. Left: capture-area mask. Right: array row gradient."
    )
    return f"""
    <section>
      <h2>Tissue Array Overview</h2>
      <p>{caption}</p>
      <div class="fig-grid">
        <div class="fig-wrap" style="max-width:900px;">
          <img src="data:image/png;base64,{b64}" alt="tissue overview">
        </div>
      </div>
    </section>
    """


def _section_spatial(adata, img_key: Optional[str] = None):
    """QC metrics plotted on the tissue — total counts and MT% per spot."""
    b64_counts = _spatial_scatter(adata, "total_counts", img_key)
    b64_mt     = _spatial_scatter(adata, "pct_counts_mt", img_key)
    if not b64_counts and not b64_mt:
        return ""
    resolved = _resolve_img_key(adata, img_key)
    title = (
        "Spatial Distribution of QC Metrics (on H&amp;E)"
        if resolved else
        "Spatial Distribution of QC Metrics"
    )
    figs = ""
    if b64_counts:
        figs += (
            '<div class="fig-wrap"><h3>Total UMI counts per spot</h3>'
            f'<img src="data:image/png;base64,{b64_counts}" alt="spatial counts"></div>'
        )
    if b64_mt:
        figs += (
            '<div class="fig-wrap"><h3>MT gene % per spot</h3>'
            f'<img src="data:image/png;base64,{b64_mt}" alt="spatial MT%"></div>'
        )
    return f"""
    <section>
      <h2>{title}</h2>
      <p>High MT% spots concentrated at tissue edges may indicate damaged cells.
         Low-count spots in the tissue interior may reflect poor RNA capture.</p>
      <div class="fig-grid">{figs}</div>
    </section>
    """


def _section_filter_breakdown(outputs):
    b64 = _threshold_bar(outputs)
    return f"""
    <section>
      <h2>Filter Breakdown</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Spots removed per filter criterion</h3>
          <img src="data:image/png;base64,{b64}" alt="threshold bar">
        </div>
      </div>
      <p class="note">Spots may fail multiple filters simultaneously &mdash; counts may overlap.</p>
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
        "    <h1>OmicSage &#8212; Spatial QC Report</h1>\n"
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

def generate_spatial_qc_report(
    adata: ad.AnnData,
    output_path: str,
    dataset_id: str = "spatial",
    img_key: Optional[str] = "hires",
) -> str:
    """Generate a self-contained HTML QC report for Visium data.

    Parameters
    ----------
    adata
        AnnData returned by :func:`spatial_qc`.
    output_path
        Path to write the ``.html`` file.
    dataset_id
        Dataset label shown in the report header.
    img_key
        H&E image key stored in ``adata.uns["spatial"][sample]["images"]``.
        Common values: ``"hires"`` (default), ``"lowres"``.
        If the image is absent all figures degrade gracefully to spots-only.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    if "omicsage_spatial_qc" not in adata.uns:
        raise ValueError(
            "adata.uns['omicsage_spatial_qc'] not found. "
            "Run spatial_qc() before generating the report."
        )
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    qc_info   = adata.uns["omicsage_spatial_qc"]
    params    = qc_info.get("params",  {})
    outputs   = qc_info.get("outputs", {})

    sections = [
        _section_summary(adata, qc_info, dataset_id, timestamp),
        _section_tissue_overview(adata, img_key),       # NEW: tissue array + H&E
        _section_distributions(adata, params),
        _section_spatial(adata, img_key),               # UPDATED: img_key support
        _section_filter_breakdown(outputs),
    ]
    html = _render_page(
        title=f"OmicSage -- Spatial QC -- {dataset_id}",
        sections=sections,
        timestamp=timestamp,
    )
    output_path = str(output_path)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    Path(output_path).write_text(html, encoding="utf-8")
    size_kb = Path(output_path).stat().st_size / 1024
    logger.info("Spatial QC report -> %s (%.1f KB)", output_path, size_kb)
    return str(Path(output_path).resolve())
