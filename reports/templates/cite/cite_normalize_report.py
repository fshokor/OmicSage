"""
OmicSage — CITE-seq ADT Normalization Report
reports/templates/cite/cite_normalize_report.py

Generated after normalize_adt step.
Output: cite_01_normalize_report.html

Usage
-----
    from reports.templates.cite.cite_normalize_report import run_cite_normalize_report
    run_cite_normalize_report(
        adt=adt_norm,
        metrics=metrics,
        report_path="reports/GSE194122/cite_01_normalize_report.html",
        dataset_name="BMMC CITE-seq (NeurIPS 2021)",
    )
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


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _get_dense(layer) -> np.ndarray:
    if hasattr(layer, "toarray"):
        return layer.toarray()
    return np.asarray(layer, dtype=float)


# ---------------------------------------------------------------------------
# Plots — CLR (always shown)
# ---------------------------------------------------------------------------

def _plot_clr_violin(adt: AnnData) -> str:
    """Violin plot of CLR values across proteins (up to 40)."""
    clr = _get_dense(adt.layers.get("adt_clr", adt.X))
    n_proteins = clr.shape[1]
    max_show = 40
    idx = np.linspace(0, n_proteins - 1, min(max_show, n_proteins), dtype=int)
    clr_sub = clr[:, idx]
    labels = [adt.var_names[i] for i in idx]

    fig, ax = plt.subplots(figsize=(max(10, len(idx) * 0.35), 4.5))
    parts = ax.violinplot(
        [clr_sub[:, j] for j in range(clr_sub.shape[1])],
        positions=range(clr_sub.shape[1]),
        showmedians=True, widths=0.7,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#4C78A8"); pc.set_alpha(0.6)
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color("#16213e"); parts[key].set_linewidth(1.0)
    ax.set_xticks(range(clr_sub.shape[1]))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("CLR value", fontsize=10)
    ax.set_title(
        f"CLR-Normalised ADT Distribution"
        + (f" (showing {len(idx)} of {n_proteins})" if n_proteins > max_show else ""),
        fontsize=12, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_clr_mean_bar(adt: AnnData) -> str:
    """Per-protein CLR mean bar chart (ranked)."""
    clr = _get_dense(adt.layers.get("adt_clr", adt.X))
    means = clr.mean(axis=0)
    order = np.argsort(means)[::-1]
    labels = [adt.var_names[i] for i in order]
    values = means[order]

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.28), 4))
    ax.bar(range(len(labels)), values, color="#e07b3a", alpha=0.8, width=0.7)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Mean CLR value", fontsize=10)
    ax.set_title("Mean CLR Expression per Protein (ranked)", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_raw_vs_clr(adt: AnnData) -> str:
    """Scatter: raw total counts vs CLR sum per cell."""
    raw = adt.layers.get("counts", None)
    clr = adt.layers.get("adt_clr", adt.X)
    if raw is None:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "Raw counts layer not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)
    raw_sum = _get_dense(raw).sum(axis=1)
    clr_sum = _get_dense(clr).sum(axis=1)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(raw_sum, clr_sum, s=2, alpha=0.4, color="#4C78A8", rasterized=True)
    ax.set_xlabel("Raw total ADT counts per cell", fontsize=10)
    ax.set_ylabel("CLR sum per cell", fontsize=10)
    ax.set_title("Raw Counts vs CLR Sum", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)



def _plot_total_counts_clr_before_after(adt: AnnData) -> str:
    """
    Two-panel histplot of total protein counts per cell before and after CLR.
    Mirrors tutorial cell 23 (raw) + cell 36 (CLR):
      Before: layers["counts"].sum(axis=1)  xlim (0, 20000)
      After:  layers["adt_clr"].sum(axis=1) xlim (0, 400)
    """
    raw = _get_dense(adt.layers["counts"]).sum(axis=1)
    clr = _get_dense(adt.layers["adt_clr"]).sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # --- Before ---
    ax = axes[0]
    ax.hist(raw, bins=50, color="#4C78A8", alpha=0.8, edgecolor="none")
    ax.set_xlim(0, 20000)
    ax.set_xlabel("Total protein counts per cell", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Total counts per cell\n(raw)", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # --- After CLR ---
    ax = axes[1]
    ax.hist(clr, bins=50, color="#e07b3a", alpha=0.8, edgecolor="none")
    ax.set_xlim(0, 400)
    ax.set_xlabel("Total CLR counts per cell", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Total counts per cell\n(CLR normalised)", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Before vs After CLR Normalisation\n"
        "Sequencing depth differences compressed; range reduced to ~0–400",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Plots — DSB (only shown when DSB was applied)
# ---------------------------------------------------------------------------

def _plot_dsb_violin(adt: AnnData) -> str:
    """
    Violin plot of DSB values per protein — the key quality check.
    The tutorial signature: bimodal shape = background left peak + signal right peak.
    """
    dsb = _get_dense(adt.layers["adt_dsb"])
    n_proteins = dsb.shape[1]
    max_show = 40
    idx = np.linspace(0, n_proteins - 1, min(max_show, n_proteins), dtype=int)
    dsb_sub = dsb[:, idx]
    labels = [adt.var_names[i] for i in idx]

    fig, ax = plt.subplots(figsize=(max(10, len(idx) * 0.35), 4.5))
    parts = ax.violinplot(
        [dsb_sub[:, j] for j in range(dsb_sub.shape[1])],
        positions=range(dsb_sub.shape[1]),
        showmedians=True, widths=0.7,
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#3cb371"); pc.set_alpha(0.65)
    for key in ("cmedians", "cbars", "cmins", "cmaxes"):
        if key in parts:
            parts[key].set_color("#1a3a1a"); parts[key].set_linewidth(1.0)
    # Draw DSB=0 reference line (background threshold)
    ax.axhline(0, color="#cc3333", linewidth=1.2, linestyle="--",
               label="DSB = 0 (background threshold)")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xticks(range(dsb_sub.shape[1]))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("DSB value", fontsize=10)
    ax.set_title(
        f"DSB-Normalised ADT Distribution"
        + (f" (showing {len(idx)} of {n_proteins})" if n_proteins > max_show else "")
        + "\n(red line = background threshold; values > 0 = truly expressed)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_dsb_vs_clr(adt: AnnData) -> str:
    """
    Scatter comparing CLR vs DSB values per protein per cell (sampled).
    Shows how DSB shifts the background towards 0 and compresses noise.
    """
    dsb = _get_dense(adt.layers["adt_dsb"])
    clr = _get_dense(adt.layers["adt_clr"])

    n_cells, n_proteins = dsb.shape
    # Sample up to 3000 cells × all proteins for the scatter
    rng = np.random.default_rng(0)
    cell_idx = rng.choice(n_cells, size=min(3000, n_cells), replace=False)
    clr_vals = clr[cell_idx, :].ravel()
    dsb_vals = dsb[cell_idx, :].ravel()

    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(clr_vals, dsb_vals, s=1, alpha=0.15, color="#5a6abf", rasterized=True)
    ax.axhline(0, color="#cc3333", linewidth=1.0, linestyle="--", label="DSB = 0")
    ax.set_xlabel("CLR value", fontsize=10)
    ax.set_ylabel("DSB value", fontsize=10)
    ax.set_title("CLR vs DSB Values\n(sampled cells × proteins)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_dsb_density_grid(adt: AnnData, n_proteins: int = 12) -> str:
    """
    Density histogram per protein showing bimodal distribution after DSB.
    Shows the top-expressed proteins (highest DSB mean).
    The tutorial Figure 1: left peak = background, right peak = expressing cells.
    """
    dsb = _get_dense(adt.layers["adt_dsb"])
    means = dsb.mean(axis=0)
    # Pick the n_proteins with highest mean DSB (most informative)
    top_idx = np.argsort(means)[::-1][:n_proteins]
    labels = [adt.var_names[i] for i in top_idx]

    ncols = 4
    nrows = int(np.ceil(n_proteins / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 2.4))
    axes = axes.ravel()

    for plot_i, (prot_i, label) in enumerate(zip(top_idx, labels)):
        ax = axes[plot_i]
        vals = dsb[:, prot_i]
        ax.hist(vals, bins=60, color="#3cb371", alpha=0.75, density=True,
                edgecolor="none")
        ax.axvline(0, color="#cc3333", linewidth=1.0, linestyle="--")
        ax.set_title(label, fontsize=8, fontweight="bold")
        ax.set_xlabel("DSB", fontsize=7)
        ax.tick_params(labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)

    # Hide unused axes
    for ax in axes[n_proteins:]:
        ax.set_visible(False)

    fig.suptitle(
        f"DSB Density per Protein (top {n_proteins} by mean DSB)\n"
        "Red line = 0 (background threshold); bimodal shape confirms DSB worked",
        fontsize=10, fontweight="bold", y=1.01,
    )
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_dsb_ambient_bar(adt: AnnData, metrics: dict) -> str:
    """
    Bar chart of mean DSB value per protein sorted by value.
    Negative mean = protein is at/below background (ambient or not expressed).
    Positive mean = protein expressed in majority of cells.
    """
    dsb = _get_dense(adt.layers["adt_dsb"])
    means = dsb.mean(axis=0)
    order = np.argsort(means)
    labels = [adt.var_names[i] for i in order]
    values = means[order]
    colors = ["#cc3333" if v < 0 else "#3cb371" for v in values]

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.28), 4))
    ax.bar(range(len(labels)), values, color=colors, alpha=0.85, width=0.8)
    ax.axhline(0, color="#333", linewidth=0.8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=7)
    ax.set_ylabel("Mean DSB value", fontsize=10)
    ax.set_title(
        "Mean DSB Value per Protein\n(red = below background / ambient; green = expressed)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML rendering
# ---------------------------------------------------------------------------

def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq ADT Normalization Report -- {dataset_name}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }}
    header {{ background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}
    section {{ background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }}
    section h2 {{ font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }}
    section h3 {{ font-size: 0.98rem; font-weight: 600; color: #16213e; margin: 18px 0 8px; }}
    section p  {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .dsb-badge {{ display: inline-block; padding: 3px 10px; border-radius: 20px;
                  font-size: 0.78rem; font-weight: 700; letter-spacing: 0.5px; }}
    .dsb-yes {{ background: #d4edda; color: #155724; }}
    .dsb-no  {{ background: #fff3cd; color: #856404; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    code {{ font-family: "SFMono-Regular", Consolas, monospace;
            background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em; }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{ background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }}
    .stat-card.dsb {{ background: #d4edda; }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-card.dsb .stat-value {{ color: #155724; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }}
    th {{ background: #f0f2ff; color: #0f3460; font-weight: 600;
          padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0; }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}
    .fig-grid {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }}
    .fig-wrap {{ flex: 1 1 300px; max-width: 520px; }}
    .fig-wrap.wide {{ flex: 1 1 100%; max-width: 100%; }}
    .fig-wrap h3 {{ font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }}
    .fig-wrap img {{ width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }}
    .callout {{ background: #e8f4fd; border-left: 4px solid #3498db;
                padding: 12px 16px; border-radius: 0 6px 6px 0;
                margin: 12px 0; font-size: 0.88rem; color: #1a3a5c; }}
    .callout.warn {{ background: #fff8e1; border-left-color: #f39c12; color: #5d3a00; }}
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage — CITE-seq ADT Normalization Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    dsb_applied = metrics.get("dsb_applied", False)
    dsb_badge = (
        '<span class="dsb-badge dsb-yes">✓ DSB applied</span>'
        if dsb_applied else
        '<span class="dsb-badge dsb-no">CLR only (no empty droplets)</span>'
    )
    active_layer = metrics.get("active_layer", "adt_clr")

    # Core stat cards (always shown)
    cards = [
        ("Cells",       f"{metrics.get('n_cells', '?'):,}",  False),
        ("Proteins",    str(metrics.get("n_proteins", "?")), False),
        ("CLR axis",    str(metrics.get("clr_axis", "?")),   False),
        ("Active layer",active_layer,                         dsb_applied),
    ]
    # DSB-specific cards (only when DSB was run)
    if dsb_applied:
        cards += [
            ("Empty droplets", f"{metrics.get('dsb_n_empty_droplets', '?'):,}", True),
            ("DSB > 0 (%)",
             f"{metrics.get('dsb_fraction_positive', 0)*100:.1f}%",             True),
            ("DSB mean",
             f"{metrics.get('dsb_mean', 0):.3f}",                               True),
        ]

    stat_html = "".join(
        f'<div class="stat-card{"  dsb" if is_dsb else ""}">'
        f'<div class="stat-value">{val}</div>'
        f'<div class="stat-label">{lbl}</div></div>'
        for lbl, val, is_dsb in cards
    )

    # Parameter table
    rows = [
        ("Normalization method",
         "DSB + CLR (DSB is active layer)" if dsb_applied else "CLR only"),
        ("CLR axis",
         "per-protein across cells (axis=0)" if metrics.get("clr_axis") == 0
         else "per-cell across proteins (axis=1)"),
        ("adata.X after normalization", f"<code>{active_layer}</code>"),
        ("Raw counts layer", f"<code>{metrics.get('raw_counts_in_layer', 'counts')}</code>"),
        ("CLR layer",        f"<code>{metrics.get('clr_in_layer', 'adt_clr')}</code>"),
        ("CLR max",          f"{metrics.get('clr_max', 0):.3f}"),
        ("CLR min",          f"{metrics.get('clr_min', 0):.3f}"),
        ("CLR mean/cell",    f"{metrics.get('clr_mean_per_cell', 0):.3f}"),
        ("Raw median total counts/cell",
         f"{metrics.get('raw_median_total_counts_per_cell', 0):.1f}"),
    ]
    if dsb_applied:
        iso = metrics.get("dsb_isotype_controls_used", [])
        rows += [
            ("DSB layer",           "<code>adt_dsb</code>"),
            ("DSB pseudocount",     str(metrics.get("dsb_pseudocount", 10))),
            ("DSB denoising",       str(metrics.get("dsb_denoise", True))),
            ("Isotype controls",    ", ".join(iso) if iso else "None"),
            ("DSB min",             f"{metrics.get('dsb_min', 0):.3f}"),
            ("DSB max",             f"{metrics.get('dsb_max', 0):.3f}"),
        ]

    rows_html = "".join(f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in rows)

    clr_only_callout = "" if dsb_applied else """
    <div class="callout warn">
      <strong>DSB not applied.</strong> No empty droplets were provided
      (<code>dsb_empty_adata=None</code>). Downstream steps will use CLR-normalised
      values. DSB is strongly recommended for CITE-seq data when empty droplets are
      available — it removes ambient protein background and produces interpretable
      bimodal distributions. See the
      <a href="https://www.sc-best-practices.org/surface_protein/normalization.html">
      sc-best-practices normalization chapter</a> for details.
    </div>"""

    dsb_callout = """
    <div class="callout">
      <strong>DSB normalization applied.</strong> DSB uses empty droplets to estimate
      and subtract ambient protein signal, then applies per-cell denoising via Gaussian
      mixture modelling. After DSB, values near 0 represent background-level expression
      and positive values represent genuine protein expression. This removes the long
      right-tailed noise typical of CLR-only normalization. Reference:
      Mulè et al. 2022 (<em>Nature Communications</em> 13:2099).
    </div>""" if dsb_applied else ""

    return f"""
    <section>
      <h2>Run Summary {dsb_badge}</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      {clr_only_callout}
      {dsb_callout}
      <div class="stat-grid">{stat_html}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows_html}</tbody></table>
    </section>"""


def _section_clr_figures(
    fig_violin: str, fig_bar: str, fig_scatter: str,
    fig_clr_before_after: str, dsb_applied: bool
) -> str:
    heading = "CLR Figures" if dsb_applied else "Figures"
    note = (
        "<p>CLR is stored in <code>adata.layers['adt_clr']</code> as a fallback. "
        "The active layer for downstream analysis is <code>adt_dsb</code>.</p>"
        if dsb_applied else ""
    )
    return f"""
    <section>
      <h2>{heading}</h2>
      {note}
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>Total Counts per Cell — Before vs After CLR</h3>
          <img src="data:image/png;base64,{fig_clr_before_after}" alt="Before after CLR histplot">
        </div>
        <div class="fig-wrap wide">
          <h3>CLR Distribution per Protein (violin)</h3>
          <img src="data:image/png;base64,{fig_violin}" alt="CLR violin">
        </div>
        <div class="fig-wrap wide">
          <h3>Mean CLR Expression per Protein (ranked)</h3>
          <img src="data:image/png;base64,{fig_bar}" alt="CLR mean bar">
        </div>
        <div class="fig-wrap">
          <h3>Raw Counts vs CLR Sum per Cell</h3>
          <img src="data:image/png;base64,{fig_scatter}" alt="Raw vs CLR scatter">
        </div>
      </div>
    </section>"""


def _plot_total_counts_before_after(adt: AnnData) -> str:
    """
    Two-panel histplot of total protein counts per cell before and after DSB.
    Mirrors tutorial cells 23 + 25:
      Before: layers["counts"].sum(axis=1)  xlim (0, 20000)
      After:  X.sum(axis=1)                 xlim (-200, 300)
    """
    raw = _get_dense(adt.layers["counts"]).sum(axis=1)
    dsb = _get_dense(adt.layers["adt_dsb"]).sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # --- Before ---
    ax = axes[0]
    ax.hist(raw, bins=50, color="#4C78A8", alpha=0.8, edgecolor="none")
    ax.set_xlim(0, 20000)
    ax.set_xlabel("Total protein counts per cell", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Total counts per cell\n(raw)", fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)

    # --- After ---
    ax = axes[1]
    ax.hist(dsb, bins=50, color="#3cb371", alpha=0.8, edgecolor="none")
    ax.axvline(0, color="#cc3333", linewidth=1.2, linestyle="--",
               label="DSB = 0")
    ax.set_xlim(-200, 300)
    ax.set_xlabel("Total DSB counts per cell", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Total counts per cell\n(DSB normalised)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Before vs After DSB Normalisation\n"
        "Large range compressed; differences become biology-driven not noise-driven",
        fontsize=10, y=1.02,
    )
    fig.tight_layout()
    return _fig_to_b64(fig)


def _section_dsb_figures(
    fig_dsb_violin: str,
    fig_dsb_vs_clr: str,
    fig_dsb_density: str,
    fig_dsb_ambient: str,
    fig_before_after: str,
) -> str:
    return f"""
    <section>
      <h2>DSB Figures</h2>
      <div class="callout">
        Key quality check: each protein should show a <strong>bimodal distribution</strong>
        after DSB — a left peak at background (DSB ≈ 0) and a right peak for
        expressing cells. Proteins where the distribution is entirely below 0
        are not expressed in this sample and represent ambient signal.
      </div>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <h3>Total Counts per Cell — Before vs After DSB</h3>
          <img src="data:image/png;base64,{fig_before_after}" alt="Before after DSB histplot">
        </div>
        <div class="fig-wrap wide">
          <h3>DSB Distribution per Protein (violin) — background threshold at red line</h3>
          <img src="data:image/png;base64,{fig_dsb_violin}" alt="DSB violin">
        </div>
        <div class="fig-wrap wide">
          <h3>Per-Protein Density Histograms — bimodal shape confirms DSB worked</h3>
          <img src="data:image/png;base64,{fig_dsb_density}" alt="DSB density grid">
        </div>
        <div class="fig-wrap wide">
          <h3>Mean DSB per Protein — negative mean = ambient / not expressed</h3>
          <img src="data:image/png;base64,{fig_dsb_ambient}" alt="DSB ambient bar">
        </div>
        <div class="fig-wrap">
          <h3>CLR vs DSB per Cell (sampled) — DSB compresses background towards 0</h3>
          <img src="data:image/png;base64,{fig_dsb_vs_clr}" alt="CLR vs DSB scatter">
        </div>
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_cite_normalize_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_01_normalize_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")
    dsb_applied = metrics.get("dsb_applied", False)

    print(f"Building CITE normalize report for '{dataset_name}' ...", flush=True)

    # CLR figures (always)
    print("  Rendering CLR violin ...", flush=True)
    fig_clr_violin  = _plot_clr_violin(adt)
    print("  Rendering CLR mean bar ...", flush=True)
    fig_clr_bar     = _plot_clr_mean_bar(adt)
    print("  Rendering raw vs CLR scatter ...", flush=True)
    fig_raw_clr          = _plot_raw_vs_clr(adt)
    print("  Rendering CLR before/after histplot ...", flush=True)
    fig_clr_before_after = _plot_total_counts_clr_before_after(adt)

    sections = [
        _section_summary(metrics, dataset_name, timestamp),
        _section_clr_figures(fig_clr_violin, fig_clr_bar, fig_raw_clr,
                              fig_clr_before_after, dsb_applied),
    ]

    # DSB figures (only when DSB was run)
    if dsb_applied and "adt_dsb" in adt.layers:
        print("  Rendering DSB violin ...", flush=True)
        fig_dsb_violin  = _plot_dsb_violin(adt)
        print("  Rendering DSB vs CLR scatter ...", flush=True)
        fig_dsb_vs_clr  = _plot_dsb_vs_clr(adt)
        print("  Rendering DSB density grid ...", flush=True)
        fig_dsb_density = _plot_dsb_density_grid(adt)
        print("  Rendering DSB ambient bar ...", flush=True)
        fig_dsb_ambient = _plot_dsb_ambient_bar(adt, metrics)
        print("  Rendering before/after histplot ...", flush=True)
        fig_before_after = _plot_total_counts_before_after(adt)
        sections.append(
            _section_dsb_figures(fig_dsb_violin, fig_dsb_vs_clr,
                                  fig_dsb_density, fig_dsb_ambient,
                                  fig_before_after)
        )

    html = _render_page(sections=sections, timestamp=timestamp, dataset_name=dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
