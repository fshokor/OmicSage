"""
OmicSage — Multiome ATAC Annotation Report
reports/templates/multiome/atac_annotate_report.py

Generated after annotate_atac step.
Output: multiome_03_annotate_report.html

Sections
--------
* Run Summary          — cells / peaks / genes in activity matrix / matched barcodes
* Gene Activity Matrix — top genes by mean activity, shown as a horizontal bar chart
* Label Transfer Table — per-cluster cell counts, majority RNA label, match rate
* Cell Type UMAP       — UMAP coloured by atac_celltype
* Leiden UMAP          — UMAP coloured by atac_leiden (for cross-reference)
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
    header { background: linear-gradient(135deg, #0d2137 0%, #1a3a5c 60%, #0f5c4a 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0d4a3a;
                  border-bottom: 2px solid #e0f0eb; padding-bottom: 10px;
                  margin-bottom: 18px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
             border-left: 3px solid #f0c040; padding: 8px 12px;
             border-radius: 4px; margin-bottom: 14px; }
    .info { font-size: 0.82rem; color: #1a4a6e; background: #e8f4fb;
             border-left: 3px solid #3a90c0; padding: 8px 12px;
             border-radius: 4px; margin-bottom: 14px; }
    code { font-family: "SFMono-Regular", Consolas, monospace;
            background: #eef6f3; padding: 1px 5px; border-radius: 3px;
            font-size: 0.85em; }
    .stat-grid { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }
    .stat-card { background: #eef6f3; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }
    .stat-value { font-size: 1.4rem; font-weight: 700; color: #0d4a3a; }
    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }
    table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }
    th { background: #eef6f3; color: #0d4a3a; font-weight: 600;
          padding: 9px 12px; text-align: left; border-bottom: 2px solid #c0ddd5; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f5fbf8; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 520px; }
    .fig-wrap.wide { flex: 1 1 100%; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #0d2137; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #d0e8e0; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa;
              padding: 24px 0 32px; }
    footer a { color: #0d4a3a; text-decoration: none; }
"""


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _placeholder(msg: str, w: float = 5, h: float = 3) -> str:
    fig, ax = plt.subplots(figsize=(w, h))
    ax.text(0.5, 0.5, msg, ha="center", va="center",
            transform=ax.transAxes, fontsize=10, color="#888")
    ax.axis("off")
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------

def _plot_top_genes_activity(adata: AnnData, top_n: int = 30) -> str:
    """
    Horizontal bar chart of the top_n genes by mean gene activity score.
    """
    if "gene_activity" not in adata.obsm:
        return _placeholder("gene_activity not found in obsm")

    gene_names: list[str] = adata.uns.get("gene_activity_var_names", [])
    mat: np.ndarray = adata.obsm["gene_activity"]

    if mat.shape[1] == 0 or not gene_names:
        return _placeholder("Gene activity matrix is empty")

    n_show = min(top_n, mat.shape[1])
    means  = mat.mean(axis=0)  # shape: (n_genes,)

    top_idx   = np.argsort(means)[-n_show:][::-1]
    top_genes = [gene_names[i] for i in top_idx]
    top_means = means[top_idx]

    fig, ax = plt.subplots(figsize=(8, max(4, n_show * 0.35)))
    cmap   = plt.get_cmap("YlGn", n_show + 2)
    colors = [cmap(n_show - i) for i in range(n_show)]

    y_pos = np.arange(n_show)
    ax.barh(y_pos, top_means[::-1], color=colors[::-1], height=0.7,
            edgecolor="white")
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_genes[::-1], fontsize=8)
    ax.set_xlabel("Mean gene activity score", fontsize=10)
    ax.set_title(
        f"Top {n_show} genes by mean gene activity",
        fontsize=12, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_celltype(adata: AnnData, obs_key: str, title: str) -> str:
    """
    UMAP coloured by any categorical obs column.
    Expects obsm["X_umap_atac"].
    """
    umap_key = "X_umap_atac"
    if umap_key not in adata.obsm or obs_key not in adata.obs.columns:
        return _placeholder(f"UMAP or {obs_key} not found")

    umap   = adata.obsm[umap_key]
    labels = adata.obs[obs_key].astype(str)
    unique = sorted(labels.unique())
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, label in enumerate(unique):
        mask = (labels == label).values
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            s=2, alpha=0.7, color=cmap(i),
            label=label, rasterized=True,
        )

    n_cols = max(1, len(unique) // 15)
    ax.legend(
        markerscale=4, frameon=False, fontsize=7,
        loc="upper right", ncol=n_cols,
        title=obs_key, title_fontsize=7,
    )
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_gene_activity_heatmap(adata: AnnData, top_n: int = 20) -> str:
    """
    Heatmap: top_n genes (rows) × Leiden clusters (cols), mean activity.
    """
    if "gene_activity" not in adata.obsm or "atac_leiden" not in adata.obs.columns:
        return _placeholder("gene_activity or atac_leiden not found")

    gene_names: list[str] = adata.uns.get("gene_activity_var_names", [])
    mat: np.ndarray = adata.obsm["gene_activity"]

    if mat.shape[1] == 0:
        return _placeholder("Gene activity matrix is empty")

    leiden = adata.obs["atac_leiden"].astype(str).values
    clusters = sorted(set(leiden))
    n_genes  = min(top_n, mat.shape[1])

    # Select top_n genes by variance across cells
    gene_var = mat.var(axis=0)
    top_idx  = np.argsort(gene_var)[-n_genes:][::-1]
    top_genes = [gene_names[i] for i in top_idx]
    mat_sub   = mat[:, top_idx]  # n_cells × n_genes

    # Compute mean per cluster
    heatmap = np.zeros((len(clusters), n_genes), dtype=np.float32)
    for ci, cid in enumerate(clusters):
        mask = leiden == cid
        heatmap[ci] = mat_sub[mask].mean(axis=0)

    # Z-score per gene (column) for visibility
    gene_std = heatmap.std(axis=0)
    gene_std[gene_std == 0] = 1
    heatmap_z = (heatmap - heatmap.mean(axis=0)) / gene_std

    fig, ax = plt.subplots(figsize=(max(8, n_genes * 0.45), max(4, len(clusters) * 0.5)))
    im = ax.imshow(heatmap_z, aspect="auto", cmap="RdYlGn", interpolation="nearest")
    ax.set_xticks(np.arange(n_genes))
    ax.set_xticklabels(top_genes, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(clusters)))
    ax.set_yticklabels([f"Cluster {c}" for c in clusters], fontsize=8)
    ax.set_title(
        f"Gene activity heatmap — top {n_genes} variable genes\n"
        "(z-scored across clusters)",
        fontsize=11, fontweight="bold",
    )
    plt.colorbar(im, ax=ax, shrink=0.6, label="z-score")
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML section builders
# ---------------------------------------------------------------------------

def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — ATAC Annotation — {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
<header>
  <h1>OmicSage &mdash; ATAC Annotation</h1>
  <p>{dataset_name} &nbsp;|&nbsp; {timestamp} &nbsp;|&nbsp; Gene activity scores + label transfer</p>
</header>
<main>
{body}
</main>
<footer>
  <p>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
     &nbsp;&bull;&nbsp; {timestamp}</p>
</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    n_cells       = metrics.get("n_cells", "—")
    n_peaks       = metrics.get("n_peaks", "—")
    n_genes       = metrics.get("n_genes_activity", "—")
    n_matched     = metrics.get("n_rna_barcodes_matched", "—")
    rna_provided  = metrics.get("rna_provided", False)
    prom_bp       = metrics.get("promoter_upstream_bp", "—")
    ann_source    = metrics.get("peak_annotation_source", "—")

    match_note = ""
    if not rna_provided:
        match_note = (
            "<p class='note'>No RNA AnnData was provided — "
            "<code>atac_celltype</code> is 'Unknown' for all cells.  "
            "Re-run with <code>rna=</code> to enable label transfer.</p>"
        )
    elif ann_source == "coordinate_fallback":
        match_note = (
            "<p class='note'>Peak annotation was not found in "
            "<code>uns['atac']['peak_annotation']</code>.  "
            "Gene names in the activity matrix are synthetic region labels, "
            "not real gene names.  For production use, load the original "
            "10x h5 file so muon populates <code>peak_annotation</code>.</p>"
        )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Generated: {timestamp} &nbsp;|&nbsp; Dataset: {dataset_name}</p>
      {match_note}
      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-value">{n_cells:,}</div>
          <div class="stat-label">Cells</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_peaks:,}</div>
          <div class="stat-label">Peaks</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_genes:,}</div>
          <div class="stat-label">Genes in activity matrix</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_matched:,}</div>
          <div class="stat-label">RNA barcodes matched</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{prom_bp:,} bp</div>
          <div class="stat-label">Promoter window</div>
        </div>
      </div>
      <p class="info">
        Peak annotation source: <code>{ann_source}</code>
        &nbsp;&bull;&nbsp;
        Gene activity: sum of peak counts overlapping gene body
        + {prom_bp:,} bp upstream promoter window.
      </p>
    </section>"""


def _section_label_transfer_table(adata: AnnData, metrics: dict) -> str:
    if "atac_leiden" not in adata.obs.columns or "atac_celltype" not in adata.obs.columns:
        return ""

    leiden   = adata.obs["atac_leiden"].astype(str)
    celltype = adata.obs["atac_celltype"].astype(str)
    total    = int(adata.n_obs)

    header = (
        "<tr>"
        "<th>Cluster</th>"
        "<th>Cells</th>"
        "<th>% of total</th>"
        "<th>Majority RNA label</th>"
        "</tr>"
    )

    rows = ""
    for cid in sorted(leiden.unique(), key=lambda x: int(x) if x.isdigit() else x):
        mask = (leiden == cid).values
        cnt  = int(mask.sum())
        pct  = 100 * cnt / total if total else 0
        # majority label is the same for all cells in a cluster (by construction)
        label = celltype[mask].iloc[0] if cnt > 0 else "Unknown"
        rows += (
            f"<tr><td>{cid}</td><td>{cnt:,}</td>"
            f"<td>{pct:.1f}%</td><td>{label}</td></tr>"
        )

    n_unknown = int((celltype == "Unknown").sum())
    unknown_note = ""
    if n_unknown > 0:
        pct_unk = 100 * n_unknown / total if total else 0
        unknown_note = (
            f"<p class='note'>"
            f"{n_unknown:,} cells ({pct_unk:.1f}%) received 'Unknown' — "
            f"their ATAC Leiden cluster had no matching RNA barcodes.  "
            f"This is expected for cells that dropped out of the RNA library.</p>"
        )

    return f"""
    <section>
      <h2>Label Transfer — Cluster Table</h2>
      <p>
        Each ATAC Leiden cluster receives the majority
        <code>cell_type_vote</code> label from RNA cells with matching barcodes.
        Clusters with no RNA barcode overlap receive 'Unknown'.
      </p>
      {unknown_note}
      <table><thead>{header}</thead><tbody>{rows}</tbody></table>
    </section>"""


def _section_gene_activity(fig_bar: str, fig_heatmap: str) -> str:
    return f"""
    <section>
      <h2>Gene Activity Matrix</h2>
      <p>
        Gene activity scores are computed by summing peak counts overlapping
        each gene body plus the upstream promoter window.  Higher score = more
        open chromatin near that gene.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Top genes by mean activity</h3>
          <img src="data:image/png;base64,{fig_bar}" alt="Top genes bar chart">
        </div>
        <div class="fig-wrap wide">
          <h3>Heatmap — top variable genes × Leiden clusters (z-scored)</h3>
          <img src="data:image/png;base64,{fig_heatmap}" alt="Gene activity heatmap">
        </div>
      </div>
    </section>"""


def _section_umaps(fig_celltype: str, fig_leiden: str) -> str:
    return f"""
    <section>
      <h2>ATAC UMAP — Cell Types and Clusters</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>UMAP — transferred cell type (atac_celltype)</h3>
          <img src="data:image/png;base64,{fig_celltype}" alt="UMAP cell type">
        </div>
        <div class="fig-wrap">
          <h3>UMAP — Leiden clusters (atac_leiden)</h3>
          <img src="data:image/png;base64,{fig_leiden}" alt="UMAP Leiden">
        </div>
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_atac_annotate_report(
    adata: AnnData,
    metrics: dict,
    report_path: str = "reports/multiome_03_annotate_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the ATAC annotation HTML report.

    Parameters
    ----------
    adata : AnnData
        ATAC AnnData returned by annotate_atac().
    metrics : dict
        Metrics dict returned by annotate_atac().
    report_path : str
        Output HTML path.
    dataset_name : str
        Dataset name shown in the report header.

    Returns
    -------
    str : Absolute path of the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building ATAC annotate report for '{dataset_name}' ...", flush=True)

    print("  • top genes bar chart ...", flush=True)
    fig_bar = _plot_top_genes_activity(adata)

    print("  • gene activity heatmap ...", flush=True)
    fig_heatmap = _plot_gene_activity_heatmap(adata)

    print("  • UMAP cell type ...", flush=True)
    fig_celltype = _plot_umap_celltype(
        adata, "atac_celltype", "ATAC UMAP — transferred cell type"
    )

    print("  • UMAP Leiden ...", flush=True)
    fig_leiden = _plot_umap_celltype(
        adata, "atac_leiden", "ATAC UMAP — Leiden clusters"
    )

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_label_transfer_table(adata, metrics),
            _section_gene_activity(fig_bar, fig_heatmap),
            _section_umaps(fig_celltype, fig_leiden),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )

    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
