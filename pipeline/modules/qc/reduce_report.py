"""
OmicSage — Dimensionality Reduction Report
reports/reduce_report.py

Generates a self-contained HTML report for the dimensionality reduction step.
Produces the same figures that would go in the notebook, but saves them
to disk and wraps them in a portable HTML file — no Quarto needed yet.

Usage
-----
    # From CLI:
    conda activate omicsage
    python reports/reduce_report.py \
        --input  data/processed/GSE194122_cite_normalized.h5ad \
        --output data/processed/GSE194122_cite_reduced.h5ad \
        --report reports/output/reduce_report.html \
        --dataset GSE194122_CITE

    # From notebook:
    from reports.reduce_report import run_reduce_report
    run_reduce_report(
        adata_reduced=adata_reduced,
        metrics=metrics,
        report_path="reports/output/reduce_report.html",
        dataset_name="GSE194122_CITE",
    )

Figures produced
----------------
1. Scree plot             — variance explained per PC + cumulative curve
                            with selected n_pcs marked
2. UMAP coloured by n_genes_by_counts  — QC metric overlay
3. UMAP coloured by total_counts       — QC metric overlay
4. UMAP coloured by pct_counts_mt      — QC metric overlay
5. PCA coloured by n_genes_by_counts   — same QC metrics on PCA axes
"""

from __future__ import annotations

import argparse
import base64
import io
import logging
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")   # non-interactive backend — safe for scripts
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Figure helpers — each returns a base64-encoded PNG string
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    """Encode a matplotlib figure as a base64 PNG for embedding in HTML."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_scree(adata_reduced: AnnData, metrics: dict) -> str:
    """
    Scree plot showing variance explained per PC and cumulative variance curve.
    A vertical dashed line marks the auto-selected (or manual) n_pcs_used.
    """
    variance_ratio = np.array(metrics["variance_explained_per_pc"])
    cumvar = np.cumsum(variance_ratio)
    n_pcs_used = metrics["n_pcs_used"]
    pc_selection_method = metrics.get("pc_selection_method", "")
    n_comps = len(variance_ratio)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # --- Left: per-PC variance bar chart ---
    ax = axes[0]
    x = np.arange(1, n_comps + 1)
    colors = ["#E07B3A" if i < n_pcs_used else "#AAAAAA" for i in range(n_comps)]
    ax.bar(x, variance_ratio * 100, color=colors, width=0.8, alpha=0.85)
    ax.axvline(n_pcs_used + 0.5, color="#C0392B", linewidth=1.5,
               linestyle="--", label=f"Selected: {n_pcs_used} PCs ({pc_selection_method})")
    ax.set_xlabel("Principal Component", fontsize=11)
    ax.set_ylabel("Variance explained (%)", fontsize=11)
    ax.set_title("Scree Plot", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))

    # --- Right: cumulative variance curve ---
    ax = axes[1]
    ax.plot(x, cumvar * 100, color="#4C78A8", linewidth=2)
    ax.axvline(n_pcs_used, color="#C0392B", linewidth=1.5, linestyle="--",
               label=f"Selected: {n_pcs_used} PCs")
    ax.axhline(
        float(cumvar[n_pcs_used - 1]) * 100,
        color="#C0392B", linewidth=0.8, linestyle=":",
    )
    ax.set_xlabel("Number of PCs", fontsize=11)
    ax.set_ylabel("Cumulative variance explained (%)", fontsize=11)
    ax.set_title("Cumulative Variance", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=10))

    # Annotate cumulative variance at selected n_pcs
    cum_at_selected = float(cumvar[n_pcs_used - 1]) * 100
    ax.annotate(
        f"{cum_at_selected:.1f}%",
        xy=(n_pcs_used, cum_at_selected),
        xytext=(n_pcs_used + max(1, n_comps // 10), cum_at_selected - 5),
        fontsize=8,
        color="#C0392B",
        arrowprops=dict(arrowstyle="->", color="#C0392B", lw=0.8),
    )

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_qc(adata_reduced: AnnData) -> str:
    """
    2×2 grid of UMAP plots coloured by QC metrics:
    n_genes_by_counts, total_counts, pct_counts_mt, doublet_score.
    Falls back gracefully when a QC column is absent (e.g. doublet_score).
    """
    if "X_umap" not in adata_reduced.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "UMAP not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    umap = adata_reduced.obsm["X_umap"]
    obs = adata_reduced.obs

    qc_cols = [
        ("n_genes_by_counts", "Genes per cell",   "YlOrRd"),
        ("total_counts",      "Total UMI counts", "Blues"),
        ("pct_counts_mt",     "MT% per cell",     "Reds"),
        ("doublet_score",     "Doublet score",    "Purples"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    axes = axes.flatten()

    for ax, (col, title, cmap) in zip(axes, qc_cols):
        if col in obs.columns:
            values = obs[col].values.astype(float)
            sc_plot = ax.scatter(
                umap[:, 0], umap[:, 1],
                c=values, cmap=cmap, s=2, alpha=0.6, rasterized=True,
            )
            plt.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02)
        else:
            ax.scatter(umap[:, 0], umap[:, 1], s=2, alpha=0.4,
                       color="#AAAAAA", rasterized=True)
            ax.text(0.5, 0.02, f"'{col}' not in obs",
                    ha="center", transform=ax.transAxes,
                    fontsize=8, color="grey")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=9)
        ax.set_ylabel("UMAP 2", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

    fig.suptitle("UMAP — QC Metric Overlays", fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_pca_qc(adata_reduced: AnnData) -> str:
    """
    PCA scatter (PC1 vs PC2) coloured by n_genes_by_counts and total_counts.
    Helps identify if major axes of variation are driven by technical factors.
    """
    if "X_pca" not in adata_reduced.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "PCA not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    pca = adata_reduced.obsm["X_pca"]
    obs = adata_reduced.obs

    qc_cols = [
        ("n_genes_by_counts", "Genes per cell",   "YlOrRd"),
        ("total_counts",      "Total UMI counts", "Blues"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, (col, title, cmap) in zip(axes, qc_cols):
        if col in obs.columns:
            values = obs[col].values.astype(float)
            sc_plot = ax.scatter(
                pca[:, 0], pca[:, 1],
                c=values, cmap=cmap, s=3, alpha=0.6, rasterized=True,
            )
            plt.colorbar(sc_plot, ax=ax, fraction=0.03, pad=0.02)
        else:
            ax.scatter(pca[:, 0], pca[:, 1], s=3, alpha=0.4,
                       color="#AAAAAA", rasterized=True)
            ax.text(0.5, 0.02, f"'{col}' not in obs",
                    ha="center", transform=ax.transAxes,
                    fontsize=8, color="grey")
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("PC 1", fontsize=10)
        ax.set_ylabel("PC 2", fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("PCA — QC Metric Overlays (PC1 vs PC2)",
                 fontsize=13, fontweight="bold")
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_batch(adata_reduced: AnnData, batch_key: Optional[str]) -> Optional[str]:
    """
    UMAP coloured by batch label.  Only rendered when batch_key is provided
    and the column exists in obs.  Returns None otherwise.
    """
    if batch_key is None:
        return None
    if "X_umap" not in adata_reduced.obsm:
        return None
    if batch_key not in adata_reduced.obs.columns:
        logger.warning("batch_key='%s' not found in obs — skipping batch UMAP", batch_key)
        return None

    umap = adata_reduced.obsm["X_umap"]
    batches = adata_reduced.obs[batch_key].astype(str)
    unique_batches = sorted(batches.unique())
    cmap = plt.get_cmap("tab20", len(unique_batches))
    color_map = {b: cmap(i) for i, b in enumerate(unique_batches)}

    fig, ax = plt.subplots(figsize=(7, 5))
    for batch in unique_batches:
        mask = batches == batch
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            s=2, alpha=0.6, color=color_map[batch],
            label=batch, rasterized=True,
        )
    ax.legend(
        markerscale=4, frameon=False, fontsize=8,
        loc="upper right", ncol=max(1, len(unique_batches) // 10),
    )
    ax.set_title(f"UMAP — coloured by {batch_key}", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML template
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage — Dimensionality Reduction Report: {dataset_name}</title>
<style>
  body       {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                sans-serif; margin: 2em auto; max-width: 1100px;
                color: #2c3e50; background: #fafafa; }}
  h1         {{ color: #1a252f; border-bottom: 2px solid #E07B3A;
                padding-bottom: 0.3em; }}
  h2         {{ color: #2c3e50; margin-top: 2em; }}
  p.meta     {{ color: #666; font-size: 0.9em; margin-top: -0.5em; }}
  table      {{ border-collapse: collapse; width: 100%; margin-bottom: 1.5em; }}
  th, td     {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
  th         {{ background: #f2f2f2; font-weight: 600; }}
  tr:nth-child(even) {{ background: #fafafa; }}
  img        {{ max-width: 100%; height: auto; border: 1px solid #eee;
                border-radius: 4px; }}
  .fig-grid  {{ display: grid;
                grid-template-columns: repeat(auto-fit, minmax(480px, 1fr));
                gap: 1.5em; margin-top: 1em; }}
  .fig-box   {{ background: white; border: 1px solid #e0e0e0; border-radius: 6px;
                padding: 1em; }}
  .fig-box.wide {{ grid-column: 1 / -1; }}
  .badge     {{ display: inline-block; padding: 2px 8px; border-radius: 3px;
                font-size:0.8em; font-weight:bold; }}
  .ok        {{ background:#d4edda; color:#155724; }}
  footer     {{ margin-top: 3em; font-size: 0.8em; color: #999;
                border-top: 1px solid #eee; padding-top: 12px; }}
</style>
</head>
<body>

<h1>🧬 OmicSage — Dimensionality Reduction Report</h1>
<p class="meta">
  Dataset: <strong>{dataset_name}</strong> &nbsp;|&nbsp;
  Generated: <strong>{timestamp}</strong> &nbsp;|&nbsp;
  OmicSage v0.1.0
</p>

<h2>Summary</h2>
<table>
  <tr><th>Parameter</th><th>Value</th></tr>
  {summary_rows}
</table>

<h2>Figures</h2>
<div class="fig-grid">
  <div class="fig-box wide">
    <strong>Scree Plot — Variance Explained per PC</strong>
    <p class="meta">
      Orange bars = PCs selected for neighbor graph.
      Grey bars = PCs computed but not used.
      Red dashed line = selected cut-off.
    </p>
    <img src="data:image/png;base64,{fig_scree}" alt="Scree plot">
  </div>
  <div class="fig-box wide">
    <strong>UMAP — QC Metric Overlays</strong>
    <p class="meta">
      Colour intensity reflects per-cell QC values.
      Cells with extreme values may indicate technical artifacts.
    </p>
    <img src="data:image/png;base64,{fig_umap_qc}" alt="UMAP QC">
  </div>
  <div class="fig-box wide">
    <strong>PCA — QC Metric Overlays (PC1 vs PC2)</strong>
    <p class="meta">
      If PC1/PC2 correlate strongly with sequencing depth or MT%, consider
      regressing out these covariates before clustering.
    </p>
    <img src="data:image/png;base64,{fig_pca_qc}" alt="PCA QC">
  </div>
  {batch_section}
</div>

<h2>Provenance</h2>
<table>
  <tr><th>Key</th><th>Value</th></tr>
  {provenance_rows}
</table>

<footer>
  Generated by OmicSage · reports/reduce_report.py ·
  <a href="https://github.com/fshokor/OmicSage">github.com/fshokor/OmicSage</a>
</footer>
</body>
</html>
"""

_BATCH_SECTION_TEMPLATE = """\
  <div class="fig-box wide">
    <strong>UMAP — Coloured by {batch_key}</strong>
    <p class="meta">Each colour represents one batch. Well-mixed batches suggest
    successful integration or that batch effects are minor.</p>
    <img src="data:image/png;base64,{fig_batch}" alt="UMAP batch">
  </div>"""


def _build_html(
    adata_reduced: AnnData,
    metrics: dict,
    dataset_name: str,
    batch_key: Optional[str] = None,
) -> str:
    """Render all figures and assemble the HTML string."""

    # --- figures ---
    print("  Rendering scree plot ...", flush=True)
    fig_scree    = _plot_scree(adata_reduced, metrics)
    print("  Rendering UMAP QC overlays ...", flush=True)
    fig_umap_qc  = _plot_umap_qc(adata_reduced)
    print("  Rendering PCA QC overlays ...", flush=True)
    fig_pca_qc   = _plot_pca_qc(adata_reduced)

    # Optional batch UMAP
    batch_section = ""
    fig_batch = _plot_umap_batch(adata_reduced, batch_key)
    if fig_batch is not None:
        print("  Rendering batch UMAP ...", flush=True)
        batch_section = _BATCH_SECTION_TEMPLATE.format(
            batch_key=batch_key,
            fig_batch=fig_batch,
        )

    # --- summary table rows ---
    cum_var_pct = metrics.get("cumulative_variance_explained", 0) * 100
    summary_items = [
        ("Cells",                           f"{metrics.get('n_cells', adata_reduced.n_obs):,}"),
        ("Genes",                           f"{metrics.get('n_genes', adata_reduced.n_vars):,}"),
        ("HVGs used for PCA",               f"{metrics.get('n_hvg_used', '—'):,}"),
        ("PCA components computed",         f"{metrics.get('n_comps_computed', '—')}"),
        ("PCs used for neighbor graph",     f"{metrics.get('n_pcs_used', '—')}"),
        ("PC selection method",             metrics.get("pc_selection_method", "—")),
        ("Cumulative variance explained",   f"{cum_var_pct:.1f}%"),
        ("Neighbors (k)",                   f"{metrics.get('n_neighbors', '—')}"),
        ("Embeddings computed",             ", ".join(metrics.get("embeddings_computed", []))),
        ("t-SNE computed",                  "✓" if metrics.get("run_tsne") else "✗"),
    ]
    summary_rows = "\n  ".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in summary_items
    )

    # --- provenance table rows (from uns) ---
    prov = adata_reduced.uns.get("omicsage_reduce", {})
    # Truncate variance_explained_per_pc — too long to display in full
    prov_display = {
        k: (f"[{v[0]:.4f}, {v[1]:.4f}, ... ] ({len(v)} values)"
            if k == "variance_explained_per_pc" and isinstance(v, list)
            else v)
        for k, v in prov.items()
    }
    provenance_rows = "\n  ".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov_display.items()
    )

    return _HTML_TEMPLATE.format(
        dataset_name=dataset_name,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        summary_rows=summary_rows,
        provenance_rows=provenance_rows,
        fig_scree=fig_scree,
        fig_umap_qc=fig_umap_qc,
        fig_pca_qc=fig_pca_qc,
        batch_section=batch_section,
    )


# ---------------------------------------------------------------------------
# Public API — callable from notebook
# ---------------------------------------------------------------------------

def run_reduce_report(
    adata_reduced: AnnData,
    metrics: dict,
    report_path: str = "reports/output/reduce_report.html",
    dataset_name: str = "dataset",
    batch_key: Optional[str] = None,
) -> str:
    """
    Generate the dimensionality reduction HTML report and write it to disk.

    Parameters
    ----------
    adata_reduced : AnnData
        Reduced AnnData returned by ``reduce()``.
        Must have obsm['X_pca'] and obsm['X_umap'].
    metrics : dict
        Metrics dict returned by ``reduce()``.
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Label shown in the report header.
    batch_key : str, optional
        obs column for the batch UMAP panel. Omitted if None.

    Returns
    -------
    str
        Absolute path to the written report file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Building dimensionality reduction report for '{dataset_name}' ...", flush=True)
    html = _build_html(adata_reduced, metrics, dataset_name, batch_key=batch_key)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved → {out.resolve()}", flush=True)
    return str(out.resolve())


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate OmicSage dimensionality reduction report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="Normalized h5ad (log1p values in .X, HVGs flagged)")
    p.add_argument("--output",  default=None,
                   help="Where to save the reduced h5ad "
                        "(if omitted, reduce runs but h5ad is not saved)")
    p.add_argument("--report",  default="reports/output/reduce_report.html",
                   help="Path for the output HTML report")
    p.add_argument("--dataset", default=None,
                   help="Dataset label shown in the report (default: input filename)")
    p.add_argument("--n-comps", type=int, default=50,
                   help="Number of PCA components to compute")
    p.add_argument("--n-pcs", type=int, default=None,
                   help="Number of PCs for neighbor graph (default: auto-select)")
    p.add_argument("--n-pcs-method", default="elbow",
                   choices=["elbow", "variance", "fixed"],
                   help="PC auto-selection strategy")
    p.add_argument("--n-neighbors", type=int, default=15,
                   help="Number of neighbors for kNN graph")
    p.add_argument("--run-tsne", action="store_true",
                   help="Also compute t-SNE (slow on large datasets)")
    p.add_argument("--batch-key", default=None,
                   help="obs column for batch UMAP panel")
    return p.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args()

    # Add OmicSage root to path if running as a script
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from pipeline.modules.qc.reduce import reduce

    dataset_name = args.dataset or Path(args.input).stem

    print(f"Loading {args.input} ...", flush=True)
    adata = sc.read_h5ad(args.input)
    print(adata, flush=True)

    print("Running dimensionality reduction ...", flush=True)
    adata_reduced, metrics = reduce(
        adata,
        n_comps=args.n_comps,
        n_pcs=args.n_pcs,
        n_pcs_method=args.n_pcs_method,
        n_neighbors=args.n_neighbors,
        run_tsne=args.run_tsne,
    )

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        adata_reduced.write_h5ad(args.output)
        print(f"Reduced h5ad saved → {args.output}", flush=True)

    run_reduce_report(
        adata_reduced=adata_reduced,
        metrics=metrics,
        report_path=args.report,
        dataset_name=dataset_name,
        batch_key=args.batch_key,
    )


if __name__ == "__main__":
    main()
