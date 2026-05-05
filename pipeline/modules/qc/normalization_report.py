"""
OmicSage — Normalization Report
reports/normalization_report.py

Generates a self-contained HTML report for the normalization step.
Produces the same figures that would go in the notebook, but saves them
to disk and wraps them in a portable HTML file — no Quarto needed yet.

Usage
-----
    # From CLI:
    conda activate omicsage
    python reports/normalization_report.py \
        --input  data/processed/GSE194122_cite_rna_qc.h5ad \
        --output data/processed/GSE194122_cite_normalized.h5ad \
        --report reports/output/normalization_report.html \
        --dataset GSE194122_CITE \
        --batch-key batch

    # From notebook (replace Cell 6 onwards):
    from reports.normalization_report import run_normalization_report
    run_normalization_report(
        adata_norm=adata_norm,
        metrics=metrics,
        report_path="reports/output/normalization_report.html",
        dataset_name="GSE194122_CITE",
    )

Figures produced
----------------
1. HVG scatter        — mean expression vs variance, HVGs highlighted
2. Library size dist  — per-cell total counts before vs after normalization
3. Top 20 HVGs        — dot plot of most variable gene names
4. Gene detection     — cumulative fraction of genes detected per cell
"""

from __future__ import annotations

import argparse
import base64
import io
import logging
import os
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


def _plot_hvg_scatter(adata_norm: AnnData) -> str:
    """
    Highly variable genes scatter plot.
    X = mean expression, Y = normalised variance (or dispersion).
    HVGs shown in orange, background genes in grey.
    """
    fig, ax = plt.subplots(figsize=(7, 4.5))

    var = adata_norm.var.copy()

    # seurat_v3 gives 'variances_norm'; seurat/cell_ranger give 'dispersions_norm'
    if "variances_norm" in var.columns:
        y_col, y_label = "variances_norm", "Normalised variance"
        x_col = "means" if "means" in var.columns else var.columns[0]
        x_label = "Mean expression"
    elif "dispersions_norm" in var.columns:
        y_col, y_label = "dispersions_norm", "Normalised dispersion"
        x_col = "means" if "means" in var.columns else var.columns[0]
        x_label = "Mean expression"
    else:
        # Fallback — just plot HVG flag as a bar count
        n_hvg = int(var["highly_variable"].sum())
        n_non = len(var) - n_hvg
        fig2, ax2 = plt.subplots(figsize=(4, 4))
        ax2.bar(["HVG", "Non-HVG"], [n_hvg, n_non], color=["#E07B3A", "#AAAAAA"])
        ax2.set_ylabel("Gene count")
        ax2.set_title("Highly Variable Genes")
        return _fig_to_b64(fig2)

    hvg_mask = var["highly_variable"]
    ax.scatter(
        var.loc[~hvg_mask, x_col], var.loc[~hvg_mask, y_col],
        s=4, alpha=0.4, color="#AAAAAA", label="Non-HVG", rasterized=True,
    )
    ax.scatter(
        var.loc[hvg_mask, x_col], var.loc[hvg_mask, y_col],
        s=6, alpha=0.7, color="#E07B3A", label=f"HVG (n={hvg_mask.sum()})",
        rasterized=True,
    )
    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel(y_label, fontsize=11)
    ax.set_title("Highly Variable Genes", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_library_size(adata_norm: AnnData) -> str:
    """
    Violin + strip of per-cell total counts before and after normalization.
    'Before' = sum of layers['counts']; 'After' = sum of layers['logcounts'].
    """
    import scipy.sparse as sp

    raw = adata_norm.layers["counts"]
    norm = adata_norm.layers["logcounts"]

    raw_arr  = raw.toarray()  if sp.issparse(raw)  else np.asarray(raw)
    norm_arr = norm.toarray() if sp.issparse(norm) else np.asarray(norm)

    raw_sums  = raw_arr.sum(axis=1)
    norm_sums = norm_arr.sum(axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))

    for ax, data, title, color in zip(
        axes,
        [raw_sums, norm_sums],
        ["Raw counts per cell", "log1p-normalised sum per cell"],
        ["#4C78A8", "#54A868"],
    ):
        ax.violinplot(data, positions=[0], showmedians=True, widths=0.6)
        # colour the violin body
        for pc in ax.collections:
            pc.set_facecolor(color)
            pc.set_alpha(0.6)
        ax.set_xticks([])
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_ylabel("Sum per cell", fontsize=10)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f"{x:,.0f}"
        ))
        ax.spines[["top", "right"]].set_visible(False)
        # Annotate median
        median = np.median(data)
        ax.axhline(median, color="black", linewidth=1, linestyle="--", alpha=0.5)
        ax.text(0.55, median, f"median={median:,.1f}", va="center",
                fontsize=8, transform=ax.get_yaxis_transform())

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_top_hvgs(adata_norm: AnnData, n: int = 20) -> str:
    """
    Horizontal bar chart of the top-n HVGs ranked by normalised variance
    (or dispersion). Shows gene names clearly.
    """
    var = adata_norm.var[adata_norm.var["highly_variable"]].copy()

    rank_col = None
    for col in ["variances_norm", "dispersions_norm", "highly_variable_rank"]:
        if col in var.columns:
            rank_col = col
            break

    if rank_col is None or rank_col == "highly_variable_rank":
        # Fall back to alphabetical if no variance metric
        top = var.head(n)
        values = np.ones(len(top))
        xlabel = "HVG (no variance metric available)"
    else:
        top = var.nlargest(n, rank_col)
        values = top[rank_col].values
        xlabel = rank_col.replace("_", " ").title()

    fig, ax = plt.subplots(figsize=(7, max(3, n * 0.32)))
    colors = ["#E07B3A"] * len(top)
    ax.barh(range(len(top)), values[::-1], color=colors, alpha=0.8)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.index.tolist()[::-1], fontsize=9)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_title(f"Top {n} Highly Variable Genes", fontsize=13, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_gene_detection(adata_norm: AnnData) -> str:
    """
    Cumulative distribution: fraction of cells in which each gene is detected
    (count > 0 in layers['counts']).  Shows how many genes are expressed in
    what fraction of cells — useful for assessing data quality after filtering.
    """
    import scipy.sparse as sp

    raw = adata_norm.layers["counts"]
    raw_arr = raw.toarray() if sp.issparse(raw) else np.asarray(raw)

    detection_rate = (raw_arr > 0).mean(axis=0)   # fraction of cells per gene
    sorted_rates = np.sort(detection_rate)[::-1]
    x = np.arange(1, len(sorted_rates) + 1)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x, sorted_rates, color="#4C78A8", linewidth=1.5)
    ax.axhline(0.01, color="grey", linestyle="--", linewidth=0.8, label="1% cells")
    ax.axhline(0.10, color="grey", linestyle=":",  linewidth=0.8, label="10% cells")
    ax.set_xlabel("Genes ranked by detection rate", fontsize=11)
    ax.set_ylabel("Fraction of cells detected", fontsize=11)
    ax.set_title("Gene Detection Rate", fontsize=13, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylim(0, 1.02)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage — Normalization Report — {dataset_name}</title>
<style>
  body       {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                max-width: 1100px; margin: 40px auto; padding: 0 24px;
                color: #222; background: #fafafa; }}
  h1         {{ color: #2a5c8a; border-bottom: 2px solid #2a5c8a; padding-bottom: 8px; }}
  h2         {{ color: #444; margin-top: 2em; }}
  table      {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td     {{ border: 1px solid #ddd; padding: 8px 14px; text-align: left; }}
  th         {{ background: #2a5c8a; color: white; }}
  tr:nth-child(even) {{ background: #f2f2f2; }}
  .fig-grid  {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px;
                margin: 1.5em 0; }}
  .fig-box   {{ background: white; border: 1px solid #ddd; border-radius: 6px;
                padding: 12px; box-shadow: 0 1px 4px rgba(0,0,0,.07); }}
  .fig-box.wide {{ grid-column: 1 / -1; }}
  img        {{ width: 100%; height: auto; display: block; }}
  .meta      {{ font-size: 0.85em; color: #666; margin-top: 0.3em; }}
  .badge     {{ display:inline-block; padding:2px 8px; border-radius:4px;
                font-size:0.8em; font-weight:bold; }}
  .ok        {{ background:#d4edda; color:#155724; }}
  footer     {{ margin-top: 3em; font-size: 0.8em; color: #999;
                border-top: 1px solid #eee; padding-top: 12px; }}
</style>
</head>
<body>

<h1>🧬 OmicSage — Normalization Report</h1>
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
    <strong>Highly Variable Genes</strong>
    <p class="meta">Orange = selected HVGs. Grey = non-HVG background genes.</p>
    <img src="data:image/png;base64,{fig_hvg}" alt="HVG scatter">
  </div>
  <div class="fig-box">
    <strong>Library Size Distribution</strong>
    <p class="meta">Per-cell total counts before and after normalization.</p>
    <img src="data:image/png;base64,{fig_libsize}" alt="Library size">
  </div>
  <div class="fig-box">
    <strong>Gene Detection Rate</strong>
    <p class="meta">Fraction of cells in which each gene is detected (raw counts > 0).</p>
    <img src="data:image/png;base64,{fig_detection}" alt="Gene detection">
  </div>
  <div class="fig-box wide">
    <strong>Top 20 Highly Variable Genes</strong>
    <p class="meta">Ranked by normalised variance / dispersion.</p>
    <img src="data:image/png;base64,{fig_top_hvg}" alt="Top HVGs">
  </div>
</div>

<h2>Provenance</h2>
<table>
  <tr><th>Key</th><th>Value</th></tr>
  {provenance_rows}
</table>

<footer>
  Generated by OmicSage · pipeline/modules/qc/normalize.py ·
  <a href="https://github.com/fshokor/OmicSage">github.com/fshokor/OmicSage</a>
</footer>
</body>
</html>
"""


def _build_html(
    adata_norm: AnnData,
    metrics: dict,
    dataset_name: str,
) -> str:
    """Render all figures and assemble the HTML string."""

    # --- figures ---
    print("  Rendering HVG scatter ...", flush=True)
    fig_hvg       = _plot_hvg_scatter(adata_norm)
    print("  Rendering library size ...", flush=True)
    fig_libsize   = _plot_library_size(adata_norm)
    print("  Rendering top HVGs ...", flush=True)
    fig_top_hvg   = _plot_top_hvgs(adata_norm)
    print("  Rendering gene detection ...", flush=True)
    fig_detection = _plot_gene_detection(adata_norm)

    # --- summary table rows ---
    summary_items = [
        ("Cells",                  f"{metrics.get('n_cells', adata_norm.n_obs):,}"),
        ("Genes",                  f"{metrics.get('n_genes', adata_norm.n_vars):,}"),
        ("HVGs selected",          f"{metrics.get('n_hvg_selected', '—'):,}"),
        ("HVG flavor",             metrics.get("hvg_flavor", "—")),
        ("Batch key",              metrics.get("batch_key") or "None"),
        ("Target sum (CP10K)",     f"{metrics.get('target_sum', 1e4):,.0f}"),
        ("log1p applied",          "✓" if metrics.get("log1p_applied") else "✗"),
        ("Raw counts layer",       metrics.get("raw_counts_in_layer", "counts")),
        ("Normalized layer",       metrics.get("normalized_in_layer", "logcounts")),
        ("Mean norm sum per cell", f"{metrics.get('mean_counts_per_cell_after_norm', 0):.3f}"),
    ]
    summary_rows = "\n  ".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in summary_items
    )

    # --- provenance table rows (from uns) ---
    prov = adata_norm.uns.get("omicsage_normalization", {})
    provenance_rows = "\n  ".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov.items()
    )

    return _HTML_TEMPLATE.format(
        dataset_name=dataset_name,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        summary_rows=summary_rows,
        provenance_rows=provenance_rows,
        fig_hvg=fig_hvg,
        fig_libsize=fig_libsize,
        fig_top_hvg=fig_top_hvg,
        fig_detection=fig_detection,
    )


# ---------------------------------------------------------------------------
# Public API — callable from notebook
# ---------------------------------------------------------------------------

def run_normalization_report(
    adata_norm: AnnData,
    metrics: dict,
    report_path: str = "reports/output/normalization_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the normalization HTML report and write it to disk.

    Parameters
    ----------
    adata_norm : AnnData
        Normalized AnnData returned by ``normalize()``.
        Must have layers['counts'] and layers['logcounts'].
    metrics : dict
        Metrics dict returned by ``normalize()``.
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Label shown in the report header.

    Returns
    -------
    str
        Absolute path to the written report file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Building normalization report for '{dataset_name}' ...", flush=True)
    html = _build_html(adata_norm, metrics, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved → {out.resolve()}", flush=True)
    return str(out.resolve())


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate OmicSage normalization report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="QC-filtered h5ad (raw counts in .X)")
    p.add_argument("--output",  default=None,
                   help="Where to save the normalized h5ad "
                        "(if omitted, normalization runs but h5ad is not saved)")
    p.add_argument("--report",  default="reports/output/normalization_report.html",
                   help="Path for the output HTML report")
    p.add_argument("--dataset", default=None,
                   help="Dataset label shown in the report (default: input filename)")
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--n-top-genes", type=int, default=2000)
    p.add_argument("--hvg-flavor", default="seurat_v3",
                   choices=["seurat_v3", "seurat", "cell_ranger"])
    p.add_argument("--batch-key", default=None,
                   help="obs column to use for per-batch HVG selection")
    return p.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args()

    # Add OmicSage root to path if running as a script
    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from pipeline.modules.qc.normalize import normalize

    dataset_name = args.dataset or Path(args.input).stem

    print(f"Loading {args.input} ...", flush=True)
    adata = sc.read_h5ad(args.input)
    print(adata, flush=True)

    print("Running normalization ...", flush=True)
    adata_norm, metrics = normalize(
        adata,
        target_sum=args.target_sum,
        n_top_genes=args.n_top_genes,
        hvg_flavor=args.hvg_flavor,
        batch_key=args.batch_key,
    )

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        adata_norm.write_h5ad(args.output)
        print(f"Normalized h5ad saved → {args.output}", flush=True)

    run_normalization_report(
        adata_norm=adata_norm,
        metrics=metrics,
        report_path=args.report,
        dataset_name=dataset_name,
    )


if __name__ == "__main__":
    main()
