"""
OmicSage — Leiden Clustering Report
reports/cluster_report.py

Generates a self-contained HTML report for the Leiden clustering step.
Produces the same figures that would go in the notebook, but saves them
to disk and wraps them in a portable HTML file — no Quarto needed yet.

Usage
-----
    # From CLI:
    conda activate omicsage
    python reports/cluster_report.py \
        --input  data/processed/GSE194122_cite_reduced.h5ad \
        --output data/processed/GSE194122_cite_clustered.h5ad \
        --report reports/output/cluster_report.html \
        --dataset GSE194122_CITE

    # From notebook:
    from reports.cluster_report import run_cluster_report
    run_cluster_report(
        adata_clustered=adata_clustered,
        metrics=metrics,
        report_path="reports/output/cluster_report.html",
        dataset_name="GSE194122_CITE",
    )

Figures produced
----------------
1. UMAP coloured by each Leiden resolution tested
2. Silhouette score vs resolution bar chart (selected resolution highlighted)
3. Cluster size distribution at the best resolution
4. UMAP coloured by ground-truth cell_type (if present in obs)
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


def _plot_umap_resolutions(adata_clustered: AnnData, metrics: dict) -> str:
    """
    Grid of UMAP plots coloured by each Leiden resolution tested.
    The best (auto-selected) resolution is marked with a gold border.
    Falls back gracefully when X_umap is absent.
    """
    if "X_umap" not in adata_clustered.obsm:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "UMAP not computed", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    resolutions = metrics["resolutions"]
    best_res = metrics["best_resolution"]
    umap = adata_clustered.obsm["X_umap"]

    n = len(resolutions)
    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 4.5 * nrows),
                              squeeze=False)
    axes_flat = axes.flatten()

    for idx, res in enumerate(resolutions):
        ax = axes_flat[idx]
        obs_key = f"leiden_{res:.4g}"

        if obs_key not in adata_clustered.obs.columns:
            ax.text(0.5, 0.5, f"leiden_{res:.4g}\nnot found",
                    ha="center", va="center", transform=ax.transAxes)
            ax.axis("off")
            continue

        labels = adata_clustered.obs[obs_key].astype("category")
        unique_labels = labels.cat.categories.tolist()
        n_clusters = len(unique_labels)
        cmap = plt.get_cmap("tab20", max(n_clusters, 1))
        color_map = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}
        colors = [color_map[lbl] for lbl in labels]

        ax.scatter(umap[:, 0], umap[:, 1], c=colors, s=1.5,
                   alpha=0.6, rasterized=True)
        ax.set_title(
            f"resolution = {res:.2g}  ({n_clusters} clusters)",
            fontsize=10, fontweight="bold",
        )
        ax.set_xlabel("UMAP 1", fontsize=8)
        ax.set_ylabel("UMAP 2", fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)

        # Gold border on the best resolution panel
        if res == best_res:
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_edgecolor("#E0A800")
                spine.set_linewidth(2.5)
            ax.set_title(
                f"★ resolution = {res:.2g}  ({n_clusters} clusters)  [selected]",
                fontsize=10, fontweight="bold", color="#B8860B",
            )

    # Hide any unused axes
    for idx in range(n, len(axes_flat)):
        axes_flat[idx].set_visible(False)

    fig.suptitle("UMAP — Leiden Clustering across Resolutions",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_silhouette_bar(metrics: dict) -> str:
    """
    Bar chart of silhouette score per resolution.
    The best resolution bar is highlighted in orange; others in steel blue.
    """
    resolutions = metrics["resolutions"]
    scores = [metrics["silhouette_scores"][r] for r in resolutions]
    best_res = metrics["best_resolution"]

    colors = ["#E07B3A" if r == best_res else "#4C78A8" for r in resolutions]

    fig, ax = plt.subplots(figsize=(max(6, len(resolutions) * 1.2), 4))
    x = np.arange(len(resolutions))
    bars = ax.bar(x, scores, color=colors, width=0.6, alpha=0.85, edgecolor="white")

    # Value labels on bars
    for bar, score in zip(bars, scores):
        ypos = bar.get_height() + 0.005 if score >= 0 else bar.get_height() - 0.02
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            ypos,
            f"{score:.3f}",
            ha="center", va="bottom", fontsize=9,
        )

    ax.set_xticks(x)
    ax.set_xticklabels([f"{r:.2g}" for r in resolutions], fontsize=10)
    ax.set_xlabel("Leiden Resolution", fontsize=11)
    ax.set_ylabel("Silhouette Score", fontsize=11)
    ax.set_title(
        "Silhouette Score vs Resolution\n"
        f"(selected: {best_res:.2g}, score: {metrics['best_silhouette']:.3f})",
        fontsize=12, fontweight="bold",
    )
    ax.axhline(0, color="#999", linewidth=0.8, linestyle="--")
    ax.spines[["top", "right"]].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#E07B3A", label=f"Selected (res={best_res:.2g})"),
        Patch(facecolor="#4C78A8", label="Other resolutions"),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=9)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_cluster_sizes(adata_clustered: AnnData, metrics: dict) -> str:
    """
    Bar chart of cluster sizes at the best resolution, sorted descending.
    """
    best_res = metrics["best_resolution"]
    obs_key = f"leiden_{best_res:.4g}"

    if obs_key not in adata_clustered.obs.columns:
        obs_key = "leiden"

    if obs_key not in adata_clustered.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "Cluster labels not found",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.axis("off")
        return _fig_to_b64(fig)

    counts = adata_clustered.obs[obs_key].value_counts().sort_values(ascending=False)
    n_clusters = len(counts)
    cmap = plt.get_cmap("tab20", n_clusters)

    fig, ax = plt.subplots(figsize=(max(8, n_clusters * 0.5), 4))
    x = np.arange(n_clusters)
    bars = ax.bar(x, counts.values,
                  color=[cmap(i) for i in range(n_clusters)],
                  width=0.7, alpha=0.85, edgecolor="white")

    ax.set_xticks(x)
    ax.set_xticklabels(counts.index.tolist(), fontsize=8, rotation=45, ha="right")
    ax.set_xlabel("Cluster", fontsize=11)
    ax.set_ylabel("Number of Cells", fontsize=11)
    ax.set_title(
        f"Cluster Size Distribution  (resolution={best_res:.2g}, "
        f"{n_clusters} clusters, {adata_clustered.n_obs:,} cells)",
        fontsize=12, fontweight="bold",
    )
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.spines[["top", "right"]].set_visible(False)

    # Median line
    median_size = int(np.median(counts.values))
    ax.axhline(median_size, color="#C0392B", linewidth=1.2, linestyle="--",
               label=f"Median: {median_size:,} cells")
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_cell_type(adata_clustered: AnnData) -> Optional[str]:
    """
    UMAP coloured by ground-truth cell_type label (if present in obs).
    Returns None when cell_type is absent — caller skips the section.
    """
    if "X_umap" not in adata_clustered.obsm:
        return None
    if "cell_type" not in adata_clustered.obs.columns:
        return None

    umap = adata_clustered.obsm["X_umap"]
    cell_types = adata_clustered.obs["cell_type"].astype(str)
    unique_types = sorted(cell_types.unique())
    cmap = plt.get_cmap("tab20", len(unique_types))
    color_map = {ct: cmap(i) for i, ct in enumerate(unique_types)}
    colors = [color_map[ct] for ct in cell_types]

    n_types = len(unique_types)
    # Scale figure width to give the legend enough room
    fig_w = 8 + max(0, (n_types - 20) // 10) * 1.5
    fig_h = max(6, n_types * 0.22)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.scatter(umap[:, 0], umap[:, 1], c=colors, s=1.5,
               alpha=0.6, rasterized=True)

    # Legend placed outside the axes to the right — never overlaps the UMAP
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color_map[ct], label=ct)
                       for ct in unique_types]
    ncol = max(1, n_types // 25)
    legend = ax.legend(
        handles=legend_elements,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=False,
        fontsize=7.5,
        ncol=ncol,
        handlelength=1.2,
        handleheight=1.2,
        labelspacing=0.4,
    )

    ax.set_title("UMAP — Ground-Truth Cell Type", fontsize=13, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    # bbox_inches="tight" in _fig_to_b64 ensures the legend is not clipped
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
<title>OmicSage — Clustering Report: {dataset_name}</title>
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

<h1>🧬 OmicSage — Leiden Clustering Report</h1>
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
    <strong>UMAP — Leiden Clustering across Resolutions</strong>
    <p class="meta">
      Each panel shows cluster assignments at one resolution.
      The gold-bordered panel is the auto-selected best resolution
      (highest silhouette score).
    </p>
    <img src="data:image/png;base64,{fig_umap_resolutions}" alt="UMAP resolutions">
  </div>
  <div class="fig-box wide">
    <strong>Silhouette Score vs Resolution</strong>
    <p class="meta">
      Higher silhouette score = better-separated clusters.
      Orange bar = auto-selected resolution.
    </p>
    <img src="data:image/png;base64,{fig_silhouette}" alt="Silhouette scores">
  </div>
  <div class="fig-box wide">
    <strong>Cluster Size Distribution — Best Resolution ({best_res:.2g})</strong>
    <p class="meta">
      Bars sorted by cluster size (descending).
      Red dashed line = median cluster size.
      Large size disparities may indicate over- or under-clustering.
    </p>
    <img src="data:image/png;base64,{fig_cluster_sizes}" alt="Cluster sizes">
  </div>
  {cell_type_section}
</div>

<h2>Provenance</h2>
<table>
  <tr><th>Key</th><th>Value</th></tr>
  {provenance_rows}
</table>

<footer>
  Generated by OmicSage · reports/cluster_report.py ·
  <a href="https://github.com/fshokor/OmicSage">github.com/fshokor/OmicSage</a>
</footer>
</body>
</html>
"""

_CELL_TYPE_SECTION_TEMPLATE = """\
  <div class="fig-box wide">
    <strong>UMAP — Ground-Truth Cell Type</strong>
    <p class="meta">Each colour is one annotated cell type from the original publication.
    Use this as a reference to assess whether Leiden clusters align with known biology.</p>
    <img src="data:image/png;base64,{fig_cell_type}" alt="UMAP cell type">
  </div>"""


def _build_html(
    adata_clustered: AnnData,
    metrics: dict,
    dataset_name: str,
) -> str:
    """Render all figures and assemble the HTML string."""

    print("  Rendering UMAP resolution panels ...", flush=True)
    fig_umap_resolutions = _plot_umap_resolutions(adata_clustered, metrics)

    print("  Rendering silhouette bar chart ...", flush=True)
    fig_silhouette = _plot_silhouette_bar(metrics)

    print("  Rendering cluster size distribution ...", flush=True)
    fig_cluster_sizes = _plot_cluster_sizes(adata_clustered, metrics)

    # Optional cell_type panel
    cell_type_section = ""
    fig_cell_type = _plot_umap_cell_type(adata_clustered)
    if fig_cell_type is not None:
        print("  Rendering ground-truth cell type UMAP ...", flush=True)
        cell_type_section = _CELL_TYPE_SECTION_TEMPLATE.format(
            fig_cell_type=fig_cell_type,
        )

    # --- summary table ---
    best_res = metrics["best_resolution"]
    silhouette_scores = metrics["silhouette_scores"]
    sil_rows = "  ".join(
        f"<tr><td>Silhouette @ res={r:.2g}</td><td>{silhouette_scores[r]:.4f}</td></tr>"
        for r in metrics["resolutions"]
    )
    summary_items = [
        ("Cells",                       f"{adata_clustered.n_obs:,}"),
        ("Resolutions tested",          ", ".join(f"{r:.2g}" for r in metrics["resolutions"])),
        ("Best resolution (selected)",  f"{best_res:.2g}"),
        ("Clusters at best resolution", str(metrics["best_n_clusters"])),
        ("Best silhouette score",       f"{metrics['best_silhouette']:.4f}"),
    ]
    summary_rows = "\n  ".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in summary_items
    )
    # Append per-resolution silhouette rows
    summary_rows += "\n  " + sil_rows

    # --- provenance table ---
    prov = adata_clustered.uns.get("omicsage_cluster", {})
    # Flatten nested dicts (n_clusters, silhouette_scores) for display
    prov_display = {}
    for k, v in prov.items():
        if isinstance(v, dict):
            for subk, subv in v.items():
                display_subv = f"{subv:.4f}" if isinstance(subv, float) else subv
                prov_display[f"{k}[{subk}]"] = display_subv
        else:
            prov_display[k] = v
    provenance_rows = "\n  ".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov_display.items()
    )

    return _HTML_TEMPLATE.format(
        dataset_name=dataset_name,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        best_res=best_res,
        summary_rows=summary_rows,
        provenance_rows=provenance_rows,
        fig_umap_resolutions=fig_umap_resolutions,
        fig_silhouette=fig_silhouette,
        fig_cluster_sizes=fig_cluster_sizes,
        cell_type_section=cell_type_section,
    )


# ---------------------------------------------------------------------------
# Public API — callable from notebook
# ---------------------------------------------------------------------------

def run_cluster_report(
    adata_clustered: AnnData,
    metrics: dict,
    report_path: str = "reports/output/cluster_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the Leiden clustering HTML report and write it to disk.

    Parameters
    ----------
    adata_clustered : AnnData
        Clustered AnnData returned by ``cluster()``.
        Must have obsm['X_umap'] and obs['leiden_*'] columns.
    metrics : dict
        Metrics dict returned by ``cluster()``.
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

    print(f"Building clustering report for '{dataset_name}' ...", flush=True)
    html = _build_html(adata_clustered, metrics, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved → {out.resolve()}", flush=True)
    return str(out.resolve())


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate OmicSage Leiden clustering report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="Reduced h5ad (obsm['X_pca'], obsp['connectivities'])")
    p.add_argument("--output",  default=None,
                   help="Where to save the clustered h5ad "
                        "(if omitted, cluster runs but h5ad is not saved)")
    p.add_argument("--report",  default="reports/output/cluster_report.html",
                   help="Path for the output HTML report")
    p.add_argument("--dataset", default=None,
                   help="Dataset label shown in the report (default: input filename)")
    p.add_argument("--resolutions", nargs="+", type=float,
                   default=[0.2, 0.4, 0.6, 0.8, 1.0],
                   help="Leiden resolutions to sweep")
    return p.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from pipeline.modules.clustering.cluster import cluster

    dataset_name = args.dataset or Path(args.input).stem

    print(f"Loading {args.input} ...", flush=True)
    adata = sc.read_h5ad(args.input)
    print(adata, flush=True)

    print("Running Leiden clustering ...", flush=True)
    adata_clustered, metrics = cluster(
        adata,
        resolution_range=args.resolutions,
    )

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        adata_clustered.write_h5ad(args.output)
        print(f"Clustered h5ad saved → {args.output}", flush=True)

    run_cluster_report(
        adata_clustered=adata_clustered,
        metrics=metrics,
        report_path=args.report,
        dataset_name=dataset_name,
    )


if __name__ == "__main__":
    main()
