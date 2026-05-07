"""
OmicSage — Cell-Type Annotation Report
reports/annotate_report.py

Generates a self-contained HTML report for the cell-type annotation step.
Mirrors the structure of cluster_report.py — same CSS, same figure helpers,
same CLI entry point pattern.

Usage
-----
    # From CLI:
    conda activate omicsage
    python reports/annotate_report.py \\
        --input  data/processed/GSE194122_cite_clustered.h5ad \\
        --output data/processed/GSE194122_cite_annotated.h5ad \\
        --report reports/annotate_report.html \\
        --dataset GSE194122_CITE

    # From notebook:
    from reports.annotate_report import run_annotate_report
    run_annotate_report(
        adata_annotated=adata_annotated,
        annotation_dict=annotation_dict,
        report_path="reports/annotate_report.html",
        dataset_name="GSE194122_CITE",
    )

Figures produced
----------------
1. UMAP coloured by CellTypist coarse labels
2. UMAP coloured by CellTypist fine labels
3. UMAP coloured by marker-score labels
4. UMAP coloured by consensus cell_type (majority vote)  [if vote was run]
5. UMAP coloured by ground-truth cell_type               [if present in obs]
6. Marker score heatmap (clusters × cell types)          [if marker_score_df present]
7. Cluster → consensus label assignment table
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
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Figure helpers — each returns a base64-encoded PNG string
# (identical helper to cluster_report.py)
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    """Encode a matplotlib figure as a base64 PNG for embedding in HTML."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_umap_labels(
    adata: AnnData,
    obs_col: str,
    title: str,
    subtitle: str = "",
) -> Optional[str]:
    """
    UMAP scatter coloured by a categorical obs column.
    Returns None gracefully when X_umap or the column is absent.
    """
    if "X_umap" not in adata.obsm:
        return None
    if obs_col not in adata.obs.columns:
        return None

    umap = adata.obsm["X_umap"]
    labels = adata.obs[obs_col].astype(str)
    unique_labels = sorted(labels.unique())
    n = len(unique_labels)

    cmap = plt.get_cmap("tab20", max(n, 1))
    color_map = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}
    colors = [color_map[lbl] for lbl in labels]

    # Scale figure height to accommodate large legends
    fig_w = 9 + max(0, (n - 20) // 10) * 1.5
    fig_h = max(6, n * 0.22)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.scatter(umap[:, 0], umap[:, 1], c=colors, s=1.5, alpha=0.6, rasterized=True)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color_map[lbl], label=lbl) for lbl in unique_labels]
    ncol = max(1, n // 25)
    ax.legend(
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

    ax.set_title(title, fontsize=13, fontweight="bold")
    if subtitle:
        ax.set_xlabel(subtitle, fontsize=9, color="#666")
    else:
        ax.set_xlabel("UMAP 1", fontsize=10)
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_side_by_side(
    adata: AnnData,
    col_left: str,
    col_right: str,
    title_left: str,
    title_right: str,
    suptitle: str,
) -> Optional[str]:
    """
    Side-by-side UMAP: predicted labels (left) vs ground-truth (right).
    Returns None if X_umap absent or both columns absent.
    """
    if "X_umap" not in adata.obsm:
        return None
    if col_left not in adata.obs.columns and col_right not in adata.obs.columns:
        return None

    umap = adata.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    for ax, col, title in [(axes[0], col_left, title_left),
                            (axes[1], col_right, title_right)]:
        if col not in adata.obs.columns:
            ax.text(0.5, 0.5, f"'{col}' not in obs",
                    ha="center", va="center", transform=ax.transAxes, fontsize=11)
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            continue

        labels = adata.obs[col].astype(str)
        unique_labels = sorted(labels.unique())
        n = len(unique_labels)
        cmap = plt.get_cmap("tab20", max(n, 1))
        color_map = {lbl: cmap(i) for i, lbl in enumerate(unique_labels)}
        colors = [color_map[lbl] for lbl in labels]

        ax.scatter(umap[:, 0], umap[:, 1], c=colors, s=0.5, alpha=0.5, rasterized=True)
        ax.set_title(f"{title}  ({n} types)", fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_xticks([]); ax.set_yticks([])

        from matplotlib.patches import Patch
        legend_el = [Patch(facecolor=color_map[lbl], label=lbl) for lbl in unique_labels]
        ax.legend(handles=legend_el, bbox_to_anchor=(1.01, 1), loc="upper left",
                  frameon=False, fontsize=6.5, handlelength=1.0)

    plt.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    return _fig_to_b64(fig)


def _plot_marker_heatmap(score_df: pd.DataFrame) -> Optional[str]:
    """
    Heatmap of marker-gene scores (clusters × cell types).
    Highest-scoring cell type per cluster is highlighted.
    Returns None if score_df is None or has no score columns.
    """
    if score_df is None:
        return None

    # Drop the best_by_score annotation column for the heatmap matrix
    score_cols = [c for c in score_df.columns if c != "best_by_score"]
    if not score_cols:
        return None

    mat = score_df[score_cols].astype(float)
    n_clusters, n_types = mat.shape

    fig_w = max(10, n_types * 0.7)
    fig_h = max(5, n_clusters * 0.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(mat.values, aspect="auto", cmap="YlOrRd", interpolation="nearest")

    # Tick labels
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(score_cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(mat.index.tolist(), fontsize=8)

    # Highlight the best-scoring cell type per cluster with a white star
    if "best_by_score" in score_df.columns:
        for row_idx, cluster_id in enumerate(mat.index):
            best_ct = score_df.loc[cluster_id, "best_by_score"]
            if best_ct in score_cols:
                col_idx = score_cols.index(best_ct)
                ax.text(col_idx, row_idx, "★", ha="center", va="center",
                        fontsize=9, color="white", fontweight="bold")

    plt.colorbar(im, ax=ax, label="Mean log-normalised expression", shrink=0.6)
    ax.set_title("Marker Gene Score Heatmap — Clusters × Cell Types\n"
                 "(★ = best-scoring type per cluster)",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Cell Type", fontsize=10)
    ax.set_ylabel("Cluster", fontsize=10)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_confidence_distribution(adata: AnnData) -> Optional[str]:
    """
    Histogram of per-cell confidence scores for the consensus vote.
    Returns None if cell_type_confidence not in obs.
    """
    if "cell_type_confidence" not in adata.obs.columns:
        return None

    conf = adata.obs["cell_type_confidence"].astype(float)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(conf, bins=20, color="#4C78A8", edgecolor="white", alpha=0.85)

    # Threshold lines
    for thresh, color, label in [
        (0.5,  "#E07B3A", "Low confidence (<0.5)"),
        (0.75, "#2CA02C", "High confidence (≥0.75)"),
    ]:
        ax.axvline(thresh, color=color, linewidth=1.5, linestyle="--", label=label)

    median_conf = conf.median()
    ax.axvline(median_conf, color="#D62728", linewidth=1.5, linestyle=":",
               label=f"Median ({median_conf:.2f})")

    low_pct  = 100 * (conf < 0.5).mean()
    high_pct = 100 * (conf >= 0.75).mean()
    ax.set_title(
        f"Consensus Vote Confidence Distribution\n"
        f"Low (<0.5): {low_pct:.1f}%  |  High (≥0.75): {high_pct:.1f}%",
        fontsize=12, fontweight="bold",
    )
    ax.set_xlabel("Confidence Score (fraction of methods agreeing)", fontsize=11)
    ax.set_ylabel("Number of Cells", fontsize=11)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
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
<title>OmicSage — Annotation Report: {dataset_name}</title>
<style>
  body       {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI",
                sans-serif; margin: 2em auto; max-width: 1100px;
                color: #2c3e50; background: #fafafa; }}
  h1         {{ color: #1a252f; border-bottom: 2px solid #4C78A8;
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
  .warn      {{ background:#fff3cd; color:#856404; }}
  footer     {{ margin-top: 3em; font-size: 0.8em; color: #999;
                border-top: 1px solid #eee; padding-top: 12px; }}
</style>
</head>
<body>

<h1>🧬 OmicSage — Cell-Type Annotation Report</h1>
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

<h2>Cluster → Cell Type Assignment</h2>
<p class="meta">
  Per-cluster label from each method run.
  The <em>Consensus</em> column is the majority-vote winner.
  Confidence = fraction of active methods agreeing.
</p>
<table>
  <tr>
    <th>Cluster</th><th>n cells</th>
    <th>CellTypist coarse</th><th>CellTypist fine</th>
    <th>Marker score</th>
    <th>Consensus</th><th>Confidence</th>
  </tr>
  {assignment_rows}
</table>

<h2>Figures</h2>
<div class="fig-grid">
  {celltypist_section}
  {vote_section}
  {groundtruth_section}
  {heatmap_section}
  {confidence_section}
</div>

<h2>Provenance</h2>
<table>
  <tr><th>Key</th><th>Value</th></tr>
  {provenance_rows}
</table>

<footer>
  Generated by OmicSage · reports/annotate_report.py ·
  <a href="https://github.com/fshokor/OmicSage">github.com/fshokor/OmicSage</a>
</footer>
</body>
</html>
"""

_FIG_BOX_WIDE = """\
  <div class="fig-box wide">
    <strong>{title}</strong>
    <p class="meta">{caption}</p>
    <img src="data:image/png;base64,{b64}" alt="{alt}">
  </div>"""

_FIG_BOX_HALF = """\
  <div class="fig-box">
    <strong>{title}</strong>
    <p class="meta">{caption}</p>
    <img src="data:image/png;base64,{b64}" alt="{alt}">
  </div>"""


def _build_html(
    adata_annotated: AnnData,
    annotation_dict: dict,
    dataset_name: str,
) -> str:
    """Render all figures and assemble the HTML string."""
    prov = adata_annotated.uns.get("omicsage_annotate", {})
    methods_run = prov.get("methods_run", [])
    leiden_col  = prov.get("leiden_col", "leiden")

    # ── Figures ───────────────────────────────────────────────────────────────

    print("  Rendering CellTypist UMAP panels ...", flush=True)
    b64_coarse = _plot_umap_labels(
        adata_annotated, "celltypist_coarse",
        title="UMAP — CellTypist Coarse Labels (Immune_All_High)",
        subtitle="Majority-voted label per Leiden cluster · coarse immune types",
    )
    b64_fine = _plot_umap_labels(
        adata_annotated, "celltypist_fine",
        title="UMAP — CellTypist Fine Labels (Immune_All_Low)",
        subtitle="Majority-voted label per Leiden cluster · fine-grained subtypes",
    )
    b64_markers = _plot_umap_labels(
        adata_annotated, "cell_type_markers",
        title="UMAP — Marker Gene Score Labels",
        subtitle="Best-scoring cell type per cluster by mean log-normalised expression",
    )

    celltypist_boxes = ""
    for b64, title, caption, alt in [
        (b64_coarse,
         "CellTypist — Coarse Labels (Immune_All_High)",
         "Cluster-level majority vote from the coarse ~30-type model.",
         "UMAP celltypist coarse"),
        (b64_fine,
         "CellTypist — Fine Labels (Immune_All_Low)",
         "Cluster-level majority vote from the fine-grained ~200-subtype model.",
         "UMAP celltypist fine"),
        (b64_markers,
         "Marker Gene Score Labels",
         "Each cluster labelled by the cell type with highest mean log-expression "
         "across its canonical marker genes.  Built-in MARKER_SETS used.",
         "UMAP marker scores"),
    ]:
        if b64 is not None:
            celltypist_boxes += _FIG_BOX_HALF.format(
                title=title, caption=caption, b64=b64, alt=alt,
            )

    # Consensus vote UMAP
    vote_section = ""
    if "vote" in methods_run and "cell_type_vote" in adata_annotated.obs.columns:
        print("  Rendering consensus vote UMAP ...", flush=True)
        b64_consensus = _plot_umap_labels(
            adata_annotated, "cell_type_vote",
            title="UMAP — Consensus Cell Type (Majority Vote)",
            subtitle="2-way vote: CellTypist fine + marker score  "
                     "(ScType and SingleR slots reserved for Session B)",
        )
        if b64_consensus is not None:
            vote_section = _FIG_BOX_WIDE.format(
                title="Consensus Cell Type (Majority Vote)",
                caption="2-way vote: CellTypist fine label + marker gene score. "
                        "ScType and SingleR will be added in Session B. "
                        "Confidence = fraction of methods agreeing (0.0–1.0).",
                b64=b64_consensus,
                alt="UMAP consensus",
            )

    # Ground-truth side-by-side: consensus vs ground-truth
    groundtruth_section = ""
    print("  Rendering ground-truth UMAP (if available) ...", flush=True)
    if "cell_type_groundtruth" in adata_annotated.obs.columns:
        # Side-by-side: consensus (left) vs ground-truth (right)
        b64_gt_side = _plot_umap_side_by_side(
            adata_annotated,
            col_left="cell_type_vote",
            col_right="cell_type_groundtruth",
            title_left="Consensus Vote",
            title_right="Ground-Truth (publication)",
            suptitle="Consensus Vote vs Ground-Truth Cell Type",
        )
        if b64_gt_side is not None:
            groundtruth_section = _FIG_BOX_WIDE.format(
                title="Consensus vs Ground-Truth Cell Type",
                caption="Left: majority-vote consensus label (our result). "
                        "Right: original publication cell_type label preserved in "
                        "obs['cell_type_groundtruth'] before annotation ran.",
                b64=b64_gt_side,
                alt="UMAP consensus vs ground truth",
            )

    # Marker heatmap
    heatmap_section = ""
    score_df = annotation_dict.get("marker_score_df")
    if score_df is not None:
        print("  Rendering marker score heatmap ...", flush=True)
        b64_hm = _plot_marker_heatmap(score_df)
        if b64_hm is not None:
            heatmap_section = _FIG_BOX_WIDE.format(
                title="Marker Gene Score Heatmap (Clusters × Cell Types)",
                caption="Mean log-normalised expression of canonical marker gene sets "
                        "per Leiden cluster.  ★ = best-scoring type for that cluster.  "
                        "Rows = clusters, columns = cell types.",
                b64=b64_hm,
                alt="Marker heatmap",
            )

    # Confidence histogram
    confidence_section = ""
    if "vote" in methods_run:
        print("  Rendering confidence distribution ...", flush=True)
        b64_conf = _plot_confidence_distribution(adata_annotated)
        if b64_conf is not None:
            confidence_section = _FIG_BOX_HALF.format(
                title="Consensus Confidence Distribution",
                caption="Per-cell confidence = fraction of active methods that agreed. "
                        "Cells below 0.5 may warrant manual review.",
                b64=b64_conf,
                alt="Confidence distribution",
            )

    # ── Summary table ─────────────────────────────────────────────────────────
    n_cells    = adata_annotated.n_obs
    n_clusters = adata_annotated.obs[leiden_col].nunique() if leiden_col in adata_annotated.obs else "?"
    n_ct_coarse  = adata_annotated.obs["celltypist_coarse"].nunique() \
                   if "celltypist_coarse" in adata_annotated.obs else "not run"
    n_ct_fine    = adata_annotated.obs["celltypist_fine"].nunique() \
                   if "celltypist_fine" in adata_annotated.obs else "not run"
    n_ct_markers = adata_annotated.obs["cell_type_markers"].nunique() \
                   if "cell_type_markers" in adata_annotated.obs else "not run"
    n_consensus  = adata_annotated.obs["cell_type"].nunique() \
                   if "cell_type" in adata_annotated.obs else "not run"

    conf_median = "—"
    if "cell_type_confidence" in adata_annotated.obs.columns:
        conf_median = f"{adata_annotated.obs['cell_type_confidence'].median():.3f}"

    summary_items = [
        ("Cells",                           f"{n_cells:,}"),
        ("Leiden clusters",                 str(n_clusters)),
        ("Methods run",                     ", ".join(methods_run) if methods_run else "—"),
        ("CellTypist coarse types",         str(n_ct_coarse)),
        ("CellTypist fine types",           str(n_ct_fine)),
        ("Marker-score types",              str(n_ct_markers)),
        ("Consensus types (vote)",          str(n_consensus)),
        ("Median consensus confidence",     conf_median),
        ("CellTypist models dir",           prov.get("celltypist_models_dir", "—")),
        ("CellTypist models",               ", ".join(prov.get("celltypist_models", []))),
        ("Marker sets used",                str(len(prov.get("marker_sets_keys", [])))),
    ]
    summary_rows = "\n  ".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in summary_items
    )

    # ── Assignment table ───────────────────────────────────────────────────────
    clusters = sorted(
        adata_annotated.obs[leiden_col].unique() if leiden_col in adata_annotated.obs else [],
        key=lambda x: (int(x) if str(x).isdigit() else x),
    )

    def _cluster_mode(col, cl):
        if col not in adata_annotated.obs.columns:
            return "—"
        mask = adata_annotated.obs[leiden_col] == cl
        vc = adata_annotated.obs.loc[mask, col].value_counts()
        return vc.index[0] if len(vc) > 0 else "—"

    def _cluster_conf(cl):
        if "cell_type_confidence" not in adata_annotated.obs.columns:
            return "—"
        mask = adata_annotated.obs[leiden_col] == cl
        vals = adata_annotated.obs.loc[mask, "cell_type_confidence"].astype(float)
        med  = vals.median()
        badge_cls = "ok" if med >= 0.75 else "warn"
        return f'<span class="badge {badge_cls}">{med:.2f}</span>'

    asgn_rows = []
    for cl in clusters:
        mask   = adata_annotated.obs[leiden_col] == cl
        n_cl   = int(mask.sum())
        coarse = _cluster_mode("celltypist_coarse", cl)
        fine   = _cluster_mode("celltypist_fine",   cl)
        marker = _cluster_mode("cell_type_markers", cl)
        consns = _cluster_mode("cell_type",          cl)
        conf   = _cluster_conf(cl)
        asgn_rows.append(
            f"<tr><td>{cl}</td><td>{n_cl:,}</td>"
            f"<td>{coarse}</td><td>{fine}</td><td>{marker}</td>"
            f"<td><strong>{consns}</strong></td><td>{conf}</td></tr>"
        )
    assignment_rows = "\n  ".join(asgn_rows)

    # ── Provenance table ───────────────────────────────────────────────────────
    prov_display = {}
    for k, v in prov.items():
        if isinstance(v, list):
            prov_display[k] = ", ".join(str(x) for x in v)
        elif isinstance(v, dict):
            for subk, subv in v.items():
                prov_display[f"{k}[{subk}]"] = subv
        else:
            prov_display[k] = v
    provenance_rows = "\n  ".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov_display.items()
    )

    return _HTML_TEMPLATE.format(
        dataset_name=dataset_name,
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M"),
        summary_rows=summary_rows,
        assignment_rows=assignment_rows,
        provenance_rows=provenance_rows,
        celltypist_section=celltypist_boxes,
        vote_section=vote_section,
        groundtruth_section=groundtruth_section,
        heatmap_section=heatmap_section,
        confidence_section=confidence_section,
    )


# ---------------------------------------------------------------------------
# Public API — callable from notebook
# ---------------------------------------------------------------------------

def run_annotate_report(
    adata_annotated: AnnData,
    annotation_dict: dict,
    report_path: str = "reports/annotate_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the annotation HTML report and write it to disk.

    Parameters
    ----------
    adata_annotated : AnnData
        Annotated AnnData returned by ``annotate()``.
        Must have obsm['X_umap'] and at least one annotation obs column.
    annotation_dict : dict
        Second return value of ``annotate()``.
        May contain 'marker_score_df' and 'vote_df' for additional figures.
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

    print(f"Building annotation report for '{dataset_name}' ...", flush=True)
    html = _build_html(adata_annotated, annotation_dict, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved → {out.resolve()}", flush=True)
    return str(out.resolve())


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate OmicSage cell-type annotation report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="Clustered h5ad (obs['leiden'] set)")
    p.add_argument("--output",  default=None,
                   help="Where to save the annotated h5ad "
                        "(if omitted, annotate runs but h5ad is not saved)")
    p.add_argument("--report",  default="reports/annotate_report.html",
                   help="Path for the output HTML report")
    p.add_argument("--dataset", default=None,
                   help="Dataset label shown in the report (default: input filename)")
    p.add_argument("--methods", nargs="+",
                   default=["celltypist", "markers", "vote"],
                   help="Annotation methods to run")
    return p.parse_args()


def main() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s")
    args = _parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    from pipeline.modules.qc.annotate import annotate

    dataset_name = args.dataset or Path(args.input).stem

    print(f"Loading {args.input} ...", flush=True)
    adata = sc.read_h5ad(args.input)
    print(adata, flush=True)

    print("Running annotation ...", flush=True)
    adata_annotated, annotation_dict = annotate(
        adata,
        methods=args.methods,
    )

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        adata_annotated.write_h5ad(args.output)
        print(f"Annotated h5ad saved → {args.output}", flush=True)

    run_annotate_report(
        adata_annotated=adata_annotated,
        annotation_dict=annotation_dict,
        report_path=args.report,
        dataset_name=dataset_name,
    )


if __name__ == "__main__":
    main()
