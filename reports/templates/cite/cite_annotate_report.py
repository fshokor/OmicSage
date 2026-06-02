"""
OmicSage — CITE-seq ADT Annotation Report
reports/templates/cite/cite_annotate_report.py

Generated after annotate_adt step.
Output: cite_05_annotate_report.html

Sections
--------
* Run Summary            — cell / cluster / resolution stats
* Cluster Table          — per-cluster cell counts and assigned cell types
* Cluster Overview       — cluster size bar chart + UMAP coloured by Leiden
* Dotplot                — sc.pl.rank_genes_groups_dotplot (top 3 markers/cluster)
* Marker UMAP grid       — small-multiple UMAPs, one per key marker protein,
                           coloured by DSB/CLR expression
* Cell Type UMAP         — UMAP coloured by adt_celltype (scoring result)
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
import scanpy as sc
from anndata import AnnData

_DPI = 130
_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px; box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
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
    td { padding: 8px 12px; border-bottom: 1px solid #eee; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9ff; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 520px; }
    .fig-wrap.wide { flex: 1 1 100%; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
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


def _get_expression_matrix(adt: AnnData) -> tuple[np.ndarray, list[str]]:
    """Return (cells × proteins ndarray, protein names).

    Prefers 'adt_dsb' layer → 'adt_clr' layer → .X in that order.
    """
    var_names = list(adt.var_names)
    for layer in ("adt_dsb", "adt_clr"):
        if layer in adt.layers:
            X = adt.layers[layer]
            if hasattr(X, "toarray"):
                X = X.toarray()
            return np.array(X, dtype=float), var_names
    # fall back to .X
    X = adt.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.array(X, dtype=float), var_names


# ---------------------------------------------------------------------------
# Existing figures (unchanged logic, same as before)
# ---------------------------------------------------------------------------

def _plot_cluster_sizes(adt: AnnData) -> str:
    if "leiden" not in adt.obs.columns:
        return _placeholder("leiden column not found")

    counts = adt.obs["leiden"].value_counts().sort_index()
    labels = [str(k) for k in counts.index]
    values = counts.values
    cmap   = plt.get_cmap("tab20", max(len(labels), 2))

    fig, ax = plt.subplots(figsize=(max(6, len(labels) * 0.7), 4))
    bars = ax.bar(labels, values,
                  color=[cmap(i) for i in range(len(labels))],
                  width=0.7, alpha=0.85, edgecolor="white")
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01,
                f"{val:,}", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("Leiden Cluster", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("Leiden Cluster Sizes", fontsize=12, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_leiden(adt: AnnData) -> str:
    """UMAP coloured by Leiden cluster ID only — no cell type labels."""
    if "X_umap_adt" not in adt.obsm or "leiden" not in adt.obs.columns:
        return _placeholder("UMAP or leiden not found")

    umap   = adt.obsm["X_umap_adt"]
    labels = adt.obs["leiden"].astype(str)
    unique = sorted(labels.unique(), key=lambda x: int(x) if x.isdigit() else x)
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, cid in enumerate(unique):
        mask = (labels == cid).values
        ax.scatter(umap[mask, 0], umap[mask, 1],
                   s=2, alpha=0.7, color=cmap(i),
                   label=f"Cluster {cid}", rasterized=True)

    ax.legend(markerscale=4, frameon=False, fontsize=8,
              loc="upper right", ncol=max(1, len(unique) // 12),
              title="Leiden", title_fontsize=8)
    ax.set_title("ADT UMAP — Leiden Clusters", fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# NEW: rank_genes_groups dotplot
# ---------------------------------------------------------------------------

def _plot_rank_genes_dotplot(adt: AnnData, n_genes: int = 3) -> str:
    """
    sc.pl.rank_genes_groups_dotplot — top marker proteins per Leiden cluster.

    Always uses groupby="leiden" because rank_genes_groups and dendrogram
    are computed on leiden in annotate_adt(). This dotplot is a diagnostic:
    inspect which proteins define each cluster, then decide on annotation.
    Cluster IDs on the y-axis can be cross-referenced against the cluster
    table to map them to cell type labels.
    """
    if "rank_genes_groups" not in adt.uns or "leiden" not in adt.obs.columns:
        return _placeholder("rank_genes_groups not found — run annotate_adt() first")

    n_clusters = adt.obs["leiden"].nunique()
    fig, ax = plt.subplots(figsize=(14, max(4, n_clusters * 0.6)))
    try:
        sc.pl.rank_genes_groups_dotplot(
            adt,
            n_genes=n_genes,
            groupby="leiden",
            ax=ax,
            show=False,
        )
        ax.set_title(
            f"Top {n_genes} marker proteins per Leiden cluster "
            f"(use this to guide annotation)",
            fontsize=11, fontweight="bold", pad=10,
        )
    except Exception as exc:
        plt.close(fig)
        return _placeholder(f"Dotplot failed: {exc}")

    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# NEW: Marker UMAP grid
# ---------------------------------------------------------------------------

def _plot_marker_umap_grid(
    adt: AnnData,
    marker_panel: Optional[dict[str, list[str]]] = None,
    n_cols: int = 4,
    max_markers: int = 16,
) -> str:
    """
    Small-multiple UMAPs, one panel per key marker protein, coloured by
    DSB-normalised (or CLR) expression.

    Selects up to ``max_markers`` representative proteins from
    ``marker_panel``.  Falls back to the most variable proteins if the
    named markers are not found in ``var_names``.
    """
    if "X_umap_adt" not in adt.obsm:
        return _placeholder("X_umap_adt not found")

    # ---- resolve marker list ------------------------------------------------
    # Flatten marker_panel to a de-duplicated ordered list of protein names
    if marker_panel:
        ordered: list[str] = []
        seen: set[str] = set()
        for markers in marker_panel.values():
            for m in markers:
                if m not in seen:
                    ordered.append(m)
                    seen.add(m)
    else:
        ordered = list(adt.var_names)

    # Keep only markers that exist in var_names
    var_set = set(adt.var_names)
    available = [m for m in ordered if m in var_set]

    # If too few matches, fall back to top variable proteins
    if len(available) < 4:
        X_mat, var_names = _get_expression_matrix(adt)
        variances = X_mat.var(axis=0)
        top_idx = np.argsort(variances)[::-1][:max_markers]
        available = [var_names[i] for i in top_idx]

    proteins = available[:max_markers]

    # ---- expression matrix --------------------------------------------------
    X_mat, var_names = _get_expression_matrix(adt)
    var_index = {v: i for i, v in enumerate(var_names)}

    umap = adt.obsm["X_umap_adt"]
    n_proteins = len(proteins)
    n_rows = (n_proteins + n_cols - 1) // n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 3.5, n_rows * 3.2),
        squeeze=False,
    )
    fig.suptitle(
        "ADT UMAP — key marker proteins (DSB / CLR expression)",
        fontsize=13, fontweight="bold", y=1.01,
    )

    for idx, protein in enumerate(proteins):
        row, col = divmod(idx, n_cols)
        ax = axes[row][col]

        if protein not in var_index:
            ax.set_visible(False)
            continue

        expr = X_mat[:, var_index[protein]]
        # Clip extreme outliers for better colour scale
        vmin, vmax = np.percentile(expr, [2, 98])
        sc_plot = ax.scatter(
            umap[:, 0], umap[:, 1],
            c=expr, cmap="RdBu_r",
            vmin=vmin, vmax=vmax,
            s=1, alpha=0.6, rasterized=True,
        )
        ax.set_title(protein, fontsize=9, fontweight="bold")
        ax.set_xticks([]); ax.set_yticks([])
        ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
        plt.colorbar(sc_plot, ax=ax, fraction=0.04, pad=0.02)

    # Hide unused axes
    for idx in range(n_proteins, n_rows * n_cols):
        row, col = divmod(idx, n_cols)
        axes[row][col].set_visible(False)

    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# NEW: UMAP coloured by assigned cell type
# ---------------------------------------------------------------------------

def _plot_umap_celltype_score(adt: AnnData) -> str:
    """
    UMAP coloured by obs["adt_celltype_score"] — the auto-scored labels
    derived from marker panel fold-change scoring, always computed
    regardless of whether an explicit annotation_map was provided.
    Allows comparison between scoring result and manual annotation.
    """
    if "X_umap_adt" not in adt.obsm:
        return _placeholder("X_umap_adt not found")
    if "adt_celltype_score" not in adt.obs.columns:
        return _placeholder(
            "adt_celltype_score not found\n"
            "Run annotate_adt() with preset='bmmc' or a custom marker_panel"
        )

    umap   = adt.obsm["X_umap_adt"]
    labels = adt.obs["adt_celltype_score"].astype(str)
    unique = sorted(labels.unique())
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(8, 6))
    for i, ct in enumerate(unique):
        mask = (labels == ct).values
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            s=2, alpha=0.7, color=cmap(i),
            label=ct, rasterized=True,
        )
    ax.legend(
        markerscale=5, frameon=False, fontsize=8,
        loc="upper right", ncol=max(1, len(unique) // 10),
        title="Cell type (score)", title_fontsize=8,
    )
    ax.set_title(
        "ADT UMAP — Auto-scored (adt_celltype_score)",
        fontsize=12, fontweight="bold",
    )
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_celltype(adt: AnnData) -> str:
    """
    UMAP coloured by obs["adt_celltype_manual"] — the manual annotation
    applied from the explicit annotation_map in config.
    Returns a placeholder if adt_celltype_manual or X_umap_adt are absent.
    """
    if "X_umap_adt" not in adt.obsm:
        return _placeholder("X_umap_adt not found")
    if "adt_celltype_manual" not in adt.obs.columns:
        return _placeholder(
            "adt_celltype_manual not found\n"
            "Provide annotation_map in config to enable manual annotation"
        )

    umap      = adt.obsm["X_umap_adt"]
    labels    = adt.obs["adt_celltype_manual"].astype(str)
    unique    = sorted(labels.unique())
    cmap      = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(8, 6))
    for i, ct in enumerate(unique):
        mask = (labels == ct).values
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            s=2, alpha=0.7, color=cmap(i),
            label=ct, rasterized=True,
        )

    ax.legend(
        markerscale=5, frameon=False, fontsize=8,
        loc="upper right", ncol=max(1, len(unique) // 10),
        title="Cell type", title_fontsize=8,
    )
    ax.set_title(
        "ADT UMAP — Manual annotation (adt_celltype_manual)",
        fontsize=12, fontweight="bold",
    )
    ax.set_xlabel("UMAP 1", fontsize=9)
    ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
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
  <title>OmicSage -- CITE-seq ADT Annotation Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq ADT Annotation Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, dataset_name: str, timestamp: str) -> str:
    cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",      f"{metrics.get('n_cells', '?'):,}"),
            ("Clusters",   str(metrics.get("n_clusters", "?"))),
            ("Resolution", str(metrics.get("resolution", "?"))),
            ("Annotated",  "Yes" if metrics.get("annotated") else "No"),
        ]
    )
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Leiden key",    metrics.get("leiden_key", "leiden")),
            ("Cell type key", metrics.get("celltype_key") or "not applied"),
            ("Random state",  str(metrics.get("random_state", 0))),
        ]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{cards}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows}</tbody></table>
    </section>"""


def _section_cluster_table(adt: AnnData, metrics: dict) -> str:
    if "leiden" not in adt.obs.columns:
        return ""
    cluster_sizes = metrics.get("cluster_sizes", {})
    has_ct = "adt_celltype" in adt.obs.columns
    header = "<tr><th>Cluster</th><th>Cells</th>"
    header += "<th>Cell Type</th>" if has_ct else ""
    header += "<th>% of total</th></tr>"
    total  = sum(cluster_sizes.values()) or 1
    rows   = ""
    for cid in sorted(cluster_sizes.keys(), key=lambda x: int(x) if x.isdigit() else x):
        cnt = cluster_sizes[cid]
        ct_cell = ""
        if has_ct:
            ct_vals = adt.obs.loc[adt.obs["leiden"] == str(cid), "adt_celltype"].unique()
            ct = ct_vals[0] if len(ct_vals) == 1 else ", ".join(ct_vals[:3])
            ct_cell = f"<td>{ct}</td>"
        pct = 100.0 * cnt / total
        rows += f"<tr><td>{cid}</td><td>{cnt:,}</td>{ct_cell}<td>{pct:.1f}%</td></tr>"

    note = ""
    if not metrics.get("annotated"):
        note = """<p class="note">No <code>annotation_map</code> was provided.
            Inspect clusters above, then add an <code>annotation_map</code>
            to the config and re-run this step with
            <code>--step annotate_adt --force</code>.</p>"""

    return f"""
    <section>
      <h2>Cluster Table</h2>
      {note}
      <table><thead>{header}</thead><tbody>{rows}</tbody></table>
    </section>"""


def _section_cluster_figures(fig_bar: str, fig_umap: str) -> str:
    return f"""
    <section>
      <h2>Cluster Overview</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Cluster Sizes</h3>
          <img src="data:image/png;base64,{fig_bar}" alt="Cluster sizes">
        </div>
        <div class="fig-wrap">
          <h3>ADT UMAP — Leiden Clusters</h3>
          <img src="data:image/png;base64,{fig_umap}" alt="UMAP Leiden">
        </div>
      </div>
    </section>"""


def _section_dotplot(fig_dot: str) -> str:
    return f"""
    <section>
      <h2>Top Marker Proteins per Cluster (Dotplot)</h2>
      <p>
        Top 3 marker proteins per Leiden cluster ranked by
        <code>sc.tl.rank_genes_groups</code> (Wilcoxon).
        Dot size = fraction of cells expressing the protein;
        dot colour = mean expression (DSB / CLR normalised).
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_dot}" alt="Dotplot rank genes">
        </div>
      </div>
    </section>"""


def _section_marker_umaps(fig_grid: str) -> str:
    return f"""
    <section>
      <h2>Key Marker Protein Expression on UMAP</h2>
      <p>
        Each panel shows the ADT UMAP coloured by DSB-normalised (or CLR)
        expression of a single lineage marker. Colour scale clipped to the
        2nd–98th percentile. Red = high expression, blue = low.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_grid}" alt="Marker UMAP grid">
        </div>
      </div>
    </section>"""


def _section_celltype_umap(
    fig_celltype_score: Optional[str],
    fig_celltype_manual: Optional[str] = None,
) -> str:
    """
    Run 1 (no annotation_map): shows score UMAP only.
    Run 2 (annotation_map provided): shows score UMAP + manual UMAP side by side.
    """
    score_panel = ""
    if fig_celltype_score:
        score_panel = (
            '<div class="fig-wrap wide">'
            '<h3>Auto-scored labels (adt_celltype_score)</h3>'
            '<p style="font-size:0.85rem;color:#555;margin-bottom:8px;">'
            'Derived automatically from marker panel fold-change scoring. '
            'Always computed — use this to verify or refine the annotation map.'
            '</p>'
            f'<img src="data:image/png;base64,{fig_celltype_score}" alt="UMAP score">'
            '</div>'
        )

    manual_panel = ""
    if fig_celltype_manual:
        manual_panel = (
            '<div class="fig-wrap wide">'
            '<h3>Manual annotation (adt_celltype_manual)</h3>'
            '<p style="font-size:0.85rem;color:#555;margin-bottom:8px;">'
            'Applied from the explicit <code>annotation_map</code> in config. '
            'Compare with the score UMAP above to validate your map.'
            '</p>'
            f'<img src="data:image/png;base64,{fig_celltype_manual}" alt="UMAP manual">'
            '</div>'
        )

    return f"""
    <section>
      <h2>ADT UMAP — Cell Type Annotation</h2>
      <div class="fig-grid">
        {score_panel}
        {manual_panel}
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_cite_annotate_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_05_annotate_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the CITE-seq ADT annotation HTML report.

    Parameters
    ----------
    adt : AnnData
        ADT AnnData returned by annotate_adt().
    metrics : dict
        Metrics dict returned by annotate_adt().
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

    print(f"Building CITE annotate report for '{dataset_name}' ...", flush=True)

    # Recover marker panel from uns if not in metrics
    marker_panel: Optional[dict] = (
        metrics.get("marker_panel")
        or adt.uns.get("omicsage_adt_marker_panel")
    )

    # -- render figures -------------------------------------------------------
    print("  • cluster size bar chart ...", flush=True)
    fig_bar         = _plot_cluster_sizes(adt)

    print("  • UMAP Leiden ...", flush=True)
    fig_umap        = _plot_umap_leiden(adt)

    print("  • rank_genes_groups dotplot ...", flush=True)
    fig_dot         = _plot_rank_genes_dotplot(adt)

    print("  • marker UMAP grid ...", flush=True)
    fig_marker_grid = _plot_marker_umap_grid(adt, marker_panel=marker_panel)

    fig_celltype_score  = None
    fig_celltype_manual = None

    if "adt_celltype_score" in adt.obs.columns:
        print("  • UMAP cell type (score) ...", flush=True)
        fig_celltype_score = _plot_umap_celltype_score(adt)

    if "adt_celltype_manual" in adt.obs.columns:
        print("  • UMAP cell type (manual) ...", flush=True)
        fig_celltype_manual = _plot_umap_celltype(adt)

    # -- assemble HTML --------------------------------------------------------
    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_cluster_table(adt, metrics),
            _section_cluster_figures(fig_bar, fig_umap),
            _section_dotplot(fig_dot),
            _section_marker_umaps(fig_marker_grid),
            _section_celltype_umap(fig_celltype_score, fig_celltype_manual),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
