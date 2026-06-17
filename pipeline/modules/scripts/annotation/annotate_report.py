"""
OmicSage -- Cell-Type Annotation Report
pipeline/modules/annotation/annotate_report.py

Usage
-----
    from pipeline.modules.annotation.annotate_report import run_annotate_report
    run_annotate_report(
        adata_annotated=adata_annotated,
        annotation_dict=annotation_dict,
        report_path="reports/GSE194122/05_annotate_report.html",
        dataset_name="GSE194122_CITE",
    )
"""

from __future__ import annotations

import base64
import io
import logging
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
# Figure helpers — plotting logic unchanged from original
# ---------------------------------------------------------------------------

def _fig_to_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_umap_labels(adata: AnnData, obs_col: str, title: str,
                      subtitle: str = "") -> Optional[str]:
    if "X_umap" not in adata.obsm or obs_col not in adata.obs.columns:
        return None
    umap   = adata.obsm["X_umap"]
    labels = adata.obs[obs_col].astype(str)
    unique = sorted(labels.unique())
    n      = len(unique)
    cmap_  = plt.get_cmap("tab20", max(n, 1))
    col_map = {l: cmap_(i) for i, l in enumerate(unique)}
    fig, ax = plt.subplots(figsize=(9 + max(0, (n - 20) // 10) * 1.5, max(6, n * 0.22)))
    ax.scatter(umap[:, 0], umap[:, 1], c=[col_map[l] for l in labels],
               s=1.5, alpha=0.6, rasterized=True)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=col_map[l], label=l) for l in unique],
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0,
              frameon=False, fontsize=7.5, ncol=max(1, n // 25),
              handlelength=1.2, labelspacing=0.4)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xlabel(subtitle if subtitle else "UMAP 1", fontsize=9, color="#666")
    ax.set_ylabel("UMAP 2", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_umap_side_by_side(adata: AnnData, col_left: str, col_right: str,
                             title_left: str, title_right: str,
                             suptitle: str) -> Optional[str]:
    if "X_umap" not in adata.obsm:
        return None
    if col_left not in adata.obs.columns and col_right not in adata.obs.columns:
        return None
    umap    = adata.obsm["X_umap"]
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))
    for ax, col, title in [(axes[0], col_left, title_left), (axes[1], col_right, title_right)]:
        if col not in adata.obs.columns:
            ax.text(0.5, 0.5, f"'{col}' not in obs", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11)
            ax.set_title(title, fontsize=12, fontweight="bold")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        labels = adata.obs[col].astype(str)
        unique = sorted(labels.unique())
        cmap_  = plt.get_cmap("tab20", max(len(unique), 1))
        col_map = {l: cmap_(i) for i, l in enumerate(unique)}
        ax.scatter(umap[:, 0], umap[:, 1], c=[col_map[l] for l in labels],
                   s=0.5, alpha=0.5, rasterized=True)
        ax.set_title(f"{title}  ({len(unique)} types)", fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
        ax.set_xticks([]); ax.set_yticks([])
        from matplotlib.patches import Patch
        ax.legend(handles=[Patch(facecolor=col_map[l], label=l) for l in unique],
                  bbox_to_anchor=(1.01, 1), loc="upper left",
                  frameon=False, fontsize=6.5, handlelength=1.0)
    plt.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    return _fig_to_b64(fig)


def _plot_marker_heatmap(score_df: pd.DataFrame) -> Optional[str]:
    if score_df is None:
        return None
    score_cols = [c for c in score_df.columns if c != "best_by_score"]
    if not score_cols:
        return None
    mat = score_df[score_cols].astype(float)
    n_clusters, n_types = mat.shape
    fig, ax = plt.subplots(figsize=(max(10, n_types * 0.7), max(5, n_clusters * 0.4)))
    im = ax.imshow(mat.values, aspect="auto", cmap="YlOrRd", interpolation="nearest")
    ax.set_xticks(range(n_types))
    ax.set_xticklabels(score_cols, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(n_clusters))
    ax.set_yticklabels(mat.index.tolist(), fontsize=8)
    if "best_by_score" in score_df.columns:
        for row_idx, cluster_id in enumerate(mat.index):
            best_ct = score_df.loc[cluster_id, "best_by_score"]
            if best_ct in score_cols:
                ax.text(score_cols.index(best_ct), row_idx, "*",
                        ha="center", va="center", fontsize=9,
                        color="white", fontweight="bold")
    plt.colorbar(im, ax=ax, label="Mean log-normalised expression", shrink=0.6)
    ax.set_title("Marker Gene Score Heatmap -- Clusters x Cell Types\n"
                 "(* = best-scoring type per cluster)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Cell Type", fontsize=10); ax.set_ylabel("Cluster", fontsize=10)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_dotplot(adata: AnnData, leiden_col: str) -> Optional[str]:
    """
    Dot plot: clusters (Y) × top marker genes (X).
    Dot size  = fraction of cells in cluster expressing the gene (> 0).
    Dot colour = mean log-normalised expression across expressing cells.

    Collects the top-3 Wilcoxon marker genes per cluster from
    adata.uns['rank_genes_groups'] when available, then falls back to
    the top-3 highest-mean genes per cluster from logcounts.
    Requires at least one cluster and one gene to render.
    """
    if leiden_col not in adata.obs.columns:
        return None

    # ── Choose expression matrix ──────────────────────────────────────────────
    if "logcounts" in adata.layers:
        import scipy.sparse as sp
        X = adata.layers["logcounts"]
    else:
        import scipy.sparse as sp
        X = adata.X

    clusters = sorted(
        adata.obs[leiden_col].unique(),
        key=lambda x: (int(x) if str(x).isdigit() else x),
    )
    if not clusters:
        return None

    # ── Collect marker genes ──────────────────────────────────────────────────
    TOP_N = 3   # genes per cluster shown in the dot plot

    selected_genes: list[str] = []

    if "rank_genes_groups" in adata.uns:
        rgg = adata.uns["rank_genes_groups"]
        gene_names = rgg.get("names")
        if gene_names is not None:
            for cl in clusters:
                cl_str = str(cl)
                if cl_str in gene_names.dtype.names:
                    top = [str(g) for g in gene_names[cl_str][:TOP_N]
                           if str(g) in adata.var_names]
                    selected_genes.extend(top)
    
    # Fallback: pick top-3 highest mean-expressing genes per cluster
    if not selected_genes:
        for cl in clusters:
            mask = (adata.obs[leiden_col] == cl).values
            sub = X[np.where(mask)[0], :]
            if sp.issparse(sub):
                sub = sub.toarray()
            means = np.asarray(sub).mean(axis=0)
            top_idx = np.argsort(means)[::-1][:TOP_N]
            selected_genes.extend([adata.var_names[i] for i in top_idx])

    # Deduplicate while preserving order
    seen: set[str] = set()
    genes: list[str] = []
    for g in selected_genes:
        if g not in seen and g in adata.var_names:
            seen.add(g)
            genes.append(g)

    if not genes:
        return None

    # ── Compute per-cluster dot statistics ────────────────────────────────────
    gene_idx = [adata.var_names.get_loc(g) for g in genes]
    n_cl = len(clusters)
    n_g  = len(genes)

    mean_expr   = np.zeros((n_cl, n_g), dtype=float)
    pct_express = np.zeros((n_cl, n_g), dtype=float)

    for ci, cl in enumerate(clusters):
        mask = (adata.obs[leiden_col] == cl).values
        sub  = X[np.where(mask)[0], :][:, gene_idx]
        if sp.issparse(sub):
            sub = sub.toarray()
        sub = np.asarray(sub, dtype=float)
        expressed = sub > 0
        pct_express[ci, :] = expressed.mean(axis=0)
        # Mean expression only over expressing cells (avoid dilution by zeros)
        with np.errstate(invalid="ignore"):
            mean_expr[ci, :] = np.where(
                expressed.sum(axis=0) > 0,
                np.where(expressed, sub, np.nan).mean(axis=0),
                0.0,
            )

    # ── Draw figure ───────────────────────────────────────────────────────────
    fig_w = max(8, n_g * 0.55 + 2)
    fig_h = max(4, n_cl * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    max_expr = mean_expr.max() if mean_expr.max() > 0 else 1.0
    MAX_RADIUS = 0.42   # fraction of grid cell

    for ci in range(n_cl):
        for gi in range(n_g):
            pct  = pct_express[ci, gi]
            expr = mean_expr[ci, gi]
            radius = MAX_RADIUS * np.sqrt(pct)            # area ∝ pct
            color  = plt.cm.YlOrRd(expr / max_expr)       # colour ∝ mean expr
            circle = plt.Circle(
                (gi, n_cl - 1 - ci), radius,
                color=color, linewidth=0.4,
                edgecolor="#aaa" if pct > 0 else "none",
            )
            ax.add_patch(circle)

    ax.set_xlim(-0.7, n_g - 0.3)
    ax.set_ylim(-0.7, n_cl - 0.3)
    ax.set_xticks(range(n_g))
    ax.set_xticklabels(genes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(n_cl))
    ax.set_yticklabels([str(c) for c in reversed(clusters)], fontsize=8)
    ax.set_xlabel("Marker genes", fontsize=10)
    ax.set_ylabel("Cluster", fontsize=10)
    ax.set_title(
        f"Dot Plot — Top {TOP_N} Marker Genes per Cluster\n"
        "Dot size = % cells expressing  ·  Colour = mean expression (expressing cells)",
        fontsize=11, fontweight="bold",
    )
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)

    # Colour bar
    sm = plt.cm.ScalarMappable(
        cmap="YlOrRd",
        norm=plt.Normalize(vmin=0, vmax=max_expr),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Mean expr.\n(expressing cells)", fontsize=8)

    # Size legend
    for pct_val, label in [(0.25, "25%"), (0.5, "50%"), (1.0, "100%")]:
        ax.scatter(
            [], [], s=(MAX_RADIUS * np.sqrt(pct_val) * 72) ** 2,
            c="gray", alpha=0.6, label=label,
        )
    ax.legend(
        title="% expressing", title_fontsize=7,
        loc="lower right", fontsize=7,
        frameon=True, framealpha=0.8,
        bbox_to_anchor=(1.18, 0),
    )

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_confidence_distribution(adata: AnnData) -> Optional[str]:
    if "cell_type_confidence" not in adata.obs.columns:
        return None
    conf = adata.obs["cell_type_confidence"].astype(float)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(conf, bins=20, color="#4C78A8", edgecolor="white", alpha=0.85)
    for thresh, color, label in [
        (0.5,  "#e07b3a", "Low confidence (<0.5)"),
        (0.75, "#2ca02c", "High confidence (>=0.75)"),
    ]:
        ax.axvline(thresh, color=color, linewidth=1.5, linestyle="--", label=label)
    median_conf = conf.median()
    ax.axvline(median_conf, color="#d62728", linewidth=1.5, linestyle=":",
               label=f"Median ({median_conf:.2f})")
    low_pct  = 100 * (conf < 0.5).mean()
    high_pct = 100 * (conf >= 0.75).mean()
    ax.set_title(f"Consensus Vote Confidence Distribution\n"
                 f"Low (<0.5): {low_pct:.1f}%  |  High (>=0.75): {high_pct:.1f}%",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Confidence Score (fraction of methods agreeing)", fontsize=11)
    ax.set_ylabel("Number of Cells", fontsize=11)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML renderer — matches deg_report._render_page exactly
# ---------------------------------------------------------------------------

def _render_page(title: str, sections: list[str], timestamp: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{title}</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc;
    }}
    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white; padding: 32px 40px 24px;
    }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}
    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}
    section {{
      background: white; border-radius: 10px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.07);
      padding: 28px 32px; margin-bottom: 24px;
    }}
    section h2 {{
      font-size: 1.15rem; font-weight: 700; color: #0f3460;
      border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px;
    }}
    section h3 {{ font-size: 1rem; font-weight: 600; color: #16213e; margin: 18px 0 10px; }}
    section p {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em;
    }}
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{
      background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
      min-width: 130px; text-align: center; flex: 1 1 130px;
    }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }}
    th {{
      background: #f0f2ff; color: #0f3460; font-weight: 600;
      padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0;
    }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}
    .ok   {{ color: #155724; font-weight: 600; }}
    .warn {{ color: #856404; font-weight: 600; }}
    .fig-grid {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }}
    .fig-wrap {{ flex: 1 1 300px; max-width: 520px; }}
    .fig-wrap.wide {{ flex: 1 1 100%; max-width: 100%; }}
    .fig-wrap h3 {{ font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }}
    .fig-wrap img {{ width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }}
    footer {{ text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage -- Cell-Type Annotation Report</h1>
    <p>Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    &middot; MIT License
  </footer>
</body>
</html>
"""


def _section_summary(adata_annotated: AnnData, prov: dict,
                     dataset_name: str, timestamp: str) -> str:
    methods_run = prov.get("methods_run", [])
    leiden_col  = prov.get("leiden_col", "leiden")
    n_clusters  = adata_annotated.obs[leiden_col].nunique() \
                  if leiden_col in adata_annotated.obs else "?"
    conf_median = "?"
    if "cell_type_confidence" in adata_annotated.obs.columns:
        conf_median = f"{adata_annotated.obs['cell_type_confidence'].median():.3f}"

    # Mean scANVI posterior probability (shown only when scANVI ran)
    scanvi_stat = ""
    if "cell_type_scanvi" in adata_annotated.obs.columns:
        scanvi_types = adata_annotated.obs["cell_type_scanvi"].nunique()
        scanvi_stat = (
            f'<div class="stat-card"><div class="stat-value">{scanvi_types}</div>'
            f'<div class="stat-label">scANVI types</div></div>'
        )

    # SingleR stat card (shown only when SingleR ran)
    singler_stat = ""
    if "cell_type_singler" in adata_annotated.obs.columns:
        n_singler = adata_annotated.obs["cell_type_singler"].nunique()
        n_assigned = (adata_annotated.obs["cell_type_singler"] != "Unassigned").sum()
        pct = 100 * n_assigned / max(adata_annotated.n_obs, 1)
        singler_stat = (
            f'<div class="stat-card"><div class="stat-value">{n_singler}</div>'
            f'<div class="stat-label">SingleR types ({pct:.0f}% assigned)</div></div>'
        )

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",           f"{adata_annotated.n_obs:,}"),
            ("Clusters",        str(n_clusters)),
            ("Methods run",     str(len(methods_run))),
            ("Consensus types", str(adata_annotated.obs["cell_type_vote"].nunique()
                                   if "cell_type_vote" in adata_annotated.obs else "?")),
            ("Med. confidence", conf_median),
        ]
    ) + scanvi_stat + singler_stat

    sctype_tissue = prov.get("sctype_tissue", "—")
    scanvi_model  = prov.get("scanvi_model_path", "—")
    singler_ref   = prov.get("singler_ref", "—")
    param_rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Methods run",           ", ".join(methods_run) if methods_run else "?"),
            ("CellTypist models",     ", ".join(prov.get("celltypist_models", []))),
            ("Marker sets used",      str(len(prov.get("marker_sets_keys", [])))),
            ("ScType tissue",         sctype_tissue),
            ("SingleR reference",     singler_ref),
            ("scANVI model",          scanvi_model),
            ("CellTypist models dir", prov.get("celltypist_models_dir", "?")),
        ]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{param_rows}</tbody>
      </table>
    </section>
    """


def _section_assignment_table(adata_annotated: AnnData, prov: dict) -> str:
    leiden_col = prov.get("leiden_col", "leiden")
    if leiden_col not in adata_annotated.obs.columns:
        return ""
    clusters = sorted(adata_annotated.obs[leiden_col].unique(),
                      key=lambda x: (int(x) if str(x).isdigit() else x))

    def _mode(col, cl):
        mask = adata_annotated.obs[leiden_col] == cl
        vc   = adata_annotated.obs.loc[mask, col].value_counts()
        return vc.index[0] if len(vc) > 0 else "?"

    def _conf_badge(cl):
        """Return (score_html, label_str) for the two split confidence columns."""
        if "cell_type_confidence" not in adata_annotated.obs.columns:
            return "?", "?"
        mask = adata_annotated.obs[leiden_col] == cl
        med  = adata_annotated.obs.loc[mask, "cell_type_confidence"].astype(float).median()
        css  = "ok" if med >= 0.75 else "warn"
        if med >= 0.75:
            label = "High"
        elif med >= 0.50:
            label = "Medium"
        else:
            label = "Low"
        return f'<span class="{css}">{med:.2f}</span>', label

    # ── Discover method columns that actually exist in obs ────────────────────
    # Each entry: (obs_col, header_label)
    # Order: CellTypist models first (sorted), then other methods, then consensus last.

    obs_cols = set(adata_annotated.obs.columns)

    # All celltypist_* columns that are present and not "not_run"
    def _col_active(col):
        return col in obs_cols and adata_annotated.obs[col].iloc[0] != "not_run"

    ct_col_headers = []
    for col in sorted(c for c in obs_cols if c.startswith("celltypist_")):
        if not _col_active(col):
            continue
        stem  = col.replace("celltypist_", "").replace("_", " ").title()
        label = f"CellTypist {stem}"
        ct_col_headers.append((col, label))

    # Other fixed method columns — only include if present and active
    other_col_headers = []
    for col, label in [
        ("cell_type_markers", "Marker score"),
        ("cell_type_sctype",  "ScType"),
        ("cell_type_singler", "SingleR"),
        ("cell_type_scanvi",  "scANVI"),
    ]:
        if _col_active(col):
            other_col_headers.append((col, label))

    # Final column order: n cells | CellTypist* | other methods | Consensus | Confidence
    method_cols = ct_col_headers + other_col_headers

    # ── Build header row ──────────────────────────────────────────────────────
    header_cells = "".join(f"<th>{label}</th>" for _, label in method_cols)

    # ── Build data rows ───────────────────────────────────────────────────────
    rows = ""
    for cl in clusters:
        n = int((adata_annotated.obs[leiden_col] == cl).sum())
        method_cells = "".join(f"<td>{_mode(col, cl)}</td>" for col, _ in method_cols)
        consensus = _mode("cell_type_vote", cl) if "cell_type_vote" in obs_cols else "?"
        score_html, label_str = _conf_badge(cl)
        rows += (
            f"<tr><td>{cl}</td><td>{n:,}</td>"
            f"{method_cells}"
            f"<td><strong>{consensus}</strong></td>"
            f"<td>{score_html}</td>"
            f"<td>{label_str}</td></tr>"
        )

    return f"""
    <section>
      <h2>Cluster to Cell Type Assignment</h2>
      <p>Per-cluster label from each method.
         <strong>Consensus score</strong> = weighted fraction of active methods agreeing (0.0–1.0).
         <strong>Confidence</strong> = qualitative label (High ≥ 0.75, Medium ≥ 0.50, Low &lt; 0.50).
         <em>These reflect method agreement only — not AI self-reported confidence.</em></p>
      <table>
        <thead>
          <tr><th>Cluster</th><th>n cells</th>
              {header_cells}
              <th>Consensus</th><th>Consensus score</th><th>Confidence</th></tr>
        </thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_figures(figs: dict, methods_run: list) -> str:
    panels = []

    # All CellTypist panels — keyed by obs column name, already filtered to
    # columns that are present and active (populated in run_annotate_report)
    for key, title in figs.get("_ct_panels", []):
        if figs.get(key):
            panels.append(
                f'<div class="fig-wrap"><h3>{title}</h3>'
                f'<img src="data:image/png;base64,{figs[key]}" alt="{title}"></div>'
            )

    # Other method panels — only rendered when the column exists
    for key, title in [
        ("markers", "UMAP — Marker Gene Score Labels"),
        ("sctype",  "UMAP — ScType Labels"),
        ("singler", "UMAP — SingleR Labels"),
        ("scanvi",  "UMAP — scANVI Transfer Labels"),
    ]:
        if figs.get(key):
            panels.append(
                f'<div class="fig-wrap"><h3>{title}</h3>'
                f'<img src="data:image/png;base64,{figs[key]}" alt="{title}"></div>'
            )

    # Consensus vote — full width
    if figs.get("vote") and "vote" in methods_run:
        panels.append(
            f'<div class="fig-wrap wide"><h3>UMAP — Consensus Cell Type (Majority Vote)</h3>'
            f'<img src="data:image/png;base64,{figs["vote"]}" alt="UMAP consensus"></div>'
        )

    # Ground-truth side-by-side — full width
    if figs.get("groundtruth"):
        panels.append(
            f'<div class="fig-wrap wide"><h3>Consensus vs Ground-Truth Cell Type</h3>'
            f'<img src="data:image/png;base64,{figs["groundtruth"]}" alt="UMAP vs gt"></div>'
        )

    # Marker dot plot — full width
    if figs.get("dotplot"):
        panels.append(
            f'<div class="fig-wrap wide"><h3>Dot Plot — Top Marker Genes per Cluster'
            f'<span style="font-weight:400;font-size:0.82rem;color:#666;margin-left:8px">'
            f'dot size = % expressing &nbsp;·&nbsp; colour = mean expression</span></h3>'
            f'<img src="data:image/png;base64,{figs["dotplot"]}" alt="Marker dot plot"></div>'
        )

    # Marker heatmap — full width
    if figs.get("heatmap"):
        panels.append(
            f'<div class="fig-wrap wide"><h3>Marker Gene Score Heatmap (Clusters × Cell Types)</h3>'
            f'<img src="data:image/png;base64,{figs["heatmap"]}" alt="Marker heatmap"></div>'
        )

    # Confidence histogram — half width
    if figs.get("confidence") and "vote" in methods_run:
        panels.append(
            f'<div class="fig-wrap"><h3>Consensus Confidence Distribution</h3>'
            f'<img src="data:image/png;base64,{figs["confidence"]}" alt="Confidence dist"></div>'
        )

    return f"""
    <section>
      <h2>Figures</h2>
      <div class="fig-grid">{"".join(panels)}</div>
    </section>
    """


def _section_groundtruth_accuracy(adata_annotated: AnnData, prov: dict) -> str:
    """
    Per-cluster comparison table of consensus vote label vs ground truth.
    Only rendered when obs['cell_type_groundtruth'] exists.

    For each cluster reports:
      - Consensus label
      - Most common ground-truth label in that cluster
      - Total cells in cluster
    """
    leiden_col = prov.get("leiden_col", "leiden")

    if (
        "cell_type_groundtruth" not in adata_annotated.obs.columns
        or "cell_type_vote" not in adata_annotated.obs.columns
        or leiden_col not in adata_annotated.obs.columns
    ):
        return ""

    clusters = sorted(
        adata_annotated.obs[leiden_col].unique(),
        key=lambda x: (int(x) if str(x).isdigit() else x),
    )

    obs = adata_annotated.obs[[leiden_col, "cell_type_vote", "cell_type_groundtruth"]]

    rows_html = ""
    for cl in clusters:
        mask      = obs[leiden_col] == cl
        sub       = obs[mask]
        n         = len(sub)
        consensus = sub["cell_type_vote"].mode().iloc[0] if n > 0 else "?"
        gt_mode   = sub["cell_type_groundtruth"].astype(str).mode().iloc[0] if n > 0 else "?"

        rows_html += (
            f"<tr><td>{cl}</td><td>{n:,}</td>"
            f"<td><strong>{consensus}</strong></td>"
            f"<td>{gt_mode}</td></tr>"
        )

    return f"""
    <section>
      <h2>Consensus vs Ground Truth — Per Cluster</h2>
      <p>
        Compares the consensus vote label against
        <code>cell_type_groundtruth</code> (publication ground truth preserved from
        <code>obs['cell_type']</code> before annotation ran).
        Ground-truth and consensus labels may use different naming conventions —
        use this table alongside the Consensus vs Ground-Truth UMAP for visual validation.
      </p>
      <table>
        <thead>
          <tr>
            <th>Cluster</th>
            <th>n cells</th>
            <th>Consensus label</th>
            <th>Most common ground-truth label</th>
          </tr>
        </thead>
        <tbody>{rows_html}</tbody>
      </table>
    </section>
    """


def _section_provenance(adata_annotated: AnnData) -> str:
    prov = adata_annotated.uns.get("omicsage_annotate", {})
    prov_display = {}
    for k, v in prov.items():
        if isinstance(v, list):
            prov_display[k] = ", ".join(str(x) for x in v)
        elif isinstance(v, dict):
            for subk, subv in v.items():
                prov_display[f"{k}[{subk}]"] = subv
        else:
            prov_display[k] = v
    rows = "".join(
        f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"
        for k, v in prov_display.items()
    )
    return f"""
    <section>
      <h2>Provenance</h2>
      <table>
        <thead><tr><th>Key</th><th>Value</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


# ---------------------------------------------------------------------------
# Public API
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
        Annotated AnnData returned by annotate().
    annotation_dict : dict
        Second return value of annotate().
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
    timestamp   = datetime.now().strftime("%Y-%m-%d %H:%M")
    prov        = adata_annotated.uns.get("omicsage_annotate", {})
    methods_run = prov.get("methods_run", [])

    print(f"Building annotation report for '{dataset_name}' ...", flush=True)

    # ── CellTypist UMAPs — only for columns that exist and are not "not_run" ──
    print("  Rendering CellTypist UMAP panels ...", flush=True)
    figs: dict = {}
    ct_panels = []   # list of (fig_key, title) for _section_figures

    for col in sorted(c for c in adata_annotated.obs.columns
                      if c.startswith("celltypist_")):
        # Skip placeholder columns written when a model wasn't run
        if adata_annotated.obs[col].iloc[0] == "not_run":
            continue
        stem  = col.replace("celltypist_", "").replace("_", " ").title()
        title = f"UMAP — CellTypist {stem}"
        b64   = _plot_umap_labels(adata_annotated, col, title,
                                  f"Per-cluster label from {stem} model")
        if b64:
            figs[col] = b64
            ct_panels.append((col, title))

    figs["_ct_panels"] = ct_panels

    # ── Other method UMAPs ────────────────────────────────────────────────────
    figs["markers"] = _plot_umap_labels(adata_annotated, "cell_type_markers",
                                        "UMAP — Marker Gene Score Labels",
                                        "Best-scoring cell type per cluster")
    figs["sctype"]  = _plot_umap_labels(adata_annotated, "cell_type_sctype",
                                        "UMAP — ScType Labels",
                                        "Best ScType label per cluster")
    figs["singler"] = _plot_umap_labels(adata_annotated, "cell_type_singler",
                                        "UMAP — SingleR Labels",
                                        "Per-cell label from SingleR")
    figs["scanvi"]  = _plot_umap_labels(adata_annotated, "cell_type_scanvi",
                                        "UMAP — scANVI Transfer Labels",
                                        "Per-cell label from scANVI model")

    if "vote" in methods_run:
        print("  Rendering consensus vote UMAP ...", flush=True)
        figs["vote"] = _plot_umap_labels(
            adata_annotated, "cell_type_vote",
            "UMAP -- Consensus Cell Type (Majority Vote)",
        )

    print("  Rendering ground-truth UMAP (if available) ...", flush=True)
    if "cell_type_groundtruth" in adata_annotated.obs.columns:
        figs["groundtruth"] = _plot_umap_side_by_side(
            adata_annotated,
            col_left="cell_type_vote", col_right="cell_type_groundtruth",
            title_left="Consensus Vote", title_right="Ground-Truth (publication)",
            suptitle="Consensus Vote vs Ground-Truth Cell Type",
        )

    score_df = annotation_dict.get("marker_score_df")
    if score_df is not None:
        print("  Rendering marker score heatmap ...", flush=True)
        figs["heatmap"] = _plot_marker_heatmap(score_df)

    leiden_col_for_dot = prov.get("leiden_col", "leiden")
    print("  Rendering marker dot plot ...", flush=True)
    figs["dotplot"] = _plot_dotplot(adata_annotated, leiden_col_for_dot)

    if "vote" in methods_run:
        print("  Rendering confidence distribution ...", flush=True)
        figs["confidence"] = _plot_confidence_distribution(adata_annotated)

    sections = [
        _section_summary(adata_annotated, prov, dataset_name, timestamp),
        _section_assignment_table(adata_annotated, prov),
        _section_groundtruth_accuracy(adata_annotated, prov),
        _section_figures(figs, methods_run),
        _section_provenance(adata_annotated),
    ]

    html = _render_page(
        title=f"OmicSage -- Annotate Report -- {dataset_name}",
        sections=sections,
        timestamp=timestamp,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
