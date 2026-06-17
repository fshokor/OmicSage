"""
deg_report.py — HTML report generator for OmicSage DEG module
Phase 1, Step 6

Generates a self-contained HTML report from deg() output.
Matches the style of annotate_report.py.

Usage:
    from reports.deg_report import generate_deg_report

    generate_deg_report(
        adata=adata_deg,
        deg_dict=deg_dict,
        output_path="reports/deg_report.html",
    )
"""

from __future__ import annotations

import base64
import io
import warnings
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_deg_report(
    adata: AnnData,
    deg_dict: dict,
    output_path: str = "reports/deg_report.html",
    top_n_volcano: int = 10,
    top_n_dotplot: int = 5,
    max_volcano_groups: int = 20,
) -> str:
    """
    Generate a self-contained HTML report for DEG results.

    Parameters
    ----------
    adata : AnnData
        AnnData returned by deg(), containing adata.uns['omicsage_deg']
        and adata.uns['rank_genes_groups'].
    deg_dict : dict
        deg_dict returned by deg().
    output_path : str
        Path to write the HTML file.
    top_n_volcano : int
        Number of top genes to label on each volcano plot.
    top_n_dotplot : int
        Top N DEGs per group to include in the dot plot.
    max_volcano_groups : int
        Maximum number of groups to render volcano plots for.
        Default raised to 20 — all groups rendered for typical datasets.
        When the limit is exceeded, a visible note is added to the report
        and the groups with the most DEGs are shown first.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    provenance  = deg_dict["provenance"]
    results     = deg_dict["results"]
    summary_df  = deg_dict["summary_df"]

    # ------------------------------------------------------------------
    # Build HTML sections
    # ------------------------------------------------------------------
    sections = []

    sections.append(_section_summary_stats(provenance, results))
    sections.append(_section_summary_table(summary_df))
    sections.append(_section_global_volcano(results, provenance, top_n=top_n_volcano))
    sections.append(_section_volcano_plots(
        adata, results, provenance,
        top_n=top_n_volcano,
        max_groups=max_volcano_groups,
    ))
    sections.append(_section_dotplot(adata, results, top_n=top_n_dotplot))

    html = _render_page(
        title="OmicSage — Differential Expression Report",
        sections=sections,
        provenance=provenance,
    )

    output_path.write_text(html, encoding="utf-8")
    return str(output_path.resolve())


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _section_summary_stats(provenance: dict, results: dict) -> str:
    n_groups            = provenance.get("n_groups", len(results))
    method              = provenance.get("method", "wilcoxon")
    groupby             = provenance.get("groupby", "—")
    min_logfc           = provenance.get("min_logfc", "—")
    max_pval_adj        = provenance.get("max_pval_adj", "—")
    n_genes             = provenance.get("n_genes", "—")
    exclude_prefixes    = provenance.get("exclude_gene_prefixes", [])
    timestamp           = provenance.get("timestamp", "—")

    total_sig = sum(len(df) for df in results.values())
    per_group_counts = {g: len(df) for g, df in results.items()}

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Groups tested",           n_groups),
            ("Method",                  method.capitalize()),
            ("Grouped by",              groupby),
            ("Min |log₂FC|",            min_logfc),
            ("Max adj. p-value",        max_pval_adj),
            ("Genes computed / group",  n_genes),
            ("Total significant DEGs",  total_sig),
        ]
    )

    per_group_rows = "".join(
        f"<tr><td>{group}</td><td>{count}</td></tr>"
        for group, count in sorted(per_group_counts.items())
    )

    # Show excluded prefixes note if any were applied
    exclude_note = ""
    if exclude_prefixes:
        prefix_str = ", ".join(exclude_prefixes)
        exclude_note = (
            f'<p class="timestamp">ℹ Gene prefix exclusion applied: '
            f'<code>{prefix_str}</code> — these genes are excluded from '
            f'results but were still used in fold-change computation.</p>'
        )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Analysis run: {timestamp}</p>
      <div class="stat-grid">{stat_cards}</div>
      {exclude_note}
      <h3>Significant DEGs per Group</h3>
      <table>
        <thead><tr><th>Group</th><th>Significant DEGs</th></tr></thead>
        <tbody>{per_group_rows}</tbody>
      </table>
    </section>
    """


def _section_summary_table(summary_df: pd.DataFrame) -> str:
    if summary_df.empty:
        return "<section><h2>Top DEGs per Group</h2><p>No significant DEGs found.</p></section>"

    rows = ""
    current_group = None
    for _, row in summary_df.iterrows():
        group_cell = ""
        if row["group"] != current_group:
            count = (summary_df["group"] == row["group"]).sum()
            group_cell = f'<td rowspan="{count}" class="group-cell">{row["group"]}</td>'
            current_group = row["group"]

        logfc_class = "pos-fc" if float(row["logfc"]) > 0 else "neg-fc"
        direction   = "▲ Up" if float(row["logfc"]) > 0 else "▼ Down"
        pval_fmt    = f'{float(row["pval_adj"]):.2e}'
        logfc_fmt   = f'{float(row["logfc"]):.3f}'

        rows += (
            f"<tr>{group_cell}"
            f"<td>{int(row['rank'])}</td>"
            f"<td><strong>{row['gene']}</strong></td>"
            f'<td class="{logfc_class}">{direction}</td>'
            f'<td class="{logfc_class}">{logfc_fmt}</td>'
            f"<td>{pval_fmt}</td></tr>"
        )

    return f"""
    <section>
      <h2>Top DEGs per Group</h2>
      <p>Top 5 significant genes per group, ranked by adjusted p-value.</p>
      <table>
        <thead>
          <tr>
            <th>Group</th><th>Rank</th><th>Gene</th>
            <th>Direction</th><th>log₂FC</th><th>Adj. p-value</th>
          </tr>
        </thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_global_volcano(
    results: dict,
    provenance: dict,
    top_n: int = 10,
) -> str:
    """Single volcano plot pooling DEGs from all groups combined."""
    frames = [df for df in results.values() if not df.empty]
    if not frames:
        return (
            "<section><h2>Global Volcano Plot — All Groups Combined</h2>"
            "<p>No significant DEGs found across any group.</p></section>"
        )

    combined = pd.concat(frames, ignore_index=True)
    min_logfc    = float(provenance.get("min_logfc", 0.25))
    max_pval_adj = float(provenance.get("max_pval_adj", 0.05))

    img_b64 = _render_volcano(
        combined, "All Groups Combined",
        min_logfc=min_logfc,
        max_pval_adj=max_pval_adj,
        top_n=top_n,
        figsize=(7, 5),
    )
    n_groups = len(results)
    n_genes  = len(combined)
    return f"""
    <section>
      <h2>Global Volcano Plot — All Groups Combined</h2>
      <p>
        All {n_genes:,} tested gene-entries across {n_groups} group(s) shown in a single plot.
        Dashed lines mark thresholds: |log₂FC| ≥ {min_logfc}, adj. p ≤ {max_pval_adj}.
        Top {top_n} genes by significance are labelled.
      </p>
      <div style="max-width:640px;">
        <img src="data:image/png;base64,{img_b64}" alt="Global volcano plot"
             style="width:100%; border-radius:6px; border:1px solid #e8eaf6;">
      </div>
    </section>
    """


def _section_volcano_plots(
    adata: AnnData,
    results: dict,
    provenance: dict,
    top_n: int = 10,
    max_groups: int = 20,
) -> str:
    all_groups = list(results.keys())
    note_html  = ""

    if len(all_groups) > max_groups:
        # Show groups with the most DEGs first so the most informative
        # groups are always included when truncation occurs
        sorted_groups = sorted(
            all_groups,
            key=lambda g: len(results[g]),
            reverse=True,
        )
        excluded = sorted(set(all_groups) - set(sorted_groups[:max_groups]))
        groups = sorted_groups[:max_groups]
        note_html = (
            f'<p class="note">⚠ Showing {max_groups} of {len(all_groups)} groups '
            f'(most DEGs first). Excluded: {", ".join(excluded)}</p>'
        )
    else:
        groups = all_groups

    min_logfc    = float(provenance.get("min_logfc", 0.25))
    max_pval_adj = float(provenance.get("max_pval_adj", 0.05))

    plots_html = ""
    for group in groups:
        df = results[group]
        if df.empty:
            plots_html += (
                f"<div class='volcano-wrap'><h3>{group}</h3>"
                "<p>No significant DEGs.</p></div>"
            )
            continue

        img_b64 = _render_volcano(
            df, group,
            min_logfc=min_logfc,
            max_pval_adj=max_pval_adj,
            top_n=top_n,
        )
        plots_html += (
            f'<div class="volcano-wrap">'
            f'<h3>{group}</h3>'
            f'<img src="data:image/png;base64,{img_b64}" alt="Volcano plot {group}">'
            f'</div>'
        )

    return f"""
    <section>
      <h2>Volcano Plots</h2>
      <p>
        Each plot shows log₂ fold-change vs. −log₁₀(adjusted p-value).
        Dashed lines mark thresholds: |log₂FC| ≥ {min_logfc}, adj. p ≤ {max_pval_adj}.
        Top {top_n} genes by significance are labelled.
      </p>
      {note_html}
      <div class="volcano-grid">{plots_html}</div>
    </section>
    """


def _section_dotplot(
    adata: AnnData,
    results: dict,
    top_n: int = 5,
) -> str:
    """
    Dot plot of top N DEGs per group using sc.pl.dotplot.
    Falls back gracefully if rank_genes_groups key is missing.
    """
    if "rank_genes_groups" not in adata.uns:
        return (
            "<section><h2>Dot Plot</h2>"
            "<p>rank_genes_groups not found in adata.uns — dot plot skipped.</p>"
            "</section>"
        )

    gene_lists: dict[str, list[str]] = {}
    for group, df in results.items():
        if not df.empty:
            gene_lists[group] = df["gene"].head(top_n).tolist()

    if not gene_lists:
        return (
            "<section><h2>Dot Plot</h2>"
            "<p>No significant DEGs to display.</p>"
            "</section>"
        )

    # Flatten to unique ordered gene list
    seen = set()
    all_genes = []
    for genes in gene_lists.values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                all_genes.append(g)

    valid_genes = [g for g in all_genes if g in adata.var_names]
    if not valid_genes:
        return (
            "<section><h2>Dot Plot</h2>"
            "<p>Top DEGs not found in adata.var_names — dot plot skipped.</p>"
            "</section>"
        )

    groupby = adata.uns.get("omicsage_deg", {}).get("groupby", None)
    if groupby is None or groupby not in adata.obs.columns:
        return (
            "<section><h2>Dot Plot</h2>"
            "<p>Group column not found in adata.obs — dot plot skipped.</p>"
            "</section>"
        )

    try:
        fig, ax = plt.subplots(figsize=(max(6, len(valid_genes) * 0.5 + 2), 4))
        sc.pl.dotplot(
            adata,
            var_names=valid_genes,
            groupby=groupby,
            ax=ax,
            show=False,
            title=f"Top {top_n} DEGs per group",
        )
        img_b64 = _fig_to_base64(fig)
        plt.close(fig)
        img_html = (
            f'<img src="data:image/png;base64,{img_b64}" '
            'alt="Dot plot top DEGs" style="max-width:100%;">'
        )
    except Exception as e:
        img_html = f"<p>Dot plot could not be rendered: {e}</p>"

    return f"""
    <section>
      <h2>Dot Plot — Top {top_n} DEGs per Group</h2>
      <p>
        Dot size encodes fraction of cells expressing the gene.
        Colour encodes mean expression in expressing cells.
      </p>
      {img_html}
    </section>
    """


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _render_volcano(
    df: pd.DataFrame,
    group: str,
    min_logfc: float,
    max_pval_adj: float,
    top_n: int = 10,
    figsize: tuple = (5, 4),
) -> str:
    """Render a single volcano plot; return base64-encoded PNG string."""
    plot_df = df.copy()
    plot_df["logfc"]    = pd.to_numeric(plot_df["logfc"],    errors="coerce").fillna(0.0)
    plot_df["pval_adj"] = pd.to_numeric(plot_df["pval_adj"], errors="coerce").fillna(1.0)

    plot_df["pval_adj"] = plot_df["pval_adj"].clip(lower=1e-300)
    plot_df["neg_log_p"] = -np.log10(plot_df["pval_adj"])

    sig_up   = (plot_df["pval_adj"] <= max_pval_adj) & (plot_df["logfc"] >=  min_logfc)
    sig_down = (plot_df["pval_adj"] <= max_pval_adj) & (plot_df["logfc"] <= -min_logfc)
    ns       = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(plot_df.loc[ns,       "logfc"], plot_df.loc[ns,       "neg_log_p"],
               s=8, color="#aaaaaa", alpha=0.5, label="NS", rasterized=True)
    ax.scatter(plot_df.loc[sig_up,   "logfc"], plot_df.loc[sig_up,   "neg_log_p"],
               s=10, color="#e05252", alpha=0.7, label="Up", rasterized=True)
    ax.scatter(plot_df.loc[sig_down, "logfc"], plot_df.loc[sig_down, "neg_log_p"],
               s=10, color="#5282e0", alpha=0.7, label="Down", rasterized=True)

    ax.axhline(-np.log10(max_pval_adj), color="#888888", linestyle="--", linewidth=0.8)
    ax.axvline( min_logfc,              color="#888888", linestyle="--", linewidth=0.8)
    ax.axvline(-min_logfc,              color="#888888", linestyle="--", linewidth=0.8)

    top_genes = plot_df.nsmallest(top_n, "pval_adj")
    for _, row in top_genes.iterrows():
        ax.annotate(
            row["gene"],
            xy=(row["logfc"], row["neg_log_p"]),
            xytext=(4, 2),
            textcoords="offset points",
            fontsize=5.5,
            color="#222222",
        )

    ax.set_xlabel("log₂ Fold-Change", fontsize=9)
    ax.set_ylabel("−log₁₀(adj. p-value)", fontsize=9)
    ax.set_title(group, fontsize=10, fontweight="bold")
    ax.legend(fontsize=7, markerscale=1.5, frameon=False)
    fig.tight_layout()

    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def _fig_to_base64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


# ---------------------------------------------------------------------------
# HTML renderer
# ---------------------------------------------------------------------------

def _render_page(title: str, sections: list[str], provenance: dict) -> str:
    body = "\n".join(sections)
    scanpy_ver = provenance.get("scanpy_version", "—")
    timestamp  = provenance.get("timestamp", "—")

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
      font-size: 14px;
      line-height: 1.6;
      color: #1a1a2e;
      background: #f7f8fc;
    }}

    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white;
      padding: 32px 40px 24px;
    }}
    header h1 {{ font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }}
    header p  {{ font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }}

    main {{ max-width: 1100px; margin: 0 auto; padding: 32px 24px; }}

    section {{
      background: white;
      border-radius: 10px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.07);
      padding: 28px 32px;
      margin-bottom: 24px;
    }}
    section h2 {{
      font-size: 1.15rem;
      font-weight: 700;
      color: #0f3460;
      border-bottom: 2px solid #e8eaf6;
      padding-bottom: 10px;
      margin-bottom: 18px;
    }}
    section h3 {{
      font-size: 1rem;
      font-weight: 600;
      color: #16213e;
      margin: 18px 0 10px;
    }}
    section p {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}

    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}

    .note {{
      font-size: 0.82rem;
      color: #7a5c00;
      background: #fffbe6;
      border-left: 3px solid #f0c040;
      padding: 8px 12px;
      border-radius: 4px;
      margin-bottom: 14px;
    }}

    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 0.85em;
    }}

    /* Stat cards */
    .stat-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      margin-bottom: 24px;
    }}
    .stat-card {{
      background: #f0f2ff;
      border-radius: 8px;
      padding: 14px 20px;
      min-width: 130px;
      text-align: center;
      flex: 1 1 130px;
    }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}

    /* Tables */
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.88rem;
      margin-top: 8px;
    }}
    th {{
      background: #f0f2ff;
      color: #0f3460;
      font-weight: 600;
      padding: 9px 12px;
      text-align: left;
      border-bottom: 2px solid #d0d4f0;
    }}
    td {{
      padding: 8px 12px;
      border-bottom: 1px solid #eee;
      vertical-align: middle;
    }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}
    .group-cell {{
      font-weight: 600;
      color: #0f3460;
      background: #f8f9ff;
      border-right: 3px solid #d0d4f0;
      vertical-align: top;
      padding-top: 10px;
    }}
    .pos-fc {{ color: #c0392b; font-weight: 600; }}
    .neg-fc {{ color: #2980b9; font-weight: 600; }}    /* Volcano grid */
    .volcano-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 18px;
      margin-top: 12px;
    }}
    .volcano-wrap {{
      flex: 1 1 300px;
      max-width: 420px;
    }}
    .volcano-wrap h3 {{
      font-size: 0.9rem;
      margin-bottom: 6px;
      color: #16213e;
    }}
    .volcano-wrap img {{
      width: 100%;
      border-radius: 6px;
      border: 1px solid #e8eaf6;
    }}

    footer {{
      text-align: center;
      font-size: 0.78rem;
      color: #aaa;
      padding: 24px 0 32px;
    }}
    footer a {{ color: #0f3460; text-decoration: none; }}
  </style>
</head>
<body>
  <header>
    <h1>OmicSage — Differential Expression Report</h1>
    <p>Generated {timestamp} · Scanpy {scanpy_ver}</p>
  </header>
  <main>
    {body}
  </main>
  <footer>
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    · MIT License
  </footer>
</body>
</html>
"""
