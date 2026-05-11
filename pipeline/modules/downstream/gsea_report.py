"""
gsea_report.py — HTML report generator for OmicSage GSEA module
Phase 1, Step 7

Generates a self-contained HTML report from gsea() output.
Matches the style of deg_report.py.

Usage:
    from reports.gsea_report import generate_gsea_report

    generate_gsea_report(
        gsea_dict=gsea_dict,
        output_path="reports/gsea_report.html",
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


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_gsea_report(
    gsea_dict: dict,
    output_path: str = "reports/gsea_report.html",
    top_n_table: int = 5,
    top_n_bar: int = 10,
    max_bubble_groups: int = 9,
) -> str:
    """
    Generate a self-contained HTML report for GSEA results.

    Parameters
    ----------
    gsea_dict : dict
        gsea_dict returned by gsea(). Must contain 'results', 'summary_df',
        'provenance', and 'skipped' keys.
    output_path : str
        Path to write the HTML file.
    top_n_table : int
        Number of top pathways per group to show in the summary table.
    top_n_bar : int
        Number of top pathways per group to show in the bar chart.
    max_bubble_groups : int
        Maximum number of groups to include in the bubble plot.
        If more groups exist the bubble plot is skipped.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    provenance = gsea_dict["provenance"]
    results    = gsea_dict["results"]
    summary_df = gsea_dict["summary_df"]
    skipped    = gsea_dict.get("skipped", [])

    sections = []
    sections.append(_section_summary_stats(provenance, results, skipped))
    sections.append(_section_top_pathways_table(results, top_n=top_n_table))
    sections.append(_section_bar_charts(results, top_n=top_n_bar))

    if len(results) >= 2:
        sections.append(_section_bubble_plot(results, max_groups=max_bubble_groups))

    html = _render_page(
        title="OmicSage — GSEA Pathway Enrichment Report",
        sections=sections,
        provenance=provenance,
    )

    output_path.write_text(html, encoding="utf-8")
    return str(output_path.resolve())


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _section_summary_stats(
    provenance: dict,
    results: dict,
    skipped: list,
) -> str:
    gene_sets      = provenance.get("gene_sets", [])
    organism       = provenance.get("organism", "—")
    direction      = provenance.get("direction", "up")
    n_tested       = provenance.get("n_groups_tested", len(results))
    n_skipped      = provenance.get("n_groups_skipped", len(skipped))
    min_logfc      = provenance.get("min_logfc", "—")
    max_pval_adj   = provenance.get("max_pval_adj", "—")
    timestamp      = provenance.get("timestamp", "—")
    gseapy_version = provenance.get("gseapy_version", "—")

    direction_label = {
        "up":   "▲ Upregulated only",
        "down": "▼ Downregulated only",
        "both": "▲▼ Both directions (separate queries)",
    }.get(direction, direction)

    # Count total significant pathways (adj. p ≤ 0.05)
    total_sig = 0
    for df in results.values():
        if df.empty or "Adjusted P-value" not in df.columns:
            continue
        total_sig += int((pd.to_numeric(df["Adjusted P-value"], errors="coerce") <= 0.05).sum())

    gene_sets_str = ", ".join(gene_sets) if gene_sets else "—"

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Groups tested",              n_tested),
            ("Groups skipped",             n_skipped),
            ("Gene sets queried",           len(gene_sets)),
            ("Organism",                   organism.capitalize()),
            ("Min |log₂FC| (query)",       min_logfc),
            ("Sig. pathways (adj.p≤0.05)", total_sig),
        ]
    )

    # Build per-group summary rows — for "both", show up and down on
    # separate rows under the same group cell
    per_group_rows = ""
    if direction == "both":
        # Group the result keys by base group name
        base_groups: dict[str, list[str]] = {}
        for key in sorted(results.keys()):
            base = key.rsplit("__", 1)[0] if "__" in key else key
            base_groups.setdefault(base, []).append(key)

        for base, keys in sorted(base_groups.items()):
            rowspan = len(keys)
            first = True
            for key in sorted(keys):   # __down before __up alphabetically
                df = results[key]
                dir_label = key.rsplit("__", 1)[1] if "__" in key else "—"
                dir_badge = (
                    '<span class="dir-up">▲ Up</span>'
                    if dir_label == "up"
                    else '<span class="dir-down">▼ Down</span>'
                )
                group_cell = (
                    f'<td rowspan="{rowspan}" class="group-cell">{base}</td>'
                    if first else ""
                )
                first = False
                per_group_rows += (
                    f"<tr>{group_cell}"
                    f"<td>{dir_badge}</td>"
                    f"<td>{_count_sig(df)}</td>"
                    f"<td>{_top_term(df)}</td></tr>"
                )
    else:
        for key, df in sorted(results.items()):
            per_group_rows += (
                f"<tr><td>{key}</td>"
                f"<td>{_count_sig(df)}</td>"
                f"<td>{_top_term(df)}</td></tr>"
            )

    # Format skipped — may be list of strings or list of (group, dir) tuples
    skipped_html = ""
    if skipped:
        skipped_items = [
            f"{g} ({d})" if isinstance(s, tuple) and len(s) == 2
            else str(s)
            for s in skipped
            for g, d in [s if isinstance(s, tuple) else (s, "")]
        ]
        skipped_list = ", ".join(skipped_items)
        skipped_html = (
            f'<p class="timestamp">⚠ Skipped (too few DEGs): {skipped_list}</p>'
        )

    dir_header = "<th>Direction</th>" if direction == "both" else ""

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Analysis run: {timestamp} · gseapy {gseapy_version}</p>
      <p class="timestamp">Gene sets: {gene_sets_str}</p>
      <p class="timestamp">Query direction: {direction_label}</p>
      <div class="stat-grid">{stat_cards}</div>
      {skipped_html}
      <h3>Significant Pathways per Group</h3>
      <table>
        <thead>
          <tr>
            <th>Group</th>
            {dir_header}
            <th>Sig. Pathways (adj.p≤0.05)</th>
            <th>Top Pathway</th>
          </tr>
        </thead>
        <tbody>{per_group_rows}</tbody>
      </table>
    </section>
    """


def _section_top_pathways_table(results: dict, top_n: int = 5) -> str:
    if not results:
        return (
            "<section><h2>Top Pathways per Group</h2>"
            "<p>No enrichment results to display.</p></section>"
        )

    # Detect whether keys are suffixed with __up / __down
    has_direction = any("__up" in k or "__down" in k for k in results)

    rows = ""
    current_group_cell = None

    for key in sorted(results.keys()):
        df = results[key]
        if df.empty:
            continue

        # Parse display group and direction badge
        if has_direction and "__" in key:
            base_group, dir_label = key.rsplit("__", 1)
            dir_badge = (
                '<span class="dir-up">▲ Up</span>'
                if dir_label == "up"
                else '<span class="dir-down">▼ Down</span>'
            )
        else:
            base_group = key
            dir_badge = ""

        top = df.head(top_n)

        for i, (_, row) in enumerate(top.iterrows()):
            # Group cell spans all rows for this key
            group_cell = ""
            if key != current_group_cell:
                count = min(top_n, len(df))
                group_cell = (
                    f'<td rowspan="{count}" class="group-cell">'
                    f'{base_group}'
                    f'{"<br>" + dir_badge if dir_badge else ""}'
                    f'</td>'
                )
                current_group_cell = key

            padj    = row.get("Adjusted P-value", np.nan)
            padj_fmt = f"{float(padj):.2e}" if pd.notna(padj) else "—"
            term     = row.get("Term", "—")
            overlap  = row.get("Overlap", "—")
            if not isinstance(overlap, str) and pd.isna(overlap):
                overlap = "—"
            genes = row.get("Genes", "—")
            if isinstance(genes, str) and len(genes) > 80:
                genes = genes[:80] + "…"

            rows += (
                f"<tr>{group_cell}"
                f"<td>{term}</td>"
                f"<td>{overlap}</td>"
                f"<td>{padj_fmt}</td>"
                f"<td class='gene-list'>{genes}</td></tr>"
            )

    return f"""
    <section>
      <h2>Top Pathways per Group</h2>
      <p>Top {top_n} pathways per group, ranked by adjusted p-value.</p>
      <table>
        <thead>
          <tr>
            <th>Group</th>
            <th>Pathway Term</th>
            <th>Genes Matched</th>
            <th>Adj. p-value</th>
            <th>Genes</th>
          </tr>
        </thead>
        <tbody>{rows}</tbody>
      </table>
    </section>
    """


def _section_bar_charts(results: dict, top_n: int = 10) -> str:
    if not results:
        return (
            "<section><h2>Pathway Bar Charts</h2>"
            "<p>No enrichment results to display.</p></section>"
        )

    plots_html = ""
    for group in sorted(results.keys()):
        df = results[group]
        if df.empty:
            plots_html += (
                f"<div class='bar-wrap'><h3>{group}</h3>"
                "<p>No significant pathways.</p></div>"
            )
            continue

        img_b64 = _render_bar_chart(df, group, top_n=top_n)
        plots_html += (
            f'<div class="bar-wrap">'
            f"<h3>{group}</h3>"
            f'<img src="data:image/png;base64,{img_b64}" '
            f'alt="Bar chart {group}">'
            f"</div>"
        )

    return f"""
    <section>
      <h2>Pathway Bar Charts</h2>
      <p>
        Top {top_n} pathways per group ranked by −log₁₀(adjusted p-value).
        Longer bars = more significant enrichment.
      </p>
      <div class="bar-grid">{plots_html}</div>
    </section>
    """


def _section_bubble_plot(results: dict, max_groups: int = 9) -> str:
    """
    Bubble plot: pathway × group matrix.
    X-axis = groups, Y-axis = top pathways (union of top 5 per group).
    Bubble size = overlap count, colour = −log10(adj. p-value).

    If there are more groups than max_groups, selects the top max_groups
    by number of significant pathways (adj.p ≤ 0.05) so the plot stays
    readable. A note is added to the section header.
    """
    note_html = ""
    plot_results = results

    if len(results) > max_groups:
        # Rank groups by number of significant pathways, take top N
        sig_counts = {
            group: _count_sig(df)
            for group, df in results.items()
        }
        top_groups = sorted(sig_counts, key=sig_counts.get, reverse=True)[:max_groups]
        plot_results = {g: results[g] for g in top_groups}
        excluded = sorted(set(results.keys()) - set(top_groups))
        note_html = (
            f'<p class="timestamp">⚠ Showing top {max_groups} groups by '
            f'significant pathway count. Excluded: {", ".join(excluded)}</p>'
        )

    try:
        img_b64 = _render_bubble_plot(plot_results)
        img_html = (
            f'<img src="data:image/png;base64,{img_b64}" '
            'alt="Bubble plot pathways across groups" '
            'style="max-width:100%;">'
        )
    except Exception as e:
        img_html = f"<p>Bubble plot could not be rendered: {e}</p>"

    return f"""
    <section>
      <h2>Bubble Plot — Pathways Across Groups</h2>
      <p>
        Each bubble represents a pathway in a group.
        Size encodes the number of matched genes; colour encodes −log₁₀(adj. p-value).
        Only the top 5 pathways per group are shown.
      </p>
      {note_html}
      {img_html}
    </section>
    """


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _render_bar_chart(df: pd.DataFrame, group: str, top_n: int = 10) -> str:
    """Render a horizontal bar chart of top pathways; return base64 PNG."""
    plot_df = df.head(top_n).copy()
    padj = pd.to_numeric(plot_df["Adjusted P-value"], errors="coerce").fillna(1.0)
    padj = padj.clip(lower=1e-300)
    plot_df["neg_log_padj"] = -np.log10(padj)

    # Truncate long term names
    plot_df["term_short"] = plot_df["Term"].apply(
        lambda t: (t[:50] + "…") if isinstance(t, str) and len(t) > 50 else t
    )

    n = len(plot_df)
    fig_height = max(3, n * 0.35 + 0.8)
    fig, ax = plt.subplots(figsize=(7, fig_height))

    colors = plt.cm.RdYlBu_r(
        np.linspace(0.3, 0.9, n)
    )
    bars = ax.barh(
        range(n),
        plot_df["neg_log_padj"].values,
        color=colors,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(range(n))
    ax.set_yticklabels(plot_df["term_short"].values, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("−log₁₀(Adjusted p-value)", fontsize=8)
    ax.set_title(group, fontsize=9, fontweight="bold")
    ax.axvline(-np.log10(0.05), color="#888888", linestyle="--",
               linewidth=0.8, label="adj.p=0.05")
    ax.legend(fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def _render_bubble_plot(results: dict) -> str:
    """
    Render a bubble plot of pathways × groups.
    Returns base64 PNG string.
    """
    # Collect top 5 unique pathway terms per group
    top_terms: dict[str, list[str]] = {}
    for group, df in results.items():
        if not df.empty and "Term" in df.columns:
            top_terms[group] = df["Term"].head(5).tolist()

    all_terms = []
    seen = set()
    for terms in top_terms.values():
        for t in terms:
            if t not in seen:
                seen.add(t)
                all_terms.append(t)

    groups = sorted(results.keys())
    n_terms  = len(all_terms)
    n_groups = len(groups)

    if n_terms == 0 or n_groups == 0:
        raise ValueError("No terms to plot")

    # Build matrices: neg_log_padj and overlap_size
    neg_log_mat = np.zeros((n_terms, n_groups))
    size_mat    = np.zeros((n_terms, n_groups))

    for j, group in enumerate(groups):
        df = results[group]
        if df.empty:
            continue
        for i, term in enumerate(all_terms):
            match = df[df["Term"] == term]
            if match.empty:
                continue
            row = match.iloc[0]
            padj = pd.to_numeric(row.get("Adjusted P-value", 1.0), errors="coerce")
            padj = max(float(padj) if pd.notna(padj) else 1.0, 1e-300)
            neg_log_mat[i, j] = -np.log10(padj)
            # Parse overlap — may be "k/n" (old gseapy) or "k" (new gseapy ≥1.0)
            overlap_str = str(row.get("Overlap", "0"))
            try:
                # Handle both "14/500" and "14"
                size_mat[i, j] = int(overlap_str.split("/")[0])
            except (ValueError, IndexError):
                size_mat[i, j] = 0

    fig_w = max(5, n_groups * 1.5 + 2)
    fig_h = max(4, n_terms * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for i in range(n_terms):
        for j in range(n_groups):
            val  = neg_log_mat[i, j]
            size = size_mat[i, j]
            if val == 0 and size == 0:
                continue
            ax.scatter(
                j, i,
                s=max(20, size * 15),
                c=[[plt.cm.RdYlBu_r(min(val / 5.0, 1.0))]],
                alpha=0.85,
                edgecolors="white",
                linewidth=0.5,
            )

    # Truncate long term names
    short_terms = [
        (t[:45] + "…") if isinstance(t, str) and len(t) > 45 else t
        for t in all_terms
    ]

    ax.set_xticks(range(n_groups))
    ax.set_xticklabels(groups, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(range(n_terms))
    ax.set_yticklabels(short_terms, fontsize=7)
    ax.set_xlim(-0.5, n_groups - 0.5)
    ax.set_ylim(-0.5, n_terms - 0.5)
    ax.invert_yaxis()
    ax.grid(True, linewidth=0.3, alpha=0.4)
    ax.set_title("Pathway enrichment across groups", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Colourbar
    sm = plt.cm.ScalarMappable(
        cmap="RdYlBu_r",
        norm=plt.Normalize(vmin=0, vmax=5),
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("−log₁₀(adj. p-value)", fontsize=7)

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
# Stat helpers
# ---------------------------------------------------------------------------

def _count_sig(df: pd.DataFrame, threshold: float = 0.05) -> int:
    if df.empty or "Adjusted P-value" not in df.columns:
        return 0
    return int((pd.to_numeric(df["Adjusted P-value"], errors="coerce") <= threshold).sum())


def _top_term(df: pd.DataFrame) -> str:
    if df.empty or "Term" not in df.columns:
        return "—"
    term = str(df["Term"].iloc[0])
    return (term[:60] + "…") if len(term) > 60 else term


# ---------------------------------------------------------------------------
# HTML renderer
# ---------------------------------------------------------------------------

def _render_page(title: str, sections: list[str], provenance: dict) -> str:
    body           = "\n".join(sections)
    gseapy_version = provenance.get("gseapy_version", "—")
    timestamp      = provenance.get("timestamp", "—")

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
      font-size: 0.85rem;
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
    .gene-list {{
      font-size: 0.78rem;
      color: #555;
      font-family: "SFMono-Regular", Consolas, monospace;
    }}

    /* Direction badges */
    .dir-up {{
      display: inline-block;
      background: #fdecea;
      color: #c0392b;
      font-weight: 600;
      font-size: 0.75rem;
      padding: 1px 7px;
      border-radius: 10px;
      white-space: nowrap;
    }}
    .dir-down {{
      display: inline-block;
      background: #e8f0fe;
      color: #2980b9;
      font-weight: 600;
      font-size: 0.75rem;
      padding: 1px 7px;
      border-radius: 10px;
      white-space: nowrap;
    }}

    /* Bar chart grid */
    .bar-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 18px;
      margin-top: 12px;
    }}
    .bar-wrap {{
      flex: 1 1 340px;
      max-width: 520px;
    }}
    .bar-wrap h3 {{
      font-size: 0.9rem;
      margin-bottom: 6px;
      color: #16213e;
    }}
    .bar-wrap img {{
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
    <h1>OmicSage — GSEA Pathway Enrichment Report</h1>
    <p>Generated {timestamp} · gseapy {gseapy_version}</p>
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
