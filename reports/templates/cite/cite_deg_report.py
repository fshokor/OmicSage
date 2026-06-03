"""
OmicSage — CITE-seq DEG Report
reports/templates/cite/cite_deg_report.py

Generated after cite_deg step (Step 7).
Output: cite_07_deg_report.html

Two sections
------------
1. DPE (Differential Protein Expression)
   Top differentially expressed proteins per ADT-defined cell type.
   Bar charts of top proteins per cluster + summary table.

2. Cross-modal RNA DEG
   Top differentially expressed genes per ADT-defined cell type.
   Volcano-style ranked plot + summary table.
   Only rendered when input was MuData (cite_06 path).
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
import pandas as pd

_DPI = 130

_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1200px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6;
                  padding-bottom: 10px; margin-bottom: 18px; }
    section h3 { font-size: 1rem; font-weight: 600; color: #16213e;
                  margin: 16px 0 8px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
             border-left: 3px solid #f0c040; padding: 8px 12px;
             border-radius: 4px; margin-bottom: 14px; }
    code { font-family: "SFMono-Regular", Consolas, monospace;
            background: #f0f2ff; padding: 1px 5px;
            border-radius: 3px; font-size: 0.85em; }
    .stat-grid { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }
    .stat-card { background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }
    .stat-value { font-size: 1.4rem; font-weight: 700; color: #0f3460; }
    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }
    table { width: 100%; border-collapse: collapse;
             font-size: 0.88rem; margin-top: 8px; }
    th { background: #f0f2ff; color: #0f3460; font-weight: 600;
          padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee;
          vertical-align: middle; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9ff; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 520px; }
    .fig-wrap.wide { flex: 1 1 100%; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }
    img { max-width: 100%; height: auto; }
    .badge { display:inline-block; padding:2px 8px; border-radius:20px;
              font-size:0.78rem; font-weight:600; margin-left:6px; }
    .badge-blue   { background:#dbeafe; color:#1e40af; }
    .badge-orange { background:#ffedd5; color:#9a3412; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa;
              padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


# Platelet / megakaryocyte markers that should not be top DPE hits
# for non-platelet cell types in BMMC / PBMC data
_PLATELET_MARKERS = {"CD41", "CD42b", "CD62P", "CD42a", "CD61", "CD9"}
# Cell types where platelet markers are expected
_PLATELET_CELLTYPES = {"Platelet", "Megakaryocyte", "MK", "Thrombocyte"}


def _check_dpe_flags(dpe_results: dict[str, pd.DataFrame]) -> str:
    """
    Scan DPE results for two artifact classes and return an HTML warning block.

    1. Platelet marker contamination: CD41/CD42b/CD62P in top 5 for non-platelet
       cell types (platelet-monocyte adhesion artifact in BMMC preparations).
    2. All-negative top 5: cell type whose top 5 DPE hits are all log2FC < 0
       (cell type defined by what it lacks, not what it expresses).
    """
    flags: list[str] = []

    for ct, df in dpe_results.items():
        if df.empty:
            continue
        top5 = df.head(5)

        # --- Platelet contamination check ---
        is_platelet_ct = any(p.lower() in ct.lower() for p in _PLATELET_CELLTYPES)
        if not is_platelet_ct:
            prot_col = "protein" if "protein" in top5.columns else top5.columns[0]
            hits = set(top5[prot_col].astype(str))
            platelet_hits = hits & _PLATELET_MARKERS
            if platelet_hits:
                flags.append(
                    f"<li><strong>{ct}</strong>: platelet marker(s) "
                    f"<code>{'</code>, <code>'.join(sorted(platelet_hits))}</code> "
                    f"in top 5 DPE hits. In BMMC/PBMC data this indicates "
                    f"platelet-cell adhesion artefact — these proteins reflect "
                    f"platelet contamination, not genuine {ct} surface expression. "
                    f"Treat with caution.</li>"
                )

        # --- All-negative top 5 check ---
        if "logfc" in top5.columns:
            if (top5["logfc"].values < 0).all():
                flags.append(
                    f"<li><strong>{ct}</strong>: all top 5 DPE hits are "
                    f"downregulated (log2FC &lt; 0). This cell type is defined "
                    f"by low surface protein expression relative to the rest — "
                    f"check whether positive markers exist further down the ranked "
                    f"list, or whether this cluster may be under-annotated.</li>"
                )

    if not flags:
        return ""

    items = "\n".join(flags)
    return (
        f"<div class='note'><strong>⚠ DPE flags detected:</strong>"
        f"<ul style='margin:6px 0 0 16px;'>{items}</ul></div>"
    )


def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq DEG Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq DEG Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
  &middot; MIT License</footer>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _plot_top_features_bar(
    results: dict[str, pd.DataFrame],
    feature_col: str,
    n_top: int = 10,
    title_prefix: str = "Top features",
    color: str = "#4C78A8",
) -> str:
    """
    One horizontal bar chart per cell type showing top n_top features by logFC.
    All charts combined into a grid figure.
    """
    cell_types = [ct for ct, df in results.items() if not df.empty]
    if not cell_types:
        fig, ax = plt.subplots(figsize=(5, 2))
        ax.text(0.5, 0.5, "No significant results to display.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    n_cols = min(3, len(cell_types))
    n_rows = (len(cell_types) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(n_cols * 5, n_rows * 3.5))
    axes = np.array(axes).flatten()

    for i, ct in enumerate(cell_types):
        ax = axes[i]
        df = results[ct].head(n_top).copy()
        if df.empty:
            ax.axis("off")
            continue
        df = df.sort_values("logfc", ascending=True)
        bar_colors = [
            "#e07b3a" if v >= 0 else "#4C78A8"
            for v in df["logfc"]
        ]
        ax.barh(df[feature_col].astype(str), df["logfc"],
                color=bar_colors, alpha=0.85, height=0.65)
        ax.axvline(0, color="#aaa", linewidth=0.8, linestyle="--")
        ax.set_title(ct, fontsize=9, fontweight="bold")
        ax.set_xlabel("log2FC", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=7)
        ax.spines[["top", "right"]].set_visible(False)

    # Hide empty subplots
    for j in range(len(cell_types), len(axes)):
        axes[j].axis("off")

    fig.suptitle(title_prefix, fontsize=12, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_ranked_genes(
    results: dict[str, pd.DataFrame],
    n_top: int = 20,
    title: str = "Top DEGs by log2FC",
) -> str:
    """
    Dot plot: x = log2FC rank, y = cell type, dot size = -log10(pval_adj).
    Shows the top n_top genes per cell type as a summary.
    """
    cell_types = [ct for ct, df in results.items() if not df.empty]
    if not cell_types:
        fig, ax = plt.subplots(figsize=(5, 2))
        ax.text(0.5, 0.5, "No significant results to display.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    fig_height = max(4, len(cell_types) * 1.2)
    fig, ax = plt.subplots(figsize=(10, fig_height))

    y_positions = {ct: i for i, ct in enumerate(cell_types)}
    cmap = plt.get_cmap("tab10", len(cell_types))

    for i, ct in enumerate(cell_types):
        df = results[ct].head(n_top).copy()
        if df.empty:
            continue
        df = df.sort_values("logfc", ascending=False).reset_index(drop=True)
        x = df["logfc"].values
        pval_adj = df["pval_adj"].clip(lower=1e-300).values
        sizes = np.clip(-np.log10(pval_adj) * 8, 10, 200)

        ax.scatter(
            x,
            np.full(len(x), y_positions[ct]),
            s=sizes,
            color=cmap(i),
            alpha=0.7,
            edgecolors="white",
            linewidths=0.3,
        )

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()), fontsize=9)
    ax.axvline(0, color="#aaa", linewidth=0.8, linestyle="--")
    ax.set_xlabel("log2 Fold Change", fontsize=10)
    ax.set_title(
        f"{title}\n(dot size = −log₁₀(adj p-value), top {n_top} per group)",
        fontsize=11, fontweight="bold",
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML sections
# ---------------------------------------------------------------------------

def _section_run_summary(
    metrics: dict,
    dataset_name: str,
    timestamp: str,
) -> str:
    prov = metrics.get("provenance", metrics)
    n_cells      = prov.get("n_cells",      prov.get("n_obs", "?"))
    n_cts        = prov.get("n_cell_types", "?")
    n_dpe        = prov.get("n_dpe_significant", "?")
    n_rna        = prov.get("n_rna_crossmodal_significant", "?")
    groupby      = prov.get("groupby", "adt_celltype_manual")
    input_type   = prov.get("input_type", "?")
    method       = prov.get("method", "wilcoxon")
    min_logfc    = prov.get("min_logfc", "?")
    max_pval     = prov.get("max_pval_adj", "?")

    crossmodal_note = (
        "<p class='note'>Cross-modal RNA DEG was skipped — "
        "cite_deg received an AnnData (cite_05) instead of a MuData (cite_06). "
        "Re-run cite_deg with the cite_06 checkpoint to enable RNA DEG and GSEA.</p>"
        if input_type == "anndata" else ""
    )

    try:
        cells_display = f"{int(n_cells):,}"
    except (TypeError, ValueError):
        cells_display = str(n_cells) if n_cells not in (None, "?", "") else "?"

    cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",             cells_display),
            ("Cell types",        n_cts),
            ("DPE sig. proteins", n_dpe),
            ("RNA sig. genes",    n_rna),
        ]
    )
    rows = "".join(
        f"<tr><td>{k}</td><td><code>{v}</code></td></tr>"
        for k, v in [
            ("Groupby",     groupby),
            ("Method",      method),
            ("Min log2FC",  min_logfc),
            ("Max adj p",   max_pval),
            ("Input type",  input_type),
        ]
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong>
         &middot; {timestamp}</p>
      {crossmodal_note}
      <div class="stat-grid">{cards}</div>
      <table>
        <thead><tr><th>Parameter</th><th>Value</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>"""


def _section_dpe(
    dpe_results: dict[str, pd.DataFrame],
    dpe_summary: pd.DataFrame,
) -> str:
    if not dpe_results:
        return """
    <section>
      <h2>Differential Protein Expression (DPE)</h2>
      <p class="note">No DPE results available.</p>
    </section>"""

    n_sig = sum(len(df) for df in dpe_results.values())
    dpe_flags_html = _check_dpe_flags(dpe_results)
    fig_bars = _plot_top_features_bar(
        dpe_results, feature_col="protein",
        n_top=10,
        title_prefix="Top Differential Proteins per Cell Type (log2FC)",
        color="#4C78A8",
    )

    # Summary table — top 5 proteins per cell type
    if not dpe_summary.empty:
        trows = ""
        for _, row in dpe_summary.iterrows():
            logfc_color = "#c0392b" if float(row["logfc"]) >= 0 else "#1e40af"
            trows += (
                f"<tr>"
                f"<td>{row['group']}</td>"
                f"<td><strong>{row['protein']}</strong></td>"
                f"<td style='color:{logfc_color};font-weight:600'>"
                f"{float(row['logfc']):.3f}</td>"
                f"<td>{float(row['pval_adj']):.2e}</td>"
                f"</tr>"
            )
        table_html = f"""
        <table>
          <thead>
            <tr>
              <th>Cell type</th><th>Protein</th>
              <th>log2FC</th><th>adj p-value</th>
            </tr>
          </thead>
          <tbody>{trows}</tbody>
        </table>"""
    else:
        table_html = "<p>No significant proteins passed thresholds.</p>"

    return f"""
    <section>
      <h2>Differential Protein Expression (DPE)
        <span class="badge badge-blue">{n_sig} significant</span>
      </h2>
      <p>
        Wilcoxon rank-sum test on CLR-normalised ADT values, grouped by
        <code>adt_celltype_manual</code> (one-vs-rest). Proteins are filtered
        by log2FC &ge; threshold and BH-adjusted p-value &le; threshold.
        Isotype controls are excluded from results.
      </p>
      {dpe_flags_html}
      <h3>Top Differential Proteins per Cell Type</h3>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_bars}"
               alt="DPE bar charts per cell type">
        </div>
      </div>
      <h3>Top 5 Proteins per Cell Type</h3>
      {table_html}
    </section>"""


def _section_rna_crossmodal(
    rna_results: dict[str, pd.DataFrame],
    rna_summary: pd.DataFrame,
    input_type: str,
) -> str:
    if input_type == "anndata":
        return """
    <section>
      <h2>Cross-Modal RNA DEG</h2>
      <p class="note">
        Cross-modal RNA DEG requires MuData from the integration step (cite_06).
        cite_deg was run on an AnnData (cite_05) — only DPE was computed.
        Re-run with the cite_06 MuData to enable cross-modal DEG and GSEA.
      </p>
    </section>"""

    if not rna_results:
        return """
    <section>
      <h2>Cross-Modal RNA DEG</h2>
      <p class="note">No cross-modal RNA DEG results — no groups passed
      significance thresholds.</p>
    </section>"""

    n_sig = sum(len(df) for df in rna_results.values())

    # Flag cell types whose top 5 RNA DEGs are all downregulated
    rna_neg_flags = []
    for ct, df in rna_results.items():
        if not df.empty and "logfc" in df.columns:
            if (df.head(5)["logfc"].values < 0).all():
                rna_neg_flags.append(ct)
    rna_neg_html = ""
    if rna_neg_flags:
        cts = ", ".join(f"<strong>{c}</strong>" for c in rna_neg_flags)
        rna_neg_html = (
            f"<div class='note'>⚠ All top 5 RNA DEGs are downregulated for: "
            f"{cts}. These cell types are transcriptionally defined by what they "
            f"lack relative to other populations. This is biologically meaningful "
            f"(e.g. erythroid cells downregulate most transcription at maturity) "
            f"but check that upregulated markers exist further down the ranked list."
            f"</div>"
        )

    fig_ranked = _plot_ranked_genes(
        rna_results,
        n_top=20,
        title="Cross-Modal RNA DEG — Top Genes per Cell Type",
    )
    fig_bars = _plot_top_features_bar(
        rna_results, feature_col="gene",
        n_top=10,
        title_prefix="Top DEGs per ADT-Defined Cell Type (log2FC)",
        color="#e07b3a",
    )

    # Summary table
    if not rna_summary.empty:
        trows = ""
        for _, row in rna_summary.iterrows():
            logfc_color = "#c0392b" if float(row["logfc"]) >= 0 else "#1e40af"
            trows += (
                f"<tr>"
                f"<td>{row['group']}</td>"
                f"<td><strong>{row['gene']}</strong></td>"
                f"<td style='color:{logfc_color};font-weight:600'>"
                f"{float(row['logfc']):.3f}</td>"
                f"<td>{float(row['pval_adj']):.2e}</td>"
                f"</tr>"
            )
        table_html = f"""
        <table>
          <thead>
            <tr>
              <th>Cell type</th><th>Gene</th>
              <th>log2FC</th><th>adj p-value</th>
            </tr>
          </thead>
          <tbody>{trows}</tbody>
        </table>"""
    else:
        table_html = "<p>No significant genes passed thresholds.</p>"

    return f"""
    <section>
      <h2>Cross-Modal RNA DEG
        <span class="badge badge-orange">{n_sig} significant</span>
      </h2>
      <p>
        RNA-layer DEG using ADT-defined cell type labels as the grouping variable.
        Cells are grouped by <code>adt_celltype_manual</code> (surface phenotype),
        and DEG is run on the RNA layer. This identifies transcriptional programs
        underlying each immunophenotypically pure population. These gene lists
        are the input to the GSEA step (cite_08).
      </p>
      {rna_neg_html}
      <h3>Ranked DEG Dot Plot</h3>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_ranked}"
               alt="Cross-modal RNA DEG ranked dot plot">
        </div>
      </div>
      <h3>Top DEGs per Cell Type (bar charts)</h3>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_bars}"
               alt="Cross-modal RNA DEG bar charts">
        </div>
      </div>
      <h3>Top 5 Genes per Cell Type</h3>
      {table_html}
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_cite_deg_report(
    deg_dict: dict,
    report_path: str = "reports/cite_07_deg_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the CITE-seq DEG HTML report.

    Parameters
    ----------
    deg_dict : dict
        Dict returned by cite_deg(). Keys: dpe, rna_crossmodal,
        dpe_summary, rna_summary, provenance, input_type.
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Dataset display name for the report header.

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    dpe_results  = deg_dict.get("dpe", {})
    rna_results  = deg_dict.get("rna_crossmodal", {})
    dpe_summary  = deg_dict.get("dpe_summary", pd.DataFrame())
    rna_summary  = deg_dict.get("rna_summary", pd.DataFrame())
    input_type   = deg_dict.get("input_type", "mudata")
    provenance   = deg_dict.get("provenance", {})

    print(f"Building CITE DEG report for '{dataset_name}' ...", flush=True)

    sections = [
        _section_run_summary(
            {"provenance": provenance},
            dataset_name=dataset_name,
            timestamp=timestamp,
        ),
        _section_dpe(dpe_results, dpe_summary),
        _section_rna_crossmodal(rna_results, rna_summary, input_type),
    ]

    html = _render_page(sections, timestamp, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
