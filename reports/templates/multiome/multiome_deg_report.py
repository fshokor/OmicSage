"""
OmicSage — Multiome DEG Report
reports/templates/multiome/multiome_deg_report.py

Generated after multiome_deg step (Session 5).
Output: multiome_05_deg_report.html

Sections
--------
1. Run Summary       — cells / cell types / significant hits / method
2. RNA DEG           — top genes per ATAC-defined cell type (bar chart + table)
                       Only rendered when input was MuData.
3. DCA               — top accessible peaks per cell type (bar chart + table)
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
    header { background: linear-gradient(135deg, #0d2137 0%, #1b3a5c 60%, #2a5298 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1200px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #1b3a5c;
                  border-bottom: 2px solid #dde6f5; padding-bottom: 10px;
                  margin-bottom: 18px; }
    section h3 { font-size: 1rem; font-weight: 600; color: #0d2137;
                  margin: 16px 0 8px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
             border-left: 3px solid #f0c040; padding: 8px 12px;
             border-radius: 4px; margin-bottom: 14px; }
    .info { font-size: 0.82rem; color: #1a3a6e; background: #e8eef8;
             border-left: 3px solid #2a5298; padding: 8px 12px;
             border-radius: 4px; margin-bottom: 14px; }
    code { font-family: "SFMono-Regular", Consolas, monospace;
            background: #eef1fa; padding: 1px 5px; border-radius: 3px;
            font-size: 0.85em; }
    .stat-grid { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }
    .stat-card { background: #eef1fa; border-radius: 8px; padding: 14px 20px;
                  min-width: 130px; text-align: center; flex: 1 1 130px; }
    .stat-value { font-size: 1.4rem; font-weight: 700; color: #1b3a5c; }
    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }
    table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }
    th { background: #eef1fa; color: #1b3a5c; font-weight: 600;
          padding: 9px 12px; text-align: left; border-bottom: 2px solid #c8d4ec; }
    td { padding: 8px 12px; border-bottom: 1px solid #eee; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f5f7fc; }
    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }
    .fig-wrap { flex: 1 1 300px; max-width: 520px; }
    .fig-wrap.wide { flex: 1 1 100%; max-width: 100%; }
    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #0d2137; }
    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #d0daea; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa;
              padding: 24px 0 32px; }
    footer a { color: #1b3a5c; text-decoration: none; }
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

def _plot_top_features_bar(
    results: dict[str, pd.DataFrame],
    feature_col: str,
    top_n: int = 10,
    title_prefix: str = "Top features",
    color: str = "#2a5298",
) -> str:
    """
    Grid of horizontal bar charts — one per cell type — showing top_n features
    by absolute logfc.
    """
    groups = [g for g, df in results.items() if not df.empty]
    if not groups:
        return _placeholder(f"No significant {feature_col}s found")

    n_groups  = len(groups)
    n_cols    = min(3, n_groups)
    n_rows    = (n_groups + n_cols - 1) // n_cols
    fig_h     = max(3, n_rows * 3.5)
    fig_w     = n_cols * 5

    fig, axes = plt.subplots(n_rows, n_cols,
                              figsize=(fig_w, fig_h),
                              squeeze=False)

    for idx, group in enumerate(groups):
        ax  = axes[idx // n_cols][idx % n_cols]
        df  = results[group].head(top_n)
        if df.empty:
            ax.axis("off")
            continue

        features = df[feature_col].tolist()
        logfcs   = df["logfc"].abs().tolist()

        y = np.arange(len(features))
        ax.barh(y, logfcs[::-1], color=color, alpha=0.8, height=0.65)
        ax.set_yticks(y)
        ax.set_yticklabels(features[::-1], fontsize=7)
        ax.set_xlabel("|log2FC|", fontsize=8)
        ax.set_title(group, fontsize=9, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)

    # Hide unused axes
    for idx in range(len(groups), n_rows * n_cols):
        axes[idx // n_cols][idx % n_cols].axis("off")

    fig.suptitle(f"{title_prefix} — top {top_n} per cell type",
                 fontsize=11, fontweight="bold", y=1.01)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _build_summary_table_html(
    summary: pd.DataFrame,
    feature_col: str,
) -> str:
    """Render the top-features summary DataFrame as an HTML table."""
    if summary.empty:
        return "<p class='note'>No significant results to display.</p>"

    cols = ["group", "rank", feature_col, "logfc", "pval_adj"]
    cols = [c for c in cols if c in summary.columns]

    header = "".join(f"<th>{c}</th>" for c in cols)
    rows   = ""
    for _, row in summary.iterrows():
        cells = ""
        for c in cols:
            v = row[c]
            if c in ("logfc", "pval_adj"):
                v = f"{float(v):.4f}"
            cells += f"<td>{v}</td>"
        rows += f"<tr>{cells}</tr>"

    return f"<table><thead><tr>{header}</tr></thead><tbody>{rows}</tbody></table>"


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
  <title>OmicSage — Multiome DEG — {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
<header>
  <h1>OmicSage &mdash; Multiome Differential Analysis</h1>
  <p>{dataset_name} &nbsp;|&nbsp; {timestamp} &nbsp;|&nbsp; RNA DEG + Differential Chromatin Accessibility</p>
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
    n_cells     = metrics.get("n_cells", "—")
    n_types     = metrics.get("n_cell_types", "—")
    n_rna       = metrics.get("n_rna_significant", "—")
    n_dca       = metrics.get("n_dca_significant", "—")
    method      = metrics.get("method", "wilcoxon")
    groupby     = metrics.get("groupby", "—")
    input_type  = metrics.get("input_type", "—")

    anndata_note = ""
    if input_type == "anndata":
        anndata_note = (
            "<p class='note'>Input was an AnnData (ATAC only). "
            "RNA DEG was skipped. Re-run with a MuData containing both "
            "'rna' and 'atac' modalities to enable RNA DEG.</p>"
        )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Generated: {timestamp} &nbsp;|&nbsp; Dataset: {dataset_name}</p>
      {anndata_note}
      <p class="info">
        Groupby: <code>{groupby}</code>
        &nbsp;&bull;&nbsp; Method: <code>{method}</code>
        &nbsp;&bull;&nbsp; Input type: <code>{input_type}</code>
      </p>
      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-value">{n_cells:,}</div>
          <div class="stat-label">Cells</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_types}</div>
          <div class="stat-label">Cell types</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_rna}</div>
          <div class="stat-label">Significant RNA DEG hits</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_dca}</div>
          <div class="stat-label">Significant DCA peaks</div>
        </div>
      </div>
    </section>"""


def _section_rna_deg(
    rna_deg: dict,
    rna_summary: pd.DataFrame,
) -> str:
    if not rna_deg:
        return ""

    fig_bar = _plot_top_features_bar(
        rna_deg, feature_col="gene",
        title_prefix="RNA DEG", color="#1b5e20",
    )
    table_html = _build_summary_table_html(rna_summary, feature_col="gene")

    return f"""
    <section>
      <h2>RNA Differential Gene Expression</h2>
      <p>
        Genes differentially expressed between ATAC-defined cell types.
        Cell type labels transferred from <code>atac_celltype</code> to RNA
        cells using shared barcodes. Test: Wilcoxon rank-sum, BH correction.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_bar}" alt="RNA DEG bar charts">
        </div>
      </div>
      <h3>Top 5 genes per cell type</h3>
      {table_html}
    </section>"""


def _section_dca(
    dca: dict,
    dca_summary: pd.DataFrame,
) -> str:
    if not dca:
        return "<section><h2>Differential Chromatin Accessibility</h2>" \
               "<p class='note'>No significant DCA results.</p></section>"

    fig_bar = _plot_top_features_bar(
        dca, feature_col="peak",
        title_prefix="DCA", color="#1565c0",
    )
    table_html = _build_summary_table_html(dca_summary, feature_col="peak")

    return f"""
    <section>
      <h2>Differential Chromatin Accessibility (DCA)</h2>
      <p>
        Genomic peaks differentially accessible between cell types.
        Raw peak counts used for testing (Wilcoxon rank-sum, BH correction).
        Peak IDs are in <code>chr:start-end</code> format.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_bar}" alt="DCA bar charts">
        </div>
      </div>
      <h3>Top 5 peaks per cell type</h3>
      {table_html}
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_multiome_deg_report(
    deg_dict: dict,
    report_path: str = "reports/multiome_05_deg_report.html",
    dataset_name: str = "dataset",
) -> str:
    """
    Generate the multiome DEG HTML report.

    Parameters
    ----------
    deg_dict : dict
        Dict returned by multiome_deg() with keys:
        rna_deg, dca, rna_summary, dca_summary, provenance, input_type.
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

    print(f"Building multiome DEG report for '{dataset_name}' ...", flush=True)

    provenance  = deg_dict.get("provenance", {})
    rna_deg     = deg_dict.get("rna_deg", {})
    dca         = deg_dict.get("dca", {})
    rna_summary = deg_dict.get("rna_summary", pd.DataFrame())
    dca_summary = deg_dict.get("dca_summary", pd.DataFrame())

    sections = [
        _section_summary(provenance, dataset_name, timestamp),
    ]

    if rna_deg:
        print("  • RNA DEG bar charts ...", flush=True)
        sections.append(_section_rna_deg(rna_deg, rna_summary))

    print("  • DCA bar charts ...", flush=True)
    sections.append(_section_dca(dca, dca_summary))

    html = _render_page(sections, timestamp, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
