"""
OmicSage — Multiome GRN Report
reports/templates/multiome/multiome_grn_report.py

Generated after multiome_grn step.
Output: multiome_06_grn_report.html

Sections
--------
* Run Summary    — cells / TFs discovered / GRN edges / method info
* TF Activity    — bar chart of top TFs by mean AUCell score (RNA + ATAC)
* GRN Network    — top edges table (TF → target gene, scores, cell type)
* Cell Type View — per-cell-type bar chart of top active TFs
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
from mudata import MuData

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
    .stat-label { font-size: 0.78rem; color: #5a6a8a; margin-top: 4px; }
    table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }
    th { background: #eef1fa; color: #1b3a5c; font-weight: 600;
          padding: 8px 12px; text-align: left; }
    td { padding: 7px 12px; border-bottom: 1px solid #eee; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f5f7fd; }
    img { max-width: 100%; border-radius: 6px; margin-top: 12px; }
"""


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _fig_top_tfs(grn_df: pd.DataFrame, n_top: int = 20) -> str:
    """Bar chart of top TFs by mean combined_score across all cell types."""
    if grn_df.empty or "tf" not in grn_df.columns:
        return ""

    tf_scores = (
        grn_df.groupby("tf")["combined_score"]
        .mean()
        .sort_values(ascending=False)
        .head(n_top)
    )
    if tf_scores.empty:
        return ""

    fig, ax = plt.subplots(figsize=(10, max(3, len(tf_scores) * 0.35)))
    colors = plt.cm.Blues(np.linspace(0.4, 0.85, len(tf_scores)))[::-1]
    ax.barh(tf_scores.index[::-1], tf_scores.values[::-1], color=colors[::-1])
    ax.set_xlabel("Mean combined score", fontsize=10)
    ax.set_title(f"Top {len(tf_scores)} TFs by GRN activity", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _fig_celltype_tfs(grn_df: pd.DataFrame, n_top: int = 10) -> str:
    """Per-cell-type stacked bar showing top TF activity."""
    if grn_df.empty or "cell_type" not in grn_df.columns:
        return ""

    cell_types = grn_df["cell_type"].unique()
    if len(cell_types) == 0:
        return ""

    # Top TFs overall
    top_tfs = (
        grn_df.groupby("tf")["combined_score"]
        .mean()
        .sort_values(ascending=False)
        .head(n_top)
        .index.tolist()
    )

    pivot = (
        grn_df[grn_df["tf"].isin(top_tfs)]
        .groupby(["cell_type", "tf"])["combined_score"]
        .mean()
        .unstack(fill_value=0.0)
    )
    pivot = pivot.reindex(columns=top_tfs, fill_value=0.0)

    if pivot.empty:
        return ""

    fig, ax = plt.subplots(figsize=(10, max(3, len(cell_types) * 0.6 + 1)))
    cmap = plt.cm.get_cmap("tab10", len(top_tfs))
    bottom = np.zeros(len(pivot))
    for i, tf in enumerate(top_tfs):
        vals = pivot[tf].values
        ax.barh(pivot.index, vals, left=bottom, color=cmap(i), label=tf, height=0.6)
        bottom += vals

    ax.set_xlabel("Mean combined score", fontsize=10)
    ax.set_title("Top TF activity per cell type", fontsize=11, fontweight="bold")
    ax.tick_params(labelsize=9)
    ax.legend(
        title="TF", fontsize=7, title_fontsize=8,
        loc="lower right", ncol=2, framealpha=0.7,
    )
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML section builders
# ---------------------------------------------------------------------------

def _section_summary(result: dict, dataset_name: str, timestamp: str) -> str:
    prov      = result.get("provenance", {})
    params    = prov.get("params", {})
    n_rna     = result.get("n_tfs_rna", 0)
    n_atac    = result.get("n_tfs_atac", 0)
    n_edges   = result.get("n_grn_edges", 0)
    motif_db  = params.get("motif_db", "jaspar")
    groupby   = params.get("groupby", "atac_celltype")
    n_top     = params.get("n_top_peaks", 500)

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Generated: {timestamp}</p>
      <p>GRN inference on <strong>{dataset_name}</strong> using
         <code>pyscenic</code> (RNA regulons) + <code>decoupler</code>
         (ATAC motif enrichment, motif DB: <code>{motif_db}</code>).
         Cell types defined by <code>{groupby}</code>.
         Top <code>{n_top}</code> DCA peaks per cell type used for
         ATAC motif enrichment.</p>
      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-value">{n_rna}</div>
          <div class="stat-label">RNA TF regulons</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_atac}</div>
          <div class="stat-label">ATAC motif TFs</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_edges:,}</div>
          <div class="stat-label">GRN edges</div>
        </div>
      </div>
      <div class="info">
        Implementation: <code>pyscenic</code> GRN (correlation fallback if
        <code>arboreto</code> not installed) + <code>decoupler</code> AUCell
        for ATAC motif enrichment using JASPAR 2022 (via CollecTRI).
        Upgrade path: full SCENIC+ (pycistarget + cisTarget database) for
        coordinate-level peak–motif matching.
      </div>
    </section>"""


def _section_top_tfs(fig_b64: str) -> str:
    if not fig_b64:
        return """
    <section>
      <h2>TF Activity — Top TFs</h2>
      <div class="note">No TF scores available (pyscenic or decoupler not installed,
      or no TFs overlapped expression data).</div>
    </section>"""

    return f"""
    <section>
      <h2>TF Activity — Top TFs</h2>
      <p>Mean combined GRN score (RNA AUCell + ATAC motif enrichment) across
         all cell types. Higher score = stronger evidence for TF regulatory activity.</p>
      <img src="data:image/png;base64,{fig_b64}" alt="Top TFs">
    </section>"""


def _section_celltype_view(fig_b64: str) -> str:
    if not fig_b64:
        return """
    <section>
      <h2>TF Activity per Cell Type</h2>
      <div class="note">No per-cell-type TF scores available.</div>
    </section>"""

    return f"""
    <section>
      <h2>TF Activity per Cell Type</h2>
      <p>Top TFs broken down by cell type. Each bar segment represents the
         mean combined GRN score for that TF in the given cell type.</p>
      <img src="data:image/png;base64,{fig_b64}" alt="TF activity per cell type">
    </section>"""


def _section_grn_table(grn_df: pd.DataFrame, n_show: int = 50) -> str:
    if grn_df.empty:
        return """
    <section>
      <h2>GRN Edge Table</h2>
      <div class="note">No GRN edges produced. Check that pyscenic and decoupler
      are installed and that TFs overlap the expression matrix.</div>
    </section>"""

    show = grn_df.head(n_show)
    rows = ""
    for _, row in show.iterrows():
        rows += (
            f"<tr>"
            f"<td><code>{row['tf']}</code></td>"
            f"<td>{row['target_gene']}</td>"
            f"<td>{row.get('cell_type', '')}</td>"
            f"<td>{row.get('rna_score', 0):.3f}</td>"
            f"<td>{row.get('atac_score', 0):.3f}</td>"
            f"<td><strong>{row.get('combined_score', 0):.3f}</strong></td>"
            f"</tr>"
        )

    note = ""
    if len(grn_df) > n_show:
        note = (
            f'<p class="note">Showing top {n_show} of {len(grn_df):,} edges '
            f'(sorted by combined score). Full network in '
            f'<code>mdata.uns["grn_network"]</code>.</p>'
        )

    return f"""
    <section>
      <h2>GRN Edge Table</h2>
      {note}
      <table>
        <thead>
          <tr>
            <th>TF</th><th>Target gene</th><th>Cell type</th>
            <th>RNA score</th><th>ATAC score</th><th>Combined score</th>
          </tr>
        </thead>
        <tbody>{rows}</tbody>
      </table>
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_multiome_grn_report(
    result: dict,
    report_path: str = "reports/multiome_06_grn_report.html",
    dataset_name: str = "dataset",
    mdata: Optional[MuData] = None,
) -> str:
    """
    Generate the multiome GRN HTML report.

    Parameters
    ----------
    result : dict
        Dict returned by ``multiome_grn()``.
    report_path : str
        Output HTML path.
    dataset_name : str
        Dataset name shown in the report header.
    mdata : MuData, optional
        Not currently used; reserved for future embedding plots.

    Returns
    -------
    str : Absolute path of the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building multiome GRN report for '{dataset_name}' ...", flush=True)

    grn_df = result.get("grn_df", pd.DataFrame())

    fig_top_tfs  = _fig_top_tfs(grn_df)
    fig_cell_ct  = _fig_celltype_tfs(grn_df)

    sections = [
        _section_summary(result, dataset_name, timestamp),
        _section_top_tfs(fig_top_tfs),
        _section_celltype_view(fig_cell_ct),
        _section_grn_table(grn_df),
    ]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — GRN Report — {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage — Gene Regulatory Network Report</h1>
    <p>{dataset_name} &nbsp;·&nbsp; {timestamp}</p>
  </header>
  <main>
    {"".join(sections)}
  </main>
</body>
</html>"""

    out.write_text(html, encoding="utf-8")
    print(f"GRN report written → {out.resolve()}", flush=True)
    return str(out.resolve())
