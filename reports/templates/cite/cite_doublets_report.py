"""
OmicSage — CITE-seq ADT Doublet Detection Report
reports/templates/cite/cite_doublets_report.py

Generated after doublets step.
Output: cite_02_doublets_report.html
"""

from __future__ import annotations

import base64
import io
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData

_DPI = 130
_CSS = """
    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; }
    header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
              color: white; padding: 32px 40px 24px; }
    header h1 { font-size: 1.8rem; font-weight: 700; letter-spacing: -0.5px; }
    header p  { font-size: 0.85rem; opacity: 0.7; margin-top: 4px; }
    main { max-width: 1100px; margin: 0 auto; padding: 32px 24px; }
    section { background: white; border-radius: 10px;
               box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }
    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }
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
    td { padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }
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


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_score_histogram(adt: AnnData) -> str:
    if "adt_doublet_score" not in adt.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "adt_doublet_score not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    scores = adt.obs["adt_doublet_score"].values.astype(float)
    n_doublets = int((scores > 0).sum())
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(scores, bins=40, color="#4C78A8", alpha=0.8, edgecolor="white", linewidth=0.4)
    ax.axvline(0, color="#c0392b", linewidth=1.5, linestyle="--",
               label=f"Threshold (score > 0)\n{n_doublets:,} doublets flagged")
    ax.set_xlabel("ADT Doublet Score", fontsize=10)
    ax.set_ylabel("Number of cells", fontsize=10)
    ax.set_title("ADT Doublet Score Distribution", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_score_by_donor(adt: AnnData) -> str:
    """Box plot of doublet scores per donor/batch."""
    score_col = "adt_doublet_score"
    batch_col = next((c for c in ["donor", "batch", "sample"] if c in adt.obs.columns), None)
    if score_col not in adt.obs.columns or batch_col is None:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "No batch column found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    groups = adt.obs.groupby(batch_col)[score_col].apply(list)
    labels = list(groups.index)
    data   = [groups[l] for l in labels]

    fig, ax = plt.subplots(figsize=(max(5, len(labels) * 0.6), 4))
    ax.boxplot(data, labels=labels, patch_artist=True,
               boxprops=dict(facecolor="#4C78A8", alpha=0.6),
               medianprops=dict(color="#16213e", linewidth=1.5),
               whiskerprops=dict(color="#555"), capprops=dict(color="#555"),
               flierprops=dict(marker=".", color="#aaa", markersize=2, alpha=0.5))
    ax.set_xlabel(batch_col, fontsize=10)
    ax.set_ylabel("ADT Doublet Score", fontsize=10)
    ax.set_title(f"Doublet Score by {batch_col}", fontsize=12, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq ADT Doublet Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq ADT Doublet Detection Report</h1>
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
            ("Cells before",    f"{metrics.get('n_cells_before', '?'):,}"),
            ("Doublets flagged", f"{metrics.get('n_doublets_detected', '?'):,}"),
            ("% doublets",      f"{metrics.get('pct_doublets', 0):.1f}%"),
            ("Cells after",     f"{metrics.get('n_cells_after', '?'):,}"),
            ("Threshold",       str(metrics.get("threshold", 2.5))),
            ("Filtered",        "Yes" if metrics.get("filter_doublets") else "No"),
        ]
    )
    pairs_eval = metrics.get("pairs_evaluated", [])
    pairs_skip = metrics.get("pairs_skipped", [])
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Pairs evaluated", ", ".join(f"{a}/{b}" for a, b in pairs_eval) or "none"),
            ("Pairs skipped",   ", ".join(f"{a}/{b}" for a, b in pairs_skip) or "none"),
            ("Filter applied",  "Yes — cells removed" if metrics.get("filter_doublets")
                                else "No — flagged only (obs[adt_predicted_doublet])"),
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


def _section_figures(fig_hist: str, fig_donor: str) -> str:
    return f"""
    <section>
      <h2>Figures</h2>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Doublet Score Distribution</h3>
          <img src="data:image/png;base64,{fig_hist}" alt="Doublet score histogram">
        </div>
        <div class="fig-wrap">
          <h3>Doublet Score by Donor/Batch</h3>
          <img src="data:image/png;base64,{fig_donor}" alt="Doublet score by donor">
        </div>
      </div>
    </section>"""


def run_cite_doublets_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_02_doublets_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building CITE doublets report for '{dataset_name}' ...", flush=True)
    fig_hist  = _plot_score_histogram(adt)
    fig_donor = _plot_score_by_donor(adt)

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_figures(fig_hist, fig_donor),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
