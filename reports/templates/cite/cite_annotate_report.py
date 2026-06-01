"""
OmicSage — CITE-seq ADT Annotation Report
reports/templates/cite/cite_annotate_report.py

Generated after annotate_adt step.
Output: cite_05_annotate_report.html
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


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def _plot_cluster_sizes(adt: AnnData) -> str:
    if "leiden" not in adt.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5, "leiden column not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

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
    if "X_umap_adt" not in adt.obsm or "leiden" not in adt.obs.columns:
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.text(0.5, 0.5, "UMAP or leiden not found", ha="center", va="center",
                transform=ax.transAxes, fontsize=11, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    umap   = adt.obsm["X_umap_adt"]
    labels = adt.obs["leiden"].astype(str)
    unique = sorted(labels.unique(), key=lambda x: int(x) if x.isdigit() else x)
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    # If cell types are annotated, use those labels in the legend
    ct_col = "adt_celltype" if "adt_celltype" in adt.obs.columns else None
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, cid in enumerate(unique):
        mask = labels == cid
        if ct_col:
            ct_vals = adt.obs.loc[mask, ct_col].unique()
            legend_label = ct_vals[0] if len(ct_vals) == 1 else f"Cluster {cid}"
        else:
            legend_label = f"Cluster {cid}"
        ax.scatter(umap[mask.values, 0], umap[mask.values, 1],
                   s=2, alpha=0.7, color=cmap(i),
                   label=legend_label, rasterized=True)

    ax.legend(markerscale=4, frameon=False, fontsize=8,
              loc="upper right", ncol=max(1, len(unique) // 12))
    title = "ADT UMAP — coloured by Leiden cluster"
    if ct_col:
        title += " (with cell type labels)"
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


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
            ("Cells",     f"{metrics.get('n_cells', '?'):,}"),
            ("Clusters",  str(metrics.get("n_clusters", "?"))),
            ("Resolution", str(metrics.get("resolution", "?"))),
            ("Annotated", "Yes" if metrics.get("annotated") else "No"),
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
            to the config and re-run this step with <code>--step annotate_adt --force</code>.</p>"""

    return f"""
    <section>
      <h2>Cluster Table</h2>
      {note}
      <table><thead>{header}</thead><tbody>{rows}</tbody></table>
    </section>"""


def _section_figures(fig_bar: str, fig_umap: str) -> str:
    return f"""
    <section>
      <h2>Figures</h2>
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


def run_cite_annotate_report(
    adt: AnnData,
    metrics: dict,
    report_path: str = "reports/cite_05_annotate_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    print(f"Building CITE annotate report for '{dataset_name}' ...", flush=True)
    fig_bar  = _plot_cluster_sizes(adt)
    fig_umap = _plot_umap_leiden(adt)

    html = _render_page(
        sections=[
            _section_summary(metrics, dataset_name, timestamp),
            _section_cluster_table(adt, metrics),
            _section_figures(fig_bar, fig_umap),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
