"""
OmicSage — CITE-seq Integration Report
reports/templates/cite/cite_integration_report.py

Generated after integration step.
Output: cite_06_integration_report.html
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


def _scatter_by_labels(coords, labels, title: str,
                        xlabel="UMAP 1", ylabel="UMAP 2") -> str:
    unique = sorted(set(labels))
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, b in enumerate(unique):
        mask = np.array(labels) == b
        ax.scatter(coords[mask, 0], coords[mask, 1], s=2, alpha=0.6,
                   color=cmap(i), label=b, rasterized=True)
    ax.legend(markerscale=4, frameon=False, fontsize=8,
              ncol=max(1, len(unique) // 12))
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=9); ax.set_ylabel(ylabel, fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _scatter_plain(coords, title: str) -> str:
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(coords[:, 0], coords[:, 1], s=2, alpha=0.5,
               color="#4C78A8", rasterized=True)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_integration_umaps(mdata, method: str, batch_key: str) -> tuple:
    """Return (fig_batch_b64, fig_celltype_b64, fig_cluster_b64)."""
    umap_key = "X_umap_mofa" if method == "mofa" else "X_umap_totalVI"
    if umap_key not in mdata.obsm:
        placeholder = _scatter_plain(np.random.randn(10, 2), "Embedding not found")
        return placeholder, placeholder, placeholder

    umap = mdata.obsm[umap_key]

    # --- batch ---
    obs_top = mdata.obs
    obs_rna = mdata["rna"].obs
    batch_src = obs_top if batch_key in obs_top.columns else \
                obs_rna if batch_key in obs_rna.columns else None
    if batch_src is not None:
        fig_batch = _scatter_by_labels(
            umap, batch_src[batch_key].astype(str).values,
            f"{method.upper()} UMAP — coloured by {batch_key}",
        )
    else:
        fig_batch = _scatter_plain(umap, f"{method.upper()} UMAP")

    # --- cell type (from RNA obs["cell_type_vote"] if available) ---
    ct_col = next((c for c in ["cell_type_vote", "cell_type", "adt_celltype", "leiden"]
                   if c in obs_rna.columns or c in mdata["adt"].obs.columns), None)
    if ct_col and ct_col in obs_rna.columns:
        labels_ct = obs_rna[ct_col].astype(str).values
    elif ct_col and ct_col in mdata["adt"].obs.columns:
        labels_ct = mdata["adt"].obs[ct_col].astype(str).values
    else:
        labels_ct = None

    if labels_ct is not None:
        fig_ct = _scatter_by_labels(
            umap, labels_ct,
            f"{method.upper()} UMAP — coloured by {ct_col}",
        )
    else:
        fig_ct = _scatter_plain(umap, f"{method.upper()} UMAP")

    # --- ADT Leiden cluster ---
    adt_leiden = mdata["adt"].obs.get("leiden")
    if adt_leiden is not None:
        fig_cluster = _scatter_by_labels(
            umap, adt_leiden.astype(str).values,
            f"{method.upper()} UMAP — coloured by ADT Leiden cluster",
        )
    else:
        fig_cluster = _scatter_plain(umap, f"{method.upper()} UMAP")

    return fig_batch, fig_ct, fig_cluster


def _render_page(sections: list[str], timestamp: str, dataset_name: str) -> str:
    body = "\n".join(sections)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage -- CITE-seq Integration Report -- {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
  <header>
    <h1>OmicSage -- CITE-seq Integration Report</h1>
    <p>Dataset: <strong>{dataset_name}</strong> &middot; Generated {timestamp}</p>
  </header>
  <main>{body}</main>
  <footer>Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a> &middot; MIT License</footer>
</body>
</html>"""


def _section_summary(metrics: dict, method: str, dataset_name: str,
                     timestamp: str) -> str:
    if method == "mofa":
        dim_label = "MOFA+ factors"
        dim_value = str(metrics.get("n_factors", "?"))
    else:
        dim_label = "Latent dims"
        dim_value = str(metrics.get("n_latent", "?"))

    cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",     f"{metrics.get('n_cells', '?'):,}"),
            ("Method",    method.upper()),
            (dim_label,   dim_value),
            ("Batches",   str(metrics.get("n_batches", "?"))),
            ("Batch key", metrics.get("batch_key", "?")),
        ]
    )
    extra_rows = []
    if method == "mofa":
        extra_rows = [("MOFA+ embedding key", "X_mofa"),
                      ("UMAP key", "X_umap_mofa")]
    else:
        extra_rows = [("totalVI embedding key", "X_totalVI"),
                      ("UMAP key", "X_umap_totalVI"),
                      ("Max epochs", str(metrics.get("max_epochs", "?")))]

    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Batch values", ", ".join(str(b) for b in metrics.get("batch_values", []))),
        ] + extra_rows
    )
    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Dataset: <strong>{dataset_name}</strong> &middot; {timestamp}</p>
      <div class="stat-grid">{cards}</div>
      <table><thead><tr><th>Parameter</th><th>Value</th></tr></thead>
      <tbody>{rows}</tbody></table>
    </section>"""


def _section_figures(fig_batch: str, fig_ct: str, fig_cluster: str,
                     method: str, batch_key: str) -> str:
    return f"""
    <section>
      <h2>Joint Embedding ({method.upper()})</h2>
      <p>UMAP computed from the joint RNA + ADT latent space.</p>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>Coloured by {batch_key}</h3>
          <img src="data:image/png;base64,{fig_batch}" alt="Integration UMAP batch">
        </div>
        <div class="fig-wrap">
          <h3>Coloured by cell type</h3>
          <img src="data:image/png;base64,{fig_ct}" alt="Integration UMAP cell type">
        </div>
        <div class="fig-wrap">
          <h3>Coloured by ADT Leiden cluster</h3>
          <img src="data:image/png;base64,{fig_cluster}" alt="Integration UMAP cluster">
        </div>
      </div>
    </section>"""


def run_cite_integration_report(
    mdata,
    metrics: dict,
    report_path: str = "reports/cite_06_integration_report.html",
    dataset_name: str = "dataset",
) -> str:
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    method    = "mofa" if "X_mofa" in mdata.obsm else "totalvi"
    batch_key = metrics.get("batch_key", "batch")

    print(f"Building CITE integration report ({method.upper()}) for '{dataset_name}' ...",
          flush=True)
    fig_batch, fig_ct, fig_cluster = _plot_integration_umaps(mdata, method, batch_key)

    html = _render_page(
        sections=[
            _section_summary(metrics, method, dataset_name, timestamp),
            _section_figures(fig_batch, fig_ct, fig_cluster, method, batch_key),
        ],
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
