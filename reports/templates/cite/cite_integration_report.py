"""
OmicSage — CITE-seq Integration Report
reports/templates/cite/cite_integration_report.py

Generated after integration step.
Output: cite_06_integration_report.html

Changes vs previous version
-----------------------------
- Before/after integration UMAP comparison (pre-integration ADT UMAP alongside
  the joint embedding UMAP, same colour keys, same layout).
- scib metrics table + bar chart (when metrics["scib"] or metrics["scib_comparison"]
  is present).
- color_keys parameter: pass any list of obs columns to colour the UMAPs
  (e.g. ["batch", "cell_type_vote", "donor"]).  Falls back to auto-detection
  when not supplied.
- run_both support: report renders MOFA+ and totalVI sections side-by-side when
  both embeddings are present.
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
    section { background: white; border-radius: 10px; box-shadow: 0 1px 4px rgba(0,0,0,0.07);
               padding: 28px 32px; margin-bottom: 24px; }
    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;
                  border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }
    section h3 { font-size: 1rem; font-weight: 600; color: #16213e; margin: 16px 0 8px; }
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
    .badge { display:inline-block; padding:2px 8px; border-radius:20px;
              font-size:0.78rem; font-weight:600; margin-left:6px; }
    .badge-good { background:#d4edda; color:#155724; }
    .badge-warn { background:#fff3cd; color:#856404; }
    footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }
    footer a { color: #0f3460; text-decoration: none; }
"""


def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


# ---------------------------------------------------------------------------
# Generic scatter helpers
# ---------------------------------------------------------------------------

def _scatter_by_labels(coords, labels, title: str,
                        xlabel="UMAP 1", ylabel="UMAP 2") -> str:
    unique = sorted(set(str(x) for x in labels))
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, b in enumerate(unique):
        mask = np.array([str(x) for x in labels]) == b
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


def _scatter_numeric(coords, values, title: str) -> str:
    v = np.asarray(values, dtype=float)
    vmin, vmax = np.nanpercentile(v, [2, 98])
    fig, ax = plt.subplots(figsize=(7, 5))
    sc = ax.scatter(coords[:, 0], coords[:, 1],
                    c=v, cmap="viridis", vmin=vmin, vmax=vmax,
                    s=2, alpha=0.6, rasterized=True)
    plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
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


def _resolve_labels(mdata, key: str) -> Optional[np.ndarray]:
    """Find obs column `key` across mdata.obs / rna.obs / adt.obs."""
    for obs in [mdata.obs, mdata["rna"].obs, mdata["adt"].obs]:
        if key in obs.columns:
            return obs[key].values
    return None


def _scatter_for_key(coords, mdata, key: str, title_prefix: str) -> str:
    """Render one scatter for a given obs key (auto numeric vs categorical)."""
    vals = _resolve_labels(mdata, key)
    if vals is None:
        return _scatter_plain(coords, f"{title_prefix}\n({key} not found)")

    import pandas as pd
    if pd.api.types.is_numeric_dtype(vals):
        return _scatter_numeric(coords, vals, f"{title_prefix} — {key}")
    else:
        return _scatter_by_labels(coords, vals, f"{title_prefix} — {key}")


# ---------------------------------------------------------------------------
# Before/after UMAP comparison
# ---------------------------------------------------------------------------

def _plot_before_after(
    mdata,
    method: str,
    color_keys: list[str],
) -> list[tuple[str, str, str]]:
    """
    Return a list of (key, before_b64, after_b64) tuples — one per color_key.
    'before' = pre-integration ADT UMAP (X_umap_adt from reduce step)
    'after'  = joint embedding UMAP (X_umap_mofa or X_umap_totalVI)
    """
    after_key = "X_umap_mofa" if method == "mofa" else "X_umap_totalVI"
    before_key = "X_umap_adt"  # written by adt_reduce.py

    panels = []
    for key in color_keys:
        # Before
        if before_key in mdata["adt"].obsm:
            before_coords = mdata["adt"].obsm[before_key]
            b64_before = _scatter_for_key(
                before_coords, mdata, key,
                f"Pre-integration ADT UMAP",
            )
        else:
            dummy = np.random.randn(mdata.n_obs, 2)
            b64_before = _scatter_plain(dummy, "Pre-integration UMAP not found")

        # After
        if after_key in mdata.obsm:
            after_coords = mdata.obsm[after_key]
            b64_after = _scatter_for_key(
                after_coords, mdata, key,
                f"Post-integration {method.upper()} UMAP",
            )
        else:
            dummy = np.random.randn(mdata.n_obs, 2)
            b64_after = _scatter_plain(dummy, f"{after_key} not found")

        panels.append((key, b64_before, b64_after))

    return panels


# ---------------------------------------------------------------------------
# scib metrics visualisation
# ---------------------------------------------------------------------------

def _plot_scib_bar(scib_dict: dict, title: str = "scib Integration Metrics") -> str:
    """
    Horizontal bar chart of scib metrics.
    Skips *_error keys.  Colour-codes by metric type:
      Batch correction (iLISI, graph_conn) — blue
      Bio conservation (cLISI, asw_label)  — orange
    Higher is better for all included metrics.
    """
    # Filter to numeric values only
    items = {k: v for k, v in scib_dict.items()
             if not k.endswith("_error") and isinstance(v, (int, float))}

    if not items:
        fig, ax = plt.subplots(figsize=(5, 2))
        ax.text(0.5, 0.5, "No numeric scib metrics to display.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    labels = list(items.keys())
    values = list(items.values())

    # Assign colours
    batch_keywords = {"ilisi", "graph_conn", "kbet"}
    colors = []
    for lbl in labels:
        if any(k in lbl.lower() for k in batch_keywords):
            colors.append("#4C78A8")   # blue = batch correction
        else:
            colors.append("#e07b3a")   # orange = bio conservation

    fig, ax = plt.subplots(figsize=(8, max(3, len(labels) * 0.5)))
    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, values, color=colors, alpha=0.8, height=0.55)

    # Value labels
    for bar, val in zip(bars, values):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", ha="left", fontsize=8)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlim(0, max(values) * 1.18 + 0.05)
    ax.set_xlabel("Score (higher is better)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.axvline(0, color="#aaa", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#4C78A8", alpha=0.8, label="Batch correction"),
        Patch(facecolor="#e07b3a", alpha=0.8, label="Bio conservation"),
    ]
    ax.legend(handles=legend_elements, frameon=False, fontsize=8,
              loc="lower right")

    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_scib_comparison_bar(comparison: dict) -> str:
    """
    Grouped bar chart comparing MOFA+ vs totalVI for each metric.
    Shows metrics that both methods have (prefix-matched: mofa_X vs totalvi_X).
    """
    mofa_metrics    = {k[5:]: v  for k, v in comparison.items()
                       if k.startswith("mofa_") and isinstance(v, (int, float))}
    totalvi_metrics = {k[8:]: v  for k, v in comparison.items()
                       if k.startswith("totalvi_") and isinstance(v, (int, float))}

    common_keys = [k for k in mofa_metrics if k in totalvi_metrics]
    if not common_keys:
        fig, ax = plt.subplots(figsize=(5, 2))
        ax.text(0.5, 0.5, "No common metrics to compare.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=10, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    x = np.arange(len(common_keys))
    width = 0.35

    fig, ax = plt.subplots(figsize=(max(6, len(common_keys) * 1.4), 4.5))
    bars_mofa    = ax.bar(x - width / 2,
                          [mofa_metrics[k]    for k in common_keys],
                          width=width, label="MOFA+",
                          color="#4C78A8", alpha=0.85)
    bars_totalvi = ax.bar(x + width / 2,
                          [totalvi_metrics[k] for k in common_keys],
                          width=width, label="totalVI",
                          color="#e07b3a", alpha=0.85)

    for bars in [bars_mofa, bars_totalvi]:
        for bar in bars:
            h = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, h + 0.005,
                    f"{h:.3f}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x)
    ax.set_xticklabels(common_keys, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("Score (higher is better)", fontsize=10)
    ax.set_title("scib Metrics: MOFA+ vs totalVI", fontsize=12, fontweight="bold")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# MOFA+ diagnostics (variance + weights)
# ---------------------------------------------------------------------------

def _plot_mofa_variance(mdata) -> str:
    mofa_uns = mdata.uns.get("mofa", {})
    variance = mofa_uns.get("variance")

    if not variance:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5,
                "Variance data not found in mdata.uns[mofa][variance].\n"
                "Run MOFA+ with muon.tl.mofa to populate this.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="#888")
        ax.axis("off")
        return _fig_to_b64(fig)

    view_arrays: dict[str, np.ndarray] = {}
    for view, val in variance.items():
        if isinstance(val, dict):
            arrays = list(val.values())
            view_arrays[view] = np.mean(arrays, axis=0)
        else:
            view_arrays[view] = np.asarray(val)

    views   = list(view_arrays.keys())
    n_facts = max(len(v) for v in view_arrays.values())
    factors = [f"F{i+1}" for i in range(n_facts)]

    view_colors = ["#4C78A8", "#e07b3a", "#3cb371", "#e15759", "#9467bd"]
    fig, ax = plt.subplots(figsize=(max(8, n_facts * 0.6), 4.5))

    bottoms = np.zeros(n_facts)
    for i, view in enumerate(views):
        vals = view_arrays[view] * 100
        ax.bar(factors, vals, bottom=bottoms,
               color=view_colors[i % len(view_colors)],
               alpha=0.85, label=view, width=0.65)
        bottoms += vals

    ax.set_xlabel("MOFA+ Factor", fontsize=10)
    ax.set_ylabel("Variance explained (%)", fontsize=10)
    ax.set_title(
        "MOFA+ Variance Explained per Factor\n(broken down by modality)",
        fontsize=12, fontweight="bold",
    )
    ax.legend(title="Modality", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_mofa_weights(mdata, n_top: int = 5, n_factors: int = 4) -> str:
    if "LFs" not in mdata.varm:
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.text(0.5, 0.5,
                "MOFA+ weights not found.\n"
                "Run MOFA+ with muon.tl.mofa to populate mdata.varm[LFs].",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=9, color="#888", wrap=True)
        ax.axis("off")
        return _fig_to_b64(fig)

    W_all     = np.asarray(mdata.varm["LFs"])
    n_rna     = mdata["rna"].n_vars
    n_adt     = mdata["adt"].n_vars
    rna_W     = W_all[:n_rna, :]
    adt_W     = W_all[n_rna:n_rna + n_adt, :]
    rna_names = list(mdata["rna"].var_names)
    adt_names = list(mdata["adt"].var_names)

    n_facts = min(n_factors, rna_W.shape[1], adt_W.shape[1])
    factor_labels = [f"F{i+1}" for i in range(n_facts)]

    def _top_weights(W, names, n_top, n_facts):
        top_idx = set()
        for f in range(n_facts):
            idx = np.argsort(np.abs(W[:, f]))[::-1][:n_top]
            top_idx.update(idx.tolist())
        top_idx = sorted(top_idx,
                         key=lambda i: np.abs(W[i, :n_facts]).max(),
                         reverse=True)
        top_idx = top_idx[:n_top * n_facts]
        labels  = [names[i] for i in top_idx]
        matrix  = W[np.array(top_idx), :n_facts]
        return labels, matrix

    rna_labels, rna_mat = _top_weights(rna_W, rna_names, n_top, n_facts)
    adt_labels, adt_mat = _top_weights(adt_W, adt_names, n_top, n_facts)

    fig, axes = plt.subplots(1, 2,
                             figsize=(5 + n_facts,
                                      max(4, len(rna_labels) * 0.35)))
    for ax, mat, labels, title in [
        (axes[0], rna_mat, rna_labels, f"Top RNA genes (n={len(rna_labels)})"),
        (axes[1], adt_mat, adt_labels, f"Top ADT proteins (n={len(adt_labels)})"),
    ]:
        vmax = np.abs(mat).max() or 1
        im = ax.imshow(mat, aspect="auto", cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(n_facts))
        ax.set_xticklabels(factor_labels, fontsize=9)
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("MOFA+ Factor", fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="Weight")

    fig.suptitle(
        f"MOFA+ Top Feature Weights "
        f"(top {n_top} per factor, top {n_facts} factors)",
        fontsize=11, fontweight="bold", y=1.02,
    )
    fig.tight_layout()
    return _fig_to_b64(fig)


# ---------------------------------------------------------------------------
# HTML sections
# ---------------------------------------------------------------------------

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
        dim_label, dim_value = "MOFA+ factors", str(metrics.get("n_factors", "?"))
    elif method == "totalvi":
        dim_label, dim_value = "Latent dims", str(metrics.get("n_latent", "?"))
    else:  # both
        dim_label = "Methods"
        dim_value = "MOFA+ + totalVI"

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
    batch_vals = metrics.get("batch_values") or (
        metrics.get("mofa", {}).get("batch_values", [])
        if method == "both" else []
    )
    rows = "".join(
        f"<tr><td>{k}</td><td>{v}</td></tr>"
        for k, v in [
            ("Batch values", ", ".join(str(b) for b in batch_vals)),
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


def _section_before_after(panels: list[tuple[str, str, str]], method: str) -> str:
    """Render before/after comparison rows — one row per color_key."""
    inner = ""
    for key, b64_before, b64_after in panels:
        inner += f"""
        <div style="margin-bottom:28px;">
          <h3 style="color:#0f3460;border-bottom:1px solid #eee;padding-bottom:6px;">
            Coloured by: <code>{key}</code>
          </h3>
          <div class="fig-grid">
            <div class="fig-wrap">
              <h3>Before integration (ADT UMAP)</h3>
              <img src="data:image/png;base64,{b64_before}" alt="Before integration UMAP">
            </div>
            <div class="fig-wrap">
              <h3>After integration ({method.upper()} UMAP)</h3>
              <img src="data:image/png;base64,{b64_after}" alt="After integration UMAP">
            </div>
          </div>
        </div>"""

    return f"""
    <section>
      <h2>Before / After Integration</h2>
      <p>
        Left column: ADT UMAP before integration (pre-Harmony, pre-joint embedding).<br>
        Right column: Joint RNA + ADT embedding after {method.upper()}.<br>
        A well-integrated embedding shows less batch separation and clearer
        cell-type clustering.
      </p>
      {inner}
    </section>"""


def _section_scib(scib_dict: dict, fig_bar: str,
                  is_comparison: bool = False) -> str:
    # Build metrics table
    rows = ""
    for k, v in scib_dict.items():
        if k.endswith("_error"):
            rows += f"<tr><td><code>{k}</code></td><td style='color:#c0392b'>{v}</td></tr>"
        elif isinstance(v, float):
            rows += f"<tr><td><code>{k}</code></td><td>{v:.4f}</td></tr>"
        else:
            rows += f"<tr><td><code>{k}</code></td><td>{v}</td></tr>"

    title = "scib Benchmark Comparison (MOFA+ vs totalVI)" if is_comparison \
            else "scib Integration Benchmark Metrics"
    desc = (
        "Comparing MOFA+ and totalVI on the same dataset using scib metrics. "
        "Higher is better for all metrics shown."
        if is_comparison else
        "scib integration benchmark. "
        "Batch correction: iLISI (higher = better batch mixing), "
        "graph_connectivity (higher = better global connectivity). "
        "Bio conservation: cLISI (higher = cleaner cell type clusters), "
        "ASW_label (higher = better silhouette)."
    )
    return f"""
    <section>
      <h2>{title}</h2>
      <p>{desc}</p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_bar}" alt="scib metrics bar chart">
        </div>
      </div>
      <table style="margin-top:16px;">
        <thead><tr><th>Metric</th><th>Value</th></tr></thead>
        <tbody>{rows}</tbody>
      </table>
    </section>"""


def _section_mofa_diagnostics(fig_variance: str, fig_weights: str) -> str:
    return (
        '<section>'
        '<h2>MOFA+ Factor Interpretation</h2>'
        '<p style="color:#444;font-size:0.9rem;margin-bottom:12px;">'
        'Variance decomposition shows which factors are driven by RNA vs ADT. '
        'The weights heatmap shows the top contributing features per factor — '
        'use this to assign biological meaning to each factor.</p>'
        '<div class="fig-grid">'
        '<div class="fig-wrap wide">'
        '<h3>Variance Explained per Factor (stacked by modality)</h3>'
        f'<img src="data:image/png;base64,{fig_variance}" alt="MOFA variance">'
        '</div>'
        '<div class="fig-wrap wide">'
        '<h3>Top Feature Weights per Factor (RNA | ADT)</h3>'
        f'<img src="data:image/png;base64,{fig_weights}" alt="MOFA weights">'
        '</div>'
        '</div>'
        '</section>'
    )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_cite_integration_report(
    mdata,
    metrics: dict,
    report_path: str = "reports/cite_06_integration_report.html",
    dataset_name: str = "dataset",
    color_keys: Optional[list[str]] = None,
) -> str:
    """
    Generate the CITE-seq integration HTML report.

    Parameters
    ----------
    mdata : MuData
        Integrated MuData (output of run_mofa / run_totalvi / run_both).
    metrics : dict
        Metrics dict returned by the integration function.
    report_path : str
        Where to write the HTML file.
    dataset_name : str
        Dataset display name shown in the report header.
    color_keys : list of str, optional
        obs columns to use for UMAP colouring in the before/after section.
        E.g. ``["batch", "cell_type_vote", "donor"]``.
        When None, auto-detects: uses batch_key + cell type column (if found).

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    method    = metrics.get("method", "mofa" if "X_mofa" in mdata.obsm else "totalvi")
    batch_key = metrics.get("batch_key", "batch")

    # Resolve color_keys
    if color_keys is None:
        color_keys = [batch_key]
        # Auto-add cell type column
        for ct in ["cell_type_vote", "cell_type", "adt_celltype"]:
            if _resolve_labels(mdata, ct) is not None:
                color_keys.append(ct)
                break

    print(f"Building CITE integration report ({method.upper()}) for '{dataset_name}' ...",
          flush=True)
    print(f"  color_keys: {color_keys}", flush=True)

    sections = [_section_summary(metrics, method, dataset_name, timestamp)]

    # ── Before / after comparison ──────────────────────────────────────────
    # For 'both', use mofa as the 'after' representative (first method run)
    after_method = "mofa" if method in ("mofa", "both") else "totalvi"
    print("  Rendering before/after UMAPs ...", flush=True)
    panels = _plot_before_after(mdata, after_method, color_keys)
    sections.append(_section_before_after(panels, after_method))

    # If 'both', also render totalVI before/after
    if method == "both" and "X_umap_totalVI" in mdata.obsm:
        print("  Rendering before/after UMAPs (totalVI) ...", flush=True)
        panels_tvi = _plot_before_after(mdata, "totalvi", color_keys)
        sections.append(_section_before_after(panels_tvi, "totalvi"))

    # ── MOFA+ diagnostics ──────────────────────────────────────────────────
    if method in ("mofa", "both"):
        print("  Rendering MOFA+ variance decomposition ...", flush=True)
        fig_variance = _plot_mofa_variance(mdata)
        print("  Rendering MOFA+ weights heatmap ...", flush=True)
        fig_weights  = _plot_mofa_weights(mdata)
        sections.append(_section_mofa_diagnostics(fig_variance, fig_weights))

    # ── scib metrics ───────────────────────────────────────────────────────
    if method == "both":
        comparison = metrics.get("scib_comparison")
        if comparison:
            print("  Rendering scib comparison bar chart ...", flush=True)
            fig_bar = _plot_scib_comparison_bar(comparison)
            sections.append(_section_scib(comparison, fig_bar, is_comparison=True))
    else:
        scib_dict = metrics.get("scib")
        if scib_dict:
            print("  Rendering scib metrics bar chart ...", flush=True)
            fig_bar = _plot_scib_bar(scib_dict)
            sections.append(_section_scib(scib_dict, fig_bar, is_comparison=False))

    html = _render_page(
        sections=sections,
        timestamp=timestamp,
        dataset_name=dataset_name,
    )
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
