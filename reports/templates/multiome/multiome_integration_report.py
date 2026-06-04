"""
OmicSage — Multiome Integration Report
reports/templates/multiome/multiome_integration_report.py

Generated after multiome_integration step.
Output: multiome_04_integration_report.html

Sections
--------
* Run Summary       — cells / modalities / method / batch info
* Joint Embedding   — UMAP coloured by batch and cell type (one panel per method)
* Pre-integration   — UMAP coloured by batch using ATAC LSI embedding (for comparison)
* Factor/Latent     — MOFA+ factor weights heatmap (if X_mofa); latent distribution (if X_multivi)
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
    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }
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

def _plot_umap_by_obs(
    mdata: MuData,
    umap_key: str,
    obs_key: str,
    title: str,
    obs_source: str = "mdata",   # "mdata" | "rna" | "atac"
) -> str:
    """
    UMAP scatter coloured by any categorical obs column.

    obs_source controls where the obs column is looked up:
      "mdata"  → mdata.obs  (top-level, after mu.update())
      "rna"    → mdata["rna"].obs
      "atac"   → mdata["atac"].obs
    """
    if umap_key not in mdata.obsm:
        return _placeholder(f"{umap_key} not found in mdata.obsm")

    umap = mdata.obsm[umap_key]

    # Resolve obs column
    if obs_source == "rna" and obs_key in mdata["rna"].obs.columns:
        labels = mdata["rna"].obs[obs_key].astype(str)
    elif obs_source == "atac" and obs_key in mdata["atac"].obs.columns:
        labels = mdata["atac"].obs[obs_key].astype(str)
    elif obs_key in mdata.obs.columns:
        labels = mdata.obs[obs_key].astype(str)
    elif obs_key in mdata["rna"].obs.columns:
        labels = mdata["rna"].obs[obs_key].astype(str)
    elif obs_key in mdata["atac"].obs.columns:
        labels = mdata["atac"].obs[obs_key].astype(str)
    else:
        return _placeholder(f"obs['{obs_key}'] not found")

    unique = sorted(labels.unique())
    cmap   = plt.get_cmap("tab20", max(len(unique), 2))

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, label in enumerate(unique):
        mask = (labels == label).values
        ax.scatter(
            umap[mask, 0], umap[mask, 1],
            s=2, alpha=0.6, color=cmap(i),
            label=label, rasterized=True,
        )

    n_cols = max(1, len(unique) // 15)
    ax.legend(markerscale=4, frameon=False, fontsize=7,
              loc="upper right", ncol=n_cols,
              title=obs_key, title_fontsize=7)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("UMAP 1", fontsize=9); ax.set_ylabel("UMAP 2", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_mofa_factor_variance(mdata: MuData) -> str:
    """
    Bar chart showing per-factor total variance explained (summed across modalities).
    Uses mdata.uns["mofa"] if available, otherwise falls back to X_mofa variance.
    """
    if "X_mofa" not in mdata.obsm:
        return _placeholder("X_mofa not found")

    factors = mdata.obsm["X_mofa"]  # (n_cells, n_factors)
    n_factors = factors.shape[1]

    # Proxy: variance of each factor across cells
    factor_var = factors.var(axis=0)
    factor_var_pct = 100 * factor_var / factor_var.sum()

    fig, ax = plt.subplots(figsize=(max(6, n_factors * 0.5), 4))
    x = np.arange(n_factors)
    ax.bar(x, factor_var_pct,
           color=plt.get_cmap("Blues")(0.65), width=0.7,
           edgecolor="white")
    ax.set_xlabel("MOFA+ Factor", fontsize=10)
    ax.set_ylabel("% variance (proxy)", fontsize=10)
    ax.set_title("MOFA+ factor variance across cells\n"
                 "(proxy: cell-level variance per factor)",
                 fontsize=11, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"F{i+1}" for i in x], fontsize=7, rotation=45)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


def _plot_multivi_latent_distribution(mdata: MuData) -> str:
    """
    Violin plot of the latent dimension distributions from MultiVI.
    Shows the first min(10, n_latent) dimensions.
    """
    if "X_multivi" not in mdata.obsm:
        return _placeholder("X_multivi not found")

    latent = mdata.obsm["X_multivi"]
    n_show = min(10, latent.shape[1])
    data   = [latent[:, i] for i in range(n_show)]

    fig, ax = plt.subplots(figsize=(max(6, n_show * 0.8), 4))
    parts = ax.violinplot(data, positions=np.arange(n_show),
                          showmedians=True, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(plt.get_cmap("Blues")(0.5))
        pc.set_alpha(0.7)
    ax.set_xticks(np.arange(n_show))
    ax.set_xticklabels([f"Z{i+1}" for i in range(n_show)], fontsize=8)
    ax.set_xlabel("MultiVI latent dimension", fontsize=10)
    ax.set_ylabel("Activation", fontsize=10)
    ax.set_title(f"MultiVI latent space — first {n_show} dimensions",
                 fontsize=11, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    return _fig_to_b64(fig)


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
  <title>OmicSage — Multiome Integration — {dataset_name}</title>
  <style>{_CSS}</style>
</head>
<body>
<header>
  <h1>OmicSage &mdash; Multiome Integration</h1>
  <p>{dataset_name} &nbsp;|&nbsp; {timestamp} &nbsp;|&nbsp; RNA + ATAC joint embedding</p>
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
    method      = metrics.get("method", "—")
    n_cells     = metrics.get("n_cells", "—")
    batch_key   = metrics.get("batch_key", "—")
    n_batches   = metrics.get("n_batches", "—")

    if method == "mofa":
        n_label = "Latent factors"
        n_val   = metrics.get("n_factors", "—")
        method_desc = "MOFA+ (linear factor model)"
    elif method == "multivi":
        n_label = "Latent dimensions"
        n_val   = metrics.get("n_latent", "—")
        method_desc = "MultiVI (deep generative model)"
    else:
        n_label = "Dimensions"
        n_val   = "—"
        method_desc = method

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Generated: {timestamp} &nbsp;|&nbsp; Dataset: {dataset_name}</p>
      <p class="info">
        Integration method: <strong>{method_desc}</strong>
        &nbsp;&bull;&nbsp;
        Batch key: <code>{batch_key}</code> ({n_batches} batches)
      </p>
      <div class="stat-grid">
        <div class="stat-card">
          <div class="stat-value">{n_cells:,}</div>
          <div class="stat-label">Cells</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">2</div>
          <div class="stat-label">Modalities (RNA + ATAC)</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_batches}</div>
          <div class="stat-label">Batches</div>
        </div>
        <div class="stat-card">
          <div class="stat-value">{n_val}</div>
          <div class="stat-label">{n_label}</div>
        </div>
      </div>
    </section>"""


def _section_joint_umaps(
    fig_batch: str,
    fig_celltype: Optional[str],
    method_label: str,
    embed_key: str,
) -> str:
    ct_panel = ""
    if fig_celltype:
        ct_panel = f"""
        <div class="fig-wrap">
          <h3>UMAP — cell type (atac_celltype / cell_type_vote)</h3>
          <img src="data:image/png;base64,{fig_celltype}" alt="UMAP cell type">
        </div>"""

    return f"""
    <section>
      <h2>Joint Embedding — {method_label} (<code>{embed_key}</code>)</h2>
      <p>
        UMAP computed from the joint RNA + ATAC latent space.
        Ideally, cells of the same biological type should cluster regardless
        of batch (no batch-driven separation visible in the batch UMAP).
      </p>
      <div class="fig-grid">
        <div class="fig-wrap">
          <h3>UMAP — batch</h3>
          <img src="data:image/png;base64,{fig_batch}" alt="UMAP batch">
        </div>
        {ct_panel}
      </div>
    </section>"""


def _section_pre_integration(fig_pre: str) -> str:
    return f"""
    <section>
      <h2>Pre-integration Reference — ATAC LSI UMAP</h2>
      <p>
        UMAP from the ATAC LSI embedding (<code>X_umap_atac</code>) coloured by
        batch.  Compare with the joint embedding above: after integration, batch
        clusters should be less pronounced.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_pre}" alt="Pre-integration UMAP">
        </div>
      </div>
    </section>"""


def _section_mofa_diagnostics(fig_var: str) -> str:
    return f"""
    <section>
      <h2>MOFA+ Diagnostics</h2>
      <p>
        Bar chart shows the proxy variance explained by each MOFA+ factor
        (variance of factor activations across cells).  Factors with low variance
        explain little variation and can be dropped from downstream analysis.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_var}" alt="MOFA+ factor variance">
        </div>
      </div>
    </section>"""


def _section_multivi_diagnostics(fig_latent: str) -> str:
    return f"""
    <section>
      <h2>MultiVI Diagnostics</h2>
      <p>
        Violin plot of the first 10 MultiVI latent dimensions across all cells.
        Well-trained latent spaces show approximately unit-normal distributions.
        Flat or degenerate dimensions may indicate training instability.
      </p>
      <div class="fig-grid">
        <div class="fig-wrap wide">
          <img src="data:image/png;base64,{fig_latent}" alt="MultiVI latent distributions">
        </div>
      </div>
    </section>"""


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_multiome_integration_report(
    mdata: MuData,
    metrics: dict,
    report_path: str = "reports/multiome_04_integration_report.html",
    dataset_name: str = "dataset",
    batch_key: str = "batch",
    celltype_key: Optional[str] = "atac_celltype",
) -> str:
    """
    Generate the multiome integration HTML report.

    Parameters
    ----------
    mdata : MuData
        MuData returned by run_mofa() or run_multivi().
    metrics : dict
        Metrics dict returned by the integration function.
    report_path : str
        Output HTML path.
    dataset_name : str
        Dataset name shown in the report header.
    batch_key : str
        obs column used to colour the batch UMAP.  Default: "batch".
    celltype_key : str or None
        obs column to colour the cell type UMAP.  Default: "atac_celltype".
        Pass None to skip the cell type panel.

    Returns
    -------
    str : Absolute path of the written HTML file.
    """
    out = Path(report_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    method = metrics.get("method", "mofa")
    print(f"Building multiome integration report ({method}) for '{dataset_name}' ...",
          flush=True)

    sections = [_section_summary(metrics, dataset_name, timestamp)]

    # ── Joint embedding UMAPs ────────────────────────────────────────────────
    if method == "mofa":
        embed_key    = "X_mofa"
        umap_key     = "X_umap_mofa"
        method_label = "MOFA+"
    else:
        embed_key    = "X_multivi"
        umap_key     = "X_umap_multivi"
        method_label = "MultiVI"

    print(f"  • joint UMAP batch ({umap_key}) ...", flush=True)
    fig_batch = _plot_umap_by_obs(
        mdata, umap_key, batch_key,
        f"{method_label} UMAP — batch",
    )

    fig_celltype = None
    if celltype_key:
        print(f"  • joint UMAP cell type ({celltype_key}) ...", flush=True)
        fig_celltype = _plot_umap_by_obs(
            mdata, umap_key, celltype_key,
            f"{method_label} UMAP — {celltype_key}",
        )

    sections.append(
        _section_joint_umaps(fig_batch, fig_celltype, method_label, embed_key)
    )

    # ── Pre-integration reference ─────────────────────────────────────────────
    if "X_umap_atac" in mdata["atac"].obsm:
        print("  • pre-integration ATAC LSI UMAP ...", flush=True)
        # Build a temporary MuData-compatible object using the ATAC obsm
        # by temporarily copying the ATAC umap to mdata.obsm
        _tmp_key = "_tmp_umap_atac_preint"
        mdata.obsm[_tmp_key] = mdata["atac"].obsm["X_umap_atac"]
        fig_pre = _plot_umap_by_obs(
            mdata, _tmp_key, batch_key,
            "Pre-integration — ATAC LSI UMAP coloured by batch",
        )
        del mdata.obsm[_tmp_key]
        sections.append(_section_pre_integration(fig_pre))

    # ── Method-specific diagnostics ──────────────────────────────────────────
    if method == "mofa":
        print("  • MOFA+ factor variance chart ...", flush=True)
        fig_var = _plot_mofa_factor_variance(mdata)
        sections.append(_section_mofa_diagnostics(fig_var))

    elif method == "multivi":
        print("  • MultiVI latent distribution violin ...", flush=True)
        fig_latent = _plot_multivi_latent_distribution(mdata)
        sections.append(_section_multivi_diagnostics(fig_latent))

    html = _render_page(sections, timestamp, dataset_name)
    out.write_text(html, encoding="utf-8")
    print(f"Report saved -> {out.resolve()}", flush=True)
    return str(out.resolve())
