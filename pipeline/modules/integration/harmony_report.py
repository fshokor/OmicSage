"""
harmony_report.py — HTML report generator for OmicSage Harmony batch correction module
Phase 1, Step 8

Generates a self-contained HTML report from harmony_correct() output.
Matches the style of gsea_report.py.

Usage:
    from reports.harmony_report import generate_harmony_report

    generate_harmony_report(
        adata=adata,
        output_path="reports/harmony_report.html",
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

def generate_harmony_report(
    adata,
    output_path: str = "reports/harmony_report.html",
    top_n_pcs: int = 20,
    max_batches_bar: int = 20,
) -> str:
    """
    Generate a self-contained HTML report for Harmony batch correction results.

    Parameters
    ----------
    adata : AnnData
        AnnData returned by harmony_correct(). Must contain
        ``obsm['X_pca_harmony']``, ``obsm['X_umap']``, and
        ``uns['omicsage_harmony']``.
    output_path : str
        Path to write the HTML file.
    top_n_pcs : int
        Number of PCs to show in the PC variance comparison plot.
    max_batches_bar : int
        Maximum number of batches to label individually in bar charts.
        Batches beyond this limit are grouped as "other".

    Returns
    -------
    str
        Absolute path to the written HTML file.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    provenance = adata.uns.get("omicsage_harmony", {})

    sections = []
    sections.append(_section_summary_stats(adata, provenance))
    sections.append(_section_batch_composition(adata, provenance, max_batches_bar))
    sections.append(_section_embedding_plots(adata, provenance))
    sections.append(_section_mixing_metrics(adata, provenance))
    sections.append(_section_pc_shift(adata, provenance, top_n_pcs))

    html = _render_page(
        title="OmicSage — Harmony Batch Correction Report",
        sections=sections,
        provenance=provenance,
    )

    output_path.write_text(html, encoding="utf-8")
    return str(output_path.resolve())


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

def _section_summary_stats(adata, provenance: dict) -> str:
    batch_key     = provenance.get("batch_key", "—")
    n_batches     = provenance.get("n_batches", "—")
    n_pcs         = provenance.get("n_pcs", "—")
    n_neighbors   = provenance.get("n_neighbors", "—")
    umap_min_dist = provenance.get("umap_min_dist", "—")
    max_iter      = provenance.get("max_iter_harmony", "—")
    theta         = provenance.get("theta", "—")
    elapsed       = provenance.get("elapsed_seconds", "—")
    embedding_key = provenance.get("embedding_key", "X_pca_harmony")
    neighbors_key = provenance.get("neighbors_key", "neighbors_harmony")
    umap_key      = provenance.get("umap_key", "X_umap_harmony")

    # Timestamp: try scanpy's logging or fall back to provenance
    import datetime
    timestamp = provenance.get(
        "timestamp",
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    )

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # Verify key outputs exist
    has_harmony_emb = embedding_key in adata.obsm
    has_umap        = umap_key in adata.obsm
    has_neighbors   = neighbors_key in adata.uns

    checks = {
        f"obsm['{embedding_key}']":  "✓" if has_harmony_emb else "✗ MISSING",
        f"obsm['{umap_key}']":       "✓" if has_umap        else "✗ MISSING",
        f"uns['{neighbors_key}']":   "✓" if has_neighbors   else "✗ MISSING",
    }
    check_rows = "".join(
        f"<tr><td><code>{k}</code></td>"
        f"<td style='color:{'#27ae60' if v == '✓' else '#c0392b'};font-weight:600'>{v}</td></tr>"
        for k, v in checks.items()
    )

    stat_cards = "".join(
        f'<div class="stat-card"><div class="stat-value">{v}</div>'
        f'<div class="stat-label">{k}</div></div>'
        for k, v in [
            ("Cells",           f"{n_cells:,}"),
            ("Genes",           f"{n_genes:,}"),
            ("Batches",         n_batches),
            ("PCs corrected",   n_pcs),
            ("Neighbours (k)",  n_neighbors),
            ("Elapsed (s)",     elapsed),
        ]
    )

    return f"""
    <section>
      <h2>Run Summary</h2>
      <p class="timestamp">Analysis run: {timestamp}</p>
      <p class="timestamp">
        batch_key=<code>{batch_key}</code> ·
        theta={theta} ·
        max_iter={max_iter} ·
        umap_min_dist={umap_min_dist}
      </p>
      <div class="stat-grid">{stat_cards}</div>
      <h3>Output Key Verification</h3>
      <table style="max-width:420px;">
        <thead>
          <tr><th>Key</th><th>Status</th></tr>
        </thead>
        <tbody>{check_rows}</tbody>
      </table>
    </section>
    """


def _section_batch_composition(
    adata,
    provenance: dict,
    max_batches_bar: int,
) -> str:
    batch_key = provenance.get("batch_key", "batch")
    if batch_key not in adata.obs.columns:
        return (
            "<section><h2>Batch Composition</h2>"
            f"<p>Column <code>{batch_key}</code> not found in obs.</p></section>"
        )

    counts = adata.obs[batch_key].astype(str).value_counts().sort_index()
    n_batches = len(counts)

    # Bar chart of cells per batch
    img_b64 = _render_batch_bar(counts, batch_key, max_batches_bar)

    # Table rows
    rows = ""
    pct_total = counts.sum()
    for batch, n in counts.items():
        pct = 100.0 * n / pct_total if pct_total > 0 else 0.0
        rows += (
            f"<tr><td>{batch}</td>"
            f"<td>{n:,}</td>"
            f"<td>{pct:.1f}%</td></tr>"
        )

    return f"""
    <section>
      <h2>Batch Composition</h2>
      <p>{n_batches} batches detected in <code>obs['{batch_key}']</code>.</p>
      <div style="display:flex;gap:32px;flex-wrap:wrap;align-items:flex-start;">
        <div style="flex:1 1 320px;max-width:540px;">
          <img src="data:image/png;base64,{img_b64}"
               alt="Batch cell counts"
               style="width:100%;border-radius:6px;border:1px solid #e8eaf6;">
        </div>
        <div style="flex:1 1 220px;">
          <table>
            <thead>
              <tr><th>Batch</th><th>Cells</th><th>%</th></tr>
            </thead>
            <tbody>{rows}</tbody>
          </table>
        </div>
      </div>
    </section>
    """


def _section_embedding_plots(adata, provenance: dict) -> str:
    """Side-by-side UMAP coloured by batch before and after correction."""
    batch_key     = provenance.get("batch_key", "batch")
    umap_key      = provenance.get("umap_key", "X_umap_harmony")
    precorr_key   = provenance.get("umap_precorrection_key", "X_umap_precorrection")
    embedding_key = provenance.get("embedding_key", "X_pca_harmony")

    if umap_key not in adata.obsm:
        return (
            "<section><h2>UMAP Embeddings</h2>"
            f"<p><code>obsm['{umap_key}']</code> not found — "
            "run harmony_correct() first.</p></section>"
        )

    # Convert to plain string Series — avoids MultiIndex/Categorical .map() bugs
    batch_col = (
        adata.obs[batch_key].astype(str)
        if batch_key in adata.obs.columns
        else None
    )

    # ------------------------------------------------------------------ #
    # Build colour palette for batches                                    #
    # ------------------------------------------------------------------ #
    if batch_col is not None:
        batches = sorted(batch_col.unique())
        cmap = plt.cm.get_cmap("tab20", len(batches))
        colour_map = {b: cmap(i) for i, b in enumerate(batches)}
        colours = [colour_map[b] for b in batch_col]
    else:
        colours = "#0f3460"
        batches = []

    umap_post  = adata.obsm[umap_key]                 # post-correction UMAP
    harmony_coords = adata.obsm[embedding_key][:, :2] # first 2 harmony PCs

    # Pre-correction UMAP — fall back to X_umap if precorrection key absent
    if precorr_key in adata.obsm:
        umap_pre = adata.obsm[precorr_key]
        pre_label = "Before Harmony\n(original UMAP, coloured by batch)"
    elif "X_umap" in adata.obsm:
        umap_pre = adata.obsm["X_umap"]
        pre_label = "Before Harmony\n(UMAP, coloured by batch)"
    else:
        # Last resort: raw PCA dims 0–1 (clearly labelled as such)
        umap_pre  = adata.obsm["X_pca"][:, :2]
        pre_label = "Before Harmony\n(PCA dim 1–2 fallback, coloured by batch)"

    # ------------------------------------------------------------------ #
    # Plot: pre-correction UMAP vs post-correction UMAP                  #
    # ------------------------------------------------------------------ #
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: pre-correction UMAP coloured by batch
    for batch in batches:
        mask = (batch_col == batch).values
        axes[0].scatter(
            umap_pre[mask, 0], umap_pre[mask, 1],
            s=3, alpha=0.5, label=str(batch),
            c=[colour_map[batch]],
        )
    axes[0].set_title(pre_label, fontsize=9, fontweight="bold")
    axes[0].set_xlabel("UMAP1", fontsize=8)
    axes[0].set_ylabel("UMAP2", fontsize=8)
    axes[0].tick_params(labelsize=7)
    axes[0].spines["top"].set_visible(False)
    axes[0].spines["right"].set_visible(False)

    # Right: post-correction UMAP coloured by batch
    for batch in batches:
        mask = (batch_col == batch).values
        axes[1].scatter(
            umap_post[mask, 0], umap_post[mask, 1],
            s=3, alpha=0.5, label=str(batch),
            c=[colour_map[batch]],
        )
    axes[1].set_title("After Harmony\n(UMAP on corrected embedding, coloured by batch)",
                      fontsize=9, fontweight="bold")
    axes[1].set_xlabel("UMAP1", fontsize=8)
    axes[1].set_ylabel("UMAP2", fontsize=8)
    axes[1].tick_params(labelsize=7)
    axes[1].spines["top"].set_visible(False)
    axes[1].spines["right"].set_visible(False)

    # Shared legend (max 20 items to keep plot readable)
    if batches and len(batches) <= 20:
        handles, labels = axes[1].get_legend_handles_labels()
        fig.legend(
            handles, labels,
            loc="lower center",
            ncol=min(len(batches), 5),
            fontsize=7,
            frameon=False,
            bbox_to_anchor=(0.5, -0.05),
            markerscale=3,
        )

    fig.tight_layout()
    b64_side = _fig_to_base64(fig)
    plt.close(fig)

    # ------------------------------------------------------------------ #
    # Second plot: post-correction UMAP coloured by harmony PC1 value    #
    # ------------------------------------------------------------------ #
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    sc2 = ax2.scatter(
        umap_post[:, 0], umap_post[:, 1],
        c=harmony_coords[:, 0],
        cmap="RdYlBu_r",
        s=3, alpha=0.6,
    )
    cbar = fig2.colorbar(sc2, ax=ax2, shrink=0.7, pad=0.02)
    cbar.set_label("Harmony PC1 value", fontsize=7)
    ax2.set_title("UMAP coloured by Harmony PC1", fontsize=9, fontweight="bold")
    ax2.set_xlabel("UMAP1", fontsize=8)
    ax2.set_ylabel("UMAP2", fontsize=8)
    ax2.tick_params(labelsize=7)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig2.tight_layout()
    b64_pc1 = _fig_to_base64(fig2)
    plt.close(fig2)

    return f"""
    <section>
      <h2>UMAP Embeddings</h2>
      <p>
        Left: UMAP computed <em>before</em> Harmony correction, coloured by
        <code>{batch_key}</code> — batch-driven separation is visible as
        distinct clouds per batch.<br>
        Right: UMAP recomputed on the Harmony-corrected embedding — batches
        should intermix while biological structure is preserved.
        Both plots use the same colour palette so batch clusters are
        directly comparable.
      </p>
      <img src="data:image/png;base64,{b64_side}"
           alt="Before/after UMAP"
           style="width:100%;border-radius:6px;border:1px solid #e8eaf6;margin-bottom:18px;">
      <h3>Harmony PC1 Distribution on Post-Correction UMAP</h3>
      <p>
        Cells with extreme Harmony PC1 values were shifted most by the correction.
        A smooth gradient (no sharp boundaries aligned with batch) suggests
        effective integration.
      </p>
      <div style="max-width:520px;">
        <img src="data:image/png;base64,{b64_pc1}"
             alt="UMAP coloured by Harmony PC1"
             style="width:100%;border-radius:6px;border:1px solid #e8eaf6;">
      </div>
    </section>
    """


def _section_mixing_metrics(adata, provenance: dict) -> str:
    """
    Compute and display simple batch-mixing metrics on the corrected UMAP.

    Metric: for each cell, compute the fraction of its k nearest neighbours
    (from the Harmony graph) that share the same batch. A well-mixed dataset
    should have a low batch-fraction (close to the expected random fraction
    = 1/n_batches).
    """
    batch_key     = provenance.get("batch_key", "batch")
    neighbors_key = provenance.get("neighbors_key", "neighbors_harmony")
    n_batches     = provenance.get("n_batches", 2)
    expected_frac = round(1.0 / max(n_batches, 1), 3)

    if neighbors_key not in adata.uns or batch_key not in adata.obs.columns:
        return (
            "<section><h2>Batch Mixing Metrics</h2>"
            "<p>Neighbor graph or batch column not available — skipping metrics.</p></section>"
        )

    try:
        # Retrieve connectivities from the harmony neighbor graph
        conn_key = f"{neighbors_key}_connectivities" if f"{neighbors_key}_connectivities" in adata.obsp else None
        # Scanpy stores connectivities under obsp with key from neighbors_key
        conn = None
        if conn_key and conn_key in adata.obsp:
            conn = adata.obsp[conn_key]
        elif "connectivities" in adata.obsp:
            conn = adata.obsp["connectivities"]

        if conn is None:
            raise ValueError("Connectivity matrix not found in obsp.")

        batch_labels = adata.obs[batch_key].astype(str).values
        n_cells = adata.n_obs

        # Compute per-cell same-batch fraction among neighbours
        same_batch_fracs = np.zeros(n_cells, dtype=float)
        cx = conn.tocsr()
        for i in range(n_cells):
            row = cx.getrow(i)
            neighbours = row.indices
            if len(neighbours) == 0:
                same_batch_fracs[i] = 1.0
                continue
            same = np.sum(batch_labels[neighbours] == batch_labels[i])
            same_batch_fracs[i] = same / len(neighbours)

        mean_frac  = float(np.mean(same_batch_fracs))
        median_frac = float(np.median(same_batch_fracs))

        # Mixing score: 1 = perfect mixing, 0 = no mixing
        # Normalised so that random = 1
        mixing_score = round(
            1.0 - max(0.0, (mean_frac - expected_frac) / (1.0 - expected_frac + 1e-9)),
            3,
        )

        img_b64 = _render_mixing_histogram(same_batch_fracs, expected_frac, batch_key)

        stat_cards = "".join(
            f'<div class="stat-card"><div class="stat-value">{v}</div>'
            f'<div class="stat-label">{k}</div></div>'
            for k, v in [
                ("Mean same-batch frac.",   f"{mean_frac:.3f}"),
                ("Median same-batch frac.", f"{median_frac:.3f}"),
                ("Expected (random)",       f"{expected_frac:.3f}"),
                ("Mixing score",            f"{mixing_score:.3f}"),
            ]
        )

        mixing_note = ""
        if mixing_score >= 0.8:
            mixing_note = '<p class="note">✓ Mixing score ≥ 0.8 — batches are well integrated.</p>'
        elif mixing_score >= 0.5:
            mixing_note = '<p class="note">⚠ Mixing score 0.5–0.8 — integration is moderate. Consider increasing theta or checking for biological confounding.</p>'
        else:
            mixing_note = '<p class="note">✗ Mixing score &lt; 0.5 — batches are poorly integrated. Consider higher theta, more iterations, or verifying batch_key is correct.</p>'

        return f"""
        <section>
          <h2>Batch Mixing Metrics</h2>
          <p>
            For each cell, the <em>same-batch fraction</em> measures what proportion of
            its k nearest neighbours (in the Harmony-corrected graph) share the same
            <code>{batch_key}</code> label. A perfectly random mixture gives
            {expected_frac:.3f} (= 1 / {n_batches} batches).
          </p>
          <div class="stat-grid">{stat_cards}</div>
          {mixing_note}
          <h3>Same-Batch Fraction Distribution</h3>
          <div style="max-width:560px;">
            <img src="data:image/png;base64,{img_b64}"
                 alt="Same-batch fraction histogram"
                 style="width:100%;border-radius:6px;border:1px solid #e8eaf6;">
          </div>
        </section>
        """

    except Exception as e:
        return (
            "<section><h2>Batch Mixing Metrics</h2>"
            f"<p>Mixing metrics could not be computed: {e}</p></section>"
        )


def _section_pc_shift(adata, provenance: dict, top_n_pcs: int) -> str:
    """
    Per-PC shift plot: mean absolute difference between X_pca and
    X_pca_harmony across cells, per PC dimension.

    Large shifts on early PCs suggest those PCs were dominated by batch.
    Ideally biological PCs (low index) shift less than noise PCs.
    """
    embedding_key = provenance.get("embedding_key", "X_pca_harmony")

    if "X_pca" not in adata.obsm or embedding_key not in adata.obsm:
        return (
            "<section><h2>Per-PC Correction Shift</h2>"
            "<p>X_pca or X_pca_harmony not found — skipping shift plot.</p></section>"
        )

    raw     = adata.obsm["X_pca"]
    harmony = adata.obsm[embedding_key]

    n_pcs_avail = min(raw.shape[1], harmony.shape[1], top_n_pcs)
    shift = np.mean(np.abs(raw[:, :n_pcs_avail] - harmony[:, :n_pcs_avail]), axis=0)

    img_b64 = _render_pc_shift(shift, n_pcs_avail)

    # Table: top 5 most-shifted PCs
    top_idx = np.argsort(shift)[::-1][:5]
    rows = "".join(
        f"<tr><td>PC{i+1}</td><td>{shift[i]:.4f}</td></tr>"
        for i in top_idx
    )

    return f"""
    <section>
      <h2>Per-PC Correction Shift</h2>
      <p>
        Mean absolute difference between raw PCA and Harmony-corrected values
        per principal component. PCs with large shift were adjusted most by
        Harmony — typically these correspond to batch-driven variance rather
        than biology.
      </p>
      <div style="display:flex;gap:32px;flex-wrap:wrap;align-items:flex-start;">
        <div style="flex:1 1 360px;max-width:580px;">
          <img src="data:image/png;base64,{img_b64}"
               alt="Per-PC correction shift"
               style="width:100%;border-radius:6px;border:1px solid #e8eaf6;">
        </div>
        <div style="flex:0 0 200px;">
          <h3>Top 5 Shifted PCs</h3>
          <table>
            <thead><tr><th>PC</th><th>Mean |shift|</th></tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>
      </div>
    </section>
    """


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _render_batch_bar(counts: pd.Series, batch_key: str, max_bars: int) -> str:
    """Horizontal bar chart of cells per batch."""
    if len(counts) > max_bars:
        top = counts.nlargest(max_bars - 1)
        other = pd.Series({"other": counts.iloc[max_bars - 1:].sum()})
        counts = pd.concat([top, other])

    n = len(counts)
    fig_h = max(3, n * 0.35 + 0.8)
    fig, ax = plt.subplots(figsize=(7, fig_h))

    colors = plt.cm.tab20(np.linspace(0, 1, n))
    ax.barh(range(n), counts.values, color=colors,
            edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(n))
    ax.set_yticklabels(counts.index.astype(str), fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cell count", fontsize=8)
    ax.set_title(f"Cells per {batch_key}", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Value labels
    for i, v in enumerate(counts.values):
        ax.text(v + counts.values.max() * 0.01, i, f"{v:,}",
                va="center", fontsize=7, color="#333")

    fig.tight_layout()
    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def _render_mixing_histogram(
    same_batch_fracs: np.ndarray,
    expected_frac: float,
    batch_key: str,
) -> str:
    """Histogram of per-cell same-batch neighbour fractions."""
    fig, ax = plt.subplots(figsize=(6, 3.5))

    ax.hist(same_batch_fracs, bins=40, color="#5c6bc0", edgecolor="white",
            linewidth=0.4, alpha=0.85)
    ax.axvline(expected_frac, color="#e74c3c", linestyle="--",
               linewidth=1.2, label=f"Expected random ({expected_frac:.3f})")
    ax.axvline(float(np.mean(same_batch_fracs)), color="#27ae60", linestyle="-",
               linewidth=1.2, label=f"Observed mean ({np.mean(same_batch_fracs):.3f})")
    ax.set_xlabel("Same-batch fraction", fontsize=8)
    ax.set_ylabel("Number of cells", fontsize=8)
    ax.set_title(f"Batch mixing — {batch_key}", fontsize=9, fontweight="bold")
    ax.legend(fontsize=7, frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    b64 = _fig_to_base64(fig)
    plt.close(fig)
    return b64


def _render_pc_shift(shift: np.ndarray, n_pcs: int) -> str:
    """Bar chart of mean absolute per-PC shift after Harmony correction."""
    fig, ax = plt.subplots(figsize=(max(6, n_pcs * 0.35 + 1.5), 3.5))

    x = np.arange(n_pcs)
    colors = plt.cm.RdYlBu_r(shift / (shift.max() + 1e-9))
    ax.bar(x, shift, color=colors, edgecolor="white", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"PC{i+1}" for i in x], rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Mean |X_pca − X_pca_harmony|", fontsize=8)
    ax.set_title("Per-PC Harmony correction shift", fontsize=9, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
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
    body      = "\n".join(sections)
    batch_key = provenance.get("batch_key", "—")
    n_batches = provenance.get("n_batches", "—")

    import datetime
    timestamp = provenance.get(
        "timestamp",
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
    )

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

    /* Note box (amber warning / success) */
    .note {{
      background: #fffbea;
      border-left: 4px solid #f0ad4e;
      padding: 10px 14px;
      border-radius: 0 6px 6px 0;
      font-size: 0.85rem;
      color: #7a5c00;
      margin: 12px 0;
    }}

    code {{
      background: #f0f2ff;
      border-radius: 3px;
      padding: 1px 5px;
      font-family: "SFMono-Regular", Consolas, monospace;
      font-size: 0.85em;
      color: #0f3460;
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
    <h1>OmicSage — Harmony Batch Correction Report</h1>
    <p>Generated {timestamp} · batch_key={batch_key} · {n_batches} batches</p>
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
