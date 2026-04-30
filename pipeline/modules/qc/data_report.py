"""
data_report.py
==============
Generates a self-contained HTML data intake report for any .h5ad file.

Two modes:
  - Public dataset  : pass --geo <accession> to fetch metadata from GEO
  - Personal dataset: omit --geo, report analyses the data only

Usage
-----
    # Public dataset with GEO context
    python pipeline/modules/qc/data_report.py \
        --input data/benchmark/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad \
        --geo GSE194122 \
        --output reports/data_intake_cite.html

    # Personal / private dataset
    python pipeline/modules/qc/data_report.py \
        --input data/mylab/sample.h5ad \
        --output reports/data_intake_sample.html

Output
------
A single self-contained .html file with:
  - Dataset summary card
  - QC metrics table
  - Distribution plots (genes/cell, counts/cell, MT%)
  - obs / var metadata inventory
  - GEO description (if --geo provided)
"""

import argparse
import base64
import io
import json
import re
import sys
import urllib.request
import urllib.error
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scanpy as sc


# ── Matplotlib style ────────────────────────────────────────────────────────

plt.rcParams.update({
    "font.family":      "DejaVu Sans",
    "font.size":        11,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "axes.grid":        True,
    "grid.alpha":       0.3,
    "grid.linewidth":   0.5,
    "figure.dpi":       130,
})

PALETTE = {
    "blue":   "#2563EB",
    "teal":   "#0D9488",
    "coral":  "#E85D04",
    "purple": "#7C3AED",
    "gray":   "#6B7280",
    "bg":     "#F8FAFC",
}


# ── GEO fetch ───────────────────────────────────────────────────────────────

def fetch_geo_metadata(accession: str) -> dict:
    """
    Fetch basic metadata from NCBI GEO for a given accession.
    Returns a dict with keys: title, organism, summary, pubmed_ids.
    Returns empty dict on any failure — caller handles gracefully.
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&term={accession}[Accession]&retmode=json"
    )
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            ids = json.loads(r.read())["esearchresult"]["idlist"]
        if not ids:
            return {}

        url2 = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&id={ids[0]}&retmode=json"
        )
        with urllib.request.urlopen(url2, timeout=10) as r:
            summary = json.loads(r.read())["result"][ids[0]]

        return {
            "title":      summary.get("title", ""),
            "organism":   summary.get("taxon", ""),
            "summary":    summary.get("summary", ""),
            "pubmed_ids": summary.get("pubmedids", []),
            "samples":    summary.get("n_samples", ""),
            "gse":        accession,
        }
    except Exception:
        return {}


# ── AnnData analysis ────────────────────────────────────────────────────────

def analyse(adata: sc.AnnData) -> dict:
    """Extract all metrics and statistics from an AnnData object."""
    results = {}

    # Basic shape
    results["n_cells"] = adata.n_obs
    results["n_genes"] = adata.n_vars

    # Modalities
    modalities = ["RNA"]
    if adata.obsm:
        keys = list(adata.obsm.keys())
        results["obsm_keys"] = keys
        if any("atac" in k.lower() for k in keys):
            modalities.append("ATAC")
        if any("adt"  in k.lower() or "protein" in k.lower() for k in keys):
            modalities.append("Protein (ADT)")
    else:
        results["obsm_keys"] = []

    # Check layers for modality hints
    if adata.layers:
        results["layers"] = list(adata.layers.keys())
        for l in adata.layers.keys():
            if "atac" in l.lower() and "ATAC" not in modalities:
                modalities.append("ATAC")
            if ("adt" in l.lower() or "protein" in l.lower()) and "Protein (ADT)" not in modalities:
                modalities.append("Protein (ADT)")
    else:
        results["layers"] = []

    results["modalities"] = modalities

    # Sparsity
    import scipy.sparse as sp
    X = adata.X
    if sp.issparse(X):
        nnz   = X.nnz
        total = X.shape[0] * X.shape[1]
    else:
        nnz   = int(np.count_nonzero(X))
        total = X.size
    results["sparsity_pct"] = round(100 * (1 - nnz / total), 2) if total > 0 else 0

    # Per-cell stats
    if sp.issparse(X):
        counts_per_cell = np.asarray(X.sum(axis=1)).flatten()
        genes_per_cell  = np.asarray((X > 0).sum(axis=1)).flatten()
    else:
        counts_per_cell = X.sum(axis=1)
        genes_per_cell  = (X > 0).sum(axis=1)

    results["counts_per_cell"] = counts_per_cell
    results["genes_per_cell"]  = genes_per_cell
    results["median_counts"]   = round(float(np.median(counts_per_cell)), 1)
    results["median_genes"]    = round(float(np.median(genes_per_cell)), 1)
    results["total_counts"]    = int(counts_per_cell.sum())

    # Mitochondrial genes
    mt_mask = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    results["has_mt"] = bool(mt_mask.any())
    if results["has_mt"]:
        if sp.issparse(X):
            mt_counts   = np.asarray(X[:, mt_mask].sum(axis=1)).flatten()
        else:
            mt_counts   = X[:, mt_mask].sum(axis=1)
        mt_pct = mt_counts / (counts_per_cell + 1e-9) * 100
        results["mt_pct"]        = mt_pct
        results["median_mt_pct"] = round(float(np.median(mt_pct)), 2)

    # obs columns inventory
    results["obs_cols"] = list(adata.obs.columns)
    results["var_cols"] = list(adata.var.columns)

    # Cell type / cluster annotation
    label_hints = ["cell_type", "celltype", "leiden", "louvain",
                   "cluster", "annotation", "label"]
    results["label_cols"] = [
        c for c in adata.obs.columns
        if any(h in c.lower() for h in label_hints)
    ]

    # Batch / donor columns
    batch_hints = ["batch", "donor", "sample", "site", "condition",
                   "patient", "replicate"]
    results["batch_cols"] = [
        c for c in adata.obs.columns
        if any(h in c.lower() for h in batch_hints)
    ]

    # Unique values for label cols (top 5)
    results["label_previews"] = {}
    for col in results["label_cols"][:3]:
        vals = adata.obs[col].astype(str).unique().tolist()
        results["label_previews"][col] = vals[:8]

    # Embeddings available
    results["embeddings"] = [
        k for k in adata.obsm.keys()
        if any(e in k.lower() for e in ["umap", "tsne", "pca", "lsi"])
    ]

    return results


# ── Plots ───────────────────────────────────────────────────────────────────

def _fig_to_b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight",
                facecolor=PALETTE["bg"])
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def make_plots(metrics: dict) -> dict:
    """Generate all distribution plots, return as base64 PNG strings."""
    plots = {}
    counts = metrics["counts_per_cell"]
    genes  = metrics["genes_per_cell"]

    n_plots = 3 if metrics.get("has_mt") else 2

    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 4),
                             facecolor=PALETTE["bg"])
    fig.patch.set_facecolor(PALETTE["bg"])

    # Counts per cell
    ax = axes[0]
    ax.hist(np.log1p(counts), bins=60, color=PALETTE["blue"],
            edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_xlabel("log1p(total counts)")
    ax.set_ylabel("cells")
    ax.set_title("Total counts per cell")
    ax.axvline(np.log1p(np.median(counts)), color=PALETTE["coral"],
               linewidth=1.5, linestyle="--", label=f"median={metrics['median_counts']:,.0f}")
    ax.legend(fontsize=9)
    ax.set_facecolor(PALETTE["bg"])

    # Genes per cell
    ax = axes[1]
    ax.hist(genes, bins=60, color=PALETTE["teal"],
            edgecolor="white", linewidth=0.3, alpha=0.85)
    ax.set_xlabel("genes detected")
    ax.set_ylabel("cells")
    ax.set_title("Genes detected per cell")
    ax.axvline(np.median(genes), color=PALETTE["coral"],
               linewidth=1.5, linestyle="--", label=f"median={metrics['median_genes']:,.0f}")
    ax.legend(fontsize=9)
    ax.set_facecolor(PALETTE["bg"])

    # MT%
    if metrics.get("has_mt"):
        ax = axes[2]
        ax.hist(metrics["mt_pct"], bins=60, color=PALETTE["purple"],
                edgecolor="white", linewidth=0.3, alpha=0.85)
        ax.set_xlabel("MT%")
        ax.set_ylabel("cells")
        ax.set_title("Mitochondrial % per cell")
        ax.axvline(metrics["median_mt_pct"], color=PALETTE["coral"],
                   linewidth=1.5, linestyle="--",
                   label=f"median={metrics['median_mt_pct']:.1f}%")
        ax.legend(fontsize=9)
        ax.set_facecolor(PALETTE["bg"])

    fig.tight_layout(pad=2)
    plots["distributions"] = _fig_to_b64(fig)

    # Counts vs genes scatter (sample 5000 cells max)
    n_sample = min(5000, metrics["n_cells"])
    idx      = np.random.choice(metrics["n_cells"], n_sample, replace=False)
    fig2, ax2 = plt.subplots(figsize=(5, 4), facecolor=PALETTE["bg"])
    fig2.patch.set_facecolor(PALETTE["bg"])
    ax2.scatter(np.log1p(counts[idx]), genes[idx],
                s=1.5, alpha=0.35, color=PALETTE["blue"], linewidths=0)
    ax2.set_xlabel("log1p(total counts)")
    ax2.set_ylabel("genes detected")
    ax2.set_title(f"Counts vs genes  (n={n_sample:,} cells sampled)")
    ax2.set_facecolor(PALETTE["bg"])
    fig2.tight_layout(pad=2)
    plots["scatter"] = _fig_to_b64(fig2)

    return plots


# ── HTML template ────────────────────────────────────────────────────────────

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>OmicSage · Data Intake Report · {title_short}</title>
<style>
  :root {{
    --blue:   #2563EB;
    --teal:   #0D9488;
    --coral:  #E85D04;
    --purple: #7C3AED;
    --gray:   #6B7280;
    --bg:     #F8FAFC;
    --card:   #FFFFFF;
    --border: #E2E8F0;
    --text:   #1E293B;
    --muted:  #64748B;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{
    font-family: 'Segoe UI', system-ui, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.6;
    padding: 2rem 1rem;
  }}
  .container {{ max-width: 960px; margin: 0 auto; }}

  /* Header */
  .header {{
    border-left: 4px solid var(--blue);
    padding: 0.5rem 1rem;
    margin-bottom: 2rem;
  }}
  .header h1 {{ font-size: 1.6rem; font-weight: 700; color: var(--text); }}
  .header .subtitle {{ color: var(--muted); font-size: 0.9rem; margin-top: 0.25rem; }}
  .badge {{
    display: inline-block;
    background: var(--blue);
    color: white;
    font-size: 0.72rem;
    font-weight: 600;
    padding: 0.2rem 0.6rem;
    border-radius: 99px;
    margin-right: 0.3rem;
    letter-spacing: 0.03em;
  }}
  .badge.teal   {{ background: var(--teal); }}
  .badge.coral  {{ background: var(--coral); }}
  .badge.purple {{ background: var(--purple); }}
  .badge.gray   {{ background: var(--gray); }}

  /* Cards */
  .card {{
    background: var(--card);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 1.5rem;
    margin-bottom: 1.5rem;
  }}
  .card h2 {{
    font-size: 1rem;
    font-weight: 700;
    color: var(--blue);
    text-transform: uppercase;
    letter-spacing: 0.06em;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid var(--border);
  }}

  /* Metrics grid */
  .metrics-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
    gap: 1rem;
  }}
  .metric {{
    background: var(--bg);
    border-radius: 8px;
    padding: 1rem;
    text-align: center;
  }}
  .metric .value {{
    font-size: 1.8rem;
    font-weight: 700;
    color: var(--blue);
    line-height: 1;
  }}
  .metric .label {{
    font-size: 0.78rem;
    color: var(--muted);
    margin-top: 0.3rem;
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}

  /* Table */
  table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.88rem;
  }}
  th {{
    background: var(--bg);
    text-align: left;
    padding: 0.5rem 0.75rem;
    font-size: 0.78rem;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    color: var(--muted);
    border-bottom: 1px solid var(--border);
  }}
  td {{
    padding: 0.5rem 0.75rem;
    border-bottom: 1px solid var(--border);
    vertical-align: top;
  }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: var(--bg); }}
  .mono {{ font-family: 'Consolas', monospace; font-size: 0.83rem; }}

  /* Plots */
  .plot-row {{
    display: grid;
    grid-template-columns: 3fr 2fr;
    gap: 1rem;
  }}
  .plot-row img, .card img {{ width: 100%; border-radius: 6px; }}

  /* GEO section */
  .geo-summary {{
    background: var(--bg);
    border-radius: 8px;
    padding: 1rem 1.25rem;
    font-size: 0.9rem;
    color: var(--text);
    line-height: 1.7;
  }}
  .personal-note {{
    background: #FFF7ED;
    border-left: 3px solid var(--coral);
    border-radius: 0 8px 8px 0;
    padding: 0.75rem 1rem;
    font-size: 0.88rem;
    color: #92400E;
  }}

  /* Footer */
  .footer {{
    text-align: center;
    color: var(--muted);
    font-size: 0.8rem;
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border);
  }}

  /* Tag list */
  .tag-list {{ display: flex; flex-wrap: wrap; gap: 0.4rem; }}
  .tag {{
    background: #EFF6FF;
    color: var(--blue);
    border: 1px solid #BFDBFE;
    border-radius: 4px;
    padding: 0.15rem 0.5rem;
    font-size: 0.8rem;
    font-family: 'Consolas', monospace;
  }}
  .tag.green {{
    background: #F0FDF4;
    color: #166534;
    border-color: #BBF7D0;
  }}
  .tag.orange {{
    background: #FFF7ED;
    color: #9A3412;
    border-color: #FED7AA;
  }}
</style>
</head>
<body>
<div class="container">

  <!-- Header -->
  <div class="header">
    <div style="margin-bottom:0.5rem">
      <span class="badge">OmicSage</span>
      <span class="badge gray">Data Intake Report</span>
      {modality_badges}
    </div>
    <h1>{title_short}</h1>
    <div class="subtitle">Generated {timestamp} · {input_path}</div>
  </div>

  <!-- Source info (GEO or personal) -->
  <div class="card">
    <h2>Dataset Source</h2>
    {source_section}
  </div>

  <!-- Key metrics -->
  <div class="card">
    <h2>Key Metrics</h2>
    <div class="metrics-grid">
      <div class="metric">
        <div class="value">{n_cells}</div>
        <div class="label">Cells</div>
      </div>
      <div class="metric">
        <div class="value">{n_genes}</div>
        <div class="label">Features</div>
      </div>
      <div class="metric">
        <div class="value">{median_genes}</div>
        <div class="label">Median genes/cell</div>
      </div>
      <div class="metric">
        <div class="value">{median_counts}</div>
        <div class="label">Median counts/cell</div>
      </div>
      <div class="metric">
        <div class="value">{sparsity}%</div>
        <div class="label">Sparsity</div>
      </div>
      <div class="metric">
        <div class="value">{mt_stat}</div>
        <div class="label">Median MT%</div>
      </div>
    </div>
  </div>

  <!-- QC distributions -->
  <div class="card">
    <h2>QC Distributions</h2>
    <div class="plot-row">
      <img src="data:image/png;base64,{plot_distributions}" alt="QC distributions">
      <img src="data:image/png;base64,{plot_scatter}" alt="Counts vs genes scatter">
    </div>
  </div>

  <!-- Metadata inventory -->
  <div class="card">
    <h2>Metadata Inventory</h2>
    <table>
      <thead>
        <tr>
          <th>Category</th>
          <th>Columns found</th>
          <th>Details</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td><strong>Cell annotations</strong></td>
          <td>{n_label_cols} columns</td>
          <td>{label_cols_html}</td>
        </tr>
        <tr>
          <td><strong>Batch / donor</strong></td>
          <td>{n_batch_cols} columns</td>
          <td>{batch_cols_html}</td>
        </tr>
        <tr>
          <td><strong>Embeddings</strong></td>
          <td>{n_embeddings} found</td>
          <td>{embeddings_html}</td>
        </tr>
        <tr>
          <td><strong>Layers</strong></td>
          <td>{n_layers} found</td>
          <td>{layers_html}</td>
        </tr>
        <tr>
          <td><strong>All obs columns</strong></td>
          <td>{n_obs_cols} total</td>
          <td>{obs_cols_html}</td>
        </tr>
        <tr>
          <td><strong>All var columns</strong></td>
          <td>{n_var_cols} total</td>
          <td>{var_cols_html}</td>
        </tr>
      </tbody>
    </table>
  </div>

  <!-- Cell type preview -->
  {label_preview_section}

  <!-- Footer -->
  <div class="footer">
    OmicSage Data Intake Report · Generated by pipeline/modules/qc/data_report.py<br>
    scanpy {scanpy_version} · numpy {numpy_version}
  </div>

</div>
</body>
</html>
"""


# ── Tag helpers ─────────────────────────────────────────────────────────────

def tags(items: list, cls: str = "") -> str:
    if not items:
        return '<span style="color:#9CA3AF;font-size:0.85rem">none detected</span>'
    return "".join(f'<span class="tag {cls}">{i}</span>' for i in items)


def _fmt_num(n: int) -> str:
    if n >= 1_000_000:
        return f"{n/1_000_000:.1f}M"
    if n >= 1_000:
        return f"{n/1_000:.1f}k"
    return str(n)


# ── Build HTML ───────────────────────────────────────────────────────────────

def build_html(input_path: Path, metrics: dict,
               plots: dict, geo: dict) -> str:

    # Modality badges
    badge_colors = {"RNA": "", "ATAC": "teal", "Protein (ADT)": "purple"}
    modality_badges = "".join(
        f'<span class="badge {badge_colors.get(m,"coral")}">{m}</span>'
        for m in metrics["modalities"]
    )

    # Title
    if geo:
        title_short = geo.get("title", input_path.stem)[:80]
    else:
        title_short = input_path.stem

    # Source section
    if geo:
        pubmed_html = ""
        if geo.get("pubmed_ids"):
            links = " · ".join(
                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" '
                f'target="_blank">PMID {p}</a>'
                for p in geo["pubmed_ids"]
            )
            pubmed_html = f'<p style="margin-top:0.75rem;font-size:0.85rem">'
            pubmed_html += f'<strong>Publications:</strong> {links}</p>'

        source_section = f"""
        <table style="margin-bottom:1rem">
          <tbody>
            <tr><td style="width:140px;color:var(--muted);padding:0.3rem 0.75rem 0.3rem 0">GEO Accession</td>
                <td><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo['gse']}"
                   target="_blank" class="mono">{geo['gse']}</a></td></tr>
            <tr><td style="color:var(--muted);padding:0.3rem 0.75rem 0.3rem 0">Organism</td>
                <td>{geo.get('organism','—')}</td></tr>
            <tr><td style="color:var(--muted);padding:0.3rem 0.75rem 0.3rem 0">Samples</td>
                <td>{geo.get('samples','—')}</td></tr>
          </tbody>
        </table>
        <div class="geo-summary">{geo.get('summary','No summary available.')}</div>
        {pubmed_html}
        """
    else:
        source_section = """
        <div class="personal-note">
          No GEO accession provided — this is a personal or private dataset.
          The report below is based entirely on data analysis.
          To add a public dataset description, re-run with <code>--geo &lt;accession&gt;</code>.
        </div>
        """

    # Cell type preview section
    if metrics["label_previews"]:
        rows = ""
        for col, vals in metrics["label_previews"].items():
            more = f" <span style='color:var(--muted)'>+more</span>" \
                   if len(vals) == 8 else ""
            rows += f"""
            <tr>
              <td class="mono">{col}</td>
              <td>{len(vals)}{"+" if len(vals)==8 else ""} unique values</td>
              <td>{tags(vals, 'green')}{more}</td>
            </tr>"""
        label_preview_section = f"""
        <div class="card">
          <h2>Cell Type / Cluster Annotations</h2>
          <table>
            <thead><tr><th>Column</th><th>Unique values</th><th>Preview</th></tr></thead>
            <tbody>{rows}</tbody>
          </table>
        </div>"""
    else:
        label_preview_section = ""

    # MT stat
    mt_stat = f"{metrics['median_mt_pct']}%" if metrics.get("has_mt") else "N/A"

    return HTML_TEMPLATE.format(
        title_short       = title_short,
        timestamp         = datetime.now().strftime("%Y-%m-%d %H:%M"),
        input_path        = str(input_path),
        modality_badges   = modality_badges,
        source_section    = source_section,
        n_cells           = _fmt_num(metrics["n_cells"]),
        n_genes           = _fmt_num(metrics["n_genes"]),
        median_genes      = f"{metrics['median_genes']:,.0f}",
        median_counts     = f"{metrics['median_counts']:,.0f}",
        sparsity          = metrics["sparsity_pct"],
        mt_stat           = mt_stat,
        plot_distributions= plots["distributions"],
        plot_scatter      = plots["scatter"],
        n_label_cols      = len(metrics["label_cols"]),
        n_batch_cols      = len(metrics["batch_cols"]),
        n_embeddings      = len(metrics["embeddings"]),
        n_layers          = len(metrics["layers"]),
        n_obs_cols        = len(metrics["obs_cols"]),
        n_var_cols        = len(metrics["var_cols"]),
        label_cols_html   = tags(metrics["label_cols"]),
        batch_cols_html   = tags(metrics["batch_cols"], "orange"),
        embeddings_html   = tags(metrics["embeddings"]),
        layers_html       = tags(metrics["layers"]),
        obs_cols_html     = tags(metrics["obs_cols"][:20]),
        var_cols_html     = tags(metrics["var_cols"][:20]),
        label_preview_section = label_preview_section,
        scanpy_version    = sc.__version__,
        numpy_version     = np.__version__,
    )


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate a data intake report for any .h5ad file."
    )
    parser.add_argument("--input",  required=True,
                        help="Path to .h5ad file")
    parser.add_argument("--output", required=True,
                        help="Path for output .html report")
    parser.add_argument("--geo",    default=None,
                        help="GEO accession (e.g. GSE194122) for public datasets")
    args = parser.parse_args()

    input_path  = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Loading {input_path.name} ...")
    adata = sc.read_h5ad(input_path)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} features")

    print("Analysing ...")
    metrics = analyse(adata)

    geo = {}
    if args.geo:
        print(f"Fetching GEO metadata for {args.geo} ...")
        geo = fetch_geo_metadata(args.geo)
        if geo:
            print(f"  Found: {geo['title'][:60]}...")
        else:
            print(f"  Could not fetch GEO metadata — continuing without it.")

    print("Generating plots ...")
    plots = make_plots(metrics)

    print("Building HTML report ...")
    html = build_html(input_path, metrics, plots, geo)

    output_path.write_text(html, encoding="utf-8")
    print(f"\nReport saved → {output_path}")
    print(f"Open in browser: file://{output_path.resolve()}")


if __name__ == "__main__":
    main()
