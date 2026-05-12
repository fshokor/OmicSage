"""
combined_report.py — OmicSage tabbed pipeline report generator
==============================================================
Reads all individual step HTML reports from reports_dir and assembles
them into a single self-contained HTML file with a tab per step.

Only tabs for reports that actually exist on disk are shown.
Zero changes required to existing report generators.

Usage (called automatically by run_pipeline.py at the end of a run):
    from reports.combined_report import generate_combined_report

    generate_combined_report(
        reports_dir=Path("data/GSE166635/reports"),
        dataset_name="GSE166635 — HCC",
        output_path=Path("data/GSE166635/reports/00_combined_report.html"),
    )

Can also be run standalone to rebuild from existing reports:
    python -m reports.combined_report --reports-dir data/GSE166635/reports
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Tab registry — order and display names
# ---------------------------------------------------------------------------

TAB_REGISTRY = [
    ("01_qc_report.html",             "QC",          "🔬"),
    ("02_normalization_report.html",   "Normalize",   "📐"),
    ("03_reduce_report.html",          "Reduce",      "🔭"),
    ("04_cluster_report.html",         "Cluster",     "🫧"),
    ("05_annotate_report.html",        "Annotate",    "🏷️"),
    ("06_deg_report.html",             "DEG",         "📊"),
    ("07_gsea_report.html",            "GSEA",        "🧬"),
    ("08_harmony_report.html",         "Harmony",     "🎵"),
    ("10_pseudobulk_deg_report.html",  "Pseudobulk",  "🔢"),
]


# ---------------------------------------------------------------------------
# HTML extraction helpers
# ---------------------------------------------------------------------------

def _extract_main(html: str) -> str:
    """
    Pull the content of <main>...</main> from a step report.
    Falls back to the full <body> if no <main> tag found.
    """
    m = re.search(r"<main[^>]*>(.*?)</main>", html, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1).strip()
    # Fallback: strip head + header + footer, return body content
    body = re.search(r"<body[^>]*>(.*?)</body>", html, re.DOTALL | re.IGNORECASE)
    if body:
        content = body.group(1)
        # Remove header and footer blocks
        content = re.sub(r"<header[^>]*>.*?</header>", "", content, flags=re.DOTALL | re.IGNORECASE)
        content = re.sub(r"<footer[^>]*>.*?</footer>", "", content, flags=re.DOTALL | re.IGNORECASE)
        return content.strip()
    return html


def _extract_step_css(html: str) -> str:
    """
    Extract any <style> blocks from a step report so step-specific
    CSS (volcano grid, dot plot sizing, etc.) is preserved.
    Strips the shared base styles to avoid duplication — keeps only
    rules that reference step-specific class names.
    """
    styles = re.findall(r"<style[^>]*>(.*?)</style>", html, re.DOTALL | re.IGNORECASE)
    if not styles:
        return ""
    combined = "\n".join(styles)
    # Remove rules that are already in the combined report base CSS
    # (identified by the shared class names defined in _render_combined)
    shared_selectors = {
        "body", "header", "main", "footer", "section", "table",
        "th", "td", "tr", "code", "h1", "h2", "h3",
        ".stat-grid", ".stat-card", ".stat-value", ".stat-label",
        ".timestamp", ".note", ".pos-fc", ".neg-fc",
    }
    # Keep everything — de-duplication via CSS cascade handles overlap cleanly
    return combined


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_combined_report(
    reports_dir: Path,
    dataset_name: str = "OmicSage Analysis",
    output_path: Path | None = None,
) -> str:
    """
    Assemble all available step reports into one tabbed HTML file.

    Parameters
    ----------
    reports_dir : Path
        Directory containing the individual step HTML reports.
    dataset_name : str
        Human-readable dataset name shown in the report header.
    output_path : Path, optional
        Where to write the combined report.
        Defaults to reports_dir / "00_combined_report.html".

    Returns
    -------
    str
        Absolute path to the written combined report.
    """
    reports_dir = Path(reports_dir)
    if output_path is None:
        output_path = reports_dir / "00_combined_report.html"
    output_path = Path(output_path)

    # ------------------------------------------------------------------
    # Collect available tabs
    # ------------------------------------------------------------------
    tabs = []
    for filename, label, icon in TAB_REGISTRY:
        path = reports_dir / filename
        if path.exists():
            html = path.read_text(encoding="utf-8", errors="replace")
            main_content = _extract_main(html)
            step_css     = _extract_step_css(html)
            tabs.append({
                "id":      filename.replace(".", "_").replace("-", "_"),
                "label":   label,
                "icon":    icon,
                "content": main_content,
                "css":     step_css,
            })

    if not tabs:
        print(f"[combined_report] No step reports found in {reports_dir} — nothing to combine.")
        return ""

    # ------------------------------------------------------------------
    # Render
    # ------------------------------------------------------------------
    html = _render_combined(tabs, dataset_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"[combined_report] {len(tabs)} tabs → {output_path}")
    return str(output_path.resolve())


# ---------------------------------------------------------------------------
# HTML renderer
# ---------------------------------------------------------------------------

def _render_combined(tabs: list[dict], dataset_name: str) -> str:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── Tab nav buttons ───────────────────────────────────────────────
    nav_buttons = ""
    for i, tab in enumerate(tabs):
        active_class = " active" if i == 0 else ""
        nav_buttons += (
            f'<button class="tab-btn{active_class}" '
            f'data-tab="{tab["id"]}" '
            f'onclick="switchTab(\'{tab["id"]}\')">'
            f'<span class="tab-icon">{tab["icon"]}</span>'
            f'<span class="tab-label">{tab["label"]}</span>'
            f'</button>\n'
        )

    # ── Tab content panels ────────────────────────────────────────────
    panels = ""
    for i, tab in enumerate(tabs):
        display = "block" if i == 0 else "none"
        panels += (
            f'<div class="tab-panel" id="panel_{tab["id"]}" '
            f'style="display:{display};">\n'
            f'{tab["content"]}\n'
            f'</div>\n'
        )

    # ── Step-specific CSS (collected from all reports) ────────────────
    step_styles = ""
    seen_css = set()
    for tab in tabs:
        css = tab["css"].strip()
        if css and css not in seen_css:
            step_styles += f"\n/* --- {tab['label']} --- */\n{css}\n"
            seen_css.add(css)

    # ── Progress bar ─────────────────────────────────────────────────
    n_steps_total = len(TAB_REGISTRY)
    n_steps_done  = len(tabs)
    pct           = int(100 * n_steps_done / n_steps_total)
    step_labels_done = " · ".join(t["label"] for t in tabs)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — {dataset_name}</title>
  <style>
    /* ── Reset & base ─────────────────────────────────────────────── */
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}

    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      font-size: 14px;
      line-height: 1.6;
      color: #1a1a2e;
      background: #f7f8fc;
      min-height: 100vh;
    }}

    /* ── Top header ───────────────────────────────────────────────── */
    .report-header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white;
      padding: 28px 40px 20px;
      position: sticky;
      top: 0;
      z-index: 100;
      box-shadow: 0 2px 12px rgba(0,0,0,0.25);
    }}
    .report-header h1 {{
      font-size: 1.5rem;
      font-weight: 700;
      letter-spacing: -0.4px;
    }}
    .report-header h1 span {{
      opacity: 0.55;
      font-weight: 400;
      font-size: 1.1rem;
      margin-left: 10px;
    }}
    .header-meta {{
      font-size: 0.78rem;
      opacity: 0.55;
      margin-top: 4px;
    }}

    /* ── Progress bar ─────────────────────────────────────────────── */
    .pipeline-progress {{
      margin-top: 14px;
      display: flex;
      align-items: center;
      gap: 12px;
    }}
    .progress-track {{
      flex: 1;
      height: 5px;
      background: rgba(255,255,255,0.15);
      border-radius: 3px;
      overflow: hidden;
    }}
    .progress-fill {{
      height: 100%;
      width: {pct}%;
      background: linear-gradient(90deg, #4fc3f7, #81d4fa);
      border-radius: 3px;
      transition: width 0.4s ease;
    }}
    .progress-label {{
      font-size: 0.72rem;
      opacity: 0.6;
      white-space: nowrap;
    }}

    /* ── Tab navigation ───────────────────────────────────────────── */
    .tab-nav {{
      background: #ffffff;
      border-bottom: 2px solid #e8eaf6;
      padding: 0 24px;
      display: flex;
      gap: 2px;
      overflow-x: auto;
      scrollbar-width: none;
      position: sticky;
      top: 97px;   /* below header */
      z-index: 90;
      box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    }}
    .tab-nav::-webkit-scrollbar {{ display: none; }}

    .tab-btn {{
      display: flex;
      align-items: center;
      gap: 6px;
      padding: 12px 18px;
      border: none;
      background: transparent;
      color: #666;
      font-size: 0.85rem;
      font-weight: 500;
      font-family: inherit;
      cursor: pointer;
      border-bottom: 3px solid transparent;
      margin-bottom: -2px;
      white-space: nowrap;
      transition: color 0.15s, border-color 0.15s, background 0.15s;
      border-radius: 6px 6px 0 0;
    }}
    .tab-btn:hover {{
      color: #0f3460;
      background: #f0f2ff;
    }}
    .tab-btn.active {{
      color: #0f3460;
      border-bottom-color: #0f3460;
      font-weight: 700;
      background: #f8f9ff;
    }}
    .tab-icon {{
      font-size: 1rem;
      line-height: 1;
    }}

    /* ── Tab panels ───────────────────────────────────────────────── */
    .tab-panel {{
      max-width: 1100px;
      margin: 0 auto;
      padding: 28px 24px 48px;
    }}

    /* ── Shared section styles (mirrors individual reports) ───────── */
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
    .note {{
      font-size: 0.82rem;
      color: #7a5c00;
      background: #fffbe6;
      border-left: 3px solid #f0c040;
      padding: 8px 12px;
      border-radius: 4px;
      margin-bottom: 14px;
    }}
    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 0.85em;
    }}

    /* ── Stat cards ───────────────────────────────────────────────── */
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

    /* ── Tables ───────────────────────────────────────────────────── */
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.88rem;
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
    .pos-fc {{ color: #c0392b; font-weight: 600; }}
    .neg-fc {{ color: #2980b9; font-weight: 600; }}

    /* ── Volcano grid ─────────────────────────────────────────────── */
    .volcano-grid {{
      display: flex;
      flex-wrap: wrap;
      gap: 18px;
      margin-top: 12px;
    }}
    .volcano-wrap {{
      flex: 1 1 300px;
      max-width: 420px;
    }}
    .volcano-wrap h3 {{
      font-size: 0.9rem;
      margin-bottom: 6px;
      color: #16213e;
    }}
    .volcano-wrap img {{
      width: 100%;
      border-radius: 6px;
      border: 1px solid #e8eaf6;
    }}

    /* ── Images (general) ─────────────────────────────────────────── */
    img {{ max-width: 100%; height: auto; }}

    /* ── Footer ───────────────────────────────────────────────────── */
    .report-footer {{
      text-align: center;
      font-size: 0.78rem;
      color: #aaa;
      padding: 24px 0 32px;
    }}
    .report-footer a {{ color: #0f3460; text-decoration: none; }}

    /* ── Step-specific styles from individual reports ─────────────── */
    {step_styles}
  </style>
</head>
<body>

  <!-- ── Header ──────────────────────────────────────────────────────── -->
  <header class="report-header">
    <h1>OmicSage <span>— {dataset_name}</span></h1>
    <p class="header-meta">Generated {timestamp} · {n_steps_done} of {n_steps_total} pipeline steps complete</p>
    <div class="pipeline-progress">
      <div class="progress-track"><div class="progress-fill"></div></div>
      <span class="progress-label">{pct}% · {step_labels_done}</span>
    </div>
  </header>

  <!-- ── Tab navigation ──────────────────────────────────────────────── -->
  <nav class="tab-nav" role="tablist" aria-label="Pipeline steps">
    {nav_buttons}
  </nav>

  <!-- ── Tab panels ──────────────────────────────────────────────────── -->
  {panels}

  <!-- ── Footer ──────────────────────────────────────────────────────── -->
  <footer class="report-footer">
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    · MIT License
  </footer>

  <!-- ── Tab switching JS ────────────────────────────────────────────── -->
  <script>
    function switchTab(tabId) {{
      // Hide all panels
      document.querySelectorAll('.tab-panel').forEach(p => p.style.display = 'none');
      // Deactivate all buttons
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      // Show selected panel
      document.getElementById('panel_' + tabId).style.display = 'block';
      // Activate selected button
      document.querySelector('[data-tab="' + tabId + '"]').classList.add('active');
      // Scroll tab into view if it's off-screen (mobile / many tabs)
      document.querySelector('[data-tab="' + tabId + '"]').scrollIntoView({{
        behavior: 'smooth', block: 'nearest', inline: 'center'
      }});
    }}

    // Keyboard navigation (left/right arrow keys)
    document.addEventListener('keydown', function(e) {{
      const btns = Array.from(document.querySelectorAll('.tab-btn'));
      const active = btns.findIndex(b => b.classList.contains('active'));
      if (e.key === 'ArrowRight' && active < btns.length - 1) {{
        btns[active + 1].click();
      }} else if (e.key === 'ArrowLeft' && active > 0) {{
        btns[active - 1].click();
      }}
    }});
  </script>

</body>
</html>
"""


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Rebuild the OmicSage combined report from existing step reports."
    )
    parser.add_argument(
        "--reports-dir", required=True,
        help="Path to the reports directory (e.g. data/GSE166635/reports)"
    )
    parser.add_argument(
        "--dataset-name", default="OmicSage Analysis",
        help="Human-readable dataset name for the report header"
    )
    parser.add_argument(
        "--output", default=None,
        help="Output path (default: <reports-dir>/00_combined_report.html)"
    )
    args = parser.parse_args()

    result = generate_combined_report(
        reports_dir=Path(args.reports_dir),
        dataset_name=args.dataset_name,
        output_path=Path(args.output) if args.output else None,
    )
    if result:
        print(f"Combined report written to: {result}")
    else:
        sys.exit(1)
