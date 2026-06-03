"""
cite_combined_report.py — OmicSage CITE-seq tabbed pipeline report
===================================================================
Reads all individual CITE-seq step HTML reports from reports_dir and
assembles them into a single self-contained tabbed HTML file.

Only tabs for reports that actually exist on disk are shown — partial
runs produce a partial combined report without errors. This means you
can stop the pipeline at any step and still get a valid combined report
covering everything that ran.

Output: cite_00_combined_report.html
(sorted to the top of the reports directory due to the "cite_00_" prefix)

Usage — called automatically by run_cite_pipeline.py at the end of a run:

    from reports.templates.cite.cite_combined_report import generate_cite_combined_report

    generate_cite_combined_report(
        reports_dir=Path("reports/GSE194122"),
        dataset_name="BMMC CITE-seq (NeurIPS 2021)",
    )

Can also be run standalone to rebuild from existing reports:

    python -m reports.templates.cite.cite_combined_report \\
        --reports-dir reports/GSE194122 \\
        --dataset-name "BMMC CITE-seq (NeurIPS 2021)"

    # Custom output path:
    python -m reports.templates.cite.cite_combined_report \\
        --reports-dir reports/GSE194122 \\
        --output reports/GSE194122/cite_combined.html
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
    ("cite_01_normalize_report.html",    "Normalize ADT", "🧪"),
    ("cite_02_doublets_report.html",     "Doublets",      "🔍"),
    ("cite_03_reduce_report.html",       "Reduce ADT",    "🔭"),
    ("cite_04_harmony_report.html",      "Harmony ADT",   "🎵"),
    ("cite_05_annotate_report.html",     "Annotate ADT",  "🏷️"),
    ("cite_06_integration_report.html",  "Integration",   "🔗"),
    ("cite_07_deg_report.html",          "DEG / DPE",     "📊"),
    ("cite_08_gsea_report.html",         "GSEA",          "🧬"),
]


# ---------------------------------------------------------------------------
# HTML extraction helpers — identical to combined_report.py
# ---------------------------------------------------------------------------

def _extract_main(html: str) -> str:
    """
    Pull the content of <main>...</main> from a step report.
    Falls back to <body> minus header/footer if no <main> tag found.
    """
    m = re.search(r"<main[^>]*>(.*?)</main>", html, re.DOTALL | re.IGNORECASE)
    if m:
        return m.group(1).strip()
    body = re.search(r"<body[^>]*>(.*?)</body>", html, re.DOTALL | re.IGNORECASE)
    if body:
        content = body.group(1)
        content = re.sub(r"<header[^>]*>.*?</header>", "", content,
                         flags=re.DOTALL | re.IGNORECASE)
        content = re.sub(r"<footer[^>]*>.*?</footer>", "", content,
                         flags=re.DOTALL | re.IGNORECASE)
        return content.strip()
    return html


def _extract_step_css(html: str) -> str:
    """
    Extract <style> blocks from a step report so step-specific CSS is
    preserved in the combined report. Keeps everything — CSS cascade
    handles overlap cleanly.
    """
    styles = re.findall(r"<style[^>]*>(.*?)</style>", html,
                        re.DOTALL | re.IGNORECASE)
    return "\n".join(styles) if styles else ""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def generate_cite_combined_report(
    reports_dir: Path,
    dataset_name: str = "CITE-seq Analysis",
    output_path: Path | None = None,
) -> str:
    """
    Assemble all available CITE-seq step reports into one tabbed HTML file.

    Parameters
    ----------
    reports_dir : Path
        Directory containing the individual step HTML reports.
    dataset_name : str
        Human-readable dataset name shown in the report header.
    output_path : Path, optional
        Where to write the combined report.
        Defaults to reports_dir / "cite_00_combined_report.html".

    Returns
    -------
    str
        Absolute path to the written combined report.
        Returns "" if no step reports were found.
    """
    reports_dir = Path(reports_dir)
    if output_path is None:
        output_path = reports_dir / "cite_00_combined_report.html"
    output_path = Path(output_path)

    # ------------------------------------------------------------------
    # Collect available tabs
    # ------------------------------------------------------------------
    tabs = []
    for filename, label, icon in TAB_REGISTRY:
        path = reports_dir / filename
        if path.exists():
            html         = path.read_text(encoding="utf-8", errors="replace")
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
        print(
            f"[cite_combined_report] No CITE-seq step reports found in {reports_dir} "
            f"— nothing to combine."
        )
        return ""

    # ------------------------------------------------------------------
    # Render and write
    # ------------------------------------------------------------------
    html = _render_combined(tabs, dataset_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"[cite_combined_report] {len(tabs)} tabs → {output_path}")
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

    # ── Step-specific CSS ─────────────────────────────────────────────
    step_styles = ""
    seen_css = set()
    for tab in tabs:
        css = tab["css"].strip()
        if css and css not in seen_css:
            step_styles += f"\n/* --- {tab['label']} --- */\n{css}\n"
            seen_css.add(css)

    # ── Progress bar ──────────────────────────────────────────────────
    n_steps_total    = len(TAB_REGISTRY)
    n_steps_done     = len(tabs)
    pct              = int(100 * n_steps_done / n_steps_total)
    step_labels_done = " · ".join(t["label"] for t in tabs)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>OmicSage — CITE-seq — {dataset_name}</title>
  <style>
    /* ── Reset & base ─────────────────────────────────────────────── */
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      font-size: 14px; line-height: 1.6; color: #1a1a2e;
      background: #f7f8fc; min-height: 100vh;
    }}

    /* ── Top header ───────────────────────────────────────────────── */
    .report-header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: white; padding: 28px 40px 20px;
      position: sticky; top: 0; z-index: 100;
      box-shadow: 0 2px 12px rgba(0,0,0,0.25);
    }}
    .report-header h1 {{
      font-size: 1.5rem; font-weight: 700; letter-spacing: -0.4px;
    }}
    .report-header h1 span {{ opacity: 0.55; }}
    .header-meta {{ font-size: 0.8rem; opacity: 0.65; margin-top: 4px; }}

    /* ── Progress bar ─────────────────────────────────────────────── */
    .pipeline-progress {{ margin-top: 10px; }}
    .progress-track {{
      height: 4px; background: rgba(255,255,255,0.2);
      border-radius: 2px; overflow: hidden; margin-bottom: 4px;
    }}
    .progress-fill {{
      height: 100%; width: {pct}%;
      background: linear-gradient(90deg, #54a868, #4C78A8);
      border-radius: 2px; transition: width 0.4s ease;
    }}
    .progress-label {{ font-size: 0.75rem; opacity: 0.65; }}

    /* ── Tab navigation ───────────────────────────────────────────── */
    .tab-nav {{
      display: flex; flex-wrap: nowrap; overflow-x: auto;
      background: white; border-bottom: 2px solid #e8eaf6;
      padding: 0 16px; gap: 2px;
      position: sticky; top: 88px; z-index: 99;
      box-shadow: 0 1px 4px rgba(0,0,0,0.06);
    }}
    .tab-nav::-webkit-scrollbar {{ height: 3px; }}
    .tab-nav::-webkit-scrollbar-thumb {{
      background: #d0d4f0; border-radius: 2px;
    }}
    .tab-btn {{
      display: flex; align-items: center; gap: 6px;
      padding: 12px 18px; border: none; background: none;
      font-size: 0.85rem; color: #666; cursor: pointer;
      border-bottom: 3px solid transparent; margin-bottom: -2px;
      white-space: nowrap; transition: color 0.15s, border-color 0.15s;
      flex-shrink: 0;
    }}
    .tab-btn:hover {{ color: #0f3460; background: #f8f9ff; }}
    .tab-btn.active {{
      color: #0f3460; border-bottom-color: #0f3460;
      font-weight: 700; background: #f8f9ff;
    }}
    .tab-icon {{ font-size: 1rem; line-height: 1; }}

    /* ── Tab panels ───────────────────────────────────────────────── */
    .tab-panel {{
      max-width: 1100px; margin: 0 auto; padding: 28px 24px 48px;
    }}

    /* ── Shared section styles ────────────────────────────────────── */
    section {{
      background: white; border-radius: 10px;
      box-shadow: 0 1px 4px rgba(0,0,0,0.07);
      padding: 28px 32px; margin-bottom: 24px;
    }}
    section h2 {{
      font-size: 1.15rem; font-weight: 700; color: #0f3460;
      border-bottom: 2px solid #e8eaf6;
      padding-bottom: 10px; margin-bottom: 18px;
    }}
    section h3 {{
      font-size: 1rem; font-weight: 600; color: #16213e;
      margin: 18px 0 10px;
    }}
    section p {{ color: #444; margin-bottom: 12px; font-size: 0.9rem; }}
    .timestamp {{ font-size: 0.8rem; color: #888; margin-bottom: 6px; }}
    .note {{
      font-size: 0.82rem; color: #7a5c00; background: #fffbe6;
      border-left: 3px solid #f0c040; padding: 8px 12px;
      border-radius: 4px; margin-bottom: 14px;
    }}
    code {{
      font-family: "SFMono-Regular", Consolas, monospace;
      background: #f0f2ff; padding: 1px 5px;
      border-radius: 3px; font-size: 0.85em;
    }}

    /* ── Stat cards ───────────────────────────────────────────────── */
    .stat-grid {{ display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }}
    .stat-card {{
      background: #f0f2ff; border-radius: 8px; padding: 14px 20px;
      min-width: 130px; text-align: center; flex: 1 1 130px;
    }}
    .stat-value {{ font-size: 1.4rem; font-weight: 700; color: #0f3460; }}
    .stat-label {{ font-size: 0.75rem; color: #666; margin-top: 2px; }}

    /* ── Tables ───────────────────────────────────────────────────── */
    table {{ width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }}
    th {{
      background: #f0f2ff; color: #0f3460; font-weight: 600;
      padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0;
    }}
    td {{ padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }}
    tr:last-child td {{ border-bottom: none; }}
    tr:hover td {{ background: #f8f9ff; }}

    /* ── Figure grids ─────────────────────────────────────────────── */
    .fig-grid {{ display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }}
    .fig-wrap {{ flex: 1 1 300px; max-width: 520px; }}
    .fig-wrap.wide {{ flex: 1 1 100%; max-width: 100%; }}
    .fig-wrap h3 {{ font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }}
    .fig-wrap img {{ width: 100%; border-radius: 6px; border: 1px solid #e8eaf6; }}
    img {{ max-width: 100%; height: auto; }}

    /* ── Footer ───────────────────────────────────────────────────── */
    .report-footer {{
      text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px;
    }}
    .report-footer a {{ color: #0f3460; text-decoration: none; }}

    /* ── Step-specific styles from individual reports ─────────────── */
    {step_styles}
  </style>
</head>
<body>

  <!-- ── Header ──────────────────────────────────────────────────────── -->
  <header class="report-header">
    <h1>OmicSage <span>— CITE-seq — {dataset_name}</span></h1>
    <p class="header-meta">
      Generated {timestamp}
      &middot; {n_steps_done} of {n_steps_total} pipeline steps complete
    </p>
    <div class="pipeline-progress">
      <div class="progress-track"><div class="progress-fill"></div></div>
      <span class="progress-label">{pct}% &middot; {step_labels_done}</span>
    </div>
  </header>

  <!-- ── Tab navigation ──────────────────────────────────────────────── -->
  <nav class="tab-nav" role="tablist" aria-label="CITE-seq pipeline steps">
    {nav_buttons}
  </nav>

  <!-- ── Tab panels ──────────────────────────────────────────────────── -->
  {panels}

  <!-- ── Footer ──────────────────────────────────────────────────────── -->
  <footer class="report-footer">
    Generated by <a href="https://github.com/fshokor/OmicSage">OmicSage</a>
    &middot; MIT License
  </footer>

  <!-- ── Tab switching JS ────────────────────────────────────────────── -->
  <script>
    function switchTab(tabId) {{
      document.querySelectorAll('.tab-panel').forEach(p => p.style.display = 'none');
      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
      document.getElementById('panel_' + tabId).style.display = 'block';
      document.querySelector('[data-tab="' + tabId + '"]').classList.add('active');
      document.querySelector('[data-tab="' + tabId + '"]').scrollIntoView({{
        behavior: 'smooth', block: 'nearest', inline: 'center'
      }});
    }}

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
        description="Rebuild the OmicSage CITE-seq combined report from existing step reports."
    )
    parser.add_argument(
        "--reports-dir", required=True,
        help="Path to the reports directory (e.g. reports/GSE194122)",
    )
    parser.add_argument(
        "--dataset-name", default="CITE-seq Analysis",
        help="Human-readable dataset name for the report header",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output path (default: <reports-dir>/cite_00_combined_report.html)",
    )
    args = parser.parse_args()

    result = generate_cite_combined_report(
        reports_dir=Path(args.reports_dir),
        dataset_name=args.dataset_name,
        output_path=Path(args.output) if args.output else None,
    )
    if result:
        print(f"Combined report written to: {result}")
    else:
        sys.exit(1)
