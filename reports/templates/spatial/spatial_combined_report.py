"""
spatial_combined_report.py — OmicSage Phase 7
Tabbed combined report for the spatial transcriptomics pipeline.
Identical tab UI to the RNA combined_report.py.
"""

from __future__ import annotations

import re
import sys
from datetime import datetime
from pathlib import Path


TAB_REGISTRY = [
    ("spatial_qc_report.html",          "QC",         "🔬"),
    ("spatial_reduce_report.html",      "Reduce",     "🔭"),
    ("spatial_cluster_report.html",     "Cluster",    "🫧"),
    ("spatial_deconvolve_report.html",  "Deconvolve", "🧬"),
    ("spatial_downstream_report.html",  "Downstream", "🔗"),
    ("spatial_impute_report.html",      "Impute",     "🧩"),
]


def _extract_main(html):
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


def _extract_step_css(html):
    styles = re.findall(r"<style[^>]*>(.*?)</style>", html, re.DOTALL | re.IGNORECASE)
    return "\n".join(styles) if styles else ""


def generate_spatial_combined_report(
    reports_dir,
    dataset_name="OmicSage Spatial Analysis",
    output_path=None,
):
    reports_dir = Path(reports_dir)
    if output_path is None:
        output_path = reports_dir / "00_spatial_combined_report.html"
    output_path = Path(output_path)

    tabs = []
    for filename, label, icon in TAB_REGISTRY:
        # Try bare filename and dataset-prefixed variants
        candidates = list(reports_dir.glob(f"*{filename}")) + [reports_dir / filename]
        path = next((c for c in candidates if c.exists()), None)
        if path is None:
            continue
        html = path.read_text(encoding="utf-8", errors="replace")
        tabs.append({
            "id":      filename.replace(".", "_").replace("-", "_"),
            "label":   label,
            "icon":    icon,
            "content": _extract_main(html),
            "css":     _extract_step_css(html),
        })

    if not tabs:
        print(f"[spatial_combined_report] No step reports found in {reports_dir}")
        return ""

    html = _render_combined(tabs, dataset_name)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html, encoding="utf-8")
    print(f"[spatial_combined_report] {len(tabs)} tabs -> {output_path}")
    return str(output_path.resolve())


def _render_combined(tabs, dataset_name):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M")

    nav_buttons = ""
    for i, tab in enumerate(tabs):
        active = " active" if i == 0 else ""
        nav_buttons += (
            f'<button class="tab-btn{active}" data-tab="{tab["id"]}" '
            f'onclick="switchTab(\'{tab["id"]}\')">'
            f'<span class="tab-icon">{tab["icon"]}</span>'
            f'<span class="tab-label">{tab["label"]}</span>'
            f'</button>\n'
        )

    panels = ""
    for i, tab in enumerate(tabs):
        display = "block" if i == 0 else "none"
        panels += (
            f'<div class="tab-panel" id="panel_{tab["id"]}" style="display:{display};">\n'
            f'{tab["content"]}\n</div>\n'
        )

    step_styles = ""
    seen: set = set()
    for tab in tabs:
        css = tab["css"].strip()
        if css and css not in seen:
            step_styles += f"\n/* --- {tab['label']} --- */\n{css}\n"
            seen.add(css)

    n_total = len(TAB_REGISTRY)
    n_done  = len(tabs)
    pct     = int(100 * n_done / n_total)
    done_labels = " · ".join(t["label"] for t in tabs)

    # Build HTML as string concatenation to avoid nested f-string issues
    progress_fill_style = f"width:{pct}%"

    html = (
        "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n"
        "  <meta charset=\"UTF-8\">\n"
        "  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n"
        f"  <title>OmicSage &mdash; {dataset_name}</title>\n"
        "  <style>\n"
        "    *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }\n"
        "    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;\n"
        "           font-size: 14px; line-height: 1.6; color: #1a1a2e; background: #f7f8fc; min-height: 100vh; }\n"
        "    .report-header { background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);\n"
        "                     color: white; padding: 28px 40px 20px;\n"
        "                     position: sticky; top: 0; z-index: 100;\n"
        "                     box-shadow: 0 2px 12px rgba(0,0,0,0.25); }\n"
        "    .report-header h1 { font-size: 1.5rem; font-weight: 700; letter-spacing: -0.4px; }\n"
        "    .report-header h1 span { opacity: 0.55; font-weight: 400; }\n"
        "    .header-meta { font-size: 0.8rem; opacity: 0.65; margin-top: 4px; }\n"
        "    .pipeline-progress { margin-top: 12px; }\n"
        "    .progress-track { height: 4px; background: rgba(255,255,255,0.2);\n"
        "                      border-radius: 2px; overflow: hidden; margin-bottom: 4px; }\n"
        "    .progress-fill { height: 100%; background: linear-gradient(90deg, #54a868, #4C78A8);\n"
        f"                     border-radius: 2px; {progress_fill_style}; " + "}\n"
        "    .progress-label { font-size: 0.75rem; opacity: 0.7; }\n"
        "    .tab-nav { display: flex; overflow-x: auto; background: white;\n"
        "               border-bottom: 2px solid #e8eaf6;\n"
        "               position: sticky; top: 88px; z-index: 99;\n"
        "               box-shadow: 0 1px 4px rgba(0,0,0,0.06); }\n"
        "    .tab-btn { display: flex; align-items: center; gap: 6px;\n"
        "               padding: 12px 20px; border: none; background: none;\n"
        "               font-size: 0.85rem; color: #666; cursor: pointer;\n"
        "               white-space: nowrap; border-bottom: 3px solid transparent;\n"
        "               margin-bottom: -2px; transition: color 0.15s, border-color 0.15s; }\n"
        "    .tab-btn:hover { color: #0f3460; background: #f8f9ff; }\n"
        "    .tab-btn.active { color: #0f3460; border-bottom-color: #0f3460;\n"
        "                      font-weight: 700; background: #f8f9ff; }\n"
        "    .tab-icon { font-size: 1rem; line-height: 1; }\n"
        "    .tab-panel { max-width: 1400px; margin: 0 auto; padding: 28px 24px 48px; }\n"
        "    section { background: white; border-radius: 10px;\n"
        "              box-shadow: 0 1px 4px rgba(0,0,0,0.07);\n"
        "              padding: 28px 32px; margin-bottom: 24px; }\n"
        "    section h2 { font-size: 1.15rem; font-weight: 700; color: #0f3460;\n"
        "                 border-bottom: 2px solid #e8eaf6; padding-bottom: 10px; margin-bottom: 18px; }\n"
        "    section h3 { font-size: 1rem; font-weight: 600; color: #16213e; margin: 18px 0 10px; }\n"
        "    section p  { color: #444; margin-bottom: 12px; font-size: 0.9rem; }\n"
        "    .timestamp { font-size: 0.8rem; color: #888; margin-bottom: 6px; }\n"
        "    .note { font-size: 0.82rem; color: #7a5c00; background: #fffbe6;\n"
        "            border-left: 3px solid #f0c040; padding: 8px 12px;\n"
        "            border-radius: 4px; margin-bottom: 14px; }\n"
        "    code { font-family: 'SFMono-Regular', Consolas, monospace;\n"
        "           background: #f0f2ff; padding: 1px 5px; border-radius: 3px; font-size: 0.85em; }\n"
        "    .stat-grid { display: flex; flex-wrap: wrap; gap: 14px; margin-bottom: 24px; }\n"
        "    .stat-card { background: #f0f2ff; border-radius: 8px; padding: 14px 20px;\n"
        "                 min-width: 130px; text-align: center; flex: 1 1 130px; }\n"
        "    .stat-value { font-size: 1.4rem; font-weight: 700; color: #0f3460; }\n"
        "    .stat-label { font-size: 0.75rem; color: #666; margin-top: 2px; }\n"
        "    table { width: 100%; border-collapse: collapse; font-size: 0.88rem; margin-top: 8px; }\n"
        "    th { background: #f0f2ff; color: #0f3460; font-weight: 600;\n"
        "         padding: 9px 12px; text-align: left; border-bottom: 2px solid #d0d4f0; }\n"
        "    td { padding: 8px 12px; border-bottom: 1px solid #eee; vertical-align: middle; }\n"
        "    tr:last-child td { border-bottom: none; }\n"
        "    tr:hover td { background: #f8f9ff; }\n"
        "    .fig-grid { display: flex; flex-wrap: wrap; gap: 18px; margin-top: 12px; }\n"
        "    .fig-wrap { flex: 1 1 420px; max-width: 100%; }\n"
        "    .fig-wrap h3 { font-size: 0.9rem; margin-bottom: 6px; color: #16213e; }\n"
        "    .fig-wrap img { width: 100%; border-radius: 6px; border: 1px solid #e8eaf6;\n"
        "                    cursor: zoom-in; transition: box-shadow 0.15s;\n"
        "                    display: block; }\n"
        "    .fig-wrap img:hover { box-shadow: 0 4px 16px rgba(15,52,96,0.18); }\n"
        "    img { max-width: 100%; height: auto; }\n"
        "    /* Lightbox overlay */\n"
        "    #lb-overlay { display:none; position:fixed; inset:0; z-index:9999;\n"
        "                  background:rgba(0,0,0,0.85); align-items:center;\n"
        "                  justify-content:center; cursor:zoom-out; padding:24px; }\n"
        "    #lb-overlay.open { display:flex; }\n"
        "    #lb-img { max-width:95vw; max-height:92vh; border-radius:6px;\n"
        "              box-shadow:0 8px 48px rgba(0,0,0,0.6); object-fit:contain; }\n"
        "    #lb-close { position:absolute; top:16px; right:24px; color:#fff;\n"
        "                font-size:2rem; line-height:1; cursor:pointer;\n"
        "                background:none; border:none; opacity:0.75; }\n"
        "    #lb-close:hover { opacity:1; }\n"
        "    #lb-caption { position:absolute; bottom:16px; left:50%; transform:translateX(-50%);\n"
        "                  color:#ddd; font-size:0.82rem; text-align:center;\n"
        "                  background:rgba(0,0,0,0.5); padding:4px 12px; border-radius:4px;\n"
        "                  max-width:80vw; white-space:nowrap; overflow:hidden;\n"
        "                  text-overflow:ellipsis; }\n"
        "    .report-footer { text-align: center; font-size: 0.78rem; color: #aaa; padding: 24px 0 32px; }\n"
        "    .report-footer a { color: #0f3460; text-decoration: none; }\n"
        f"    {step_styles}\n"
        "  </style>\n"
        "</head>\n<body>\n"
        "  <header class=\"report-header\">\n"
        f"    <h1>OmicSage <span>&mdash; {dataset_name}</span></h1>\n"
        f"    <p class=\"header-meta\">Generated {timestamp} &middot; {n_done} of {n_total} pipeline steps complete</p>\n"
        "    <div class=\"pipeline-progress\">\n"
        "      <div class=\"progress-track\"><div class=\"progress-fill\"></div></div>\n"
        f"      <span class=\"progress-label\">{pct}% &middot; {done_labels}</span>\n"
        "    </div>\n"
        "  </header>\n"
        "  <nav class=\"tab-nav\" role=\"tablist\" aria-label=\"Spatial pipeline steps\">\n"
        f"    {nav_buttons}"
        "  </nav>\n"
        f"  {panels}\n"
        "  <div id=\"lb-overlay\" role=\"dialog\" aria-modal=\"true\" aria-label=\"Image viewer\"\n"
        "       onclick=\"closeLightbox(event)\">\n"
        "    <button id=\"lb-close\" onclick=\"closeLightbox()\" aria-label=\"Close\">&times;</button>\n"
        "    <img id=\"lb-img\" src=\"\" alt=\"\">\n"
        "    <div id=\"lb-caption\"></div>\n"
        "  </div>\n"
        "  <footer class=\"report-footer\">\n"
        "    Generated by <a href=\"https://github.com/fshokor/OmicSage\">OmicSage</a> &middot; MIT License\n"
        "  </footer>\n"
        "  <script>\n"
        "    function switchTab(tabId) {\n"
        "      document.querySelectorAll('.tab-panel').forEach(p => p.style.display = 'none');\n"
        "      document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));\n"
        "      document.getElementById('panel_' + tabId).style.display = 'block';\n"
        "      document.querySelector('[data-tab=\"' + tabId + '\"]').classList.add('active');\n"
        "      document.querySelector('[data-tab=\"' + tabId + '\"]').scrollIntoView({behavior:'smooth',block:'nearest',inline:'center'});\n"
        "    }\n"
        "    document.addEventListener('keydown', function(e) {\n"
        "      const btns = Array.from(document.querySelectorAll('.tab-btn'));\n"
        "      const active = btns.findIndex(b => b.classList.contains('active'));\n"
        "      if (e.key === 'ArrowRight' && active < btns.length - 1) btns[active+1].click();\n"
        "      else if (e.key === 'ArrowLeft' && active > 0) btns[active-1].click();\n"
        "      else if (e.key === 'Escape') closeLightbox();\n"
        "    });\n"
        "    // Lightbox\n"
        "    function openLightbox(img) {\n"
        "      const ov = document.getElementById('lb-overlay');\n"
        "      document.getElementById('lb-img').src = img.src;\n"
        "      const cap = img.closest('.fig-wrap');\n"
        "      const h3  = cap ? cap.querySelector('h3') : null;\n"
        "      document.getElementById('lb-caption').textContent = h3 ? h3.textContent : '';\n"
        "      ov.classList.add('open');\n"
        "      document.body.style.overflow = 'hidden';\n"
        "    }\n"
        "    function closeLightbox(e) {\n"
        "      if (e && e.target === document.getElementById('lb-img')) return;\n"
        "      document.getElementById('lb-overlay').classList.remove('open');\n"
        "      document.body.style.overflow = '';\n"
        "    }\n"
        "    // Wire up all current + future images via delegation\n"
        "    document.addEventListener('click', function(e) {\n"
        "      const img = e.target.closest('.fig-wrap img');\n"
        "      if (img) { e.preventDefault(); openLightbox(img); }\n"
        "    });\n"
        "  </script>\n"
        "</body>\n</html>"
    )
    return html


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--reports-dir", required=True)
    parser.add_argument("--dataset-name", default="OmicSage Spatial Analysis")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    result = generate_spatial_combined_report(
        reports_dir=Path(args.reports_dir),
        dataset_name=args.dataset_name,
        output_path=Path(args.output) if args.output else None,
    )
    if result:
        print(f"Combined report: {result}")
    else:
        sys.exit(1)
