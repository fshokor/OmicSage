"""
OmicSage UI — Page 4: Report
==============================
Fixes:
  - Always resolve reports_dir to an absolute path before checking
  - Show the path being searched so the user can diagnose mismatches
  - Allow manual override of the reports directory
  - Scan common locations if the configured path is not found
"""
import logging
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from ui.state import *

logging.getLogger("streamlit").setLevel(logging.ERROR)

# Combined report filenames — each modality uses a different prefix
COMBINED_REPORT_NAMES = [
    "00_combined_report.html",           # scRNA-seq
    "multiome_00_combined_report.html",  # Multiome
    "cite_00_combined_report.html",      # CITE-seq
    "00_spatial_combined_report.html",   # Spatial (actual filename from pipeline)
]

STEP_REPORTS = [
    # ── scRNA-seq ──────────────────────────────────────────────────────────────
    ("00_data_report.html",            "Data Intake"),
    ("01_qc_report.html",              "QC"),
    ("02_normalization_report.html",   "Normalization"),
    ("03_reduce_report.html",          "Dimensionality Reduction"),
    ("04_cluster_report.html",         "Clustering"),
    ("05_annotate_report.html",        "Cell Type Annotation"),
    ("06_deg_report.html",             "Differential Expression"),
    ("07_gsea_report.html",            "Pathway Analysis"),
    ("08_harmony_report.html",         "Harmony Batch Correction"),
    ("10_pseudobulk_deg_report.html",  "Pseudobulk DEG"),
    # ── CITE-seq ───────────────────────────────────────────────────────────────
    ("cite_01_normalize_adt_report.html",    "ADT Normalization"),
    ("cite_02_doublets_report.html",         "ADT Doublets"),
    ("cite_03_reduce_adt_report.html",       "ADT Reduction"),
    ("cite_04_harmony_adt_report.html",      "ADT Harmony"),
    ("cite_05_annotate_adt_report.html",     "ADT Annotation"),
    ("cite_06_integration_report.html",      "RNA+ADT Integration"),
    ("cite_07_deg_report.html",              "CITE-seq DEG"),
    ("cite_08_gsea_report.html",             "CITE-seq GSEA"),
    # ── Multiome ───────────────────────────────────────────────────────────────
    ("multiome_01_qc_report.html",           "ATAC QC"),
    ("multiome_02_reduce_report.html",       "ATAC Reduction"),
    ("multiome_03_annotate_report.html",     "ATAC Annotation"),
    ("multiome_04_integration_report.html",  "Multiome Integration"),
    ("multiome_05_deg_report.html",          "Multiome DEG"),
    ("multiome_06_grn_report.html",          "GRN Inference"),
    # ── Spatial (no number prefix — actual filenames from run_spatial_pipeline.py)
    ("spatial_ingest_report.html",       "Spatial Ingest"),
    ("spatial_qc_report.html",           "Spatial QC"),
    ("spatial_reduce_report.html",       "Spatial Reduction"),
    ("spatial_cluster_report.html",      "Spatial Clustering"),
    ("spatial_deconvolve_report.html",   "Deconvolution"),
    ("spatial_downstream_report.html",   "Downstream"),
    ("spatial_impute_report.html",       "Spatial Imputation"),
]


def _scan_all_reports(rdir: Path) -> list[tuple[str, str]]:
    """
    Dynamically scan rdir for any HTML files that look like step reports.
    Used as a fallback when the static STEP_REPORTS list misses files.
    Returns list of (filename, label) sorted by filename.
    """
    known = {f for f, _ in STEP_REPORTS} | set(COMBINED_REPORT_NAMES)
    extra = []
    for f in sorted(rdir.glob("*.html")):
        if f.name not in known and "combined" not in f.name.lower():
            # Convert filename to readable label
            label = f.stem.replace("_", " ").replace("-", " ").title()
            extra.append((f.name, label))
    return extra


def _embed(html_path: Path, height: int):
    html = html_path.read_text(encoding="utf-8", errors="replace")
    components.html(html, height=height, scrolling=True)


def _resolve_dir(raw: str) -> Path | None:
    """
    Try to resolve a reports_dir path to an existing directory.
    Tries: as-is, relative to cwd, relative to repo root.
    Returns Path if found, None otherwise.
    """
    if not raw:
        return None
    candidates = [
        Path(raw),                              # absolute or relative to cwd
        Path.cwd() / raw,                       # explicit relative to cwd
        Path(__file__).resolve().parents[2] / raw,  # relative to repo root
    ]
    for c in candidates:
        if c.exists() and c.is_dir():
            return c.resolve()
    return None


def _scan_reports_root() -> list[Path]:
    """Scan the reports/ directory at repo root for any subdirectories."""
    repo_root = Path(__file__).resolve().parents[2]
    reports_root = repo_root / "reports"
    if reports_root.exists():
        return sorted(
            [d for d in reports_root.iterdir() if d.is_dir()],
            key=lambda p: p.stat().st_mtime,
            reverse=True,
        )
    return []


def render():
    init_state()
    dataset_name = st.session_state.get(KEY_DATASET_NAME, "Dataset")
    reports_dir_raw = st.session_state.get(KEY_REPORTS_DIR, "")

    st.markdown(f"## Report  ·  *{dataset_name}*")

    # ── Reports dir resolution ────────────────────────────────────────────────
    rdir = _resolve_dir(reports_dir_raw) if reports_dir_raw else None

    # If not found, show diagnostic + let user override or pick from list
    if rdir is None:
        st.warning("Reports directory not found.")

        if reports_dir_raw:
            st.caption(f"Configured path: `{reports_dir_raw}`")
            st.caption(f"Looked in: `{Path.cwd() / reports_dir_raw}`")

        # ── Manual override ───────────────────────────────────────────────────
        st.markdown("### Set reports directory manually")
        manual = st.text_input(
            "Full path to your reports folder",
            value=reports_dir_raw or "",
            placeholder="/home/shoko/OmicSage/reports/GSE166635",
            key="manual_reports_dir",
        )
        if manual:
            p = Path(manual)
            if p.exists() and p.is_dir():
                st.success(f"✓ Found: `{manual}`")
                st.session_state[KEY_REPORTS_DIR] = manual
                rdir = p.resolve()
            else:
                st.error("Directory not found — check the path.")

        # ── Auto-scan reports/ ────────────────────────────────────────────────
        existing = _scan_reports_root()
        if existing and rdir is None:
            st.divider()
            st.markdown("### Or pick from existing report folders")
            for d in existing[:10]:
                has_combined = (d / "00_combined_report.html").exists()
                any_reports  = any((d / f).exists() for f, _ in STEP_REPORTS)
                icon = "📊" if has_combined else ("📄" if any_reports else "📁")
                label = f"{icon} {d.name}"
                if st.button(label, key=f"pick_{d.name}", use_container_width=True):
                    st.session_state[KEY_REPORTS_DIR] = str(d)
                    st.rerun()

        col_back = st.columns([1, 3])[0]
        with col_back:
            if st.button("← Back to Run"):
                st.session_state[KEY_PAGE] = 2
                st.rerun()
        return

    # ── Found the directory — render reports ──────────────────────────────────
    # Try all known combined report filenames
    combined = None
    for fname in COMBINED_REPORT_NAMES:
        candidate = rdir / fname
        if candidate.exists():
            combined = candidate
            break

    # Show resolved path so user can verify it's the right one
    st.caption(f"Reports directory: `{rdir}`")

    # ── Nav + download ────────────────────────────────────────────────────────
    col_back, col_change, col_dl = st.columns([1, 1, 2])
    with col_back:
        if st.button("← Back to Run"):
            st.session_state[KEY_PAGE] = 2
            st.rerun()
    with col_change:
        if st.button("📂 Change reports folder"):
            st.session_state[KEY_REPORTS_DIR] = ""
            st.rerun()
    with col_dl:
        if combined and combined.exists():
            with open(combined, "rb") as f:
                st.download_button(
                    "⬇ Download combined report",
                    data=f.read(),
                    file_name=f"{dataset_name.replace(' ', '_')}_report.html",
                    mime="text/html",
                )

    st.divider()

    # ── Combined report inline ────────────────────────────────────────────────
    if combined and combined.exists():
        st.markdown("### Combined report")
        st.caption(f"File: `{combined.name}`")
        _embed(combined, height=850)
    else:
        tried = ", ".join(f"`{n}`" for n in COMBINED_REPORT_NAMES)
        st.warning(f"Combined report not found. Tried: {tried}")
        found = [lbl for f, lbl in STEP_REPORTS if (rdir / f).exists()]
        if found:
            st.info(
                f"Found **{len(found)}** individual step report(s): "
                f"{', '.join(found)}. Expand them below."
            )
        else:
            st.error(
                "No reports found in this directory at all. "
                "The pipeline may not have completed successfully. "
                "Check the Run log."
            )

    # ── Individual step reports ───────────────────────────────────────────────
    # Static list + dynamic scan for any reports not in the static list
    available = [(f, lbl) for f, lbl in STEP_REPORTS if (rdir / f).exists()]
    extra     = _scan_all_reports(rdir)
    available = available + [(f, lbl) for f, lbl in extra
                             if not any(f == af for af, _ in available)]
    if available:
        st.divider()
        st.markdown(f"### Individual step reports  ({len(available)} found)")
        for fname, label in available:
            fpath = rdir / fname
            with st.expander(f"📄 {label}  —  `{fname}`"):
                c1, c2 = st.columns([4, 1])
                with c1:
                    _embed(fpath, height=600)
                with c2:
                    with open(fpath, "rb") as f:
                        st.download_button(
                            "⬇ Download",
                            data=f.read(),
                            file_name=fname,
                            mime="text/html",
                            key=f"dl_{fname}",
                        )

    # ── AI outputs ────────────────────────────────────────────────────────────
    for ai_file, ai_label in [
        ("NEXT_STEPS.md",   "🤖 AI Recommended Next Steps"),
        ("ai_narrative.md", "🤖 AI Narrative"),
    ]:
        p = rdir / ai_file
        if p.exists():
            st.divider()
            st.markdown(f"### {ai_label}")
            st.markdown(p.read_text(encoding="utf-8"))
