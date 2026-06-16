"""
OmicSage UI — Page 4: Report
==============================
Embeds HTML reports inline using st.components.v1.html() which is still
the correct way to render arbitrary HTML in Streamlit as of 2026-06.

The deprecation warning referred to an OLD signature of components.v1.html
that accepted a URL. The current correct usage — passing an HTML string —
is still fully supported and is what we use here.

st.iframe() only accepts URLs (http/https), not file:// or raw HTML strings,
so it cannot be used to render local report files.
"""
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components
from ui.state import *

STEP_REPORTS = [
    ("00_data_report.html",           "Data Intake"),
    ("01_qc_report.html",             "QC"),
    ("02_normalization_report.html",  "Normalization"),
    ("03_reduce_report.html",         "Dimensionality Reduction"),
    ("04_cluster_report.html",        "Clustering"),
    ("05_annotate_report.html",       "Cell Type Annotation"),
    ("06_deg_report.html",            "Differential Expression"),
    ("07_gsea_report.html",           "Pathway Analysis"),
    ("08_harmony_report.html",        "Harmony Batch Correction"),
    ("10_pseudobulk_deg_report.html", "Pseudobulk DEG"),
    # CITE-seq
    ("cite_01_normalize_adt_report.html",    "ADT Normalization"),
    ("cite_03_reduce_adt_report.html",       "ADT Reduction"),
    ("cite_04_harmony_adt_report.html",      "ADT Harmony"),
    ("cite_05_annotate_adt_report.html",     "ADT Annotation"),
    ("cite_06_integration_report.html",      "RNA+ADT Integration"),
    ("cite_07_deg_report.html",              "CITE-seq DEG"),
    # Multiome
    ("multiome_01_atac_qc_report.html",      "ATAC QC"),
    ("multiome_02_atac_reduce_report.html",  "ATAC Reduction"),
    ("multiome_03_atac_annotate_report.html","ATAC Annotation"),
    ("multiome_04_integration_report.html",  "Multiome Integration"),
    # Spatial
    ("spatial_02_qc_report.html",            "Spatial QC"),
    ("spatial_03_reduce_report.html",        "Spatial Reduction"),
    ("spatial_04_cluster_report.html",       "Spatial Clustering"),
    ("spatial_05_deconvolve_report.html",    "Deconvolution"),
    ("spatial_06_downstream_report.html",    "Downstream"),
]


def _embed(html_path: Path, height: int, key: str):
    """Read HTML file and embed it inline as a self-contained component."""
    html = html_path.read_text(encoding="utf-8", errors="replace")
    components.html(html, height=height, scrolling=True)


def render():
    init_state()
    dataset_name = st.session_state.get(KEY_DATASET_NAME, "Dataset")
    reports_dir  = st.session_state.get(KEY_REPORTS_DIR, "")

    st.markdown(f"## Report  ·  *{dataset_name}*")

    if not reports_dir:
        st.info("No report yet — run the pipeline first.")
        if st.button("← Go to Run"):
            st.session_state[KEY_PAGE] = 2
            st.rerun()
        return

    rdir     = Path(reports_dir)
    combined = rdir / "00_combined_report.html"

    # ── Nav + download ────────────────────────────────────────────────────────
    col_back, col_dl = st.columns([1, 2])
    with col_back:
        if st.button("← Back to Run"):
            st.session_state[KEY_PAGE] = 2
            st.rerun()
    with col_dl:
        if combined.exists():
            with open(combined, "rb") as f:
                st.download_button(
                    "⬇ Download combined report",
                    data=f.read(),
                    file_name=f"{dataset_name.replace(' ', '_')}_report.html",
                    mime="text/html",
                )

    st.divider()

    # ── Combined report inline ────────────────────────────────────────────────
    if combined.exists():
        st.markdown("### Combined report")
        st.caption(
            "All pipeline steps in one tabbed document. "
            "Use the Download button above for full-screen viewing."
        )
        _embed(combined, height=850, key="combined")
    else:
        st.warning(
            "Combined report not found — the pipeline may still be running "
            "or an earlier step failed."
        )
        # Show which individual reports DO exist as a diagnostic
        found = [lbl for f, lbl in STEP_REPORTS if (rdir / f).exists()]
        if found:
            st.info(f"Found individual reports for: {', '.join(found)}")

    # ── Individual step reports ───────────────────────────────────────────────
    available = [(f, lbl) for f, lbl in STEP_REPORTS if (rdir / f).exists()]
    if available:
        st.divider()
        st.markdown("### Individual step reports")
        for fname, label in available:
            fpath = rdir / fname
            with st.expander(f"📄 {label}  —  `{fname}`"):
                c1, c2 = st.columns([4, 1])
                with c1:
                    _embed(fpath, height=600, key=f"step_{fname}")
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
