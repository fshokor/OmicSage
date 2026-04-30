"""OmicSage — Streamlit Web Interface (Phase 7 stub)."""
import streamlit as st

st.set_page_config(page_title="OmicSage", page_icon="🧬", layout="wide")

st.markdown("# 🧬 OmicSage\n### AI-Assisted Single-Cell Multi-Omics Analysis")
st.info("Status: Phase 0 complete. Full UI coming in Phase 7. Use the CLI for now: `python cli/omicsage.py --help`")

cols = st.columns(4)
with cols[0]: st.metric("scRNA-seq",  "Phase 1", "Week 2-6")
with cols[1]: st.metric("scATAC-seq", "Phase 4", "Week 13-16")
with cols[2]: st.metric("Spatial",    "Phase 5", "Week 16-19")
with cols[3]: st.metric("Multiome",   "Phase 6", "Week 19-22")
