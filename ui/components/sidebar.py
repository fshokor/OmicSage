"""OmicSage UI — Sidebar"""
import streamlit as st
from ui.state import *
from ui.defaults import STEP_LABELS

STATUS_ICON = {"idle": "⏸", "running": "⚙️", "done": "✅", "error": "❌"}


def render():
    init_state()
    with st.sidebar:
        st.markdown("""
        <div style="display:flex;align-items:center;gap:10px;padding:4px 0 16px 0;">
            <span style="font-size:1.5rem;">🧬</span>
            <span style="font-size:1.2rem;font-weight:700;letter-spacing:-0.5px;color:#E2F0F9;">
                OmicSage</span>
        </div>""", unsafe_allow_html=True)

        status = st.session_state.get(KEY_RUN_STATUS, "idle")
        st.markdown("**Current project**")
        st.markdown(f"""
        <div style="background:rgba(255,255,255,0.05);border-radius:8px;
                    padding:10px 12px;font-size:0.82rem;line-height:1.9;">
            <b>Dataset</b>: {st.session_state.get(KEY_DATASET_NAME) or '—'}<br>
            <b>Modality</b>: {st.session_state.get(KEY_MODALITY,'—')}<br>
            <b>Data</b>: <span style="color:#64748B;font-size:0.75rem;">
                {_truncate(st.session_state.get(KEY_DATA_PATH,'—'))}</span><br>
            <b>Status</b>: {STATUS_ICON.get(status,'⏸')} {status.capitalize()}
        </div>""", unsafe_allow_html=True)

        # Selected steps summary
        sel = st.session_state.get(KEY_SELECTED_STEPS, [])
        if sel:
            st.divider()
            st.markdown("**Configured steps**")
            for s in sel:
                st.markdown(f"&nbsp;&nbsp;`{s}`  {STEP_LABELS.get(s,'')}")

        st.divider()
        if st.button("↩ Start new dataset", use_container_width=True, type="secondary"):
            reset_all()
            st.rerun()


def _truncate(s: str, n: int = 35) -> str:
    if not s:
        return "—"
    return ("…" + s[-(n-1):]) if len(s) > n else s
