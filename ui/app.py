
"""
OmicSage — Streamlit UI 
==============================
Run from repo root:
    streamlit run ui/app.py
"""
import sys, os
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components
 
_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))
os.chdir(_root)
 
from ui.state import init_state, KEY_PAGE, KEY_CONFIG, KEY_RUN_STATUS, KEY_DATA_PATH
init_state()
 
# ── Sidebar open/close state ──────────────────────────────────────────────────
if "sidebar_open" not in st.session_state:
    st.session_state["sidebar_open"] = True
 
sidebar_state = "expanded" if st.session_state["sidebar_open"] else "collapsed"
 
st.set_page_config(
    page_title="OmicSage",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state=sidebar_state,
    menu_items={},
)
 
# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');
 
[data-testid="stToolbar"],[data-testid="stDecoration"],
[data-testid="stStatusWidget"],[data-testid="stSidebarNav"],
[data-testid="collapsedControl"],
.stDeployButton,#MainMenu,footer
{ display:none !important; visibility:hidden !important; }
 
[data-testid="stMain"] .block-container { padding-top:3.5rem !important; }
 
html,body,[data-testid="stAppViewContainer"] {
    background:#0D1117 !important; color:#E2F0F9 !important;
    font-family:'Inter',system-ui,sans-serif !important;
}
[data-testid="stSidebar"] {
    background:#161B22 !important;
    border-right:1px solid rgba(255,255,255,0.07) !important;
}
h1,h2,h3 { font-family:'Inter',sans-serif !important; font-weight:700 !important;
            color:#E2F0F9 !important; letter-spacing:-0.4px; }
h2 { font-size:1.5rem !important; }
h3 { font-size:1.05rem !important; color:#CBD5E1 !important; }
[data-testid="stCaption"] p { color:#64748B !important; font-size:0.79rem !important; }
code,pre,[data-testid="stCode"] {
    font-family:'JetBrains Mono',monospace !important; font-size:0.77rem !important;
    background:#0D1117 !important; border:1px solid rgba(255,255,255,0.08) !important;
    border-radius:6px !important; color:#A5F3FC !important;
}
[data-testid="stButton"] button[kind="primary"] {
    background:#3B82F6 !important; border:none !important; color:#fff !important;
    font-weight:600 !important; border-radius:8px !important;
}
[data-testid="stButton"] button[kind="primary"]:hover { background:#2563EB !important; }
[data-testid="stButton"] button[kind="secondary"] {
    background:rgba(255,255,255,0.05) !important;
    border:1px solid rgba(255,255,255,0.12) !important;
    color:#CBD5E1 !important; border-radius:8px !important;
}
[data-testid="stSlider"] > div > div > div { background:#3B82F6 !important; }
[data-testid="stTextInput"] input,
[data-testid="stSelectbox"] > div,
[data-testid="stNumberInput"] input {
    background:#161B22 !important; border:1px solid rgba(255,255,255,0.12) !important;
    color:#E2F0F9 !important; border-radius:8px !important;
}
hr { border-color:rgba(255,255,255,0.08) !important; margin:1rem 0 !important; }
[data-testid="stExpander"] {
    border:1px solid rgba(255,255,255,0.08) !important;
    border-radius:8px !important; background:#161B22 !important;
}
[data-testid="stTabs"] [data-baseweb="tab-list"] {
    background:transparent !important;
    border-bottom:1px solid rgba(255,255,255,0.08) !important;
}
 
/* Sidebar toggle button — fixed top-left */
div[data-testid="stButton"].sidebar-toggle-wrapper > button {
    position:fixed !important;
    top:0.5rem !important;
    left:0.5rem !important;
    z-index:999999 !important;
    width:2rem !important;
    height:2rem !important;
    min-height:0 !important;
    padding:0 !important;
    font-size:1.1rem !important;
    background:#1E2530 !important;
    border:1px solid rgba(255,255,255,0.2) !important;
    border-radius:6px !important;
    color:#94A3B8 !important;
    display:flex !important;
    align-items:center !important;
    justify-content:center !important;
    line-height:1 !important;
}
div[data-testid="stButton"].sidebar-toggle-wrapper > button:hover {
    background:#2D3748 !important;
    color:#E2F0F9 !important;
    border-color:rgba(255,255,255,0.35) !important;
}
</style>
""", unsafe_allow_html=True)
 
# ── Sidebar toggle button — real Streamlit button, CSS-positioned ─────────────
# We use a real st.button so clicking it triggers st.rerun() and flips state.
# CSS positions it fixed top-left so it floats above all content.
st.markdown('<div class="sidebar-toggle-wrapper">', unsafe_allow_html=True)
icon = "✕" if st.session_state["sidebar_open"] else "☰"
if st.button(icon, key="sidebar_toggle"):
    st.session_state["sidebar_open"] = not st.session_state["sidebar_open"]
    st.rerun()
st.markdown('</div>', unsafe_allow_html=True)
 
# ── Sidebar ───────────────────────────────────────────────────────────────────
from ui.components.sidebar import render as render_sidebar
render_sidebar()
 
# ── Navigation bar ────────────────────────────────────────────────────────────
def _nav():
    current    = st.session_state.get(KEY_PAGE, 0)
    data_ok    = bool(st.session_state.get(KEY_DATA_PATH))
    config_ok  = bool(st.session_state.get(KEY_CONFIG))
    run_status = st.session_state.get(KEY_RUN_STATUS, "idle")
 
    steps = [
        ("1 · Dataset",   True),
        ("2 · Configure", data_ok or config_ok),
        ("3 · Run",       config_ok),
        ("4 · Report",    run_status == "done"),
    ]
    done_flags = [
        data_ok and current != 0,
        config_ok and current != 1,
        run_status in ("done","error") and current != 2,
        run_status == "done" and current != 3,
    ]
 
    cols = st.columns(4)
    for i, ((label, enabled), done) in enumerate(zip(steps, done_flags)):
        with cols[i]:
            prefix = "✓ " if done else ""
            st.button(
                f"{prefix}{label}",
                key=f"nav_{i}",
                type="primary" if i == current else "secondary",
                use_container_width=True,
                disabled=not enabled,
                on_click=lambda idx=i: st.session_state.update({KEY_PAGE: idx}),
            )
 
_nav()
st.divider()
 
# ── Page routing ──────────────────────────────────────────────────────────────
page = st.session_state.get(KEY_PAGE, 0)
 
if page == 0:
    from ui._pages.p1_dataset import render; render()
elif page == 1:
    from ui._pages.p2_configure import render; render()
elif page == 2:
    from ui._pages.p3_run import render; render()
elif page == 3:
    from ui._pages.p4_report import render; render()