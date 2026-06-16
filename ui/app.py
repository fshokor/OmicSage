"""
OmicSage — Streamlit UI 
==============================
Run from repo root:
    streamlit run ui/app.py

Pages live in ui/_pages/ (underscore = not auto-discovered by Streamlit).
"""
import sys, os
from pathlib import Path

import streamlit as st

# Repo root on sys.path
_root = Path(__file__).resolve().parent.parent
if str(_root) not in sys.path:
    sys.path.insert(0, str(_root))
os.chdir(_root)

st.set_page_config(
    page_title="OmicSage",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={},
)

# ── CSS ───────────────────────────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

/* Kill all Streamlit chrome */
[data-testid="stToolbar"],[data-testid="stDecoration"],
[data-testid="stStatusWidget"],[data-testid="stSidebarNav"],
[data-testid="collapsedControl"],.stDeployButton,#MainMenu,footer
{ display:none !important; visibility:hidden !important; }

/* Remove toolbar padding */
[data-testid="stMain"] .block-container { padding-top:3.5rem !important; }

/* Root */
html,body,[data-testid="stAppViewContainer"] {
    background:#0D1117 !important; color:#E2F0F9 !important;
    font-family:'Inter',system-ui,sans-serif !important;
}

/* Sidebar */
[data-testid="stSidebar"] {
    background:#161B22 !important;
    border-right:1px solid rgba(255,255,255,0.07) !important;
}

/* Headings */
h1,h2,h3 { font-family:'Inter',sans-serif !important; font-weight:700 !important;
            color:#E2F0F9 !important; letter-spacing:-0.4px; }
h2 { font-size:1.5rem !important; }
h3 { font-size:1.05rem !important; color:#CBD5E1 !important; }

/* Captions */
[data-testid="stCaption"] p { color:#64748B !important; font-size:0.79rem !important; }

/* Code */
code,pre,[data-testid="stCode"] {
    font-family:'JetBrains Mono',monospace !important; font-size:0.77rem !important;
    background:#0D1117 !important; border:1px solid rgba(255,255,255,0.08) !important;
    border-radius:6px !important; color:#A5F3FC !important;
}

/* Buttons */
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

/* Inputs */
[data-testid="stSlider"] > div > div > div { background:#3B82F6 !important; }
[data-testid="stTextInput"] input,
[data-testid="stSelectbox"] > div,
[data-testid="stNumberInput"] input {
    background:#161B22 !important; border:1px solid rgba(255,255,255,0.12) !important;
    color:#E2F0F9 !important; border-radius:8px !important;
}

/* Misc */
hr { border-color:rgba(255,255,255,0.08) !important; margin:1rem 0 !important; }
[data-testid="stExpander"] {
    border:1px solid rgba(255,255,255,0.08) !important;
    border-radius:8px !important; background:#161B22 !important;
}
[data-testid="stTabs"] [data-baseweb="tab-list"] {
    background:transparent !important;
    border-bottom:1px solid rgba(255,255,255,0.08) !important;
}
</style>
""", unsafe_allow_html=True)

# ── State & sidebar ───────────────────────────────────────────────────────────
from ui.state import init_state, KEY_PAGE, KEY_CONFIG, KEY_RUN_STATUS, KEY_DATA_PATH
from ui.components.sidebar import render as render_sidebar

init_state()
render_sidebar()

# ── Navigation bar ────────────────────────────────────────────────────────────
def _nav():
    current     = st.session_state.get(KEY_PAGE, 0)
    data_ok     = bool(st.session_state.get(KEY_DATA_PATH))
    config_ok   = bool(st.session_state.get(KEY_CONFIG))
    run_status  = st.session_state.get(KEY_RUN_STATUS, "idle")

    steps = [
        ("1 · Dataset",   True),
        ("2 · Configure", data_ok or config_ok),   # unlocked when config loaded
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
