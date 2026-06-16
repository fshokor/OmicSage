"""OmicSage UI — Sidebar with project history"""
import streamlit as st
from ui.state import *
from ui.defaults import STEP_LABELS
from ui.history import get_recent, status_icon
from ui.config_io import load_config, parse_config_into_state

STATUS_ICON = {"idle": "⏸", "running": "⚙️", "done": "✅", "error": "❌"}


def _apply_parsed(parsed: dict) -> None:
    st.session_state[KEY_MODALITY]       = parsed.get("modality", "scRNA-seq")
    st.session_state[KEY_DATASET_NAME]   = parsed.get("dataset_name", "")
    st.session_state[KEY_DATASET_ID]     = parsed.get("dataset_id", "")
    st.session_state[KEY_DATA_PATH]      = parsed.get("data_path", "")
    st.session_state[KEY_RNA_PATH]       = parsed.get("rna_path", "")
    st.session_state[KEY_REPORTS_DIR]    = parsed.get("reports_dir", "")
    st.session_state[KEY_SELECTED_STEPS] = parsed.get("selected_steps", [])
    st.session_state[KEY_STEP_PARAMS]    = parsed.get("step_params", {})
    st.session_state[KEY_CONFIG]         = parsed.get("config")
    st.session_state[KEY_CONFIG_PATH]    = parsed.get("config_path", "")
    st.session_state["organism"]         = parsed.get("organism", "human")


def render():
    init_state()
    with st.sidebar:
        # ── Logo ──────────────────────────────────────────────────────────────
        st.markdown("""
        <div style="display:flex;align-items:center;gap:10px;padding:4px 0 12px 0;">
            <span style="font-size:1.5rem;">🧬</span>
            <span style="font-size:1.2rem;font-weight:700;letter-spacing:-0.5px;
                         color:#E2F0F9;">OmicSage</span>
        </div>""", unsafe_allow_html=True)

        # ── Current project ───────────────────────────────────────────────────
        status = st.session_state.get(KEY_RUN_STATUS, "idle")
        st.markdown("**Current project**")
        st.markdown(f"""
        <div style="background:rgba(255,255,255,0.05);border-radius:8px;
                    padding:10px 12px;font-size:0.82rem;line-height:1.9;">
            <b>Dataset</b>: {st.session_state.get(KEY_DATASET_NAME) or '—'}<br>
            <b>ID</b>: {st.session_state.get(KEY_DATASET_ID) or '—'}<br>
            <b>Modality</b>: {st.session_state.get(KEY_MODALITY,'—')}<br>
            <b>Status</b>: {STATUS_ICON.get(status,'⏸')} {status.capitalize()}
        </div>""", unsafe_allow_html=True)

        # ── Configured steps ──────────────────────────────────────────────────
        sel = st.session_state.get(KEY_SELECTED_STEPS, [])
        if sel:
            st.divider()
            st.markdown("**Configured steps**")
            for s in sel:
                st.markdown(
                    f"<span style='font-size:0.78rem;color:#94A3B8;'>"
                    f"<code>{s}</code>  {STEP_LABELS.get(s,'')}</span>",
                    unsafe_allow_html=True,
                )

        st.divider()

        # ── Recent runs ───────────────────────────────────────────────────────
        history = get_recent(6)
        if history:
            st.markdown("**Recent runs**")
            for entry in history:
                icon  = status_icon(entry.get("status", "idle"))
                date  = entry.get("timestamp", "")[:10]
                name  = entry.get("dataset_name", "Unknown")
                mod   = entry.get("modality", "")
                label = f"{icon} {name[:22]}{'…' if len(name)>22 else ''}"
                st.caption(f"{date}  ·  `{mod}`")
                if st.button(label, key=f"sb_hist_{entry['id']}",
                             use_container_width=True, type="secondary"):
                    cfg_path = entry.get("config_path", "")
                    from pathlib import Path
                    if cfg_path and Path(cfg_path).exists():
                        try:
                            cfg    = load_config(cfg_path)
                            parsed = parse_config_into_state(cfg, cfg_path)
                            _apply_parsed(parsed)
                            st.session_state[KEY_PAGE] = 1   # Configure first
                            st.rerun()
                        except Exception as e:
                            st.error(f"Could not load: {e}")
                    else:
                        st.warning("Config file no longer exists.")

        st.divider()

        # ── New dataset ───────────────────────────────────────────────────────
        if st.button("↩ Start new dataset", use_container_width=True, type="secondary"):
            reset_all()
            st.rerun()


def _truncate(s: str, n: int = 35) -> str:
    if not s:
        return "—"
    return ("…" + s[-(n-1):]) if len(s) > n else s
