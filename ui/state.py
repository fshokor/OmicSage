"""
OmicSage UI — Session State
============================
Single source of truth for all st.session_state keys and defaults.
Call init_state() at the top of every page.
"""
import streamlit as st

# ── Keys ──────────────────────────────────────────────────────────────────────
KEY_PAGE         = "page"           # int 0-3
KEY_MODALITY     = "modality"       # "scRNA-seq" | "CITE-seq" | "Multiome" | "Spatial"
KEY_DATASET_NAME = "dataset_name"   # str
KEY_DATASET_ID   = "dataset_id"     # slugified str
KEY_DATA_PATH    = "data_path"      # str — primary data file/dir
KEY_RNA_PATH     = "rna_path"       # str — ref RNA for CITE/Multiome/Spatial
KEY_SELECTED_STEPS = "selected_steps"  # list[str]
KEY_STEP_PARAMS  = "step_params"    # dict[step -> dict of param values]
KEY_CONFIG       = "config"         # dict — final YAML dict
KEY_CONFIG_PATH  = "config_path"    # str — path to written temp YAML
KEY_REPORTS_DIR  = "reports_dir"    # str
KEY_RUN_STATUS   = "run_status"     # "idle"|"running"|"done"|"error"
KEY_RUN_LOG      = "run_log"        # list[str]
KEY_RUN_PROCESS  = "run_process"    # subprocess.Popen | None

DEFAULTS = {
    KEY_PAGE:           0,
    KEY_MODALITY:       "scRNA-seq",
    KEY_DATASET_NAME:   "",
    KEY_DATASET_ID:     "",
    KEY_DATA_PATH:      "",
    KEY_RNA_PATH:       "",
    KEY_SELECTED_STEPS: [],
    KEY_STEP_PARAMS:    {},
    KEY_CONFIG:         None,
    KEY_CONFIG_PATH:    None,
    KEY_REPORTS_DIR:    "",
    KEY_RUN_STATUS:     "idle",
    KEY_RUN_LOG:        [],
    KEY_RUN_PROCESS:    None,
}

def init_state():
    for k, v in DEFAULTS.items():
        if k not in st.session_state:
            st.session_state[k] = v

def reset_run():
    st.session_state[KEY_RUN_STATUS]  = "idle"
    st.session_state[KEY_RUN_LOG]     = []
    st.session_state[KEY_RUN_PROCESS] = None

def reset_all():
    for k, v in DEFAULTS.items():
        st.session_state[k] = v
