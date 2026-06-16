"""
OmicSage UI — Page 1: Dataset
Options:
  A) New dataset — fill in modality, name, paths
  B) Load existing config — point at a YAML file to restore everything
"""
from pathlib import Path
import streamlit as st
import yaml

from ui.state import *
from ui.defaults import MODALITY_STEPS, STEP_DEFAULTS_ON
from ui.config_io import load_config, parse_config_into_state
from ui.history import get_recent, load_config_from_history, status_icon

MODALITY_DESCRIPTIONS = {
    "scRNA-seq": "Gene expression profiling per cell",
    "CITE-seq":  "RNA + protein surface markers (ADT)",
    "Multiome":  "Joint RNA + ATAC-seq",
    "Spatial":   "Spatially resolved transcriptomics (Visium, Xenium)",
}

NEEDS_RNA_PATH = {
    "CITE-seq": ("Annotated RNA path",
                 "Annotated RNA AnnData from the scRNA pipeline — e.g. data/processed/GSE194122/05_annotated.h5ad"),
    "Multiome": ("Annotated RNA path",
                 "Annotated RNA AnnData for label transfer — e.g. data/processed/GSE194122/05_annotated.h5ad"),
    "Spatial":  ("Reference scRNA path",
                 "Annotated scRNA-seq reference for deconvolution — e.g. data/processed/kuppe_snRNA/05_annotated.h5ad"),
}

DATA_PATH_LABELS = {
    "scRNA-seq": ("Primary data path",  "data/raw/GSE166635/  or  data/processed/GSE166635/01_qc.h5ad"),
    "CITE-seq":  ("ADT data path",      "data/processed/GSE194122/01_qc_adt.h5ad"),
    "Multiome":  ("ATAC data path",     "data/processed/GSE194122/01_qc_atac.h5ad"),
    "Spatial":   ("Spatial data path",  "data/benchmark/kuppe_visium_human_heart_2022_control.h5ad"),
}


def _apply_parsed(parsed: dict) -> None:
    """Write parsed config state into session state."""
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
    modality = st.session_state[KEY_MODALITY]

    st.markdown("## Dataset")
    st.divider()

    # ── Mode selector ─────────────────────────────────────────────────────────
    mode = st.radio(
        "How do you want to start?",
        ["🆕 New dataset", "📂 Load existing config"],
        horizontal=True,
        key="dataset_mode",
        label_visibility="collapsed",
    )

    st.divider()

    # ══════════════════════════════════════════════════════════════════════════
    # MODE A — Load existing config
    # ══════════════════════════════════════════════════════════════════════════
    if mode == "📂 Load existing config":
        st.markdown("### Load config")
        st.caption(
            "Point at any OmicSage YAML config file. "
            "The UI will extract modality, paths, steps and parameters "
            "and take you straight to Run."
        )

        config_path_input = st.text_input(
            "Config file path",
            value="",
            placeholder="config/runs/GSE166635.yaml",
            label_visibility="collapsed",
            key="load_config_path",
        )

        if config_path_input:
            p = Path(config_path_input)
            if p.exists() and p.suffix in (".yaml", ".yml"):
                st.success(f"✓ Found: `{config_path_input}`")

                # Parse the config to show a preview
                try:
                    cfg    = load_config(config_path_input)
                    parsed = parse_config_into_state(cfg, config_path_input)
                except Exception as e:
                    st.error(f"Failed to parse config: {e}")
                    parsed = None

                if parsed:
                    # ── Preview what was detected ──────────────────────────
                    st.divider()
                    st.markdown("### Confirm before loading")

                    col_info, col_mod = st.columns([1, 1])
                    with col_info:
                        st.markdown(
                            f"**Dataset**: {parsed['dataset_name'] or '—'}  \n"
                            f"**ID**: {parsed['dataset_id'] or '—'}  \n"
                            f"**Organism**: {parsed['organism']}  \n"
                            f"**Steps enabled**: {len(parsed['selected_steps'])}"
                        )
                    with col_mod:
                        # Let user confirm or override the detected modality
                        modality_opts = ["scRNA-seq", "CITE-seq", "Multiome", "Spatial"]
                        detected = parsed.get("modality", "scRNA-seq")
                        detected_idx = modality_opts.index(detected) if detected in modality_opts else 0
                        chosen_modality = st.selectbox(
                            "Modality",
                            options=modality_opts,
                            index=detected_idx,
                            key="load_modality_override",
                            help=f"Auto-detected as **{detected}** from config. Change if incorrect.",
                        )
                        if chosen_modality != detected:
                            st.caption(f"⚠ Overriding detected modality ({detected})")

                    if st.button("Load & go to Configure →", type="primary",
                                 use_container_width=True):
                        parsed["modality"] = chosen_modality
                        _apply_parsed(parsed)
                        st.session_state[KEY_PAGE] = 1
                        st.rerun()

            elif config_path_input:
                st.error("File not found or not a .yaml file.")

        col_load_spacer, col_recent = st.columns([1, 1])
        with col_load_spacer:
            pass  # config path input + preview fills this column

        # ── Recent runs ───────────────────────────────────────────────────────
        with col_recent:
            history = get_recent(8)
            if history:
                st.markdown("**Or pick a recent run:**")
                for entry in history:
                    icon = status_icon(entry.get("status", "idle"))
                    label = (
                        f"{icon} {entry['dataset_name']}  ·  "
                        f"`{entry['modality']}`  ·  "
                        f"{entry['timestamp'][:10]}"
                    )
                    if st.button(label, key=f"hist_{entry['id']}",
                                 use_container_width=True):
                        cfg_path = entry.get("config_path", "")
                        if cfg_path and Path(cfg_path).exists():
                            try:
                                cfg    = load_config(cfg_path)
                                parsed = parse_config_into_state(cfg, cfg_path)
                                _apply_parsed(parsed)
                                st.session_state[KEY_PAGE] = 2
                                st.rerun()
                            except Exception as e:
                                st.error(f"Could not reload: {e}")
                        else:
                            st.warning("Config file no longer exists at original path.")
            else:
                st.caption("No run history yet.")

        return  # don't show the new-dataset form below

    # ══════════════════════════════════════════════════════════════════════════
    # MODE B — New dataset
    # ══════════════════════════════════════════════════════════════════════════

    # ── 1. Modality ───────────────────────────────────────────────────────────
    st.markdown("### 1 · Modality")
    cols = st.columns(4)
    for i, mod in enumerate(["scRNA-seq", "CITE-seq", "Multiome", "Spatial"]):
        with cols[i]:
            is_sel = modality == mod
            border = "#3B82F6" if is_sel else "rgba(255,255,255,0.10)"
            bg     = "rgba(59,130,246,0.10)" if is_sel else "rgba(255,255,255,0.03)"
            st.markdown(
                f"""<div style="border:1.5px solid {border};border-radius:10px;
                    padding:12px 12px 10px;background:{bg};min-height:72px;">
                    <b style="font-size:0.9rem;">{mod}</b>
                    <div style="font-size:0.72rem;color:#94A3B8;margin-top:4px;
                         line-height:1.35;">{MODALITY_DESCRIPTIONS[mod]}</div>
                </div>""",
                unsafe_allow_html=True,
            )
            if st.button("✓ Selected" if is_sel else "Select",
                         key=f"mod_{mod}", use_container_width=True,
                         type="primary" if is_sel else "secondary"):
                st.session_state[KEY_MODALITY]       = mod
                st.session_state[KEY_SELECTED_STEPS] = list(STEP_DEFAULTS_ON[mod])
                st.session_state[KEY_STEP_PARAMS]    = {}
                st.rerun()

    modality = st.session_state[KEY_MODALITY]
    st.divider()

    # ── 2. Dataset name + ID ──────────────────────────────────────────────────
    st.markdown("### 2 · Dataset identity")
    col1, col2 = st.columns(2)
    with col1:
        name = st.text_input(
            "Display name",
            value=st.session_state.get(KEY_DATASET_NAME, ""),
            placeholder="e.g. HCC scRNA-seq (Wang et al. 2025)",
            help="Human-readable name used in reports and the UI.",
        )
        st.session_state[KEY_DATASET_NAME] = name
    with col2:
        dataset_id = st.text_input(
            "Dataset ID",
            value=st.session_state.get(KEY_DATASET_ID, ""),
            placeholder="e.g. GSE166635",
            help="GEO accession or project ID — used for folder names and config filename.",
        )
        st.session_state[KEY_DATASET_ID] = dataset_id

    st.divider()

    # ── 3. Organism ───────────────────────────────────────────────────────────
    st.markdown("### 3 · Organism")
    org_choice = st.radio(
        "Organism",
        options=["human", "mouse", "other"],
        index=["human", "mouse", "other"].index(
            st.session_state.get("organism", "human")
            if st.session_state.get("organism", "human") in ["human", "mouse"]
            else "other"
        ),
        horizontal=True,
        label_visibility="collapsed",
        key="organism_radio_p1",
    )
    if org_choice == "other":
        org_other = st.text_input(
            "Specify organism",
            value=st.session_state.get("organism", "")
                  if st.session_state.get("organism", "") not in ["human", "mouse"] else "",
            placeholder="e.g. zebrafish, rat, macaque",
            key="organism_other_input",
        )
        st.session_state["organism"] = org_other if org_other else "other"
    else:
        st.session_state["organism"] = org_choice

    st.divider()

    # ── 4. Data path ──────────────────────────────────────────────────────────
    path_label, path_placeholder = DATA_PATH_LABELS[modality]
    st.markdown(f"### 4 · {path_label}")
    st.caption("Enter the full path on this machine. The pipeline reads it directly — no copying.")

    data_path = st.text_input(
        path_label,
        value=st.session_state.get(KEY_DATA_PATH, ""),
        placeholder=path_placeholder,
        label_visibility="collapsed",
        key="data_path_input",
    )

    if data_path:
        p = Path(data_path)
        if p.exists():
            kind = "directory" if p.is_dir() else f"file  ·  {p.stat().st_size/1e9:.2f} GB"
            st.success(f"✓  Found: `{data_path}` — {kind}")
            st.session_state[KEY_DATA_PATH] = data_path
        else:
            st.error("Path not found — check for typos.")
    elif st.session_state.get(KEY_DATA_PATH):
        st.info(f"Currently set: `{st.session_state[KEY_DATA_PATH]}`")

    # ── 5. RNA reference path ─────────────────────────────────────────────────
    if modality in NEEDS_RNA_PATH:
        st.divider()
        rna_label, rna_caption = NEEDS_RNA_PATH[modality]
        st.markdown(f"### 5 · {rna_label}")
        st.caption(rna_caption)
        st.caption("Optional — leave blank if skipping integration/deconvolution steps.")

        rna_path = st.text_input(
            rna_label,
            value=st.session_state.get(KEY_RNA_PATH, ""),
            placeholder="data/processed/.../05_annotated.h5ad",
            label_visibility="collapsed",
            key="rna_path_input",
        )
        if rna_path:
            if Path(rna_path).exists():
                st.success(f"✓  Found: `{rna_path}`")
            else:
                st.warning("Path not found — you can still proceed.")
            st.session_state[KEY_RNA_PATH] = rna_path

    st.divider()

    # ── Continue ──────────────────────────────────────────────────────────────
    data_ok = bool(st.session_state.get(KEY_DATA_PATH))
    name_ok = bool(st.session_state.get(KEY_DATASET_NAME, "").strip())

    if not name_ok:
        st.caption("✏  Enter a dataset name to continue.")
    elif not data_ok:
        st.caption("📂  Enter a data path to continue.")

    if st.button("Continue to Configure →", type="primary",
                 disabled=not (data_ok and name_ok)):
        if not st.session_state[KEY_SELECTED_STEPS]:
            st.session_state[KEY_SELECTED_STEPS] = list(STEP_DEFAULTS_ON[modality])
        st.session_state[KEY_PAGE] = 1
        st.rerun()
