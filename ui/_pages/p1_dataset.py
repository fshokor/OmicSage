"""
OmicSage UI — Page 1: Dataset
Changes:
  - Added dataset ID field (GEO accession) separate from display name
  - Added organism selector (human / mouse / other)
  - Organism stored in session state for Configure page to pick up
"""
from pathlib import Path
import streamlit as st
from ui.state import *
from ui.defaults import MODALITY_STEPS, STEP_DEFAULTS_ON

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
    "scRNA-seq": ("Primary data path",
                  "data/raw/GSE166635/  or  data/processed/GSE166635/01_qc.h5ad"),
    "CITE-seq":  ("ADT data path",
                  "data/processed/GSE194122/01_qc_adt.h5ad"),
    "Multiome":  ("ATAC data path",
                  "data/processed/GSE194122/01_qc_atac.h5ad"),
    "Spatial":   ("Spatial data path",
                  "data/benchmark/kuppe_visium_human_heart_2022_control.h5ad"),
}


def render():
    init_state()
    modality = st.session_state[KEY_MODALITY]

    st.markdown("## Dataset")
    st.divider()

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
            help="Human-readable name shown in reports and the UI.",
        )
        st.session_state[KEY_DATASET_NAME] = name
    with col2:
        dataset_id = st.text_input(
            "Dataset ID",
            value=st.session_state.get(KEY_DATASET_ID, ""),
            placeholder="e.g. GSE166635",
            help="GEO accession or unique project ID. Used for folder names and "
                 "GEO metadata lookup in reports. Leave blank to auto-generate from name.",
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

    # ── 5. RNA reference path (CITE / Multiome / Spatial) ────────────────────
    if modality in NEEDS_RNA_PATH:
        st.divider()
        rna_label, rna_caption = NEEDS_RNA_PATH[modality]
        st.markdown(f"### 5 · {rna_label}")
        st.caption(rna_caption)
        st.caption("Optional — you can leave this blank if skipping integration/deconvolution steps.")

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
