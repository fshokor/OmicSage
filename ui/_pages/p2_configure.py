"""
OmicSage UI — Page 2: Configure
Fixes:
  - Unique key= on every widget (no more DuplicateElementId errors)
  - Organism moved to Dataset page; read from session state here
  - Sliders paired with number_input for direct value entry
  - Clustering resolution: free-text input for custom values
  - dataset_id threaded into config builder
"""
import copy
import streamlit as st
import yaml

from ui.state import *
from ui.defaults import (
    MODALITY_STEPS, STEP_LABELS, STEP_PARAM_DEFAULTS,
    SPATIAL_QC_DEFAULTS, SPATIAL_REDUCE_DEFAULTS, SPATIAL_CLUSTER_DEFAULTS,
)
from ui.config_builder import (
    build_scrna, build_cite, build_multiome, build_spatial,
    write_temp_yaml, get_reports_dir, slugify,
)


def render():
    init_state()
    modality     = st.session_state[KEY_MODALITY]
    dataset_name = st.session_state[KEY_DATASET_NAME]
    organism     = st.session_state.get("organism", "human")

    st.markdown(f"## Configure  ·  *{dataset_name}*  ·  `{modality}`")
    st.caption(f"Organism: **{organism}**  ·  Change on the Dataset page if needed.")
    st.divider()

    # ── Step selector ─────────────────────────────────────────────────────────
    st.markdown("### Steps to run")
    st.caption("Select the steps to include. Parameters appear only for selected steps.")

    all_steps    = MODALITY_STEPS[modality]
    selected_prev = set(st.session_state.get(KEY_SELECTED_STEPS, []))
    selected_now  = set()

    n_cols = 3
    rows = [all_steps[i:i+n_cols] for i in range(0, len(all_steps), n_cols)]
    for row in rows:
        cols = st.columns(n_cols)
        for j, step in enumerate(row):
            with cols[j]:
                if st.checkbox(STEP_LABELS.get(step, step),
                               value=(step in selected_prev),
                               key=f"chk_{step}"):
                    selected_now.add(step)

    st.session_state[KEY_SELECTED_STEPS] = list(selected_now)

    if not selected_now:
        st.warning("Select at least one step to continue.")
        return

    st.divider()

    # ── Per-step parameters ───────────────────────────────────────────────────
    st.markdown("### Parameters")
    st.caption("Defaults match your existing config files. Sliders also accept typed values.")

    if KEY_STEP_PARAMS not in st.session_state:
        st.session_state[KEY_STEP_PARAMS] = {}
    params = st.session_state[KEY_STEP_PARAMS]

    for step in [s for s in all_steps if s in selected_now]:
        defaults = copy.deepcopy(STEP_PARAM_DEFAULTS.get(step, {}))
        if modality == "Spatial":
            if step == "qc":     defaults = copy.deepcopy(SPATIAL_QC_DEFAULTS)
            elif step == "reduce":  defaults = copy.deepcopy(SPATIAL_REDUCE_DEFAULTS)
            elif step == "cluster": defaults = copy.deepcopy(SPATIAL_CLUSTER_DEFAULTS)

        current = params.get(step, defaults)
        with st.expander(f"**{STEP_LABELS.get(step, step)}**  —  `{step}`", expanded=False):
            params[step] = _render_step(step, current, modality, organism)

    st.session_state[KEY_STEP_PARAMS] = params
    st.divider()

    # ── YAML preview ──────────────────────────────────────────────────────────
    cfg = _build(modality, dataset_name, organism, selected_now, params)
    with st.expander("Preview YAML config", expanded=False):
        st.code(yaml.dump(cfg, default_flow_style=False, sort_keys=False,
                          allow_unicode=True), language="yaml")

    st.divider()

    col_back, col_fwd = st.columns([1, 3])
    with col_back:
        if st.button("← Back", use_container_width=True):
            st.session_state[KEY_PAGE] = 0
            st.rerun()
    with col_fwd:
        if st.button("Save & go to Run →", type="primary", use_container_width=True):
            config_path = write_temp_yaml(cfg)
            st.session_state[KEY_CONFIG]      = cfg
            st.session_state[KEY_CONFIG_PATH] = config_path
            st.session_state[KEY_REPORTS_DIR] = get_reports_dir(cfg, modality)
            st.session_state[KEY_PAGE]        = 2
            st.rerun()


# ── Number input helper (replaces slider+number combo) ───────────────────────
def _slider_num(label, min_val, max_val, default, step_val, key, fmt="%.0f"):
    """
    Simple number input — type a value directly or use arrow buttons.
    Dropped the slider entirely: two widgets sharing state caused sync bugs
    because Streamlit rerenders top-to-bottom and the second widget always
    overwrote the first. A single number_input has no sync issue.
    """
    is_float = isinstance(step_val, float) or isinstance(default, float)
    return st.number_input(
        label,
        min_value=float(min_val) if is_float else int(min_val),
        max_value=float(max_val) if is_float else int(max_val),
        value=float(default) if is_float else int(default),
        step=float(step_val) if is_float else int(step_val),
        key=f"{key}_num",
    )


# ── Step parameter renderer ───────────────────────────────────────────────────
def _render_step(step: str, cur: dict, modality: str, organism: str) -> dict:
    p = dict(cur)
    k = step  # prefix for unique widget keys

    # ── QC (scRNA / CITE-seq / Multiome) ─────────────────────────────────────
    if step == "qc" and modality != "Spatial":
        col1, col2 = st.columns(2)
        with col1:
            p["min_genes"] = int(_slider_num("Min genes per cell", 50, 2000,
                                             p.get("min_genes", 200), 50, f"{k}_ming"))
            p["max_genes"] = int(_slider_num("Max genes per cell", 1000, 20000,
                                             p.get("max_genes", 6000), 500, f"{k}_maxg"))
        with col2:
            p["max_mt_pct"] = _slider_num("Max mitochondrial %", 1.0, 60.0,
                                           p.get("max_mt_pct", 20.0), 0.5,
                                           f"{k}_mt", fmt="%.1f")
            p["remove_doublets"] = st.toggle("Remove doublets (Scrublet)",
                                             value=bool(p.get("remove_doublets", True)),
                                             key=f"{k}_dbl")

    # ── QC Spatial ────────────────────────────────────────────────────────────
    elif step == "qc" and modality == "Spatial":
        col1, col2 = st.columns(2)
        with col1:
            p["min_counts"] = int(_slider_num("Min UMI counts per spot", 100, 10000,
                                              p.get("min_counts", 500), 100, f"{k}_minc"))
            p["min_genes"]  = int(_slider_num("Min genes per spot", 50, 1000,
                                              p.get("min_genes", 200), 50, f"{k}_ming"))
        with col2:
            p["max_counts"] = st.number_input("Max UMI counts",
                                               value=int(p.get("max_counts", 100000)),
                                               step=1000, key=f"{k}_maxc")
            p["max_mt_pct"] = _slider_num("Max MT%", 1.0, 60.0,
                                           p.get("max_mt_pct", 20.0), 0.5,
                                           f"{k}_mt", fmt="%.1f")
        p["mt_prefix"] = st.text_input("MT gene prefix", value=p.get("mt_prefix","MT-"),
                                        key=f"{k}_mtp")

    # ── Normalize ─────────────────────────────────────────────────────────────
    elif step == "normalize":
        col1, col2 = st.columns(2)
        with col1:
            p["batch_key"]   = st.text_input("Batch key (obs column, blank = none)",
                                              value=p.get("batch_key","") or "",
                                              placeholder="e.g. sample, batch, donor",
                                              key=f"{k}_bk")
            p["n_top_genes"] = int(_slider_num("HVG count", 500, 6000,
                                               p.get("n_top_genes",2000), 250, f"{k}_hvg"))
        with col2:
            p["target_sum"] = st.number_input("Normalization target sum",
                                               value=int(p.get("target_sum",10000)),
                                               step=1000, key=f"{k}_ts")
            p["hvg_flavor"] = st.selectbox("HVG flavour",
                                            ["seurat","seurat_v3","cell_ranger"],
                                            index=["seurat","seurat_v3","cell_ranger"]
                                            .index(p.get("hvg_flavor","seurat")),
                                            key=f"{k}_hf")

    # ── Reduce ────────────────────────────────────────────────────────────────
    elif step == "reduce":
        col1, col2 = st.columns(2)
        with col1:
            p["n_comps"]     = int(_slider_num("PCA components to compute", 10, 100,
                                               p.get("n_comps",50), 5, f"{k}_nc"))
            p["n_neighbors"] = int(_slider_num("kNN neighbours", 5, 50,
                                               p.get("n_neighbors",15), 5, f"{k}_nn"))
        with col2:
            p["n_pcs_method"] = st.selectbox("PC selection method",
                                              ["elbow","cumvar","fixed"],
                                              index=["elbow","cumvar","fixed"]
                                              .index(p.get("n_pcs_method","elbow")),
                                              key=f"{k}_npm")

    # ── Cluster ───────────────────────────────────────────────────────────────
    elif step == "cluster" and modality != "Spatial":
        st.caption("Select preset values and/or type custom ones below.")

        preset_opts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
        current_range = p.get("resolution_range", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5])
        selected_presets = st.multiselect(
            "Resolution sweep — preset values",
            options=preset_opts,
            default=[r for r in current_range if r in preset_opts],
            key=f"{k}_rr",
        )
        custom_str = st.text_input(
            "Add custom resolution values (comma-separated)",
            value=", ".join(str(r) for r in current_range if r not in preset_opts),
            placeholder="e.g. 0.35, 0.75, 1.8",
            key=f"{k}_custom",
        )
        custom_vals = []
        for v in custom_str.split(","):
            v = v.strip()
            if v:
                try:
                    custom_vals.append(float(v))
                except ValueError:
                    st.warning(f"Ignored non-numeric value: `{v}`")

        combined = sorted(set(selected_presets + custom_vals))
        p["resolution_range"] = combined if combined else [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5]
        st.caption(f"Will sweep: `{p['resolution_range']}`")

        col1, col2 = st.columns(2)
        with col1:
            use_override = st.toggle("Pin a specific resolution (override auto-select)",
                                     value=p.get("resolution_override") is not None,
                                     key=f"{k}_uo")
        with col2:
            if use_override:
                p["resolution_override"] = st.number_input(
                    "Resolution override", min_value=0.05, max_value=5.0,
                    value=float(p.get("resolution_override") or 0.6),
                    step=0.05, key=f"{k}_ro",
                )
            else:
                p["resolution_override"] = None

    # ── Cluster Spatial ───────────────────────────────────────────────────────
    elif step == "cluster" and modality == "Spatial":
        col1, col2 = st.columns(2)
        with col1:
            p["resolution"]  = _slider_num("Clustering resolution", 0.1, 3.0,
                                            p.get("resolution",0.5), 0.05,
                                            f"{k}_res", fmt="%.2f")
            p["n_neighbors"] = int(_slider_num("kNN neighbours", 3, 30,
                                               p.get("n_neighbors",15), 1, f"{k}_nn"))
        with col2:
            p["n_pcs"]        = int(_slider_num("PCs used", 5, 50, p.get("n_pcs",30), 5, f"{k}_np"))
            p["random_state"] = st.number_input("Random state",
                                                 value=int(p.get("random_state",0)),
                                                 key=f"{k}_rs")

    # ── Annotate ──────────────────────────────────────────────────────────────
    elif step == "annotate":
        method_opts = ["celltypist","markers","sctype","singler","vote"]
        cur_methods = p.get("methods", ["celltypist","markers","vote"])
        p["methods"] = st.multiselect("Annotation methods", method_opts,
                                       default=[m for m in cur_methods if m in method_opts],
                                       key=f"{k}_meth")
        st.caption("`vote` = majority vote across selected methods. Always include it.")
        p["leiden_col"] = st.text_input("Leiden column name",
                                         value=p.get("leiden_col","leiden"),
                                         key=f"{k}_lc")
        if "singler" in p.get("methods",[]):
            singler_opts = ["hpca","blueprint_encode","dice","monaco_immune",
                            "novershtern_hematopoietic","mouse_rnaseq","hca"]
            cur_ref = p.get("singler_ref","novershtern_hematopoietic")
            idx = singler_opts.index(cur_ref) if cur_ref in singler_opts else 0
            p["singler_ref"] = st.selectbox("SingleR reference", singler_opts,
                                             index=idx, key=f"{k}_sr")

    # ── DEG / Multiome DEG ────────────────────────────────────────────────────
    elif step in ("deg","multiome_deg"):
        col1, col2 = st.columns(2)
        with col1:
            p["groupby"] = st.text_input("Group by (obs column)",
                                          value=p.get("groupby","cell_type_vote"),
                                          key=f"{k}_gb")
            p["method"]  = st.selectbox("DEG method",
                                         ["wilcoxon","t-test","logreg"],
                                         index=["wilcoxon","t-test","logreg"]
                                         .index(p.get("method","wilcoxon")),
                                         key=f"{k}_dm")
            p["n_genes"] = int(_slider_num("Max genes per group", 50, 2000,
                                           p.get("n_genes",500), 50, f"{k}_ng"))
        with col2:
            p["min_logfc"]    = _slider_num("Min log₂FC", 0.05, 3.0,
                                             p.get("min_logfc",0.25), 0.05,
                                             f"{k}_lfc", fmt="%.2f")
            p["max_pval_adj"] = _slider_num("Max adj p-value", 0.001, 0.2,
                                             p.get("max_pval_adj",0.05), 0.005,
                                             f"{k}_pv", fmt="%.3f")
        excl = p.get("exclude_gene_prefixes",[])
        excl_str = st.text_input("Exclude gene prefixes (comma-separated)",
                                  value=", ".join(excl) if excl else "",
                                  placeholder="RPL, RPS, MT-",
                                  key=f"{k}_excl")
        p["exclude_gene_prefixes"] = [x.strip() for x in excl_str.split(",") if x.strip()]

    # ── GSEA ─────────────────────────────────────────────────────────────────
    elif step in ("gsea","gsea_cite"):
        gs_opts = ["GO_Biological_Process_2023","KEGG_2021_Human","Reactome_2022",
                   "MSigDB_Hallmark_2020","WikiPathway_2023_Human"]
        cur_gs = p.get("gene_sets", ["GO_Biological_Process_2023","KEGG_2021_Human","Reactome_2022"])
        p["gene_sets"] = st.multiselect("Gene sets", gs_opts,
                                         default=[g for g in cur_gs if g in gs_opts],
                                         key=f"{k}_gs")
        col1, col2 = st.columns(2)
        with col1:
            # Key fix: unique key so it doesn't clash with organism radio on configure page
            p["organism"] = st.radio("Organism", ["human","mouse"],
                                      index=0 if p.get("organism","human")=="human" else 1,
                                      horizontal=True, key=f"{k}_org")
            p["min_logfc"] = _slider_num("Min log₂FC for gene list", 0.05, 3.0,
                                          p.get("min_logfc",0.25), 0.05,
                                          f"{k}_lfc", fmt="%.2f")
        with col2:
            p["min_genes"] = int(_slider_num("Min query genes to run", 2, 30,
                                             p.get("min_genes",5), 1, f"{k}_mg"))
            if step == "gsea_cite":
                p["direction"] = st.selectbox("Direction",["up","down","both"],
                                               index=["up","down","both"]
                                               .index(p.get("direction","up")),
                                               key=f"{k}_dir")
        excl = p.get("exclude_gene_prefixes",[])
        excl_str = st.text_input("Exclude gene prefixes",
                                  value=", ".join(excl) if excl else "",
                                  placeholder="RPL, RPS, MT-", key=f"{k}_excl")
        p["exclude_gene_prefixes"] = [x.strip() for x in excl_str.split(",") if x.strip()]

    # ── Harmony ───────────────────────────────────────────────────────────────
    elif step in ("harmony","harmony_adt"):
        col1, col2 = st.columns(2)
        with col1:
            p["batch_key"]   = st.text_input("Batch key (obs column)",
                                              value=p.get("batch_key","batch"),
                                              key=f"{k}_bk")
            p["n_pcs"]       = int(_slider_num("PCs for Harmony", 10, 100,
                                               p.get("n_pcs",50), 5, f"{k}_np"))
        with col2:
            p["n_neighbors"] = int(_slider_num("kNN after correction", 5, 50,
                                               p.get("n_neighbors",15), 5, f"{k}_nn"))
            p["random_state"]= st.number_input("Random state",
                                                value=int(p.get("random_state",0)),
                                                key=f"{k}_rs")

    # ── Cluster Harmony ───────────────────────────────────────────────────────
    elif step == "cluster_harmony":
        preset_opts = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5]
        cur_range   = p.get("resolution_range",[0.2,0.4,0.6,0.8,1.0])
        p["resolution_range"] = st.multiselect(
            "Resolution sweep (Harmony embedding)", preset_opts,
            default=[r for r in cur_range if r in preset_opts],
            key=f"{k}_rr",
        )

    # ── Pseudobulk ────────────────────────────────────────────────────────────
    elif step == "pseudobulk":
        col1, col2 = st.columns(2)
        with col1:
            p["groupby"]      = st.text_input("Group by", value=p.get("groupby","cell_type_vote"), key=f"{k}_gb")
            p["donor_key"]    = st.text_input("Donor key", value=p.get("donor_key","batch"), key=f"{k}_dk")
            p["counts_layer"] = st.text_input("Raw counts layer", value=p.get("counts_layer","counts"), key=f"{k}_cl")
        with col2:
            p["min_cells"]   = int(_slider_num("Min cells per pseudobulk", 5, 100,
                                               p.get("min_cells",10), 5, f"{k}_mc"))
            p["min_samples"] = int(_slider_num("Min donors per group", 2, 10,
                                               p.get("min_samples",3), 1, f"{k}_ms"))
        excl = p.get("exclude_gene_prefixes",[])
        excl_str = st.text_input("Exclude gene prefixes", value=", ".join(excl) if excl else "",
                                  placeholder="RPL, RPS, MT-", key=f"{k}_excl")
        p["exclude_gene_prefixes"] = [x.strip() for x in excl_str.split(",") if x.strip()]

    # ── Normalize ADT ─────────────────────────────────────────────────────────
    elif step == "normalize_adt":
        p["clr_axis"] = st.radio(
            "CLR axis",
            options=[0,1],
            format_func=lambda x: "0 — per-protein across cells (recommended)" if x==0
                                  else "1 — per-cell across proteins",
            index=int(p.get("clr_axis",0)),
            horizontal=True, key=f"{k}_ca",
        )

    # ── Doublets ──────────────────────────────────────────────────────────────
    elif step == "doublets":
        col1, col2 = st.columns(2)
        with col1:
            p["threshold"] = _slider_num("CLR expression threshold", 0.5, 8.0,
                                          p.get("threshold",2.5), 0.1,
                                          f"{k}_th", fmt="%.1f")
        with col2:
            p["filter_doublets"] = st.toggle("Remove doublets from AnnData",
                                             value=bool(p.get("filter_doublets",False)),
                                             key=f"{k}_fd")

    # ── Reduce ADT ────────────────────────────────────────────────────────────
    elif step == "reduce_adt":
        col1, col2 = st.columns(2)
        with col1:
            p["n_comps"]     = int(_slider_num("PCA components", 10, 100,
                                               p.get("n_comps",50), 5, f"{k}_nc"))
            p["n_pcs"]       = int(_slider_num("PCs used for graph", 5, 50,
                                               p.get("n_pcs",20), 5, f"{k}_np"))
        with col2:
            p["n_neighbors"] = int(_slider_num("kNN neighbours", 5, 50,
                                               p.get("n_neighbors",15), 5, f"{k}_nn"))

    # ── Annotate ADT ──────────────────────────────────────────────────────────
    elif step == "annotate_adt":
        col1, col2 = st.columns(2)
        with col1:
            p["resolution"]   = _slider_num("Leiden resolution", 0.05, 3.0,
                                             p.get("resolution",0.1), 0.05,
                                             f"{k}_res", fmt="%.2f")
        with col2:
            p["n_iterations"] = int(_slider_num("Leiden iterations", 1, 10,
                                                p.get("n_iterations",2), 1, f"{k}_ni"))
        st.caption("After running, inspect clusters then re-run with an annotation_map to assign cell types.")

    # ── Integration ───────────────────────────────────────────────────────────
    elif step == "integration":
        col1, col2 = st.columns(2)
        with col1:
            p["method"]    = st.selectbox("Integration method",["mofa","totalvi","both"],
                                           index=["mofa","totalvi","both"]
                                           .index(p.get("method","both")),
                                           key=f"{k}_meth")
            p["batch_key"] = st.text_input("Batch key",
                                            value=p.get("batch_key","batch"),
                                            key=f"{k}_bk")
        with col2:
            p["n_factors"]  = int(_slider_num("MOFA+ latent factors", 5, 50,
                                              p.get("n_factors",15), 5, f"{k}_nf"))
            p["max_epochs"] = int(_slider_num("totalVI max epochs", 5, 500,
                                              p.get("max_epochs",10), 5, f"{k}_me"))
        p["random_state"] = st.number_input("Random state",
                                             value=int(p.get("random_state",0)),
                                             key=f"{k}_rs")

    # ── ATAC QC ───────────────────────────────────────────────────────────────
    elif step == "atac_qc":
        col1, col2 = st.columns(2)
        with col1:
            p["min_peaks"]   = st.number_input("Min peaks per cell",
                                                value=int(p.get("min_peaks",750)),
                                                step=100, key=f"{k}_mnp")
            p["min_peak_counts"] = st.number_input("Min peak counts",
                                                    value=int(p.get("min_peak_counts",1500)),
                                                    step=100, key=f"{k}_mpc")
            p["max_nucleosome_signal"] = _slider_num("Max nucleosome signal", 0.5, 8.0,
                                                      p.get("max_nucleosome_signal",2.0),
                                                      0.1, f"{k}_ns", fmt="%.1f")
        with col2:
            p["max_peaks"]       = st.number_input("Max peaks per cell",
                                                    value=int(p.get("max_peaks",500000)),
                                                    step=10000, key=f"{k}_mxp")
            p["max_peak_counts"] = st.number_input("Max peak counts",
                                                    value=int(p.get("max_peak_counts",100000)),
                                                    step=1000, key=f"{k}_mxpc")
            p["min_cells"]       = int(_slider_num("Min cells per peak", 5, 50,
                                                   p.get("min_cells",15), 5, f"{k}_mc"))

    # ── ATAC Reduce ───────────────────────────────────────────────────────────
    elif step == "atac_reduce":
        col1, col2 = st.columns(2)
        with col1:
            p["n_components"]      = int(_slider_num("LSI components", 10, 100,
                                                     p.get("n_components",50), 5, f"{k}_nc"))
            p["n_neighbors"]       = int(_slider_num("kNN neighbours", 5, 50,
                                                     p.get("n_neighbors",15), 5, f"{k}_nn"))
        with col2:
            p["leiden_resolution"] = _slider_num("Leiden resolution", 0.1, 3.0,
                                                  p.get("leiden_resolution",0.5),
                                                  0.05, f"{k}_lr", fmt="%.2f")
            p["random_state"]      = st.number_input("Random state",
                                                      value=int(p.get("random_state",0)),
                                                      key=f"{k}_rs")

    # ── Multiome Integration ──────────────────────────────────────────────────
    elif step == "multiome_integration":
        col1, col2 = st.columns(2)
        with col1:
            p["method"]     = st.selectbox("Method",["multivi","mofa"],
                                            index=["multivi","mofa"]
                                            .index(p.get("method","multivi")),
                                            key=f"{k}_meth")
            p["batch_key"]  = st.text_input("Batch key",
                                             value=p.get("batch_key","batch"),
                                             key=f"{k}_bk")
            p["n_latent"]   = int(_slider_num("MultiVI latent dims", 10, 50,
                                              p.get("n_latent",20), 5, f"{k}_nl"))
        with col2:
            p["max_epochs"] = int(_slider_num("Max epochs", 5, 500,
                                              p.get("max_epochs",10), 5, f"{k}_me"))
            p["n_top_peaks"]= st.number_input("Top peaks", value=int(p.get("n_top_peaks",20000)),
                                               step=1000, key=f"{k}_ntp")

    # ── Deconvolve ────────────────────────────────────────────────────────────
    elif step == "deconvolve":
        col1, col2 = st.columns(2)
        with col1:
            p["method"]        = st.selectbox("Method",["nnls","cell2location"],
                                               index=["nnls","cell2location"]
                                               .index(p.get("method","nnls")),
                                               key=f"{k}_meth")
            p["cell_type_key"] = st.text_input("Cell type column (in reference)",
                                                value=p.get("cell_type_key","cell_type_original"),
                                                key=f"{k}_ctk")
        with col2:
            p["layer_ref"] = st.text_input("Raw counts layer (in reference)",
                                            value=p.get("layer_ref","counts"),
                                            key=f"{k}_lr")
            p["n_jobs"]    = int(_slider_num("Parallel jobs", 1, 16,
                                             p.get("n_jobs",4), 1, f"{k}_nj"))

    # ── Downstream (Spatial) ──────────────────────────────────────────────────
    elif step == "downstream":
        col1, col2 = st.columns(2)
        with col1:
            p["run_region_clustering"]   = st.toggle("Region clustering", value=bool(p.get("run_region_clustering",True)), key=f"{k}_rc")
            p["run_celltype_expression"] = st.toggle("Cell-type expression", value=bool(p.get("run_celltype_expression",True)), key=f"{k}_ce")
            p["run_co_occurrence"]       = st.toggle("Spatial co-occurrence", value=bool(p.get("run_co_occurrence",True)), key=f"{k}_co")
            p["run_svg_gsea"]            = st.toggle("SVG pathway analysis", value=bool(p.get("run_svg_gsea",True)), key=f"{k}_sg")
        with col2:
            p["run_nhood_enrichment"] = st.toggle("Neighbourhood enrichment", value=bool(p.get("run_nhood_enrichment",True)), key=f"{k}_ne")
            p["run_ligrec"]           = st.toggle("Ligand-receptor (CCC)", value=bool(p.get("run_ligrec",True)), key=f"{k}_lr2")
            p["ligrec_organism"]      = st.radio("LigRec organism", ["human","mouse"],
                                                  horizontal=True,
                                                  index=0 if p.get("ligrec_organism","human")=="human" else 1,
                                                  key=f"{k}_lo")
            p["n_jobs"] = int(_slider_num("Parallel jobs", 1, 16, p.get("n_jobs",4), 1, f"{k}_nj"))

    # ── Impute ────────────────────────────────────────────────────────────────
    elif step == "impute":
        col1, col2 = st.columns(2)
        with col1:
            p["method"]      = st.selectbox("Method",["tangram","gimvi"],
                                             index=["tangram","gimvi"]
                                             .index(p.get("method","tangram")),
                                             key=f"{k}_meth")
            p["n_top_genes"] = int(_slider_num("Top genes", 500, 5000,
                                               p.get("n_top_genes",2000), 250, f"{k}_ntg"))
        with col2:
            p["tangram_mode"] = st.selectbox("Tangram mode",["clusters","cells"],
                                              index=["clusters","cells"]
                                              .index(p.get("tangram_mode","clusters")),
                                              key=f"{k}_tm")
            p["device"]       = st.radio("Device",["cpu","cuda"], horizontal=True,
                                         index=0 if p.get("device","cpu")=="cpu" else 1,
                                         key=f"{k}_dev")

    # ── Ingest (Spatial) ──────────────────────────────────────────────────────
    elif step == "ingest":
        col1, col2 = st.columns(2)
        with col1:
            p["spatial_type"] = st.selectbox("Format",
                                              ["h5ad","visium_dir","xenium","merfish"],
                                              index=["h5ad","visium_dir","xenium","merfish"]
                                              .index(p.get("spatial_type","h5ad")),
                                              key=f"{k}_st")
            p["library_key"]  = st.text_input("Library/sample key",
                                               value=p.get("library_key","sample_name"),
                                               key=f"{k}_lk")
        with col2:
            p["load_images"] = st.toggle("Load tissue images",
                                         value=bool(p.get("load_images",True)),
                                         key=f"{k}_li")

    # ── Spatial Reduce ────────────────────────────────────────────────────────
    elif step == "reduce" and modality == "Spatial":
        col1, col2 = st.columns(2)
        with col1:
            p["n_top_genes"] = int(_slider_num("HVGs for PCA", 500, 6000,
                                               p.get("n_top_genes",3000), 250, f"{k}_ntg"))
            p["n_comps"]     = int(_slider_num("PCA components", 10, 100,
                                               p.get("n_comps",50), 5, f"{k}_nc"))
        with col2:
            p["n_neighbors"] = int(_slider_num("kNN neighbours", 3, 30,
                                               p.get("n_neighbors",6), 1, f"{k}_nn"))
            p["target_sum"]  = st.number_input("Normalization target sum",
                                                value=int(p.get("target_sum",10000)),
                                                step=1000, key=f"{k}_ts")

    # ── ATAC Annotate ─────────────────────────────────────────────────────────
    elif step == "atac_annotate":
        col1, col2 = st.columns(2)
        with col1:
            p["promoter_upstream_bp"] = st.number_input("Promoter window (bp)",
                                                         value=int(p.get("promoter_upstream_bp",2000)),
                                                         step=500, key=f"{k}_pub")
        with col2:
            p["rna_label_key"] = st.text_input("RNA label column",
                                                value=p.get("rna_label_key","cell_type_vote"),
                                                key=f"{k}_rlk")

    # ── GRN ───────────────────────────────────────────────────────────────────
    elif step == "multiome_grn":
        col1, col2 = st.columns(2)
        with col1:
            p["motif_db"]  = st.selectbox("Motif DB",["jaspar"],
                                           index=0, key=f"{k}_mdb")
            p["groupby"]   = st.text_input("Group by", value=p.get("groupby","atac_celltype"),
                                            key=f"{k}_gb")
        with col2:
            p["n_top_peaks"] = st.number_input("Top peaks per group",
                                                value=int(p.get("n_top_peaks",500)),
                                                step=100, key=f"{k}_ntp")
            p["min_cells"]   = int(_slider_num("Min cells", 5, 50,
                                               p.get("min_cells",10), 5, f"{k}_mc"))

    else:
        st.caption(f"No parameter UI for `{step}` — using defaults.")
        st.json(p)

    return p


# ── Config builder router ─────────────────────────────────────────────────────
def _build(modality, dataset_name, organism, selected_steps, params):
    ss  = st.session_state
    did = ss.get(KEY_DATASET_ID, "")
    if modality == "scRNA-seq":
        return build_scrna(dataset_name, ss[KEY_DATA_PATH], organism,
                           selected_steps, params, dataset_id=did)
    elif modality == "CITE-seq":
        return build_cite(dataset_name, ss[KEY_DATA_PATH], ss.get(KEY_RNA_PATH,""),
                          organism, selected_steps, params, dataset_id=did)
    elif modality == "Multiome":
        return build_multiome(dataset_name, ss[KEY_DATA_PATH], ss.get(KEY_RNA_PATH,""),
                              organism, selected_steps, params, dataset_id=did)
    elif modality == "Spatial":
        return build_spatial(dataset_name, ss[KEY_DATA_PATH], ss.get(KEY_RNA_PATH,""),
                             organism, selected_steps, params, dataset_id=did)
