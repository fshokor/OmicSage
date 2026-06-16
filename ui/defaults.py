"""
OmicSage UI — Step Definitions & Default Parameters
=====================================================
Single source of truth for:
  - Which steps exist per modality (in order)
  - Default parameter values for each step
  - Human-readable labels and descriptions

Defaults are sourced directly from real config files:
  GSE166635.yaml, GSE194122.yaml, GSE194122_cite.yaml,
  GSE194122_atac_multiome.yaml, kuppe_heart.yaml
"""

# ── Step ordering per modality ─────────────────────────────────────────────────

MODALITY_STEPS = {
    "scRNA-seq": [
        "qc", "normalize", "reduce", "cluster",
        "annotate", "deg", "gsea", "harmony", "cluster_harmony", "pseudobulk",
    ],
    "CITE-seq": [
        "normalize_adt", "doublets", "reduce_adt", "harmony_adt",
        "annotate_adt", "integration", "deg_cite", "gsea_cite",
        "protein_rna_corr", "epitope_characterisation",
    ],
    "Multiome": [
        "atac_qc", "atac_reduce", "atac_annotate",
        "multiome_integration", "multiome_deg", "multiome_grn",
    ],
    "Spatial": [
        "ingest", "qc", "reduce", "cluster",
        "deconvolve", "downstream", "impute",
    ],
}

# Steps on by default vs off by default
STEP_DEFAULTS_ON = {
    "scRNA-seq":  {"qc", "normalize", "reduce", "cluster", "annotate", "deg", "gsea"},
    "CITE-seq":   {"normalize_adt", "doublets", "reduce_adt", "annotate_adt", "integration"},
    "Multiome":   {"atac_qc", "atac_reduce", "atac_annotate", "multiome_integration", "multiome_deg"},
    "Spatial":    {"ingest", "qc", "reduce", "cluster", "deconvolve", "downstream"},
}

# Human-readable step labels
STEP_LABELS = {
    # scRNA
    "qc":                    "Quality Control",
    "normalize":             "Normalization + HVG selection",
    "reduce":                "Dimensionality reduction (PCA + UMAP)",
    "cluster":               "Clustering (Leiden)",
    "annotate":              "Cell type annotation",
    "deg":                   "Differential expression (Wilcoxon DEG)",
    "gsea":                  "Pathway analysis (GSEA/ORA)",
    "harmony":               "Batch correction (Harmony)",
    "cluster_harmony":       "Re-cluster on Harmony embedding",
    "pseudobulk":            "Pseudobulk DEG (DESeq2)",
    # CITE-seq
    "normalize_adt":         "ADT normalization (CLR)",
    "doublets":              "ADT doublet detection",
    "reduce_adt":            "ADT dimensionality reduction",
    "harmony_adt":           "ADT batch correction (Harmony)",
    "annotate_adt":          "ADT cluster annotation",
    "integration":           "RNA + ADT integration (WNN / totalVI)",
    "deg_cite":              "Differential protein expression (DPE)",
    "gsea_cite":             "Pathway analysis on CITE-seq DEGs",
    "protein_rna_corr":      "Protein–RNA correlation (Spearman)",
    "epitope_characterisation": "Epitope characterisation",
    # Multiome
    "atac_qc":               "ATAC quality control",
    "atac_reduce":           "ATAC dimensionality reduction (LSI)",
    "atac_annotate":         "ATAC cluster annotation (label transfer)",
    "multiome_integration":  "RNA + ATAC joint integration (MultiVI)",
    "multiome_deg":          "Multiome differential expression",
    "multiome_grn":          "Gene regulatory network (GRN)",
    # Spatial
    "ingest":                "Data ingestion",
    "deconvolve":            "Cell type deconvolution",
    "downstream":            "Downstream analyses (SVG, niche, CCC)",
    "impute":                "Spatial imputation (Tangram)",
}

# Steps that require --force if you want to re-run them
# (all steps support --force, this is just for UI hint)
FORCE_HINT_STEPS = {"cluster", "annotate", "atac_annotate", "deconvolve"}


# ── Default parameters per step, sourced from real configs ────────────────────

STEP_PARAM_DEFAULTS = {

    # ── scRNA-seq ─────────────────────────────────────────────────────────────
    "qc": {
        "min_genes":       200,
        "max_genes":       6000,
        "max_mt_pct":      20.0,
        "remove_doublets": True,
    },
    "normalize": {
        "batch_key":    "",           # empty = None (no per-batch HVG)
        "target_sum":   10000,
        "n_top_genes":  2000,
        "hvg_flavor":   "seurat",
    },
    "reduce": {
        "n_comps":       50,
        "n_pcs_method":  "elbow",
        "n_neighbors":   15,
    },
    "cluster": {
        "resolution_range":    [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5],
        "resolution_override": None,
    },
    "annotate": {
        "methods":    ["celltypist", "markers", "vote"],
        "leiden_col": "leiden",
        "celltypist_models": ["Immune_All_High.pkl", "Immune_All_Low.pkl"],
    },
    "deg": {
        "groupby":              "cell_type_vote",
        "method":               "wilcoxon",
        "min_logfc":            0.25,
        "max_pval_adj":         0.05,
        "n_genes":              500,
        "exclude_gene_prefixes": [],
    },
    "gsea": {
        "gene_sets":            ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"],
        "organism":             "human",
        "min_logfc":            0.25,
        "max_pval_adj":         0.05,
        "min_genes":            5,
        "exclude_gene_prefixes": [],
    },
    "harmony": {
        "batch_key":     "batch",
        "n_pcs":         50,
        "n_neighbors":   15,
        "umap_min_dist": 0.3,
        "random_state":  0,
    },
    "cluster_harmony": {
        "resolution_range": [0.2, 0.4, 0.6, 0.8, 1.0],
    },
    "pseudobulk": {
        "groupby":      "cell_type_vote",
        "donor_key":    "batch",
        "counts_layer": "counts",
        "min_cells":    10,
        "min_samples":  3,
        "min_logfc":    0.25,
        "max_pval_adj": 0.05,
        "exclude_gene_prefixes": [],
    },

    # ── CITE-seq ──────────────────────────────────────────────────────────────
    "normalize_adt": {
        "clr_axis": 0,
    },
    "doublets": {
        "threshold":        2.5,
        "filter_doublets":  False,
    },
    "reduce_adt": {
        "n_comps":      50,
        "n_pcs":        20,
        "n_neighbors":  15,
    },
    "harmony_adt": {
        "batch_key":    "batch",
        "n_pcs":        20,
        "n_neighbors":  15,
        "random_state": 0,
    },
    "annotate_adt": {
        "resolution":    0.1,
        "n_iterations":  2,
    },
    "integration": {
        "method":       "both",
        "batch_key":    "batch",
        "n_factors":    15,
        "max_epochs":   10,
        "random_state": 0,
    },
    "deg_cite": {
        "groupby":      "adt_celltype_manual",
        "method":       "wilcoxon",
        "min_logfc":    0.25,
        "max_pval_adj": 0.05,
        "n_genes":      200,
    },
    "gsea_cite": {
        "organism":     "human",
        "min_logfc":    0.25,
        "max_pval_adj": 0.05,
        "min_genes":    5,
        "direction":    "up",
    },
    "protein_rna_corr": {
        "method":    "spearman",
        "min_cells": 30,
    },
    "epitope_characterisation": {
        "groupby": "adt_celltype_manual",
        "preset":  "bmmc",
    },

    # ── Multiome ──────────────────────────────────────────────────────────────
    "atac_qc": {
        "min_peaks":             750,
        "max_peaks":             500000,
        "min_peak_counts":       1500,
        "max_peak_counts":       100000,
        "max_nucleosome_signal": 2.0,
        "min_cells":             15,
    },
    "atac_reduce": {
        "n_components":    50,
        "n_neighbors":     15,
        "leiden_resolution": 0.5,
        "random_state":    0,
    },
    "atac_annotate": {
        "promoter_upstream_bp": 2000,
        "rna_label_key":        "cell_type_vote",
    },
    "multiome_integration": {
        "method":       "multivi",
        "batch_key":    "batch",
        "n_latent":     20,
        "max_epochs":   10,
        "random_state": 0,
        "n_top_peaks":  20000,
    },
    "multiome_deg": {
        "groupby":    "atac_celltype",
        "method":     "wilcoxon",
        "min_logfc":  0.25,
        "max_pval_adj": 0.05,
        "n_genes":    200,
        "exclude_gene_prefixes": ["MT-", "RPL", "RPS"],
    },
    "multiome_grn": {
        "motif_db":    "jaspar",
        "groupby":     "atac_celltype",
        "n_top_peaks": 500,
        "min_cells":   10,
        "random_state": 0,
    },

    # ── Spatial ───────────────────────────────────────────────────────────────
    "ingest": {
        "spatial_type": "h5ad",
        "load_images":  True,
        "library_key":  "sample_name",
    },
    # qc, reduce, cluster reuse scRNA defaults with spatial-appropriate values
    "deconvolve": {
        "method":         "nnls",
        "cell_type_key":  "cell_type_original",
        "layer_ref":      "counts",
        "n_jobs":         4,
        "target_sum":     10000,
    },
    "downstream": {
        "run_region_clustering":  True,
        "run_celltype_expression": True,
        "run_co_occurrence":      True,
        "run_nhood_enrichment":   True,
        "run_ligrec":             True,
        "run_svg_gsea":           True,
        "ligrec_organism":        "human",
        "n_jobs":                 4,
    },
    "impute": {
        "method":         "tangram",
        "n_top_genes":    2000,
        "tangram_mode":   "clusters",
        "device":         "cpu",
    },
}

# Spatial QC defaults differ slightly from scRNA
SPATIAL_QC_DEFAULTS = {
    "min_counts":  500,
    "max_counts":  100000,
    "min_genes":   200,
    "max_genes":   10000,
    "max_mt_pct":  20.0,
    "mt_prefix":   "MT-",
}

SPATIAL_REDUCE_DEFAULTS = {
    "n_top_genes":  3000,
    "n_comps":      50,
    "n_neighbors":  6,
    "target_sum":   10000,
}

SPATIAL_CLUSTER_DEFAULTS = {
    "resolution":  0.5,
    "n_neighbors": 15,
    "n_pcs":       30,
    "random_state": 0,
}


# ── Runner script per modality ─────────────────────────────────────────────────
MODALITY_RUNNER = {
    "scRNA-seq": "run_scrna_pipeline.py",
    "CITE-seq":  "run_cite_pipeline.py",
    "Multiome":  "run_multiome_pipeline.py",
    "Spatial":   "run_spatial_pipeline.py",
}

# ── Organism options ───────────────────────────────────────────────────────────
ORGANISMS = ["human", "mouse"]

# ── Nextflow modality names (passed as --modality to main.nf) ─────────────────
MODALITY_NF_NAME = {
    "scRNA-seq": "scrna",
    "CITE-seq":  "cite",
    "Multiome":  "multiome",
    "Spatial":   "spatial",
}

# Modalities with Nextflow workflows currently implemented
NF_IMPLEMENTED = {"scRNA-seq"}
