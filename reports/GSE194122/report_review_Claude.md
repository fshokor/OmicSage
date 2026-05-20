# AI Interpretation Report
**Dataset:** BMMC CITE-seq (NeurIPS 2021)
**Date:** 2026-05-20
**Reviewer:** Manual Review Mode (D1)
**Model used:** Claude Sonnet 4.6
**Tissue:** Bone marrow (BMMC — bone marrow mononuclear cells)
**Disease:** Healthy (no disease; dataset designed for multimodal integration benchmarking)

---

## 1. Study Overview

- **Tissue:** Bone marrow mononuclear cells (BMMC)
- **Disease:** None — healthy donors
- **Species:** Homo sapiens
- **Conditions:** Single condition (healthy); nested batch structure across 4 collection sites
- **N cells (post-QC):** 21,778
- **N donors:** 12 (donors confirmed in report metadata; batches labelled s1d1 … s4d9)
- **Batch key:** `batch` (12 batch labels; site × donor combinations)
- **Biological question:** Inferred by reviewer — characterise the transcriptional landscape of healthy human bone marrow hematopoiesis across immune compartments using single-cell multimodal data; validate and extend cell type annotations across the erythroid, myeloid, and lymphoid lineages.

---

## 2. QC Assessment

- **Thresholds:** min_genes = 200, max_genes = 2,500, max_mt_pct = 5.0%, doublet removal = True (Scrublet)
- **Cells removed:** 75.9% — **above typical** for bone marrow CITE-seq (90,261 input → 21,778 kept; 68,483 removed)
- **Doublet detection:** Run — Scrublet used; 187 doublets removed; Scrublet score distribution reported
- **Overall QC assessment:** **Review recommended** — see flags below

| Severity | Flag | Suggestion |
|----------|------|------------|
| warning | 75.9% cells removed is substantially above the typical 20–40% range for 10X bone marrow CITE-seq | Inspect per-filter cell loss breakdown: max_genes = 2,500 removed 8,710 cells; max_mt_pct = 5% removed 64,699 cells. The 5% MT ceiling is the primary driver. |
| critical | max_mt_pct = 5% is highly stringent for bone marrow, where erythroid precursors and metabolically active HSCs routinely exceed 5% MT% legitimately. This threshold likely discards a large fraction of biologically real cells, particularly from the erythroid compartment. Published BMMC benchmarks (e.g. the NeurIPS 2021 competition dataset itself) often use 10–15% as the ceiling. | Consider re-running QC with max_mt_pct = 10% or 15% and compare cell yield and cluster composition. The ground-truth MT% Pearson r = 0.9998 confirms MT% calculation is technically correct; the issue is the threshold, not the metric. |
| info | min_genes = 200 is at the low end of the typical 200–500 range for 10X droplet data (removed only 1,124 cells); appropriate given dataset complexity | No action required |
| info | max_genes = 2,500 removed 8,710 cells; this may exclude doublets but could also remove large, transcriptionally active cells (e.g. monocytes, plasma cells) | Cross-reference with doublet scores; if Scrublet scores are low for removed high-gene cells, raising max_genes to 4,000–5,000 may recover informative cells |
| info | 187 doublets removed by Scrublet — a very small fraction (0.2% of kept cells), plausible for a well-processed sample | No action required |
| info | Pre-QC median MT% = 6.34%; cutting at 5% therefore removes cells near the population median for some lineages | Corroborates the critical flag above |

**Note on cell loss:** The report confirms that Scrublet ground-truth MT% validation passed with Pearson r = 0.9998, confirming pipeline integrity. The large cell loss is driven overwhelmingly by the max_mt_pct = 5% filter (64,699 cells) rather than gene count filters.

---

## 3. Clustering Assessment

- **Resolution:** 0.6 — **appropriate** for 21,778 cells (recommended range 0.5–0.8 for this cell count)
- **Recommended range:** 0.5–0.8 for medium datasets (5,000–30,000 cells)
- **N clusters:** 16 — **reasonable** for healthy bone marrow of this size; the tissue is expected to contain 8–15 major lineages with sub-populations
- **Silhouette score:** 0.1580 — **poor** (below 0.2 threshold for adequacy); note silhouette scores drop substantially with increasing resolution: 0.3833 at res=0.2, 0.3579 at res=0.4, 0.1580 at res=0.6 (selected). The steep drop between 0.4 and 0.6 (Δ = +6 clusters, silhouette drops 0.20) warrants attention.
- **Selection reason:** stability_plateau — the pipeline chose res=0.6 because the cluster count stabilised relative to res=0.8 (similar Δ pattern). The biological plausibility of 16 clusters is consistent with the expected diversity of BMMC.
- **Sub-clustering candidates:**
  - **Cluster 10 (n=3,763, CD4+ T cells):** Largest cluster. CD4+ T cell populations in bone marrow include naive, central memory, effector memory, and Treg subsets. Sub-clustering is strongly recommended.
  - **Cluster 0 (n=2,275, Late erythroid):** Erythroid differentiation is a continuum; late erythroid may span reticulocytes, normoblasts, and mature erythroblasts. Sub-clustering or trajectory analysis recommended.
  - **Cluster 3/4 (n=1,085 and n=1,753, Mid erythroid):** Two clusters with identical annotation (Mid erythroid); may represent the same population artificially split or two distinct developmental stages. Marker genes should be compared directly.
  - **Cluster 7/8 (n=1,660 and n=1,717, Classical monocytes):** Both annotated as Classical monocytes with confidence 1.00; likely represent transcriptional sub-states (inflammatory vs patrolling-like) worth distinguishing.
- **Literature context:** Published scRNA-seq studies of healthy BMMC at comparable cell counts typically use Leiden resolutions of 0.4–0.8 and report 12–20 clusters. The Human Cell Atlas bone marrow reference (Bandyopadhyay et al. 2024, PMID:38559168 — bioRxiv preprint) profiled 29,325 enriched non-hematopoietic cells and identified 9 transcriptionally distinct non-hematopoietic subtypes alone, underscoring that 16 clusters is conservative for the full hematopoietic compartment.

---

## 4. Cell Type Annotation

### 4a. Cluster-level predictions

| Cluster | N cells | Predicted type | Confidence | Supporting markers | Alternatives | SingleR agreement |
|---------|---------|----------------|------------|-------------------|--------------|-------------------|
| 0 | 2,275 | Late erythroid | **High** | HBB (log2FC=8.91 Up), RPLP1 (Down), TPT1 (Down) — terminal erythroid pattern | Reticulocyte (very similar) | Yes — SingleR: Erythroid cells; CellTypist: Late erythroid (1.00 confidence) |
| 1 | 997 | Late erythroid | **High** | Same marker profile as Cluster 0 | Reticulocyte | Yes — all 4 methods agree; confidence 1.00 |
| 2 | 1,238 | Late erythroid | **High** | Same marker profile | Reticulocyte | Yes — confidence 1.00 |
| 3 | 1,085 | Mid erythroid | **High** | PRDX2 (Up), AHSP (Up), HBM (Up), HEMGN (Up), CA2 (Up) | Early erythroid precursor | Yes — confidence 1.00 |
| 4 | 1,753 | Mid erythroid | **High** | Same as Cluster 3 | Early erythroid precursor | Yes — confidence 1.00 |
| 5 | 484 | Plasma cells | **Medium** | MZB1 (log2FC=6.37 Up), TXNDC5 (Up), SSR4 (Up), FKBP11 (Up), SEC11C (Up) — ER secretory pathway | Plasmablast | Partial: CellTypist = Plasma cells, SingleR = B cells; ScType = Naive B cells. Marker evidence strongly favours plasma cells; SingleR and ScType calls are inconsistent. Confidence 0.50. |
| 6 | 1,262 | Small pre-B cells | **High** | CD79B (log2FC=7.74 Up), VPREB3 (log2FC=8.83 Up), IGHM (log2FC=6.97 Up) | Pro-B cell, immature B cell | Yes — CellTypist Fine, SingleR, and consensus agree; ScType (Megakaryocyte) is discordant and should be disregarded |
| 7 | 1,660 | Classical monocytes | **High** | LYZ (log2FC=8.68 Up), S100A9 (log2FC=9.49 Up), FCN1 (log2FC=8.47 Up) | Non-classical monocyte, DC | Yes — all 4 methods agree; confidence 1.00 |
| 8 | 1,717 | Classical monocytes | **High** | Same as Cluster 7 (report does not provide cluster-specific marker breakdowns beyond DEG top 5) | Intermediate monocyte | Yes — confidence 1.00. Two classical monocyte clusters likely reflect sub-states. |
| 9 | 1,099 | pDC | **Medium** | IRF8 (log2FC=6.89 Up), GZMB (log2FC=7.27 Up), JCHAIN (log2FC=7.02 Up), CCDC50 (Up), PPP1R14B (Up) | Conventional DC (cDC2), B cell precursor | Partial: CellTypist = pDC, Markers = pDC, but ScType = Megakaryocyte and SingleR = Dendritic cells (non-specific). GZMB is characteristic of activated pDCs. IRF8 is a master regulator of pDC development. Confidence 0.50. |
| 10 | 3,763 | CD4+ T cells | **Medium** | LTB (log2FC=5.48 Up); ribosomal genes RPL30, RPL13, RPS12, RPL34 dominate top DEGs | Naive T cell, Treg, NKT | Partial: CellTypist = CD4+ T cells, SingleR = CD4+ T cells, Markers = T cell; ScType = Megakaryocyte (discordant). LTB is a naive/central memory T cell marker. Ribosomal gene dominance in top DEGs limits marker interpretation. Confidence 0.50. |
| 11 | 994 | CD8+ T cells | **Medium** | CCL5 (log2FC=6.76 Up), NKG7 (log2FC=5.87 Up), CST7 (log2FC=4.89 Up), GZMA (log2FC=5.12 Up) | NK cell, NKT cell | Partial: CellTypist = Tem/Temra cytotoxic T cells, SingleR = CD8+ T cells; Markers = T cell. CCL5+NKG7+GZMA is a canonical effector/memory CD8+ T cell signature. Confidence 0.50. |
| 12 | 1,480 | CD4+ T cells | **Medium** | Same pattern as Cluster 10 (ribosomal gene-dominated) | Naive T helper, Treg | Partial: same as Cluster 10. Two CD4+ T clusters likely reflect activation states. Confidence 0.50. |
| 13 | 332 | Naive B cells | **High** | MS4A1 (log2FC=7.97 Up), CD74 (Up), HLA-DRA (Up), CD79A (log2FC=5.78 Up), HLA-DPB1 (Up) | Memory B cell, transitional B cell | Yes — CellTypist, Markers, SingleR, consensus all agree; confidence 0.75 |
| 14 | 506 | HSC/MPP | **Medium** | SNHG29 (log2FC=3.52 Up); ribosomal: RPS24, RPLP0, RPLP1, RPL12 | CMP, GMP, MEP, CLP | Partial: CellTypist = HSC/MPP, Markers = Progenitor, SingleR = HSCs; ScType = Progenitor cells. SNHG29 is a novel lncRNA marker (see DEG notes). Ribosomal dominance is expected in highly proliferative progenitors. Confidence 0.50. |
| 15 | 1,133 | CD16+ NK cells | **Medium** | NKG7 (log2FC=8.06 Up), GNLY (log2FC=9.45 Up), GZMA (log2FC=7.00 Up), CST7 (log2FC=6.49 Up) | NKT cell, CD8+ T cell, ILC1 | Partial: CellTypist = CD16+ NK cells / ILC, Markers = NK_ILC, SingleR = NK cells, ScType = Natural killer cells. GNLY+NKG7+GZMA is the canonical NK cytotoxic signature. Confidence 0.50 reflects CellTypist's "ILC" broadness. |

**ScType systematic discordance note:** ScType called "Megakaryocyte" for clusters 6, 9, 10, 12 — cells annotated as B cells, pDCs, and T cells by all other methods. This is a systematic mis-classification by ScType (likely due to database mismatch or marker set conflict) and should be disregarded for these clusters. The consensus vote correctly overrides ScType for these cases.

### 4b. Tissue-specific marker sets

| Cell type | Canonical markers | BMMC-specific | Distinguishing markers | Counter-markers | Confidence |
|-----------|------------------|---------------|------------------------|-----------------|------------|
| Late erythroid (Cls 0–2) | HBB, HBA1, HBA2, GYPC, ALAS2 | HBB (log2FC=8.91 in Cluster 0 DEGs), HBM, HBD, BLVRB, PRDX2 — consistent with published BMMC CITE-seq data | SLC4A1, CA1 (distinguish from mid erythroid); GYPC, GYPB (distinguish from nucleated erythroblasts) | PTMA, TPT1, RPLP1 (high expression argues against terminal erythroid and suggests earlier stage) | High |
| Mid erythroid (Cls 3–4) | GATA1, KLF1, AHSP, HEMGN | HBM (log2FC=6.59), AHSP (log2FC=5.92), HEMGN (log2FC=5.20), CA2 (log2FC=6.26), PRDX2 (log2FC=5.10) — all in top 5 DEGs | AHSP (distinguishes from HSC/erythroblast), high cell cycle gene expression (TOP2A in GSEA hits) relative to late erythroid | Mature globin dominance (HBB >> HBM) would suggest late erythroid rather than mid | High |
| Plasma cells (Cls 5) | MZB1, JCHAIN, CD38, SDC1, IRF4, BLIMP1 (PRDM1) | MZB1 (log2FC=6.37), TXNDC5 (log2FC=6.70), SSR4 (log2FC=4.88), FKBP11 (log2FC=5.56) — secretory/ER stress markers typical of active antibody-secreting cells in BM | JCHAIN (distinguishes IgA/IgM-secreting plasma cells from IgG); MZB1 (plasma cell vs plasmablast) | MS4A1 (CD20; high expression would suggest B cell rather than plasma cell); CD79A (normally downregulated in terminally differentiated plasma cells) | Medium |
| Small pre-B cells (Cls 6) | CD79B, VPREB3, RAG1, VPREB1, IGHM, CD10 (MME) | VPREB3 (log2FC=8.83), CD79B (log2FC=7.74), IGHM (log2FC=6.97), RAG1 (pseudobulk log2FC=6.64 — #1 DESeq2 hit); HLA-DRA expression confirms antigen-presenting capability | RAG1 (unique to B-cell development stage; distinguishes pre-B from mature B); VPREB3 (surrogate light chain, absent in mature B cells) | MS4A1 (CD20; should be absent or low in pre-B cells) | High |
| Classical monocytes (Cls 7–8) | LYZ, CD14, S100A8, S100A9, FCN1, VCAN, SELL | S100A9 (log2FC=9.49), LYZ (log2FC=8.68), FCN1 (log2FC=8.47), S100A6 (log2FC=4.28) — top Wilcoxon DEGs; PILRA (log2FC=6.64), FPR1 (log2FC=7.65), CD163 (log2FC=6.41) in pseudobulk | SELL (CD62L, distinguishes classical from non-classical monocytes); ITGAM (CD11b); FCGR3A absence distinguishes classical from CD16+ non-classical | FCGR3A/CD16 (high expression would indicate non-classical monocyte); LILRA4 (pDC marker) | High |
| pDC (Cls 9) | LILRA4, BST2, CLEC4C, IRF7, SPIB, GZMB (activated) | IRF8 (log2FC=6.89), GZMB (log2FC=7.27), JCHAIN (log2FC=7.02), CCDC50 (log2FC=6.90) — IRF8 is master TF for pDC commitment; GZMB is characteristic of activated/mature BM pDCs; JCHAIN suggests immunoglobulin-class production capacity | LILRA4 and BST2 (distinguish pDC from cDC1/cDC2); SPIB (distinguishes pDC from plasmablast despite shared JCHAIN) | CD3D, CD3E (T cell markers; would argue against pDC), CD19 (B cell marker) | Medium |
| CD4+ T cells (Cls 10, 12) | CD4, IL7R, LTB, LEF1, TCF7 | LTB (log2FC=5.48) is a naive/central memory CD4+ T cell marker in BM; NOSIP (pseudobulk #1 DEG, log2FC=2.41); ZEB2 downregulated (log2FC=−5.33 pseudobulk) — ZEB2 repression is characteristic of naive T cells | IL7R, CCR7 (distinguish naive from effector memory CD4+); FOXP3 (would indicate Treg) | FOXP3 (Treg), CD8A (cytotoxic), GZMB (NK/CTL) | Medium |
| CD8+ T cells (Cls 11) | CD8A, CD8B, NKG7, CCL5, GZMA, GZMK | CCL5 (log2FC=6.76), NKG7 (log2FC=5.87), CST7 (log2FC=4.89), GZMA (log2FC=5.12); GZMK (pseudobulk #3, log2FC=5.57) — effector/memory CD8+ signature | GZMK (distinguishes effector memory from terminally exhausted Temra); CCL5 (canonical BM-resident memory CD8+) | CD4, FOXP3, CD19 | Medium |
| Naive B cells (Cls 13) | MS4A1, CD79A, CD74, HLA-DRA, HLA-DPB1, IGHD | MS4A1 (log2FC=7.97), CD79A (log2FC=5.78), HLA-DRA (log2FC=4.81), CD74 (log2FC=4.54), HLA-DPB1 (log2FC=4.36) — canonical mature naive B cell | IGHD co-expression (distinguishes naive from IgG-switched memory B cells); TCL1A (naive B marker in BM) | PRDM1/BLIMP1 (plasma cell TF), CD38 high (activated state) | High |
| HSC/MPP (Cls 14) | CD34, SPINK2, MLLT3, AVP, MECOM, HOXA9 | SNHG29 (log2FC=3.52 — #1 Wilcoxon DEG; log2FC=2.67 pseudobulk); PRTN3 (pseudobulk #2, log2FC=8.67) — neutrophil elastase-family serine protease, known in early myeloid progenitors; CLEC11A (pseudobulk #3, log2FC=4.57) — bone marrow stroma/osteoblast-associated growth factor expressed in HSCs | SPINK2 (distinguishes HSC from early committed progenitors — literature-supported); MLLT3 (HSC self-renewal, distinguishes from MPP) | HBB (erythroid commitment), CD3D (T cell), CD19 (B cell), LYZ (myeloid commitment) | Medium |
| CD16+ NK cells (Cls 15) | NKG7, GNLY, GZMA, GZMB, FCER1G, KLRD1, FCGR3A | NKG7 (log2FC=8.06), GNLY (log2FC=9.45), GZMA (log2FC=7.00), CST7 (log2FC=6.49), B2M (log2FC=2.96); MVD (pseudobulk #1, log2FC=2.53) | KLRD1/CD94 and NKG2 family (distinguish CD16+ NK from ILC1/ILC3); FCGR3A/CD16 expression (distinguishes CD16+ from CD56bright NK) | CD3D, CD3E, CD8A (T cell markers; their expression would challenge NK annotation) | Medium |

---

## 5. DEG Validation

This section covers the Wilcoxon one-vs-rest DEG results (grouped by `cell_type_vote`). Top 5 genes per group are provided in the report. For pseudobulk DESeq2 results, see notes below each group.

**Comparison: All groups (one cell type vs. rest)**

| Gene | log2FC | adj p-value | Group | Classification | Rationale |
|------|--------|-------------|-------|---------------|-----------|
| NKG7 | 8.059 | 0.00e+00 | CD16+ NK | Expected | Canonical NK cytotoxic granule protein; top marker in multiple NK scRNA-seq studies |
| GNLY | 9.448 | 0.00e+00 | CD16+ NK | Expected | Granulysin; NK and CTL effector molecule, highest log2FC among top NK DEGs |
| GZMA | 6.997 | 0.00e+00 | CD16+ NK | Expected | Granzyme A; canonical cytotoxic lymphocyte marker |
| CST7 | 6.487 | 0.00e+00 | CD16+ NK | Expected | Cystatin F; marks cytotoxic NK and CD8+ T cells |
| B2M | 2.955 | 0.00e+00 | CD16+ NK | Expected | MHC-I component; expected upregulation in antigen-presenting/killing cells; modest FC |
| RPL30 | 2.858 | 0.00e+00 | CD4+ T | Artefact flag | Ribosomal protein; dominance of RPL/RPS genes in top 5 CD4+ T DEGs (RPL30, RPL13, RPS12, RPL34) suggests ribosomal signal is not corrected for. CD4+ T cells in bone marrow are relatively quiescent. See cross-module note. |
| LTB | 5.476 | 0.00e+00 | CD4+ T | Expected | Lymphotoxin B; canonical naive/resting T cell and B cell marker in BM |
| CCL5 | 6.758 | 0.00e+00 | CD8+ T | Expected | RANTES; secreted by effector memory CD8+ T cells; prominent BM-resident T cell marker |
| NKG7 | 5.869 | 0.00e+00 | CD8+ T | Expected | Shared with NK cells; consistent with cytotoxic CD8+ T cell identity |
| GZMA | 5.116 | 0.00e+00 | CD8+ T | Expected | Granzyme A; effector CD8+ marker |
| LYZ | 8.683 | 0.00e+00 | Classical mono | Expected | Lysozyme; canonical myeloid/monocyte marker |
| S100A9 | 9.486 | 0.00e+00 | Classical mono | Expected | S100 calcium-binding protein; highest log2FC monocyte marker; inflammatory calprotectin complex with S100A8 |
| FCN1 | 8.471 | 0.00e+00 | Classical mono | Expected | Ficolin-1; pattern-recognition receptor; distinguishes classical monocytes in BM |
| SNHG29 | 3.523 | 5.78e-236 | HSC/MPP | **Unexpected — discovery** | Small nucleolar RNA host gene 29; lncRNA; poorly characterised in HSC biology. Top Wilcoxon and pseudobulk DESeq2 hit (padj=4.73e-14). See Discovery Highlights. |
| HBB | 8.913 | 0.00e+00 | Late erythroid | Expected | Beta-globin; hallmark of terminal erythroid differentiation |
| RPLP1 | −5.739 | 0.00e+00 | Late erythroid (down) | Expected | Ribosomal protein; downregulation in terminal erythroid consistent with nuclear condensation and translational shutdown |
| PRDX2 | 5.096 | 0.00e+00 | Mid erythroid | Expected | Peroxiredoxin 2; antioxidant; highly expressed in erythroid precursors handling haem-associated oxidative stress |
| AHSP | 5.921 | 0.00e+00 | Mid erythroid | Expected | Alpha-haemoglobin stabilising protein; canonical mid-erythroid marker |
| HBM | 6.586 | 0.00e+00 | Mid erythroid | Expected | Mu-globin chain; expressed in foetal/adult erythroblasts |
| MS4A1 | 7.968 | 7.68e-167 | Naive B | Expected | CD20; canonical B cell surface antigen |
| CD79A | 5.776 | 3.82e-153 | Naive B | Expected | Igα; B cell receptor signalling component |
| MZB1 | 6.366 | 1.91e-276 | Plasma cells | Expected | Marginal zone B and B1 cell-specific protein; plasma cell/secretory B cell marker |
| TXNDC5 | 6.703 | 2.71e-267 | Plasma cells | Expected | Thioredoxin domain protein; ER stress/protein folding; canonical secretory plasma cell marker |
| CD79B | 7.743 | 0.00e+00 | Small pre-B | Expected | Igβ; B cell signalling; marks pre-B and immature B stages |
| VPREB3 | 8.826 | 0.00e+00 | Small pre-B | Expected | Surrogate light chain component; specific to pre-B cell checkpoint |
| IGHM | 6.968 | 0.00e+00 | Small pre-B | Expected | IgM heavy chain; marks IgM-expressing immature B cells |
| PPP1R14B | 5.726 | 0.00e+00 | pDC | Unexpected | Protein phosphatase 1 regulatory subunit 14B; poorly characterised in pDC context. See Discovery Highlights. |
| GZMB | 7.269 | 0.00e+00 | pDC | Expected | Granzyme B; established marker of activated/mature pDCs; literature-supported role in innate antiviral response |
| IRF8 | 6.890 | 0.00e+00 | pDC | Expected | Master TF for pDC lineage commitment; essential for pDC development |
| JCHAIN | 7.023 | 0.00e+00 | pDC | Unexpected | Joining chain for IgA/IgM; not a canonical pDC marker. Also top DEG for Plasma cells (Wilcoxon top 5). May indicate technical bleed-through from adjacent plasma cell cluster or genuine pDC-Ig chain co-expression. See cross-module note. |

**Pseudobulk DESeq2 notable flags:**
- **HBB downregulated** in CD16+ NK (log2FC=−10.10, padj=1.13e-14), CD8+ T (log2FC=−9.75, padj=1.88e-14), and Naive B (log2FC=−10.20, padj=1.45e-16): HBB appearing as a top downregulated DEG in immune cell types is expected in a one-vs-rest pseudobulk test where erythroid cells dominate the "rest" group. This is a statistical artefact of the test design, not a biological finding.
- **ERBB2 upregulated in CD16+ NK** (pseudobulk log2FC=3.90, padj=4.18e-15): HER2/ERBB2 upregulation in NK cells is unexpected. Could reflect a small contaminating cell population, or genuine ERBB2 expression in a subset of BM NK cells. Warrants validation (see Discovery Highlights).
- **LINC00926 upregulated in Naive B** (pseudobulk log2FC=6.45, padj=2.83e-30): Long non-coding RNA, known B cell-enriched lncRNA in blood; consistent with B cell identity.
- **RAG1 upregulated in Small pre-B** (pseudobulk log2FC=6.64, padj=9.51e-69): Recombination activating gene 1; highly specific to the pre-B cell V(D)J recombination stage. Strongest possible support for the pre-B annotation.
- **ZFAT upregulated in pDC** (pseudobulk log2FC=5.26, padj=6.45e-78): Zinc finger and AT hook domain; transcription factor with roles in lymphocyte development; not well characterised in pDCs specifically — see Discovery Highlights.

**Discovery highlights:**
- **SNHG29** (HSC/MPP, Wilcoxon log2FC=3.52, adj-p=5.78e-236; pseudobulk log2FC=2.67, padj=4.73e-14): This lncRNA is the top DEG for HSC/MPP by both methods. SNHG29 is poorly characterised in haematopoietic stem cell biology. Its consistent top ranking across two orthogonal statistical tests in a well-powered healthy BM dataset makes it a high-priority candidate for functional validation as a novel HSC/MPP marker or regulatory lncRNA.
- **ERBB2** (CD16+ NK, pseudobulk log2FC=3.90, padj=4.18e-15): HER2 upregulation in bone marrow NK cells is notable. If not artefactual, this could represent an NK cell subset with relevance to HER2+ tumour surveillance. Requires validation by flow cytometry and/or CITE-seq antibody signal.
- **PPP1R14B** (pDC, Wilcoxon log2FC=5.73, adj-p=0.00): Protein phosphatase 1 inhibitor with poorly defined immune function. Top pDC marker by statistical significance. May represent a novel pDC-specific regulatory gene.
- **ZFAT** (pDC, pseudobulk log2FC=5.26, padj=6.45e-78): Highest-ranked pseudobulk pDC marker. ZFAT has known roles in T cell apoptosis and lymphocyte survival but has not been established as a pDC marker. Its top ranking here may identify a new regulatory axis in pDC biology.
- **PARK7** (pDC, pseudobulk log2FC=1.65, padj=1.09e-31): DJ-1/PARK7 is linked to Parkinson's disease and oxidative stress response. Expression in pDCs is not previously described and may reflect a role in pDC survival under reactive oxygen species conditions.

---

## 6. Literature Links

PubMed search was performed via the PubMed MCP tool. Attribution required by PubMed tool terms. Note that many specific gene–cell type combinations in healthy BM did not return PubMed hits (PubMed queries returned 0 results for several gene/context combinations). Manual search strings are provided where the tool returned no hits.

| Gene | PMID / Source | Title | Context (one sentence) |
|------|--------------|-------|------------------------|
| NKG7, GNLY (NK cells) | PMID:39074350 | "Single-Cell Analysis Reveals That CD47 mRNA Expression Correlates with Immune Cell Activation, Antiviral ISGs, and Cytotoxicity" (Cham et al. 2024, doi:10.33594/000000715) | NK and CD8+ T cells express high levels of NKG7, granzyme B, perforin, and granulysin as canonical cytotoxic markers in scRNA-seq data from healthy controls. |
| HBB, HBM (erythroid) | PMID:38559168 | "Mapping the Cellular Biogeography of Human Bone Marrow Niches Using Single-Cell Transcriptomics and Proteomic Imaging" (Bandyopadhyay et al. 2024, doi:10.1101/2024.03.14.585083) | A spatially-resolved scRNA-seq atlas of human bone marrow identifies erythroid lineage cells and their niche associations; globin gene expression (HBB, HBM) is expected in this compartment. |
| BMDB reference atlas | PMID:41334179 | "BMDB: An integrated database and web platform for single-cell transcriptomic profiling of bone marrow microenvironment" (Chen et al. 2025, doi:10.1016/j.csbj.2025.11.028) | BMDB integrates 435,682 single cells from 224 human/mouse BM samples and provides reference-guided annotation and comparative analysis — directly applicable to this dataset for cross-study validation. |

**Manual PubMed search strings for genes where the tool returned no hits:**

- SNHG29, HSC/MPP: `SNHG29 AND ("hematopoietic stem" OR HSC OR MPP) AND ("single cell" OR scRNA-seq)`
- ERBB2, NK cells: `ERBB2 AND "natural killer" AND ("bone marrow" OR hematopoiesis)`
- ZFAT, pDC: `ZFAT AND ("plasmacytoid dendritic" OR pDC) AND ("bone marrow" OR hematopoiesis)`
- PPP1R14B, pDC: `PPP1R14B AND "plasmacytoid dendritic" AND ("single cell" OR scRNA-seq)`
- PARK7, pDC: `PARK7 AND "dendritic cell" AND "bone marrow"`
- RAG1, pre-B: `RAG1 AND "pre-B" AND "bone marrow" AND scRNA-seq`
- JCHAIN, pDC: `JCHAIN AND "plasmacytoid dendritic" AND ("innate" OR "single cell")`
- FCN1, monocyte: `FCN1 AND monocyte AND "bone marrow" AND scRNA-seq`
- CLEC11A, HSC: `CLEC11A AND "hematopoietic stem" AND "bone marrow"`
- PRTN3, MPP: `PRTN3 AND (MPP OR "myeloid progenitor") AND ("single cell" OR scRNA-seq)`

*PMIDs and titles only. No abstract text. Based on PubMed articles retrieved via PubMed MCP tool (PubMed attribution per tool requirements).*

---

## 7. GSEA Coherence

Gene sets used: GO_Biological_Process_2023, KEGG_2021_Human, Reactome_2022. Analysis run with gseapy 1.1.3. Upregulated genes only. Total significant pathways: 4,048 (adj.p≤0.05).

**Consistent pathways (no action required):**
- **CD16+ NK cells — Immune System R-HSA-168256** (121 genes, adj.p=1.48e-29): consistent with NK cell immune effector identity
- **CD16+ NK cells — Natural killer cell mediated cytotoxicity** (25 genes, adj.p=1.93e-17): directly confirms NK annotation
- **CD16+ NK cells — Innate Immune System R-HSA-168249** (75 genes, adj.p=1.27e-20): consistent with innate lymphocyte compartment
- **Classical monocytes — Immune System R-HSA-168256** (196 genes, adj.p=1.22e-70): consistent with myeloid innate immune identity
- **Classical monocytes — Innate Immune System R-HSA-168249** (146 genes, adj.p=3.84e-69): consistent with monocyte annotation
- **Classical monocytes — Neutrophil Degranulation R-HSA-6798695** (103 genes, adj.p=1.48e-62): consistent with classical monocyte degranulation capacity; FCN1, S100A9, LYZ genes all included
- **Classical monocytes — Phagosome** (37 genes, adj.p=1.11e-25): consistent with monocyte phagocytic function
- **Classical monocytes — Cytokine Signaling R-HSA-1280215** (66 genes, adj.p=3.52e-18): consistent with monocyte inflammatory signalling
- **Plasma cells — Protein processing in endoplasmic reticulum** (42 genes, adj.p=2.20e-36): expected for high-secretory plasma cells; SSR4, SEC61A1, HSP90B1 included
- **Plasma cells — Asparagine N-linked Glycosylation R-HSA-446203** (30 genes, adj.p=1.99e-14): consistent with immunoglobulin glycosylation
- **Plasma cells — SRP-dependent Cotranslational Protein Targeting R-HSA-1799339** (20 genes, adj.p=3.76e-13): consistent with heavy chain co-translational ER import
- **Late erythroid — Hydrogen Peroxide Catabolic Process** (6 genes, adj.p=3.89e-07): consistent — PRDX2, HBB, HBA1 handle haem-associated ROS in maturing erythrocytes
- **Late erythroid — Porphyrin and chlorophyll metabolism** (4 genes, adj.p=3.36e-04): consistent with haem biosynthesis (ALAS2, FECH, HMBS)
- **Late erythroid — Erythrocytes Take Up Oxygen And Release Carbon Dioxide R-HSA-1247673** (3 genes, adj.p=1.56e-03): directly confirms erythroid annotation
- **Small pre-B — Translation/ribosome pathways**: consistent with highly proliferative B cell precursors

**Unexpected pathways:**

| Pathway | Comparison / group | Adj p-value | Why unexpected | Severity |
|---------|--------------------|-------------|----------------|----------|
| Eukaryotic Translation Elongation R-HSA-156842 (79 genes) | CD4+ T cells — **top pathway** | 5.81e-118 | Ribosomal/translation pathways should not be the #1 enriched pathway in a T cell type; this is the same pattern seen in the DEG list (RPL30, RPL13, RPS12 as top DEGs). Suggests ribosomal gene bias in the differential expression, possibly from cell cycle effects or library size confounding. | **warning** |
| Eukaryotic Translation Elongation (79 genes) | Naive B cells — **top pathway** | 1.50e-117 | Same ribosomal signal as CD4+ T cells; ribosomal bias across multiple non-proliferative cell types is suspicious | **warning** |
| Formation of Pool of Free 40S Subunits R-HSA-72689 (82 genes) | HSC/MPP — **top pathway** | 5.12e-108 | Ribosome biogenesis pathways dominate all top 5 HSC/MPP pathway hits. While HSCs are actively translating, ribosomal pathway dominance across CD4+ T, Naive B, and HSC/MPP suggests a systematic technical effect rather than lineage-specific biology. | **warning** |
| Cell Cycle, Mitotic R-HSA-69278 (63 genes) | Mid erythroid — **top pathway** | 3.34e-29 | Mid erythroid cells are known to undergo terminal cell divisions but should not have cell cycle as their dominant pathway enrichment over erythroid differentiation pathways. However, this may be biologically real — erythroblast proliferation is well documented. Consistent with TOP2A, CDC20, CCNB2, MCM7 in gene list. | **info** |
| SRP-dependent Cotranslational Protein Targeting R-HSA-1799339 | pDC — **top pathway** | 1.36e-67 | Protein secretion targeting pathway as #1 pDC hit is unexpected for a cell type whose primary function is IFN-α secretion (not protein ER targeting). May reflect a real biological programme for efficient cytokine processing, or may reflect JCHAIN/secretory gene contamination from adjacent plasma cells. | **info** |
| Cytoplasmic Translation (GO:0002181) | CD8+ T cells — **top pathway** | 1.02e-70 | Cytoplasmic translation as the #1 enriched process in CD8+ T cells is unexpected; effector/immune function pathways (TCR signalling, granzyme pathways) would be expected to dominate | **warning** |

**Global GSEA concern:** The dominance of ribosomal protein gene (RPG) pathways as top hits across CD4+ T cells, Naive B cells, HSC/MPP, CD8+ T cells, and Small pre-B cells is a systematic signal. This may arise from: (1) incomplete RPG regression during normalisation, (2) library-size correlated expression of RPGs in highly represented cell types, or (3) genuine high translational activity. The pattern across biologically diverse cell types argues for (1) or (2). Filtering out RPL/RPS genes prior to GSEA is recommended.

---

## 8. Cross-Module Coherence

| Severity | Modules involved | Flag | Suggestion |
|----------|-----------------|------|------------|
| critical | QC + Annotation | 75.9% cell loss driven by max_mt_pct=5.0% likely removed a disproportionate fraction of erythroid precursors and metabolically active cells. Post-QC erythroid clusters (0–4, n=8,348 total) represent 38.3% of kept cells — the true proportion before filtering may have been higher or differently distributed across erythroid stages. | Re-run QC with max_mt_pct=10–15%; compare cluster proportions and erythroid staging before/after |
| warning | Annotation + DEG | JCHAIN (log2FC=7.02, adj-p=0.00) is a top-5 DEG for pDC (Cluster 9) AND is a known marker of plasma cells (Cluster 5, n=484). The proximity of pDC and plasma cell clusters in a bone marrow setting creates risk of gene bleed-through. The pDC cluster has only 0.50 consensus confidence. | Check JCHAIN expression per cell in both clusters via dot plot; if expression is bimodal within the pDC cluster, consider sub-clustering or reassignment |
| warning | Clustering + Annotation | Two classical monocyte clusters (7 and 8, both n~1,700, both confidence 1.00 by all methods) and two CD4+ T cell clusters (10 and 12, both >1,400 cells) and three late erythroid clusters (0, 1, 2) exist with identical consensus annotations. Either the resolution is too fine for these compartments or there are genuine transcriptional sub-states not captured by the top-5 DEG summary. | Compare full marker gene lists for Clusters 7 vs 8, 10 vs 12, and 0–2 explicitly; run sub-cluster UMAP on monocyte and T cell compartments |
| warning | DEG + GSEA | Ribosomal protein genes (RPL/RPS) dominate the top Wilcoxon DEGs for CD4+ T cells (#1–5: RPL30, RPL13, RPS12, RPL34, RPLP0), HSC/MPP (#2–5: RPS24, RPLP0, RPLP1, RPL12), and top GSEA pathways for CD4+ T (Translation Elongation adj.p=5.81e-118), Naive B (Translation adj.p=1.50e-117), HSC/MPP (40S pool adj.p=5.12e-108), and CD8+ T (Cytoplasmic Translation adj.p=1.02e-70). This cross-module pattern strongly suggests RPG bias in the Wilcoxon test. | Remove RPL/RPS/RPLP genes from DEG gene universe before re-running Wilcoxon; use SCTransform or regress out nCount_RNA; re-run GSEA on RPG-filtered DEG list |
| warning | Annotation + QC | Cluster 14 (HSC/MPP, n=506) has the lowest cell count of all clusters and moderate confidence (0.50). HSCs are rare cells. After 75.9% cell removal, the HSC cluster may be depleted of the most transcriptionally pure HSCs (which may have had higher MT%). | Compare HSC cluster size before vs after QC threshold change; include MT% as covariate in downstream trajectory analysis |
| info | Normalization + Clustering | Batch correction was applied post-clustering (Harmony tab) with theta=2.0, max_iter=50. The clustering in the report appears to have been run on uncorrected PCA (X_pca, 7 PCs used). The Harmony-corrected embedding (X_pca_harmony) was computed but may not have been used for the primary Leiden clustering. Mixing score post-Harmony=0.833 (≥0.8 threshold met). | Confirm whether clustering was performed on X_pca or X_pca_harmony; if on X_pca, re-run clustering on Harmony embedding and compare cluster assignments |
| info | Dimensionality reduction + Clustering | Only 7 PCs were used for the neighbor graph (46.3% cumulative variance explained). For 21,778 cells and 2,000 HVGs, 7 PCs is conservative; many studies use 20–30 PCs for BMMC. Low PC count may cause nearby cell types (e.g. pDC and plasma cells sharing secretory genes) to be incorrectly placed in the same neighbourhood. | Re-run with n_pcs=15–20 and compare cluster stability and silhouette scores |

**Overall coherence: review recommended**

The dataset is internally coherent in the sense that major cell type annotations are biologically plausible and well-supported by canonical markers. However, three systemic issues require attention before publication: (1) the QC MT% threshold is too stringent and drives >75% cell loss; (2) ribosomal gene dominance in DEG/GSEA outputs is a pipeline-level confound; and (3) the clustering was likely performed on uncorrected PCA before Harmony integration, which should be validated.

---

## 9. Downstream Suggestions

| Priority | Step | Rationale | Recommended tool | Expected output |
|----------|------|-----------|-----------------|----------------|
| 1 | Re-run QC with max_mt_pct=10–15% | The 5% MT threshold removed 64,699 cells (71.7% of input) and almost certainly over-filtered the erythroid and HSC compartments. This is the highest-priority correction because all downstream results are conditioned on the filtered cell set. | Scanpy (sc.pp.filter_cells with updated threshold); re-run full pipeline | Recovery of additional erythroid sub-stages, potential new rare cell types, improved cluster stability |
| 2 | Re-run clustering on Harmony-corrected embedding (X_pca_harmony) | The Harmony batch mixing score is 0.833 (adequate), but if primary clustering used uncorrected PCA, batch effects could have driven cluster splits (e.g. two classical monocyte clusters, two CD4+ T clusters). | Scanpy Leiden on `obsm['X_pca_harmony']`; compare cluster labels | Consolidated monocyte and T cell clusters; batch-robust annotations |
| 3 | RPG-filtered DEG and GSEA re-analysis | Ribosomal protein genes dominate top DEGs for CD4+ T, HSC/MPP, Naive B, and CD8+ T cells, and cascade into GSEA outputs. Translation elongation is the #1 pathway for 4+ cell types — an implausible result biologically. | Scanpy/Wilcoxon with RPL/RPS gene exclusion mask; gseapy re-run on filtered gene set; alternatively use SCTransform in Seurat for variance-stabilised normalisation | Biologically interpretable cell type markers; immune signalling pathways replacing translation pathways in enrichment |
| 4 | Sub-clustering of CD4+ T cell compartment (Clusters 10 + 12, n=5,243 combined) | Largest cell compartment in the dataset; BM CD4+ T cells span naive, central memory, effector memory, Tregs, and follicular helper subsets. Current top DEGs are RPL/RPS genes — no CD4+ sub-type resolution. | Scanpy sub-clustering (leiden) at res=0.3–0.5 on CD4+ T cell subset; or Seurat FindSubClusters | Identification of naive (IL7R+CCR7+), Treg (FOXP3+), and memory (CD44+SELL-) CD4+ sub-populations |
| 5 | Trajectory / pseudotime analysis of erythroid compartment | Clusters 0–4 (n=8,348 total, 38% of dataset) span mid-to-late erythroid stages. Continuous differentiation trajectory analysis will resolve the staging and identify regulatory genes at each transition point. Connects to HSC/MPP cluster 14 as the origin. | Scanpy PAGA + diffusion pseudotime (dpt); or Monocle3 trajectory; RNA velocity with scVelo | Ordered erythroid differentiation axis from HSC/MPP → mid erythroid → late erythroid; branch points identifying erythroid vs myeloid commitment |
| 6 | Validation of SNHG29 as HSC/MPP marker | SNHG29 is the top DEG for HSC/MPP by both Wilcoxon (log2FC=3.52, adj-p=5.78e-236) and DESeq2 (log2FC=2.67, padj=4.73e-14) — robust across two statistical frameworks. This lncRNA has no established HSC function in the literature. | Cross-reference SNHG29 expression in existing BM atlases (BMDB — PMID:41334179); ATAC-seq peak calling at SNHG29 locus in progenitor cells | Nomination of SNHG29 as a novel HSC/MPP-enriched regulatory lncRNA for wet-lab functional validation |
| 7 | ERBB2 expression validation in NK cells | ERBB2 upregulation in CD16+ NK cells (pseudobulk log2FC=3.90, padj=4.18e-15) is an unexpected pseudobulk hit. If real, this may represent a tumour-reactive NK subset of clinical relevance. | CITE-seq antibody panel for HER2 protein; cross-reference existing NK scRNA-seq datasets (e.g. Human Cell Atlas peripheral blood NK) | Confirmation or refutation of ERBB2 surface expression on BM NK cells; if confirmed, test functional relevance in HER2+ tumour killing assays |

---

## 10. Summary

**Key findings:**

1. **QC over-filtering is the primary pipeline concern:** max_mt_pct=5.0% removed 64,699 cells (71.7% of input), far exceeding the typical 20–40% range for 10X BMMC. Pre-QC median MT%=6.34% confirms many biologically real cells were discarded. All downstream results should be interpreted in this light.

2. **Cell type landscape is biologically well-supported for major lineages:** 11 of 16 clusters received high or near-unanimous annotation support across four methods (CellTypist, Markers, SingleR, ScType). Canonical erythroid markers (HBB log2FC=8.91, AHSP, HEMGN), monocyte markers (S100A9 log2FC=9.49, LYZ, FCN1), NK markers (GNLY log2FC=9.45, NKG7), and B cell markers (MS4A1 log2FC=7.97, CD79A, VPREB3) are all confirmed in the DEG outputs. ScType's systematic Megakaryocyte calls for non-megakaryocyte clusters (Clusters 5, 6, 9, 10, 12) represent a database-level artefact and should be excluded from downstream reporting.

3. **SNHG29 is a high-priority novel HSC/MPP marker:** Top DEG by both Wilcoxon (log2FC=3.52, adj-p=5.78e-236) and DESeq2 pseudobulk (log2FC=2.67, padj=4.73e-14) for the HSC/MPP cluster (n=506), with no established haematopoietic function in the literature. Its consistent top-ranking across two orthogonal tests in healthy human BM makes it a strong candidate for follow-up as a novel regulatory lncRNA in human haematopoietic stem cell biology.

**Open questions:**

- Does ERBB2 upregulation in CD16+ NK cells (pseudobulk log2FC=3.90, padj=4.18e-15) represent a genuine NK cell sub-population, a contaminating cell type, or a statistical artefact of the pseudobulk design?
- Are the two classical monocyte clusters (7 and 8) and two CD4+ T cell clusters (10 and 12) truly distinct sub-populations, or artefactual splits driven by batch or cell cycle effects before Harmony correction?
- Is the JCHAIN signal in pDCs (log2FC=7.02) reflecting a genuine immunoglobulin chain programme in pDCs, or transcript contamination from adjacent plasma cells?
- What is the functional role of SNHG29 in HSC self-renewal or differentiation?
- If QC is repeated with max_mt_pct=10–15%, will additional erythroid sub-stages (reticulocytes, pro-erythroblasts, megakaryocyte-erythroid progenitors) become visible?

**Suggested validation experiments:**

- **Flow cytometry / CITE-seq:** Validate ERBB2 (HER2) protein surface expression on BM CD16+ NK cells using existing CITE-seq antibody panel data (the dataset is CITE-seq; HER2 antibody may be in the TotalSeq-B panel — check panel metadata).
- **SNHG29 RNA-FISH or single-molecule FISH:** Co-localise SNHG29 RNA with CD34+CD38- HSC surface markers in human BM sections to confirm HSC-specific expression at the single-cell spatial level.
- **Re-run full pipeline with relaxed MT% threshold (10–15%):** Computational experiment requiring no new data; directly tests whether current conclusions are robust to QC parameter choice.
- **Trajectory analysis (Monocle3 or PAGA):** Connect the HSC/MPP cluster (n=506) through mid erythroid (n=2,838) to late erythroid (n=4,510) to generate a differentiation pseudotime axis and identify regulatory genes at each transition — directly testable against existing chromatin accessibility data if ATAC-seq modality is analysed.

---

*OmicSage D1 Manual Review Mode — report_review.md*
*Reviewed by: Claude Sonnet 4.6 via MASTER_PROMPT.md v1.0*
*Dataset: BMMC CITE-seq (NeurIPS 2021) — GEO not directly reported; competition dataset available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122*
*Literature retrieved via PubMed MCP tool (PubMed attribution per tool requirements); DOIs included where available.*
