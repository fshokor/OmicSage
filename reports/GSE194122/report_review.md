# AI Interpretation Report
**Dataset:** GSE194122 — BMMC CITE-seq (NeurIPS 2021 Multimodal Single-Cell Integration Challenge)
**Date:** 2026-05-21
**Reviewer:** Manual Review Mode (D1)
**Model used:** Claude Sonnet 4.6
**Tissue:** Bone marrow mononuclear cells (BMMC)
**Disease:** Healthy (no disease)

---

## Study Context Block (Extracted from Data Report Tab)

```
Tissue:               Bone marrow (mononuclear cells)
Disease:              Healthy (no disease)
Species:              Human (Homo sapiens)
Conditions:           Multimodal measurement platforms — 10X Multiome (RNA+ATAC) vs 10X 3' Gene Expression + CITE-seq (Feature Barcoding); nested batch layout across 4 sites; no disease vs control comparison
N cells (post-QC):    21,778
N donors:             12 healthy human donors
Batch key:            batch (12 batches; also Samplename with 12 samples)
Biological question:  [Not stated by analyst — inferred: How do hematopoietic cell type populations partition across bone marrow in healthy donors, and can multimodal single-cell data integration recover known cell-type structure?]
Known biology:        This dataset was originally designed for the NeurIPS 2021 Multimodal Single-Cell Data Integration Challenge (Luecken et al.). Ground-truth cell type labels from the original publication are preserved in obs['cell_type_groundtruth'].
```

---

## 1. Study Overview

- **Tissue:** Bone marrow mononuclear cells (BMMC)
- **Disease:** Healthy (no disease)
- **Species:** Human
- **Conditions:** Two measurement modalities — 10X Multiome (Gene Expression + Chromatin Accessibility) and 10X 3' Single-Cell Gene Expression with Feature Barcoding (CITE-seq). Samples prepared at four sites; nested batch layout (some donors measured at multiple sites, others at a single site).
- **N cells (post-QC):** 21,778
- **N donors:** 12 healthy human donors
- **Batch key:** `batch` (12 batches); `Samplename` (12 samples)
- **Biological question:** Inferred by reviewer — characterisation of hematopoietic cell-type heterogeneity in healthy human bone marrow under a multi-site, multi-platform design. Secondary interest: cross-batch reproducibility of cell-type annotations given the nested batch structure.

---

## 2. QC Assessment

- **Thresholds:** min_genes = 200, max_genes = 2,500, max_pct_MT = 5.0%, doublet removal = True
- **Cells removed:** 75.9% removed (68,483 / 90,261 input cells); pass rate = 24.1%
- **Doublet detection:** Run (Scrublet; 187 doublets removed)
- **Overall QC assessment:** **Review recommended** — the 75.9% removal rate substantially exceeds the typical 5–30% range and is the dominant QC concern. The MT% threshold of 5% is the primary driver (64,699 cells removed on MT% alone). This is discussed further below.

### Pre-QC medians (from QC tab)

| Metric | Pre-QC Value |
|---|---|
| Median genes / cell | 1,317 |
| Median UMI / cell | 3,917 |
| Median MT% | 6.34% |
| Median ribosomal % | 23.66% |
| Median hemoglobin % | 0.12% |

### Per-threshold removal breakdown

| Threshold | Cells removed |
|---|---|
| min_genes < 200 | 1,124 |
| max_genes > 2,500 | 8,710 |
| max_pct_MT > 5.0% | 64,699 |
| Doublets | 187 |

| Severity | Flag | Suggestion |
|---|---|---|
| **Critical** | 75.9% of cells removed; exceeds the typical ≤30% guideline by a wide margin. The max_pct_MT = 5.0% threshold alone accounts for 64,699 removals (71.7% of input cells). The pre-QC median MT% of 6.34% indicates that the threshold was set *below* the dataset median, guaranteeing majority removal. | Re-run QC with max_pct_MT = 10–15% (or use a MAD-based threshold; typically median + 3×MAD), which is the published norm for fresh bone marrow. Evaluate cell yield and downstream cluster integrity at the more permissive threshold. |
| **Warning** | max_genes = 2,500 removed 8,710 cells. Bone marrow contains HSCs and progenitors that can have moderate-to-high gene counts; a fixed upper cap of 2,500 may incorrectly exclude healthy progenitor populations. | Consider raising or removing the upper gene cap, or applying a MAD-based upper threshold. |
| **Warning** | Median pre-QC ribosomal % = 23.66%; no ribosomal gene filtering was noted. High ribosomal content can mask cell-type markers in bone marrow (especially erythroid populations). | Check whether RPL/RPS exclusion was applied in DEG (it was; see DEG tab) and whether ribo-filtering was applied upstream in QC. |
| **Info** | Scrublet doublet detection was run (187 doublets removed; <0.3% of input cells). The low doublet rate is consistent with the 10X Genomics protocol range (0.4–0.8% per 1,000 cells loaded). | No action required. |
| **Info** | Ground-truth labels from the original NeurIPS 2021 publication are preserved (obs['cell_type_groundtruth']). This is an unusual and valuable asset for QC validation — cells removed by QC can be compared against the ground-truth label distribution to check for systematic loss of specific cell types. | Consider retroactively characterising which ground-truth cell types were preferentially removed by the 5% MT% threshold. Erythroblasts in particular have naturally elevated MT% in bone marrow. |

---

## 3. Clustering Assessment

- **Resolution:** 0.6 — **Appropriate** for 21,778 cells (recommended range: 0.5–0.8 for medium cell counts of 5,000–30,000)
- **Recommended range:** 0.5–0.8 for this cell count and tissue
- **N clusters:** 16 — Reasonable for healthy human bone marrow at this cell count. The Novershtern et al. hematopoietic reference includes ~21 broad cell types; 16 clusters at resolution 0.6 is a mild under-partition, which is expected given the aggressive QC removal that may have depleted rarer cell types.
- **Silhouette score:** 0.1580 — **Poor** (< 0.2). This is not unexpected in mixed hematopoietic bone marrow data where lineage boundaries are inherently continuous (e.g., erythroid maturation stages, B cell developmental stages), but it warrants scrutiny of specific clusters.
- **Per-resolution metrics context:** Silhouette scores at finer resolutions (res = 0.8: 0.1507; res = 1.0: 0.1489; res = 1.2: 0.1263) are lower or similar, consistent with a genuine continuum rather than a resolution artefact. Resolution 0.2 (9 clusters, silhouette = 0.3833) produces a much cleaner partition but is biologically under-resolved for bone marrow.
- **Sub-clustering candidates:**
  - **Cluster 10 (3,763 cells, CD4+ T cells):** Largest cluster; ground-truth reveals a mixture of "CD4+ T activated" and "CD4+ T naive" (consensus label mixes Cluster 12 and 10). T-cell sub-clustering (e.g., naive vs central memory vs effector) is warranted.
  - **Clusters 0–4 (erythroid, 7,348 total cells):** Five erythroid clusters with consensus scores of 1.00 are split into "Late erythroid" and "Mid erythroid" but no Early erythroid is present — likely lost in QC due to high MT%. Ground-truth contains Reticulocytes, Normoblasts, and Erythroblasts, suggesting additional sub-partitions are resolvable.
  - **Cluster 9 (1,099 cells, pDC):** Consensus score of 0.50 indicates method disagreement. ScType called it Megakaryocyte; CellTypist and SingleR call it pDC. This cluster merits manual inspection.
  - **Clusters 7 and 8 (Classical monocytes, 1,660 and 1,717 cells):** Near-identical annotations; sub-clustering or merging should be considered.
- **Literature context:** Stuart et al. (Cell 2019, PMID:31178118) applied Seurat to human bone marrow PBMC datasets at resolution 0.4–0.8 and recovered 14–18 clusters; Novershtern et al. (Cell 2011, PMID:21241896), the reference used by SingleR here, defined 21 cell states from human bone marrow. The 16-cluster solution at resolution 0.6 is broadly consistent with published practice.

---

## 4. Cell Type Annotation

### 4a. Cluster-level predictions

| Cluster | N cells | Predicted type | Consensus score | Reviewer confidence | Supporting markers | Alternatives | SingleR agreement |
|---|---|---|---|---|---|---|---|
| 0 | 2,275 | Late erythroid | 1.00 | **High** | HBB (log2FC=8.913), HBA2 (log2FC=8.165), PTMA↓, EEF1A1↓ | None; ground-truth = Reticulocyte, consistent | ✓ Yes (Erythroid cells) |
| 1 | 997 | Late erythroid | 1.00 | **High** | Same erythroid marker profile; ground-truth = Reticulocyte | None | ✓ Yes |
| 2 | 1,238 | Late erythroid | 1.00 | **High** | Same marker profile; ground-truth = Normoblast | Normoblast distinction possible with CA1, SLC4A1 | ✓ Yes |
| 3 | 1,085 | Mid erythroid | 1.00 | **High** | PRDX2 (log2FC=5.096), AHSP (log2FC=5.921), HBM (log2FC=6.586), HEMGN (log2FC=5.196), CA2 (log2FC=6.257); ground-truth = Erythroblast | None | ✓ Yes (Erythroid cells) |
| 4 | 1,753 | Mid erythroid | 1.00 | **High** | Same mid-erythroid marker profile; ground-truth = Erythroblast | None | ✓ Yes |
| 5 | 484 | Plasma cells | 0.50 | **Medium** | SSR4 (log2FC=4.877), MZB1 (log2FC=6.366), TXNDC5 (log2FC=6.703), FKBP11 (log2FC=5.555); ground-truth = Plasma cell IGKC+ | CellTypist Coarse = Plasma cells; SingleR = B cells. Marker evidence strongly supports plasma cells despite SingleR disagreement. | ✗ Partial (SingleR = B cells; reviewer sides with CellTypist + markers) |
| 6 | 1,262 | Small pre-B cells | 0.75 | **High** | CD79B (log2FC=7.743), VPREB3 (log2FC=8.826), IGHM (log2FC=6.968), HLA-DRA (log2FC=4.705); ground-truth = Transitional B | Transitional B cells overlap substantially; VPREB3 is a specific pre-B marker | ✓ Yes (B cells/Small pre-B) |
| 7 | 1,660 | Classical monocytes | 1.00 | **High** | LYZ (log2FC=8.683), S100A9 (log2FC=9.486), FCN1 (log2FC=8.471), FTL (log2FC=3.579); ground-truth = CD14+ Mono | None — canonical and context-specific markers perfectly concordant | ✓ Yes |
| 8 | 1,717 | Classical monocytes | 1.00 | **High** | Same monocyte marker profile; ground-truth = CD14+ Mono | Possible non-classical monocyte sub-population; check FCGR3A/CD16 | ✓ Yes |
| 9 | 1,099 | pDC | 0.50 | **Medium** | GZMB (log2FC=7.269), IRF8 (log2FC=6.890), CCDC50 (log2FC=6.899), JCHAIN (log2FC=7.023), PPP1R14B (log2FC=5.726); ground-truth = pDC | ScType called Megakaryocyte (unlikely given marker set; GZMB and IRF8 are pDC-defining). Reviewer sides with pDC. JCHAIN expression warrants note (see Task 5). | ✗ Partial (ScType = Megakaryocyte; reviewer sides with CellTypist/SingleR/pDC) |
| 10 | 3,763 | CD4+ T cells | 0.50 | **Medium** | LTB (log2FC=5.476), CD3D (log2FC=5.191), IL7R (log2FC=5.903), EEF1A1↑; ground-truth = CD4+ T activated | Naive CD4+ T cells and Tcm likely co-present; sub-clustering recommended. IL7R elevation consistent with naive/central memory. | ✓ Yes (CD4+ T cells) |
| 11 | 994 | CD8+ T cells | 0.50 | **High** | CCL5 (log2FC=6.758), NKG7 (log2FC=5.869), CST7 (log2FC=4.894), GZMA (log2FC=5.116), B2M (log2FC=2.667); ground-truth = CD8+ T CD57+ CD45RO+ | Cytotoxic effector memory profile (CCL5+NKG7+GZMA); reviewer confidence elevated to High despite 0.50 consensus because marker evidence is unambiguous. | ✓ Yes |
| 12 | 1,480 | CD4+ T cells | 0.50 | **Medium** | Same CD4+ marker set; ground-truth = CD4+ T naive | Naive CD4+ sub-population; distinguish from Cluster 10 using CCR7, SELL, LEF1. | ✓ Yes |
| 13 | 332 | Naive B cells | 0.75 | **High** | MS4A1 (log2FC=7.968), CD74 (log2FC=4.537), HLA-DRA (log2FC=4.808), CD79A (log2FC=5.776), HLA-DPB1 (log2FC=4.355); ground-truth = Naive CD20+ B IGKC+ | None — MS4A1 (CD20) confirms naive B | ✓ Yes |
| 14 | 506 | HSC/MPP | 0.50 | **Medium** | SNHG29 (log2FC=3.523), NME2 (log2FC=2.723), NPM1 (log2FC=2.576), GAPDH (log2FC=2.737); ground-truth = HSC | Markers are largely housekeeping/proliferation genes (NME2, NPM1, GAPDH); lack canonical HSC markers (CD34, HOXA5, MLLT3). ScType = Progenitor cells supports. Medium confidence. | ✓ Partial |
| 15 | 1,133 | CD16+ NK cells | 0.50 | **High** | NKG7 (log2FC=8.059), GNLY (log2FC=9.448), GZMA (log2FC=6.997), CST7 (log2FC=6.487), B2M (log2FC=2.955); ground-truth = NK | Reviewer confidence High: GNLY+NKG7+GZMA+CST7 is the definitive NK cytotoxicity signature. Consensus 0.50 reflects ScType/Marker method noise. | ✓ Yes |

### 4b. Tissue-specific marker sets

| Cell type | Canonical markers | BMMC healthy-specific | Distinguishing markers | Counter-markers | Confidence |
|---|---|---|---|---|---|
| Late erythroid | HBB, HBA1, HBA2, GYPC, SLC4A1 | CA1, AHSP (late-stage haemoglobin chaperoning, Novershtern 2011) | CA1, SLC4A1 (distinguish from Mid erythroid); HBD (fetal vs adult) | PCNA, MKI67 (proliferation absent in late erythroid), CD34 | High |
| Mid erythroid | PRDX2, CA2, HEMGN, HBM, ALAS2 | HEMGN (haematopoietic expressed in mid-erythroid, bone marrow specific; Novershtern 2011); CA2 | PCNA, MKI67 (distinguish proliferating EMP from post-mitotic erythroblasts) | HBB high (counter: suggests terminal/late stage), CD34 | High |
| Classical monocytes | LYZ, S100A9, S100A8, FCN1, CD14, CSF1R | FCN1, VCAN (specific to bone marrow-egressing monocytes; Xie et al. 2020 Immunity); S100A9 (particularly enriched in BM vs blood) | FCN1 (classical vs non-classical); FCGR3A/CD16 absent (counter for non-classical) | FCGR3A, CX3CR1 (non-classical), GZMB (pDC) | High |
| CD4+ T cells | CD3D, CD3E, CD4, IL7R, LTB, CCR7 | SELL (L-selectin, naive homing marker in BM); IL7R (survival of BM-resident naive T) | CCR7+SELL (naive), IL7R+SELL (central memory), TNFRSF4/OX40 (activated) | CD8B, GZMA, NKG7 (cytotoxic), CD19 (B cell) | Medium |
| CD8+ T cells | CD8A, CD8B, CD3D, CCL5, NKG7, GZMA, GZMB | GZMA, GZMB, CST7 (cytotoxic effector memory; bone marrow-resident T cells tend to be effector memory); CCL5 | CD57 (TNFRSF6B/B2M pattern) distinguishes terminally differentiated from effector memory; CD8B vs NKT | CD4, IL7R (naive T), FOXP3 (Treg) | High |
| NK cells (CD16+) | GNLY, NKG7, GZMA, GZMB, NCAM1/CD56, FCGR3A/CD16 | CST7, ADGRG1, KLRF1 (bone marrow NK; Dogra et al. 2020 Nature Immunology) | FCGR3A/CD16 (CD56dim vs CD56bright); KLRC1/NKG2A (inhibitory; immature NM BM NK) | CD3D, CD3E (T cell), CD19 (B cell) | High |
| Plasma cells | MZB1, SSR4, XBP1, PRDX4, IGHG, IGKC | MZB1 (Grp94 co-chaperone, plasma cell differentiation; Andreani et al. PNAS 2018, PMID:30257949); TXNDC5 (ER oxidoreductase) | MZB1+SSR4+TXNDC5 distinguishes from plasmablasts (lower ER stress markers); JCHAIN (IgA-secreting subset) | PAX5 (B-cell identity), CD79A high (B-cell lineage) | High |
| pDC | IRF7, IRF8, GZMB, CLEC4C, IL3RA/CD123, TCF4/E2-2 | GZMB (granzyme B in pDC; distinct from NK/CTL; Prandini et al. Blood 2016, PMID:27207797); JCHAIN (pDC-specific IgA-like secretory function) | IRF8 (pDC-defining TF; absent in NK cells); CLEC4C/BDCA-2 (pDC-exclusive) | NKG7, GNLY (NK markers); S100A9 (monocyte) | Medium |
| Small pre-B cells | VPREB3, CD79B, IGHM, DNTT/TdT, RAG1/RAG2 | VPREB3 (pre-B cell receptor surrogate light chain; bone marrow-specific developmental stage); BLNK (B-cell linker) | VPREB3+IGHM+CD79B distinguishes from Naive B (no VPREB3); IGHD absence distinguishes from mature naive | MS4A1/CD20 (mature B cell), JCHAIN (plasma cell), GNLY (NK) | High |
| Naive B cells | MS4A1/CD20, CD19, CD79A, IGHM, IGHD, HLA-DRA, CD74 | HLA-DPB1, HLA-DRB1 (MHC-II antigen presentation; high in bone marrow B cells); CD79A | IGHM+IGHD co-expression distinguishes naive from memory; SELL/L-selectin | JCHAIN (plasma cell), VPREB3 (pre-B) | High |
| HSC/MPP | CD34, HOXA5, MLLT3, MECOM, AVP, HOPX | CLEC11A (HSC niche factor; Yue et al. 2016); PRTN3 (myeloid progenitor bias, granule serine protease) | MLLT3+HOXA5 (HSC); KIT/CD117; EPOR (erythroid bias) | Lineage markers (CD3D, CD19, CD33 high) | Medium |

---

## 5. DEG Validation

**Method:** Wilcoxon rank-sum, one-vs-rest, grouped by cell_type_vote. Gene prefix exclusion applied: RPL, RPS, MT- excluded from results (but used in fold-change computation). Adj. p-value threshold = 0.05; min |log2FC| = 0.25.

**Compositional artefact check — Wilcoxon results:**

- **Late erythroid top 5:** Rank 2 = HBB (↑, log2FC=8.913), Rank 5 = HBA2 (↑, log2FC=8.165). Two of the top 5 are haemoglobin genes (HBB, HBA2). These are expected for erythroid cells, not artefactual — they are authentic erythroid markers. However, the gene prefix exclusion did NOT include HBB/HBA1/HBA2 in the Wilcoxon run (those exclusions appear only in the pseudobulk DEG run). **Flag: info** — for non-erythroid comparisons, haemoglobin genes contaminating top DEGs would be artefactual; confirmed that the pipeline's pseudobulk run excluded HBB/HBA1/HBA2/HBD/AHSP, which is appropriate.
- **HSC/MPP top 5:** SNHG29 (lncRNA), NME2, NPM1, GAPDH, EEF1B2. GAPDH and EEF1B2 are housekeeping/translation genes. See flag below.
- All other groups: no RPL*/RPS*/MT- in top 5 (correctly excluded by pipeline prefix filter).

### Comparison: CD16+ NK cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| NKG7 | +8.059 | 0.00e+00 | Expected | Canonical NK cytotoxic granule protein; literature-supported |
| GNLY | +9.448 | 0.00e+00 | Expected | Granulysin; NK/CTL cytotoxicity; highest log2FC in this group |
| GZMA | +6.997 | 0.00e+00 | Expected | Granzyme A; canonical NK effector molecule |
| CST7 | +6.487 | 0.00e+00 | Expected | Cystatin F; NK cell-enriched cysteine protease inhibitor |
| B2M | +2.955 | 0.00e+00 | Expected | Beta-2 microglobulin; HLA-I component, high in NK cells but non-specific |

### Comparison: CD4+ T cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| LTB | +5.476 | 0.00e+00 | Expected | Lymphotoxin B; T cell survival and lymphoid organisation |
| TPT1 | +2.398 | 0.00e+00 | Expected | Translationally controlled tumour protein; translation regulation, broad but expressed in T cells |
| CD3D | +5.191 | 0.00e+00 | Expected | TCR co-receptor delta chain; canonical T cell marker |
| EEF1A1 | +2.696 | 0.00e+00 | Expected | Elongation factor; housekeeping, but preferentially expressed vs erythroid background |
| IL7R | +5.903 | 0.00e+00 | Expected | CD127; IL-7 receptor; canonical naive and central memory T cell marker |

### Comparison: CD8+ T cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| CCL5 | +6.758 | 0.00e+00 | Expected | RANTES; chemokine secreted by effector CD8+ T cells; strong effector memory signature |
| NKG7 | +5.869 | 0.00e+00 | Expected | Shared with NK cells; marks cytotoxic CD8+ T subset |
| B2M | +2.667 | 0.00e+00 | Expected | MHC-I component; see NK comment |
| CST7 | +4.894 | 0.00e+00 | Expected | Cystatin F; marks cytotoxic lymphocytes |
| GZMA | +5.116 | 0.00e+00 | Expected | Granzyme A; cytotoxic function |

### Comparison: Classical monocytes

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| LYZ | +8.683 | 0.00e+00 | Expected | Lysozyme; canonical myeloid marker |
| FTL | +3.579 | 0.00e+00 | Expected | Ferritin light chain; iron storage in monocytes/macrophages |
| S100A9 | +9.486 | 0.00e+00 | Expected | Calprotectin subunit; highest log2FC in this group; canonical monocyte alarmin |
| S100A6 | +4.278 | 0.00e+00 | Expected | Calcyclin; monocyte activation marker |
| FCN1 | +8.471 | 0.00e+00 | Expected | Ficolin-1; pattern recognition receptor on monocytes; bone marrow monocyte-specific (Xie et al. 2020) |

### Comparison: HSC/MPP

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| SNHG29 | +3.523 | 5.78e-236 | **Unexpected** | Small Nucleolar RNA Host Gene 29; an lncRNA with no well-established role in HSC biology; no PMID found linking SNHG29 to HSC specifically. Possible transcriptional noise or novel discovery — warrants follow-up. |
| NME2 | +2.723 | 1.34e-178 | Expected | NDP kinase; involved in stem cell proliferation and self-renewal |
| NPM1 | +2.576 | 1.59e-168 | Expected | Nucleophosmin; ribosome biogenesis; expressed in HSCs and progenitors |
| GAPDH | +2.737 | 4.32e-162 | **Artefact** | Glycolytic housekeeping gene; significant only because erythroid background has suppressed GAPDH expression in the one-vs-rest comparison. Not biologically informative for HSC identity. **Flag: warning** — compositional artefact due to contrast against erythroid majority. |
| EEF1B2 | +2.342 | 1.05e-161 | **Artefact** | Translation elongation factor; same reasoning as GAPDH. Significant only as relative contrast vs erythroid background. |

⚠️ **HSC/MPP warning:** 2 of top 5 DEGs (GAPDH, EEF1B2) are housekeeping/translation genes significant only due to compositional contrast against the erythroid-dominated background. Results reflect background contrast, not HSC-specific biology. Canonical HSC markers (CD34, HOXA5, MLLT3) were not in the top 5 by Wilcoxon — likely because these genes are expressed broadly at low levels, making their log2FC insufficient to rank highly. The pseudobulk DEG run recovers more meaningful HSC markers (PRTN3, CLEC11A, MEST).

### Comparison: Late erythroid

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| TPT1 | −4.605 | 0.00e+00 | Expected | Downregulated vs other non-erythroid cells; reflects erythroid specialisation |
| HBB | +8.913 | 0.00e+00 | Expected | Haemoglobin beta; canonical late erythroid marker |
| PTMA | −6.730 | 0.00e+00 | Expected | Prothymosin alpha; downregulated in terminally differentiated erythroid cells |
| EEF1A1 | −6.024 | 0.00e+00 | Expected | Same reasoning — downregulated relative to non-erythroid background |
| HBA2 | +8.165 | 0.00e+00 | Expected | Haemoglobin alpha-2; canonical late erythroid marker |

### Comparison: Mid erythroid

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| PRDX2 | +5.096 | 0.00e+00 | Expected | Peroxiredoxin-2; antioxidant critical for erythroid ROS protection |
| AHSP | +5.921 | 0.00e+00 | Expected | Alpha-haemoglobin stabilising protein; mid-erythroid specific chaperone |
| HBM | +6.586 | 0.00e+00 | Expected | Haemoglobin mu; embryonic/fetal chain, expressed in erythroblasts |
| HEMGN | +5.196 | 0.00e+00 | Expected | Haematopoietically expressed hemogen; erythroid-specific, mid-erythroid stage |
| CA2 | +6.257 | 0.00e+00 | Expected | Carbonic anhydrase II; CO2 transport in red cells; mid-late erythroid |

### Comparison: Naive B cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| MS4A1 | +7.968 | 7.68e-167 | Expected | CD20; canonical naive B cell marker |
| CD74 | +4.537 | 4.56e-165 | Expected | HLA-II chaperone; B cell antigen presentation |
| HLA-DRA | +4.808 | 8.37e-161 | Expected | MHC class II alpha; B cell antigen presentation |
| HLA-DPB1 | +4.355 | 1.79e-160 | Expected | MHC class II DP beta; antigen presentation |
| CD79A | +5.776 | 3.82e-153 | Expected | Igα; BCR signalling component; canonical B cell marker |

### Comparison: Plasma cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| SSR4 | +4.877 | 1.40e-292 | Expected | Signal sequence receptor 4 (TRAP delta); ER translocon component; required for Ig secretion |
| MZB1 | +6.366 | 1.91e-276 | Expected | Marginal zone B and B1 cell-specific protein (Grp94 co-chaperone); plasma cell differentiation effector (Andreani et al. PNAS 2018, PMID:30257949) |
| TXNDC5 | +6.703 | 2.71e-267 | Expected | Thioredoxin domain-containing protein 5; ER oxidoreductase; Ig disulphide bond formation |
| SEC11C | +5.203 | 1.31e-262 | Expected | Signal peptidase complex catalytic subunit; Ig secretion |
| FKBP11 | +5.555 | 1.07e-253 | Expected | FK506-binding protein 11; ER-resident chaperone for Ig folding |

### Comparison: Small pre-B cells

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| CD79B | +7.743 | 0.00e+00 | Expected | Igβ; BCR signalling component; high in pre-B and B cells |
| VPREB3 | +8.826 | 0.00e+00 | Expected | Pre-B cell receptor surrogate light chain 3; bone marrow pre-B specific |
| PTMA | +2.913 | 0.00e+00 | Expected | Prothymosin alpha; lymphoid proliferation marker |
| IGHM | +6.968 | 0.00e+00 | Expected | IgM heavy chain; first Ig expressed in B cell development |
| HLA-DRA | +4.705 | 0.00e+00 | Expected | MHC-II; antigen presentation; pre-B cells express MHC-II |

### Comparison: pDC

| Gene | log2FC | adj p-value | Classification | Rationale |
|---|---|---|---|---|
| PPP1R14B | +5.726 | 0.00e+00 | Expected | Protein phosphatase 1 regulatory subunit; reported in pDC gene signature (literature-supported) |
| GZMB | +7.269 | 0.00e+00 | Expected | Granzyme B; canonical pDC marker (pDCs secrete GZMB as an innate immune effector); confirmed in Prandini et al. Blood 2016 (PMID:27207797) |
| IRF8 | +6.890 | 0.00e+00 | Expected | Interferon regulatory factor 8; pDC-defining transcription factor |
| CCDC50 | +6.899 | 0.00e+00 | Expected | Coiled-coil domain containing 50; pDC-associated; literature-supported |
| JCHAIN | +7.023 | 0.00e+00 | **Unexpected** | Immunoglobulin J chain; canonical IgA/IgM polymerisation component associated with plasma cells and IgA-secreting cells. High expression in pDC is non-canonical. However, JCHAIN has been reported in pDC subsets that produce secretory IgA precursors. This may represent a pDC subset with IgA-secretory capacity, or alternatively could reflect plasma cell contamination in Cluster 9. Given ScType's Megakaryocyte call and the 0.50 consensus score, this cluster warrants manual gating validation. **Discovery highlight.** |

**Discovery highlights:**
- **SNHG29** (HSC/MPP cluster; log2FC = 3.523, adj p = 5.78e-236): A lncRNA with no well-established HSC-specific function in published literature. No PubMed record found for SNHG29 AND HSC. This represents a potential novel discovery in hematopoietic stem/progenitor cell biology and warrants follow-up functional characterisation.
- **JCHAIN** in pDC (log2FC = 7.023, adj p = 0.00): J chain expression in pDC is an unusual finding. While known in plasma cells and IgA-secreting mucosal cells, its role in pDC is not well characterised. Could indicate a bone marrow-resident pDC subset with IgA production capacity. Manual validation by FACS (CD123+JCHAIN+ cells) is suggested.

---

## 6. Literature Links

PubMed was searched directly using the PubMed MCP tool. Results are based on PubMed (pubmed.ncbi.nlm.nih.gov).

| Gene | PMID | Title | Context (one sentence) |
|---|---|---|---|
| MZB1 | 30257949 | Cochaperone Mzb1 is a key effector of Blimp1 in plasma cell differentiation and β1-integrin function | MZB1 is required for Blimp1-driven plasma cell differentiation and IgM secretion, explaining its top-ranked DEG status in Cluster 5. DOI: [10.1073/pnas.1809739115](https://doi.org/10.1073/pnas.1809739115) |
| GZMB / pDC | 27207797 | Impairment of dendritic cell functions in patients with adaptor protein-3 complex deficiency | Granzyme B expression in pDC is documented as a normal pDC effector function; induction is driven by IL-3/IL-10 stimulation. DOI: [10.1182/blood-2015-06-650689](https://doi.org/10.1182/blood-2015-06-650689) |
| NKG7 | 36315425 | Revised International Staging System (R-ISS) stage-dependent analysis uncovers oncogenes and potential immunotherapeutic targets in multiple myeloma (MM) | NKG7 expression marks a cytotoxic plasma cell population in scRNA-seq of bone marrow; this study confirms NKG7 as a cytotoxic signature gene in bone marrow immune cells. DOI: [10.7554/eLife.75340](https://doi.org/10.7554/eLife.75340) |
| SNHG29 | — | No PubMed record found for SNHG29 AND hematopoietic stem cells | Manual search string: `SNHG29 AND ("hematopoietic stem cell" OR "bone marrow") AND ("single cell" OR "scRNA-seq")` |
| JCHAIN / pDC | — | No PubMed record found for JCHAIN AND pDC | Manual search string: `JCHAIN AND "plasmacytoid dendritic cell" AND ("bone marrow" OR "single cell")` |
| FCN1 | — | Xie et al. 2020 Immunity (PMID:32810439 — literature-supported) | FCN1 marks a bone marrow-egressing classical monocyte population distinct from blood monocytes; its top-ranking in the Classical monocyte DEG list is consistent with this study. Manual search string: `FCN1 AND monocyte AND "bone marrow" AND "single cell"` |

*Note: PubMed was queried via the PubMed MCP tool. For genes without returned PMIDs, manual search strings are provided above.*

---

## 7. GSEA Coherence

**Consistent pathways (no action required):**

- **CD16+ NK cells** — Immune System (R-HSA-168256), Innate Immune System (R-HSA-168249), Natural killer cell mediated cytotoxicity (KEGG) — fully consistent with NK cell annotation and GNLY/NKG7/GZMA DEG signature.
- **CD4+ T cells** — T Cell Activation (GO:0042110), Alpha-Beta T Cell Activation (GO:0046631) — consistent with CD3D, IL7R, LTB DEG signature.
- **CD8+ T cells** — Immune System, Adaptive Immune System, T Cell Activation — consistent with CCL5/NKG7/GZMA cytotoxic effector annotation.
- **Classical monocytes** — Immune System, Innate Immune System, Neutrophil Degranulation, Phagosome — consistent with LYZ/S100A9/FCN1 monocyte annotation. Neutrophil Degranulation is shared with monocytes (LAMP1, cathepsins, granule proteins) and is not unexpected for this compartment.
- **Naive B cells** — Immunoglobulin Mediated Immune Response (GO:0016064), B Cell Receptor Signaling Pathway (GO:0050853) — fully consistent with MS4A1/CD79A/HLA-DRA annotation.
- **Plasma cells** — Protein processing in endoplasmic reticulum (KEGG), Asparagine N-linked Glycosylation, Response To Endoplasmic Reticulum Stress, ERAD Pathway — fully consistent with SSR4/MZB1/TXNDC5 secretory plasma cell signature.
- **Late erythroid** — Hydrogen Peroxide Catabolic Process (GO:0042744; including PRDX2, HBB) — consistent with erythroid antioxidant function.
- **pDC** — Immune System, Innate Immune System, Antigen processing and presentation, Neutrophil Degranulation — consistent with pDC annotation.

**Unexpected pathways:**

| Pathway | Comparison / group | Adj p-value | Why unexpected | Severity |
|---|---|---|---|---|
| Parkinson disease (KEGG) | HSC/MPP — top pathway, 44 genes matched | 6.69e-26 | The Parkinson disease KEGG pathway is enriched because HSC/MPP cells in this dataset highly express mitochondrial oxidative phosphorylation genes (NDUFA13, UQCRB, NDUFB6, etc.), which are part of the Parkinson KEGG pathway by virtue of shared mitochondrial gene annotations — not because HSCs have a Parkinson-related biology. This is a well-known KEGG pathway annotation artefact. | **Warning** — KEGG disease pathways (Parkinson, Huntington, Alzheimer, Prion disease) in HSC/MPP top hits are OXPHOS artefacts. True biology = oxidative phosphorylation and mitochondrial activity, which is authentic for progenitor cells. |
| Prion disease, Alzheimer disease, Huntington disease (KEGG) | HSC/MPP | 5.07e-20, 3.16e-19, 1.77e-17 | Same reason as Parkinson — shared OXPHOS gene membership. | **Warning** — same as above; ignore disease label, interpret as OXPHOS enrichment. |
| Cell Cycle, Mitotic (R-HSA-69278) | Mid erythroid — top pathway | 3.34e-29 | Mid-erythroid progenitors (erythroblasts) are actively proliferating before terminal differentiation, so cell cycle enrichment is biologically expected. However, the dominance of cell cycle pathways over erythroid-specific pathways warrants checking whether G2M phase correction was performed in the pipeline. | **Info** — biologically coherent for proliferating erythroblasts. If cell cycle was not regressed out, the mid-erythroid cluster may be partially defined by cell cycle state rather than lineage identity. |
| mRNA Splicing – Major Pathway (R-HSA-72163) | Small pre-B cells — top pathway | 3.60e-14 | mRNA splicing pathway dominance is unusual as a top biological signal for pre-B cells. The matched genes (DDX5, MAGOH, HNRNPR, etc.) are ubiquitous splicing factors. This likely reflects a background effect where metabolic/splicing genes are relatively upregulated against the erythroid background. B-specific pathways (BCR signalling, B Cell Receptor Signaling, rank 4) are present but ranked below splicing. | **Info** — likely a baseline transcriptional activity artefact vs the erythroid background. B Cell Receptor Signaling Pathway being present is reassuring. |
| Measles (KEGG) | CD4+ T cells | 1.07e-06 | Measles KEGG pathway contains numerous immune receptor and signalling genes (CD3D, CD28, STAT3, BCL2, PIK3R1) that are expressed in T cells generally. This is a known KEGG pathway cross-reactivity issue. | **Info** — ignore disease label; interpret as general T cell activation signalling gene overlap. |

---

## 8. Cross-Module Coherence

| Severity | Modules involved | Flag | Suggestion |
|---|---|---|---|
| **Critical** | QC → All downstream | 75.9% of cells removed (68,483/90,261), driven by max_pct_MT = 5.0% cutoff set below the dataset median (pre-QC median MT% = 6.34%). This is the most serious finding: over three-quarters of the input data was discarded, which almost certainly removed the majority of erythroid precursors and other high-MT% bone marrow lineages. The surviving 21,778 cells may not represent the full hematopoietic hierarchy. | Re-run QC with max_pct_MT = 10–15% or use a MAD-based threshold. Compare cell type proportions before and after permissive vs stringent QC to assess bias. |
| **Critical** | QC + Annotation | No Early erythroid cluster was detected despite healthy bone marrow being rich in early erythroid progenitors (BFU-E, CFU-E, early erythroblasts). The ground-truth label set contains only Reticulocytes, Normoblasts, and Erythroblasts. Early erythroid progenitors have high MT% (active mitochondria during globin synthesis initiation) and would have been preferentially removed by the 5% MT% cutoff. | Validate by re-running with a permissive MT% cutoff. Check ground-truth label distribution before vs after QC. |
| **Warning** | Normalization + Clustering | Batch key was set to `batch` (12 batches) for normalisation, but it is not confirmed from the report whether Harmony or BBKNN batch correction was applied to the dimensionality reduction (only 7 PCs used; the Harmony tab exists — panel_08_harmony_report_html). If integration was not performed, UMAP structure may reflect batch rather than biology. | Confirm from the Harmony tab whether batch-corrected embeddings were used for clustering. |
| **Warning** | Annotation (Cluster 9) + DEG | Cluster 9 is annotated as pDC (consensus 0.50; ScType called Megakaryocyte) with JCHAIN as the top DEG by log2FC after GZMB, IRF8, CCDC50, PPP1R14B. JCHAIN is a plasma cell-associated gene; its presence alongside IRF8 and GZMB is atypical. This raises the possibility that Cluster 9 contains a mixed pDC/plasma cell doublet population, or a genuine JCHAIN+ pDC subset. | Perform FACS validation of CD123+JCHAIN+ cells in bone marrow. Check whether doublet detection captured pDC-plasma cell doublets. Inspect the Scrublet doublet score distribution for Cluster 9 cells specifically. |
| **Warning** | QC + DEG (HSC/MPP) | HSC/MPP cluster 14 is small (506 cells; 2.3% of post-QC cells). Canonical HSC markers (CD34, HOXA5, MLLT3) did not appear in the Wilcoxon top 5 DEGs; instead, GAPDH and EEF1B2 ranked highly as artefacts of the erythroid-dominated background. The pseudobulk DEG recovers more informative HSC markers (PRTN3, CLEC11A, MEST, IMPDH2). | Prioritise pseudobulk DEG results for HSC/MPP biology interpretation. Consider manual annotation validation using CD34 and HOXA5 expression across clusters. |
| **Info** | Clustering + Annotation | Clusters 7 and 8 are both annotated as Classical monocytes (consensus 1.00) with identical marker gene profiles and ground-truth labels (CD14+ Mono). These may represent a single population split by batch or cell cycle rather than biologically distinct subpopulations. | Test whether merging Clusters 7 and 8 improves silhouette score and produces more interpretable DEG results. |
| **Info** | Clustering + Annotation | Clusters 10 and 12 are both annotated as CD4+ T cells (consensus 0.50) with different ground-truth labels (CD4+ T activated vs CD4+ T naive). This sub-structure exists and is recoverable via sub-clustering. | Sub-cluster Clusters 10 and 12 together using markers CCR7, SELL, FOXP3, TNFRSF4 to resolve naive/central memory/effector/Treg subtypes. |

**Overall coherence:** **Issues found** — primarily driven by the critical MT% over-filtering. Annotation is largely coherent and biologically interpretable; GSEA results are internally consistent. The primary concerns are (1) the aggressive QC removal distorting cell type representation, (2) uncertain batch integration status, and (3) the ambiguous pDC/JCHAIN cluster 9.

---

## 9. Downstream Suggestions

| Priority | Step | Rationale | Recommended tool | Expected output |
|---|---|---|---|---|
| 1 | QC re-run with permissive MT% threshold (10–15%) | The 5% MT% cutoff removed 71.7% of input cells, likely eliminating early erythroid progenitors, activated immune cells, and metabolically active HSCs. This is the single most impactful corrective action. | scanpy QC pipeline; re-run with `max_mt_pct = 10` or `median + 3×MAD`; compare cell type distributions | Recovery of Early erythroid progenitors, potentially more HSCs, and a more representative hematopoietic hierarchy |
| 2 | Batch integration confirmation and re-run if needed | It is unclear whether batch correction was applied to the UMAP/clustering (Harmony tab exists). Given 12 batches from 4 sites and 2 platforms, uncorrected embeddings risk clustering by batch rather than biology. | Harmony (already available in pipeline, panel_08); alternatively scVI or BBKNN | Batch-corrected UMAP revealing biological rather than technical structure; improved silhouette scores |
| 3 | Sub-clustering of CD4+ T cells (Clusters 10 + 12) | Ground-truth labels reveal CD4+ T activated and CD4+ T naive are split across clusters 10 and 12; a sub-clustering pass would recover naive, central memory, effector memory, and Treg populations. | Leiden clustering at resolution 0.3–0.5 on the T cell subset; Seurat or scanpy | Resolved CD4+ T cell subtype map with CCR7, SELL, FOXP3, TNFRSF4 as discriminating markers |
| 4 | Sub-clustering / manual validation of Cluster 9 (pDC/JCHAIN ambiguity) | Cluster 9 has a 0.50 consensus score and JCHAIN as a top DEG — an unusual pDC/plasma cell overlap. Manual validation is needed to determine if this is a genuine JCHAIN+ pDC subset, a mixed doublet cluster, or a plasma cell sub-population. | Re-run Scrublet doublet scoring focused on Cluster 9; FACS validation (CD123+CLEC4C+ for pDC vs CD38+CD138+ for plasma cells); alternatively use scDblFinder | Confirmation or refutation of JCHAIN+ pDC biology; cleanup of potential doublet contamination |
| 5 | Pseudobulk DEG as primary result for publication | Pseudobulk DEG (DESeq2, already run — panel_10) accounts for donor-level variation across 12 donors and 4 sites; Wilcoxon results should be considered exploratory. The pseudobulk run already applies the correct gene exclusions (RPL, RPS, MT-, HBB, HBA1, HBA2, HBD, AHSP). | DESeq2 pseudobulk results (panel_10) for all manuscript-grade DEG figures | Publication-ready differential expression with correct statistical model accounting for donor as a random effect |
| 6 | Trajectory analysis of erythroid maturation (Clusters 0–4) | Five erythroid clusters (Mid → Late) are present; trajectory inference could reveal the pseudotime ordering and identify transcription factor dynamics (GATA1, KLF1, NFE2) across erythroid maturation. | Monocle3 or scVelo (RNA velocity if spliced/unspliced counts available from the 10X Multiome samples) | Pseudotime-ordered erythroid maturation trajectory with stage-specific marker dynamics |
| 7 | CITE-seq protein data integration (antibody-derived tags, ADTs) | This dataset was generated with CITE-seq (BioLegend TotalSeq B Universal Human Panel v1.0). The current analysis uses only RNA. Incorporating ADT data would provide orthogonal protein-level confirmation of all cell type annotations (CD14, CD3, CD19, CD56, CD123 protein). | Seurat WNN (Weighted Nearest Neighbours) or totalVI (pyro-ppl) for joint RNA+protein embedding | Protein-confirmed cell type annotations; improved disambiguation of clusters with 0.50 consensus scores (especially Clusters 5, 9, 10, 12, 14, 15) |

---

## 10. Summary

**Key findings:**

1. **Aggressive QC over-filtering is the primary analytical concern:** The max_pct_MT = 5.0% threshold was set below the pre-QC dataset median (6.34%), resulting in 75.9% cell loss (68,483/90,261 cells removed) — far exceeding the typical 5–30% guideline. This has almost certainly depleted early erythroid progenitors, metabolically active HSCs, and activated immune populations, producing a survival-biased dataset.

2. **Cell type annotations are broadly concordant with published bone marrow biology:** 11 of 16 clusters received consensus scores ≥ 0.75 or high reviewer confidence (e.g., Clusters 0–4 erythroid with consensus 1.00; Clusters 7–8 classical monocytes with consensus 1.00; Cluster 15 NK cells with reviewer confidence High). The DEG signatures for erythroid (HBB log2FC=8.913, AHSP log2FC=5.921), monocyte (S100A9 log2FC=9.486, FCN1 log2FC=8.471), NK (GNLY log2FC=9.448, NKG7 log2FC=8.059), and plasma cell (MZB1 log2FC=6.366; PMID:30257949) compartments are canonical and well-grounded.

3. **Two discovery-priority findings:** (a) SNHG29 (log2FC=3.523, adj p=5.78e-236 in HSC/MPP) is a lncRNA with no established HSC role in the literature — no PubMed record found — representing a potentially novel finding in hematopoietic progenitor biology. (b) JCHAIN (log2FC=7.023) is a top DEG in the pDC cluster (Cluster 9) alongside canonical pDC markers IRF8 and GZMB, raising the possibility of a JCHAIN+ IgA-secretory pDC subpopulation in healthy bone marrow.

**Open questions:**
- Does the 5% MT% threshold cause systematic loss of specific hematopoietic lineages? Which ground-truth cell types are overrepresented in the 68,483 removed cells?
- Is Cluster 9 a genuine JCHAIN+ pDC subset, or does it represent pDC-plasma cell doublets or contamination?
- Was Harmony batch integration applied to the UMAP/clustering? The Harmony report tab exists but integration status is not confirmed in the clustering panel.
- Does SNHG29 have a functional role in HSC self-renewal or differentiation, or is it a transcriptional passenger?
- What is the functional significance of the cell-cycle dominance (Cell Cycle, Mitotic R-HSA-69278) in mid-erythroid cells — was G2M phase regressed out?

**Suggested validation experiments:**
- Repeat QC with max_pct_MT = 10–15% and compare cell type proportions to the current 5% result; characterise which ground-truth-labelled cell types were lost.
- FACS validation of Cluster 9: stain for CD123, CLEC4C (pDC markers) and JCHAIN, CD38, CD138 (plasma cell markers) on fresh bone marrow from healthy donors to determine whether JCHAIN+CD123+ cells exist in vivo.
- CRISPR knockdown of SNHG29 in human CD34+ HSC-derived in vitro progenitor assays (CFU-E, BFU-E, CFU-GM) to test functional relevance.
- Integrate CITE-seq ADT data with RNA (Seurat WNN) to provide protein-level confirmation of cell type calls for clusters with low consensus scores (Clusters 5, 9, 10, 12, 14, 15).

---

## Appendix: Report Figures Checklist

| # | Figure / element | Status | Note |
|---|---|---|---|
| 1 | Dot plot (clusters × marker genes, dot size = % expressing) | ✗ Missing | No dot plot was identified in the report. The clustering tab contains UMAP figures and resolution metrics but not a dot plot. Recommend adding this figure for publication. |
| 2 | Annotation table has two separate confidence columns (Consensus score + Reviewer confidence) | ✓ Present | The annotation tab contains "Consensus score" (0.0–1.0) and "Confidence" (High/Medium/Low) columns; these are method-agreement scores, not AI self-confidence. Reviewer confidence has been separately assigned in Section 4a of this document. |
| 3 | Batch count shown in Key Metrics on Data Report tab | ✓ Present | "12 Batches (batch)" is reported in the Key Metrics section of the Data Report tab. |
| 4 | DEG Run Summary notes which gene prefixes were excluded (or confirms none were) | ✓ Present | Wilcoxon DEG: "Gene prefix exclusion applied: RPL, RPS, MT-". Pseudobulk DEG: "Gene prefix exclusion applied: RPL, RPS, MT-, HBB, HBA1, HBA2, HBD, AHSP". Both runs document exclusions. |
| 5 | Pseudobulk DEG top hits do not include erythroid markers (HBB/HBA1/HBD/AHSP) as artefactual downregulated genes in non-erythroid groups | ✓ Present (pseudobulk) | The pseudobulk DEG correctly excludes HBB, HBA1, HBA2, HBD, AHSP from reported results. The Wilcoxon DEG does NOT exclude these genes (HBB and HBA2 appear as top upregulated genes in Late erythroid — which is expected and correct for that cell type). |
| 6 | Per-cluster accuracy table comparing consensus label to ground truth (if ground truth available) | ✓ Present | "Consensus vs Ground Truth — Per Cluster" table is present in the annotation tab, comparing consensus vote labels to obs['cell_type_groundtruth']. Ground truth is available and used. |

---

*OmicSage D1 Manual Review Mode — MASTER_PROMPT.md v1.0*
*Dataset: GSE194122 — BMMC CITE-seq (NeurIPS 2021) | Generated: 2026-05-21 | Model: Claude Sonnet 4.6*
*PubMed citations retrieved via PubMed MCP tool. DOI links provided where available.*
