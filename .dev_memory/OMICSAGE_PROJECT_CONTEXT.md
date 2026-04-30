# OmicSage — Complete Project Context
> Paste this file at the start of every new conversation in this project.
> Last updated: April 2026
> GitHub: https://github.com/[your-username]/OmicSage
> Domain: omicsage.io (confirmed available on Namecheap — do not buy yet, wait until Month 2-3)

---

## 1. What Is OmicSage

OmicSage is an **open-source, AI-assisted modular platform for single-cell multi-omics analysis**. It covers four data modalities:
- Single-cell RNA sequencing (scRNA-seq)
- Single-cell ATAC sequencing (scATAC-seq)
- Spatial transcriptomics (Visium, MERFISH, Xenium)
- Multiome integration (RNA + ATAC jointly)

The platform serves two users simultaneously:
- **Biologists**: no-code interface, guided workflow, biological interpretation in plain language
- **Bioinformaticians**: project templates, multi-project management, automated reports, no repeated code

**Name rationale**: "Omics" covers the full multi-omics scope. "Sage" reflects AI wisdom and guidance. Together: a wise guide through your omics data. Short, memorable, professional, easy to say in any language. Domain omicsage.io confirmed available.

---

## 2. The Developer

- **Background**: Bioinformatician, comfortable in Python (primary) and R (secondary)
- **Setup**: Windows PC + VSCode + WSL2 + Docker Desktop
- **Schedule**: 5 days/week, 6 hours/day = 30 hours/week, fully dedicated
- **Status**: Currently unemployed — OmicSage is the sole project
- **Goals** (in priority order):
  1. **Month 1-2**: Portfolio piece for job interviews (working tool + clean GitHub + demo)
  2. **Month 3-4**: Freelance/consulting ready (charge labs $500-2000 per dataset)
  3. **Month 6+**: Product/startup (SaaS, core facility licensing, or company sale)

---

## 3. Why We Are Building This — The Problem

The single-cell analysis landscape has three broken things:

1. **Fragmentation**: No single tool covers scRNA + scATAC + spatial + multiome in one pipeline
2. **Inaccessibility**: Biologists cannot analyze their own data without a bioinformatician
3. **Repetition**: Bioinformaticians rewrite the same QC/clustering code for every project

Existing tools each solve one piece but not the whole:

| Tool | What It Does | What It Misses |
|------|-------------|----------------|
| nf-core/scdownstream | scRNA pipeline, reproducible | No AI, no reports, no spatial/ATAC |
| CellWhisperer (Nature Biotech 2025) | Chat with scRNA data | Human only, no pipeline, no reports |
| Operon | Claude-powered HPC IDE | macOS only, no workflow engine, no reports |
| SCassist | LLM parameter suggestions | scRNA only, no UI, no reports |
| SCAPE (bioRxiv Nov 2025) | Multi-omics platform attempt | No reports, no AI-free mode, no multi-project |
| BioChatter (Nature Biotech 2025) | LLM framework for biomedicine | No pipeline, no single-cell focus, no reports |
| OmicsPilot | 2 niche tools (metabolic scoring + raw QC) | No pipeline, no AI, no reports, no clustering, no annotation, no spatial/ATAC |

**The gap OmicSage fills**: end-to-end pipeline + automated reports + AI interpretation + multi-modal + works without AI (graceful degradation) + multi-project management. No existing tool delivers all of these.

---

## 4. Competitor Note — OmicsPilot (omicspilot.com)

We checked OmicsPilot carefully. It is **not a real competitor** — it is validation that this market exists.

**What they actually have:**
- CellMetPro: scRNA-seq → metabolic profiling via COMPASS algorithm (v0.1.0 on PyPI)
- BioSeqFlow: raw FASTQ QC and preprocessing (Illumina/PacBio/ONT)
- EDA and DEG tools: listed as "Coming Soon" — not built

**What they don't have:** clustering, annotation, trajectory, spatial, ATAC, multiome, AI layer, report generation, slides, CLI project management, workflow engine, no-code UI.

**Strategic implication**: OmicsPilot proves demand for this category while leaving the actual integrated platform space completely open. CellMetPro (metabolic profiling) is a potential future integration module for OmicSage, not a competitor. Their website design is polished — use it as a visual benchmark for our own documentation site.

---

## 5. Architecture — Three Layers

```
┌─────────────────────────────────────────────────────┐
│  LAYER 3: AI Intelligence (OPTIONAL)                │
│  BioChatter + LLM → threshold suggestions,          │
│  cluster interpretation, PubMed RAG, narratives     │
├─────────────────────────────────────────────────────┤
│  LAYER 2: Report Engine (ALWAYS ON)                 │
│  Quarto + python-pptx → HTML/PDF reports +          │
│  PowerPoint slides after every analysis step        │
├─────────────────────────────────────────────────────┤
│  LAYER 1: Core Pipeline (ALWAYS ON)                 │
│  Nextflow DSL2 → QC → normalize → integrate →       │
│  cluster → annotate → downstream analysis           │
└─────────────────────────────────────────────────────┘
```

**Key design principle**: Setting `ai_features: false` in config.yaml disables the AI layer completely. The full pipeline still runs. Reports still generate. No paid API needed. This makes OmicSage genuinely free for any researcher globally.

---

## 6. Technology Stack

| Component | Technology |
|-----------|------------|
| Workflow engine | Nextflow DSL2 (reuse nf-core modules where possible) |
| Containers | Docker + Singularity (HPC compatibility) |
| Python analysis | Scanpy, scvi-tools, Muon, squidpy, BioChatter |
| R analysis | Seurat v5, Signac, ArchR, BayesSpace |
| AI middleware | BioChatter (RAG, PubMed, LLM routing) — do NOT build from scratch |
| Report generation | Quarto (HTML/PDF) + python-pptx (slides) |
| No-code UI | Streamlit web interface |
| CLI | Python Click-based project manager |
| Reproducibility | Git + DVC + structured JSON parameter logs |
| CI/CD | GitHub Actions |

---

## 7. Repository Structure

```
omicsage/
├── .dev_memory/                  ← SESSION MEMORY (update every day)
│   ├── NEXT_SESSION.md           ← Paste this at start of every chat
│   ├── CURRENT_STATUS.md         ← What is built and tested right now
│   ├── DECISIONS.md              ← Why we chose X over Y
│   ├── PROGRESS.md               ← Completed checklist
│   └── BLOCKERS.md               ← Broken or unclear things
│
├── config/
│   └── schema.yaml               ← Master config schema for all modalities
│
├── pipeline/
│   ├── modules/                  ← Individual Nextflow modules
│   │   ├── qc/
│   │   ├── processing/
│   │   ├── clustering/
│   │   ├── annotation/
│   │   └── downstream/
│   ├── workflows/
│   │   ├── scrna.nf              ← scRNA-seq workflow
│   │   ├── scatac.nf             ← scATAC-seq workflow
│   │   ├── spatial.nf            ← Spatial workflow
│   │   └── integration.nf        ← Multiome workflow
│   └── main.nf                   ← Entry point
│
├── ai/
│   ├── biochatter_client.py      ← BioChatter integration
│   ├── threshold_suggester.py    ← QC threshold AI suggestions
│   ├── cluster_interpreter.py    ← Marker → cell type via LLM
│   └── narrative_generator.py    ← Report text generation
│
├── reports/
│   ├── templates/                ← Quarto templates per analysis type
│   └── slides/                   ← python-pptx slide templates
│
├── ui/
│   └── app.py                    ← Streamlit web interface
│
├── cli/
│   └── omicsage.py               ← Click CLI for bioinformaticians
│
├── data/
│   ├── benchmark/                ← Test datasets (GSE166635 HCC etc.)
│   └── references/               ← Cell type reference datasets
│
├── tests/                        ← pytest + nf-test
├── docs/                         ← Documentation site source
├── docker/
│   ├── Dockerfile.python         ← Python analysis environment
│   └── Dockerfile.r              ← R analysis environment
│
├── .github/
│   └── workflows/                ← CI/CD GitHub Actions
│
├── environment.yml               ← Conda environment
├── nextflow.config               ← Nextflow configuration
├── docker-compose.yml            ← Local development setup
└── README.md
```

---

## 8. Modality Pipelines

### scRNA-seq
QC (genes/cell, MT%, doublets via Scrublet, ambient RNA via SoupX) → normalization (scran) → HVG selection (top 2000) → PCA → batch correction (Harmony/scVI) → Leiden clustering → cell type annotation (SingleR + optional LLM) → DEG (Wilcoxon + pseudobulk DESeq2) → GSEA (GO/KEGG/Reactome) → trajectory (Slingshot/PAGA) → survival analysis → drug targets (DGIdb)

### scATAC-seq
Fragment QC (TSS enrichment, fragment size distribution) → peak calling (MACS3) → TF-IDF/LSI → ArchR/Signac clustering → motif enrichment (chromVAR) → TF activity (JASPAR) → gene activity scores → co-accessibility → integration with scRNA

### Spatial Transcriptomics
Spot QC → spatially variable genes → spatial normalization → BayesSpace clustering → deconvolution (RCTD/SPOTlight) → niche analysis → cell-cell communication (CellChat/LIANA) → ligand-receptor mapping

### Multiome Integration
Cross-modal QC → WNN (Seurat v5) / MOFA+ / MultiVI → joint embedding → cross-modal label transfer → SCENIC+ GRN inference → multi-modal trajectory

---

## 9. Development Roadmap

### Phase 0 — Foundation (Week 1-2) ← WE ARE HERE
- [x] GitHub repo created
- [x] Name confirmed: OmicSage
- [x] Domain confirmed available: omicsage.io
- [ ] Full repo structure scaffolded
- [ ] Docker base images (Python + R)
- [ ] CI/CD via GitHub Actions
- [ ] YAML config schema defined
- [ ] Benchmark dataset downloaded (GSE166635 HCC scRNA-seq)
- [ ] .dev_memory/ system initialized and filled
- [ ] README with project description and badges

### Phase 1 — Core scRNA Pipeline (Week 2-6)
- [ ] Data ingestion: 10x MEX, H5, AnnData auto-detection
- [ ] QC module: MT%, genes/cell, SoupX ambient RNA, Scrublet doublets
- [ ] Normalization + HVG selection
- [ ] PCA + UMAP + t-SNE
- [ ] Harmony + scVI integration
- [ ] Leiden clustering with resolution sweep
- [ ] SingleR annotation
- [ ] DEG: Wilcoxon + pseudobulk
- [ ] GSEA: GO/KEGG/Reactome
- [ ] **MILESTONE**: Reproduce key findings of Wang et al. 2025 HCC paper from raw counts

### Phase 2 — Report Engine (Week 6-9)
- [ ] Quarto QC report template
- [ ] Quarto analysis report template (clustering, DEG, pathway)
- [ ] python-pptx slide deck generator
- [ ] Auto-figure captioning from metadata
- [ ] Auto-methods text from config + software versions
- [ ] **MILESTONE**: Biologist receives complete PDF + slides from one command, no coding

### Phase 3 — AI Layer (Week 9-13)
- [ ] BioChatter integration
- [ ] QC threshold suggestion via augmented prompts
- [ ] Cluster interpretation (markers → LLM → cell type + evidence)
- [ ] PubMed RAG tied to DEG results
- [ ] Biological narrative generator for reports
- [ ] AI audit log (all LLM calls logged with input/output/model/timestamp)
- [ ] Multi-LLM support (Claude, GPT-4o, local Ollama)
- [ ] **MILESTONE**: AI report narrative groundedness score > 0.85

### Phase 4 — scATAC Module (Week 13-16)
- [ ] Fragment file ingestion and QC
- [ ] Peak calling (MACS3)
- [ ] LSI + ArchR/Signac clustering
- [ ] Motif enrichment + chromVAR
- [ ] Gene activity scores
- [ ] scATAC report template

### Phase 5 — Spatial Module (Week 16-19)
- [ ] Visium data ingestion
- [ ] Spatially variable genes
- [ ] BayesSpace clustering
- [ ] RCTD deconvolution
- [ ] Spatial report template

### Phase 6 — Multiome Integration (Week 19-22)
- [ ] WNN joint embedding
- [ ] MOFA+ integration
- [ ] SCENIC+ GRN inference
- [ ] Joint report template

### Phase 7 — User Interfaces (Week 22-25)
- [ ] Streamlit web UI (drag-drop upload, guided config, progress tracking)
- [ ] Click CLI (create-project, run, list, compare commands)
- [ ] Project dashboard
- [ ] Shared protocol library (Markdown files per biological question)

### Phase 8 — Benchmarking + Paper (Week 25-30)
- [ ] Benchmark on 5 published datasets
- [ ] User study (biologists + bioinformaticians)
- [ ] Paper written and submitted (target: Nature Methods or Bioinformatics)
- [ ] bioRxiv preprint posted
- [ ] v1.0 public release
- [ ] omicsage.io website launched

---

## 10. Session Management System

### The .dev_memory/ System
At the END of every session, update these files (10 minutes max):
- `NEXT_SESSION.md` — what to paste at the start of the next chat
- `CURRENT_STATUS.md` — what is currently built and tested
- `PROGRESS.md` — tick off completed items

At the START of every session, paste `NEXT_SESSION.md` into the chat before anything else.

### NEXT_SESSION.md Template
```markdown
## Session Context
Date:
Phase:
Last thing we completed:
File we were working on last:

## Today's Goal
[one clear, specific deliverable — not more than one]

## Known Issues From Last Session
[anything broken or unclear]

## Files Modified Last Session
[list them with their paths]

## Verify Last Session Works
[one command to run to confirm previous work still runs]

## Relevant Context
[anything Claude needs to know that is not covered above]
```

### Daily Session Structure
```
[Paste NEXT_SESSION.md — 2 min]
[Confirm today's single goal — 5 min]
[Build — 4-5 hours]
[Test what we built — 30-45 min]
[Update NEXT_SESSION.md + CURRENT_STATUS.md — 10 min]
```

### Weekly Rhythm
- **Monday**: Plan the week, review last week's output, update PROGRESS.md
- **Tuesday-Thursday**: Build — one module or feature per day, one goal per session
- **Friday**: Test, fix, document — no new features on Fridays

---

## 11. Key Decisions Already Made

| Decision | Choice | Reason |
|----------|--------|--------|
| Project name | OmicSage | Available on GitHub and omicsage.io, professional, memorable |
| Workflow engine | Nextflow DSL2 | Industry standard, HPC/cloud ready, nf-core module reuse |
| AI middleware | BioChatter | Don't build from scratch — already handles RAG/KG/LLM routing |
| Report format | Quarto + python-pptx | Quarto for HTML/PDF, pptx for editable slides |
| No-code UI | Streamlit | Fastest path to a working biologist-facing interface |
| Containerization | Docker + Singularity | Docker for development, Singularity for HPC |
| Language mix | Python primary, R secondary | Python for orchestration + AI; R for Seurat/Signac/ArchR |
| License | MIT | Maximum adoption, commercial use allowed |
| AI-free mode | Mandatory by design | Core pipeline must work without any API key |
| Primary test dataset | GSE166635 (HCC) | Wang et al. 2025 paper — we know the expected results |
| Upstream alignment | Not included | Defer to nf-core/scrnaseq; OmicSage starts from count matrices |
| OmicsPilot | Not a competitor | Different scope — potential future ecosystem partner |

---

## 12. Competitive Advantage

OmicSage is differentiated from all existing tools on five dimensions that nothing else currently delivers together:

1. **True multi-modal**: scRNA + scATAC + spatial + multiome in one pipeline with joint AI interpretation
2. **Automated reports + slides**: first-class output after every step, not an afterthought
3. **Bioinformatician project management**: multi-project CLI, templates, shared configs
4. **Graceful AI-free degradation**: full functionality without any paid subscription or API key
5. **Literature-grounded narrative**: BioChatter RAG links your actual DEGs to live PubMed evidence in the report

---

## 13. Commercial Strategy

### Phase 1: Open Source + Visibility (Month 1-2)
- Public GitHub, MIT license
- Used for job interviews and as portfolio evidence
- Post on bioRxiv as preprint as soon as Phase 1-2 are done
- Target: first GitHub stars and watchers

### Phase 2: Freelance (Month 3-4)
- Offer single-cell analysis services to wet-lab groups and small biotech
- OmicSage runs the pipeline, you provide scientific interpretation
- Pricing: $500-2000 per dataset depending on complexity and turnaround
- omicsage.io website launched as service landing page

### Phase 3: Product (Month 6+)
Options to evaluate based on traction:
- **SaaS**: hosted version with subscription ($99-499/month per lab)
- **Core facility licensing**: annual license for institutional deployment
- **Consulting + tool bundle**: analysis package deals with tool access
- **Acquisition target**: pharma or genomics company buys the platform

---

## 14. Important Reference Papers

- **Wang et al. 2025** (npj Precision Oncology, DOI: 10.1038/s41698-025-00952-3): HCC single-cell paper — primary benchmark. Dataset: GSE166635. We know the expected outputs.
- **CellWhisperer** (Schaefer et al., Nature Biotechnology, Nov 2025): nearest competitor for AI interpretation layer — study their CLIP + LLaVA architecture
- **BioChatter** (Lobentanzer et al., Nature Biotechnology, Jan 2025): the framework we use for AI middleware — read before implementing Phase 3
- **SCAPE** (bioRxiv Nov 2025): closest competitor to the full platform — monitor actively for publication
- **SCassist** (Bioinformatics, Aug 2025): inspiration for augmented-prompt threshold suggestion
- **nf-core/scdownstream**: Nextflow pipeline we extend and build on top of

---

## 15. Claude Project Instructions
> This is the text to paste into the Claude Project instructions field when creating the project.

```
PROJECT: OmicSage
An open-source, AI-assisted modular platform for single-cell
multi-omics analysis (scRNA-seq, scATAC-seq, spatial transcriptomics,
multiome integration).

GITHUB: https://github.com/[your-username]/OmicSage
DOMAIN: omicsage.io (available, not yet purchased)

DEVELOPER: Bioinformatician
SETUP: Windows PC + VSCode + WSL2 + Docker Desktop
LANGUAGES: Python (primary) + R (secondary)
STACK: Scanpy, Seurat v5, Nextflow DSL2, BioChatter,
       Streamlit, Quarto, python-pptx, Docker
SCHEDULE: 30h/week, fully dedicated, sole project

STRATEGY:
- Month 1-2: Working pipeline → job interviews + GitHub portfolio
- Month 3-4: Freelance ready → charge labs to analyze their data
- Month 6+:  Product/startup → SaaS or licensing

AT THE START OF EACH SESSION:
The developer will paste the contents of .dev_memory/NEXT_SESSION.md
Always read it fully before doing anything else.
Always confirm the single goal for the session before writing any code.
Never work on more than one feature or module per session.

CURRENT PHASE: 0 — Foundation
```

---

## 16. How To Start Every New Chat Session

Copy-paste this block at the very top of every new conversation, then paste NEXT_SESSION.md contents below it:

```
I am building OmicSage — an open-source AI-assisted single-cell
multi-omics analysis platform. Full context is in the project
instructions. Here is where we are right now:

[PASTE NEXT_SESSION.md CONTENTS HERE]
```

---

*Generated from the initial planning conversation, April 2026.*
*This is the source of truth for all project context.*
*Update Section 9 checkboxes as phases are completed.*
*Replace [your-username] with your actual GitHub username throughout.*
