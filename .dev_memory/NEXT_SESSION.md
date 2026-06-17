# OmicSage — Next Session
> Written: 2026-06-17
> Phase: 11 — Portfolio & Outreach (Tutorial + Demo + Colab)

---

## Session Context

**Last thing completed:**
Phase 10 — Nextflow DSL2 orchestration complete.
All 4 modalities wrapped in DSL2. 33 modules with skip guards.
Streamlit UI reverted to Python-only with GPU toggle.
Tested: scRNA pipeline runs end-to-end via Nextflow on GSE166635.

**Files added/modified this session:**
```
main.nf                              ← updated (all 4 modalities, if/else routing)
nextflow.config                      ← updated (GPU flag, entrypoint bypass, memory)
pipeline/workflows/scrna.nf          ← updated (val channels)
pipeline/workflows/cite.nf           ← new
pipeline/workflows/multiome.nf       ← new
pipeline/workflows/spatial.nf        ← new
pipeline/modules/scrna/  (10 files)  ← updated (skip guard, full Python path)
pipeline/modules/cite/   (10 files)  ← new (skip guard)
pipeline/modules/multiome/(6 files)  ← new (skip guard)
pipeline/modules/spatial/ (7 files)  ← new (skip guard)
ui/_pages/p3_run.py                  ← updated (Python-only, GPU toggle)
ui/defaults.py                       ← updated (MODALITY_NF_NAME, NF_IMPLEMENTED)
```

**Current test count:**
```bash
conda activate omicsage
python -m pytest tests/ -q --tb=short
```
Expected: 1472 passing, 58 skipped. Confirm before writing any new code.

**Verify Nextflow still works:**
```bash
nextflow run main.nf \
  --config config/runs/GSE166635.yaml \
  --modality scrna \
  -profile local \
  -ansi-log false
```
Expected: all steps submit, harmony + pseudobulk print "disabled in config -- skipping".

---

## Today's Goals

Three deliverables this session (in order of priority):

### 1. Tutorial notebook
A Jupyter notebook `tutorials/01_scrna_quickstart.ipynb` that walks through
the full scRNA pipeline on GSE166635 (Wang et al. 2025 HCC dataset).

**Structure:**
- Section 0: Setup and installation check
- Section 1: Data overview (what GSE166635 is, expected outputs)
- Section 2: Run via Python runner (step by step with explanations)
- Section 3: Run via Streamlit (screenshots + instructions)
- Section 4: Run via Nextflow CLI (for HPC users)
- Section 5: Understanding the outputs (reports, checkpoints)
- Section 6: Customizing parameters (how to edit the YAML)

**Style:** Biologist-friendly. Explain what each step does biologically,
not just technically. Show expected outputs and what they mean.

### 2. Google Colab notebook
A notebook `tutorials/02_colab_quickstart.ipynb` that lets anyone run
OmicSage on Google Colab with free T4 GPU.

**Structure:**
- Cell 1: Clone repo + install dependencies (pip install from requirements)
- Cell 2: Mount Google Drive (for data + output persistence)
- Cell 3: Download GSE166635 test data (or upload from Drive)
- Cell 4: Run pipeline steps with OMICSAGE_GPU=1
- Cell 5: View reports inline

**Key constraints:**
- No Docker on free Colab — use Python runner directly
- Must install into Colab's default Python (not conda)
- Data must be on Google Drive or downloaded fresh each session
- GPU steps: set OMICSAGE_GPU=1 for totalVI/Tangram/MultiVI
- Keep install cell under 5 minutes

**Requirements file needed:** `requirements-colab.txt`
(lightweight subset of environment.yml that pip can install on Colab)

### 3. Demo video script
Not the actual recording — a structured script + checklist for the demo video
so you can record it yourself after the session.

**Script outline:**
- 0:00 — Landing on GitHub repo, README overview (30s)
- 0:30 — streamlit run ui/app.py, load GSE166635 (30s)
- 1:00 — Configure page walkthrough (45s)
- 1:45 — Run page: click run, show live log (60s)
- 2:45 — Report page: scroll through combined report tabs (60s)
- 3:45 — Terminal: show nextflow run command for HPC users (30s)
- 4:15 — Colab notebook: show it running in browser (30s)
- 4:45 — Wrap up: GitHub stars CTA (15s)

---

## Build Order

```
1. requirements-colab.txt        ← needed by Colab notebook
2. tutorials/01_scrna_quickstart.ipynb
3. tutorials/02_colab_quickstart.ipynb
4. docs/demo_script.md           ← video recording guide
```

---

## Key Design Decisions for Colab

- Use `pip install` not conda (Colab doesn't support conda well)
- Pin versions to match environment.yml exactly
- Store outputs in `/content/drive/MyDrive/OmicSage/` when Drive is mounted
- Fall back to `/content/OmicSage/` if Drive not mounted
- OMICSAGE_GPU auto-detected: check `torch.cuda.is_available()`
- Skip Docker-dependent features (Nextflow) in Colab notebook

---

## Files to Update at End of Session

```
.dev_memory/NEXT_SESSION.md    ← Phase 12 (README rewrite) or next priority
.dev_memory/CURRENT_STATUS.md  ← add tutorial + colab to completed
.dev_memory/PROGRESS.md        ← tick Phase 11 items
```
