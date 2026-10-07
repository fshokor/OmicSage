# Single-section cell2location analysis

Select one section through `spatial.ingest.sample_id`, with an explicit
`spatial.ingest.library_key` (for Kuppe, `patient_region_id`). Ingestion retains
only matching observations and their spatial image library. An unknown sample
or incompatible image-library ID raises an error rather than silently selecting
the wrong tissue. Omit `sample_id` to retain the existing all-section behavior.

## Kuppe control_P1 CPU pilot

`config/runs/kuppe_heart_P1_pilot.yaml` selects 4,279 spatial observations while
retaining all 41,663 nuclei from the four reference donors. QC metrics are
measured without filtering spots. The pilot uses 10 reference epochs, 50 spatial
epochs, 512-observation batches, and 50 posterior draws. It tests execution and
resource use; its abundances are not a converged or validated biological result.
Downstream analyses and imputation are disabled.

```bash
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
python run_spatial_pipeline.py --config config/runs/kuppe_heart_P1_pilot.yaml --to-step deconvolve
```

Outputs and HTML reports are under `data/processed/kuppe_heart_P1_pilot` and
`reports/kuppe_heart_P1_pilot`. Runtime and peak resident memory for the measured
run are in `logs/kuppe_heart_P1_pilot_resources.txt`; execution details are in
`logs/kuppe_heart_P1_pilot.log`.

## Full-training configuration

`config/runs/kuppe_heart_P1_cell2location.yaml` supplies 250 reference epochs,
30,000 spatial epochs and 1,000 posterior draws in a separate output directory.
This configuration retains CPU-friendly batches of 512 and therefore does not
exactly reproduce the tutorial's full-batch training. It is not run automatically
as part of the pilot. Review convergence and reconstruction before interpreting
results; use measured pilot timings to estimate CPU feasibility.

## Implementation repairs

- The joint cell2location path accepts the logging interval passed by the runner.
- Reference counts are restricted to shared non-mitochondrial genes before gene
  filtering and reference fitting. Unused reference embeddings/expression layers
  are omitted from the fitting copy to reduce memory.
- Reference signatures and reference/spatial training histories are retained in
  checkpoint `uns` as `cell2location_reference_signatures`,
  `cell2location_reference_history` and `cell2location_spatial_history`.
- Region clustering no longer caps the requested neighbor count by the number
  of cell types.

The existing posterior cell-type-specific expression layers and full model/QC
exports from the best-practices notebook are not implemented by these changes.
The pilot does not run composition clustering, neighborhood or communication
analyses. The existing `fit_validated: false` status remains in provenance.

## Measured pilot outcome

The control_P1 pilot completed with exit status 0 in 19 min 39 sec on CPU,
with maximum resident memory 9,042,188 KiB (8.62 GiB). Reference fitting plus
export took about 11.9 min; spatial training took about 5.5 min. The checkpoint
contains 4,279 spots, one image library, 11 cell types and 14,364 fitting genes.
All abundances are finite and nonnegative; no cell type has entirely zero weight.
The spatial loss is still decreasing after 50 epochs, and mean abundances are
similar across types, so this execution pilot is not suitable for biological
interpretation. `pilot_summary.json` and `pilot_training_history.png` in its
report directory preserve the diagnostics.

Linear extrapolation from these short runs suggests about five hours for 250
reference epochs and roughly 55 hours for 30,000 spatial epochs on this CPU.
These are rough estimates, not guarantees; batches, convergence and machine
load affect runtime. Full training was not launched.

Validation: the selected spatial test suite passed 81 tests with one optional
test skipped; the final four focused regression checks also passed, including
checkpoint history serialization and region-neighbor selection. The diff
whitespace check passed.
