# Exploratory TabPFN comparison

## Scope

`tabpfn_targets.R` adds a separate branch to the existing `_targets.R` pipeline
and store. It fits AD, FTD and LBD/PD versus control, each subtype versus the
other dementias, and LBD/PD versus FTD (seven comparisons), using the nine primary
biomarkers plus age and sex from `roc_data_prepared` (`df_with_csf`). Existing
SuperLearner fitting functions and the original report are unchanged. The
comparison uses **normal dependencies** on the seven historical primary ROC
results. Their target declarations now freeze successful results and fail
before fitting if a baseline needs execution while frozen.

This is an exploratory, non-publication comparison, not a strictly paired
validation. It intentionally preserves existing upstream biomarker truncation,
complete-case selection and whole-comparison standardization before outer CV.
Evaluation rows therefore influence preprocessing. New TabPFN folds need not
match historical SuperLearner folds; historical cohorts/software may also differ.
Do not interpret differences as controlled learner superiority or use a paired
AUC test without a separately designed matched evaluation.

## Environment (user only)

For new TabPFN fits, the pinned GPU shell requires x86_64 Linux and an AMD
`gfx1030` GPU (verified hardware: RX 6900 XT). A different GPU architecture
needs a compatible runtime before it can fit models. To enter the GPU shell:

```sh
nix develop .#tabpfn-gpu
R
```

The shell preserves the original ADDF Nix/R package pin. It additionally sets
`ADDF_TABPFN_RSCRIPT` to an isolated, pinned R/Python GPU wrapper used only in
prediction subprocesses. The ordinary `nix develop` shell is still available for
non-TabPFN work. Do not install Python/R packages manually.

GPU packaging and its `uv.lock` were adapted from
`~/collab/tabpfn-test/environments/tabpfn-gpu/`; the R wrapper was adapted from
that project's `environments/tabpfn/r.nix`. The separate GPU nixpkgs and
pyproject tooling pins were copied from that project's `flake.lock`, without
upgrading the existing analysis pin. This does not depend on the other checkout
at runtime. Building the opt-in environment may acquire sizeable dependencies;
no weights or inference are invoked by shell entry or ordinary checks.

The existing external checkpoint must be present at:

```text
~/.local/share/tabpfn/models/tabpfn-v3.5-20260909.safetensors
SHA-256: ece4d67eadfea42eb0e610df5189bea60cb7f31073d81e9c7a019b76eacf0be3
```

The code checks the hash before loading. No model acquisition, hosted prediction
API, CPU fallback or automatic Python environment installation is configured.
The GPU child clears inherited R/Python/ROCm paths and ignores user profiles.

## Run (user only)

Use the existing configured `.env`, `DATA_DIR` and `CACHE_DIR`. First run the
synthetic offline checks if desired; they do not import TabPFN or load data:

```r
source("tests/tabpfn_analysis.R")
source("tests/targets_baseline_policy.R")
```

For one initial contrast and a resource check:

```r
targets::tar_make(names = tabpfn_roc_result_AD_vs_Control)
```

Build the original control comparisons and their PNG/CSV exports:

```r
targets::tar_make(names = tabpfn_comparison_files)
```

To export **all three comparison PNGs**, covering all seven contrasts, run:

```r
targets::tar_make(names = tabpfn_comparison_pngs)
targets::tar_read(tabpfn_comparison_pngs)
```

Files under `output/tabpfn/`:

- `subtypes_vs_control_comparison.png`
- `subtypes_vs_others_comparison.png`
- `lbd_vs_ftd_comparison.png`

Each shows cached SuperLearner and TabPFN side by side. Existing curve CSVs
are also retained. The Quarto report targets and sources have been removed:
no HTML/PDF rendering or TeX is involved. Existing cached predictions are reused;
missing or outdated prediction targets still run under normal targets rules.

New result targets include `tabpfn_roc_result_AD_vs_Others`,
`tabpfn_roc_result_FTD_vs_Others`, `tabpfn_roc_result_LBD_vs_Others`, and
`tabpfn_roc_result_LBD_vs_FTD`. Plots are in `tabpfn_comparison_others$plot`
and `tabpfn_comparison_lbd_ftd$plot`. New exports use `subtypes_vs_others` and
`lbd_vs_ftd` filename prefixes under `output/tabpfn/`.

This selected branch may rebuild outdated **input/preparation** dependencies.
It references `roc_result_AD_vs_Control`, `roc_result_FTD_vs_Control`, and
`roc_result_LBD_vs_Control` directly. There are no store reads inside targets.
Those three historical results use `tar_cue(mode = "never")` by default,
plus a guard that stops before fitting if targets attempts to execute them
(for example, after missing metadata or a previous error). Missing stored
objects fail to load rather than triggering replacement models. The four
primary differential-diagnosis baselines are also frozen, with a separate
`ADDF_REFIT_SL_SUBTYPES=true` opt-in for deliberate recomputation. All seven
cached baselines must be available to generate all comparison PNGs. Adjusted,
sex-stratified and subset-search SuperLearner targets remain unchanged.
A bare `tar_make()` can still rebuild unrelated
SuperLearner analyses; use the selected command above.

Only when deliberately authorizing recomputation of these three baselines:

```r
Sys.setenv(ADDF_REFIT_SL_CONTROLS = "true")
# Run the deliberately selected baseline targets with normal invalidation.
# Afterwards, restore the default frozen policy:
Sys.unsetenv("ADDF_REFIT_SL_CONTROLS")
```

The opt-in restores normal cues; it does not force an otherwise up-to-date
result to rebuild. Keep the existing store and metadata intact. Frozen
historical values need not match today's cohort, preprocessing or software.
This dependency repair leaves the per-fold GPU prediction implementation
unchanged and does not intentionally invalidate completed TabPFN predictions.

Outputs:

```text
output/tabpfn/subtypes_vs_control_comparison.png
output/tabpfn/subtypes_vs_control_curves.csv
```

Two matching panels use common subtype colours and axes, with per-panel AUC
labels. Curves average held-out TPR/FPR at the same 100 probability thresholds
used in the existing plots. AUC labels are arithmetic means of fold AUCs, not
pooled out-of-fold AUCs. TabPFN explicitly treats higher `.pred_1` probabilities
as the target diagnosis (fixed ROC orientation); cached historical SuperLearner
AUCs retain the original code's automatic pROC orientation. Plotting does not
call the existing plot helpers, which have shared output filenames.

## Execution and caching

- Crew is enabled, with separate CPU and GPU pools. The GPU pool defaults to
  one worker. Set `ADDF_TABPFN_GPU_WORKERS=2` before `tar_make()` only if the
  host has capacity for two concurrent jobs. Both share the same GPU, so
  monitor VRAM, host RAM and throughput before increasing concurrency.
- Each contrast's fit target executes ten stratified folds sequentially. The
  fold assignment and learner seeds are explicit and reproducible for an
  unchanged ordered input. Four ensemble members, temperature 0.9, automatic
  GPU precision, no probability balancing and no hyperparameter search.
- Each fold launches a fresh GPU R subprocess. It returns only ordinary R
  probabilities and diagnostics; exiting releases Python/model memory. Private
  temporary files carry the fold inputs/results and are removed on normal
  return or an R error. An abrupt worker/process kill can leave temporary files;
  handle the store and temporary directory as sensitive local analysis data.
- Three host threads and one preprocessing job per GPU fit leave CPU capacity
  for normal work. Default per-fold timeout is 600 seconds. Settings live in
  `tabpfn_default_config()` and are exposed by `tabpfn_config`.
- Per-fold cached results contain labels, probabilities, training/held-out row
  positions, seeds, timing, checkpoint/software identity, effective model
  settings and GPU allocation peaks. Row IDs are positions in
  `roc_data_prepared`, not persistent participant IDs. Fitted GPU objects are
  never serialized. There is no permanent participant-level prediction export.
- Changing plot code or deleting exported files reuses prediction targets.
  Changing fitting settings or input data invalidates the relevant fit targets.
  A failed contrast must rerun its folds; completed contrast targets stay cached.
- Memory peaks are PyTorch allocations/reservations, not whole-device VRAM.
  The reference project's n=250/four-member fixture used about 1 GiB reserved
  GPU memory, but that is not an ADDF capacity guarantee. Review the first
  contrast's diagnostics before expanding. CPU and GPU work can contend for
  host RAM despite separate scheduling pools.

## AD vs Control: SuperLearner weight pilot

This opt-in pilot makes a **new** fit of the original nine-learner library; it
cannot recover the historical weights. It does not depend on the historical
ROC fit targets. To get the refit's per-fold weights without running TabPFN,
run in the analysis R environment:

```r
Sys.setenv(ADDF_SL_WEIGHT_PILOT = "true")
targets::tar_make(names = ad_sl_baseline_weights_file)
targets::tar_read(ad_sl_baseline_weights)
```

To compare against the augmented library, run from the GPU shell described
above (set the opt-in again if this is a new R session):

```r
Sys.setenv(ADDF_SL_WEIGHT_PILOT = "true")
targets::tar_make(names = ad_sl_paired_weights_files)
targets::tar_read(ad_sl_paired_weights_files)
```

Outputs in `output/tabpfn/` use the `ad_vs_control_sl_` prefix: baseline and
paired weight CSVs, a side-by-side weight PNG, a fold-level paired AUROC CSV,
and a TabPFN weight summary CSV. Each outer fold has two vertical stacked bars:
the original library on the left and the library with TabPFN on the right.
Colours and stacking order are fixed by model. Segment height is the model's
weight in that fold. The paired CSV distinguishes failed fits from zero weights.
If you already have the paired CSV, regenerate just the plot without refitting:

```r
source("R/tabpfn_sl_weights.R")
stem <- "output/tabpfn/ad_vs_control_sl_paired_weights"
weights <- read.csv(paste0(stem, ".csv"))
ggplot2::ggsave(paste0(stem, ".png"), tabpfn_sl_plot_weights(weights),
                width = 11, height = 6, dpi = 150)
```

Baseline fits are reused. The augmented branch
adds TabPFN's probabilities on the same 20 inner splits per outer fold, then
scores the 10 outer holdouts (210 GPU fits). Completed outer fold branches can
be reused. Failed original learners are flagged separately from zero weights;
a TabPFN error stops its branch rather than assigning a zero weight.

Both ensembles use `method.AUC`, a rank-based objective that can have multiple
weight vectors with the same score. Read the weights alongside paired AUROC;
they are not feature importance. This pilot retains the exploratory branch's
whole-comparison scaling before cross-validation, so its AUROC is not a fully
held-out estimate. The synthetic CPU-only check is
`source("tests/tabpfn_sl_weights.R")` from the analysis R environment.

## Moving to another host

- Clone the repository and let Nix recreate the environment from `flake.lock`.
  Recreate any required local setup; do not copy `.direnv/` or assume `.Rlib/`
  is portable across hosts.
- Configure `.env` on the destination (see `env-sample`), setting `DATA_DIR`
  and `CACHE_DIR` to the destination paths. Move input files and the complete
  targets store, **including metadata**, separately through approved storage.
  Stop writes to the store before copying it. The store contains held-out
  predictions and other analysis data; keep it out of Git.
- Make the checkpoint available separately at the path above and verify its
  SHA-256. The pinned `gfx1030` runtime is required for new fits, not for
  opening existing PNGs. Copy any exports you want to keep without rebuilding.
- Use `targets::tar_outdated(names = tabpfn_comparison_pngs)` before building to
  check what would rerun. Moving paths or changing the runtime can invalidate
  cached targets. Leave the SuperLearner refit opt-ins unset to protect the
  historical results.

## Verification status

Assistant-side verification is static only: no R, R-backed language server,
Quarto, model inference, or raw-data inspection is permitted. Air parses/formats
the changed R sources without executing R. Nix checks are limited to parsing and
evaluating derivations, not building R packages (whose build phases execute R).
Synthetic JavaScript tests cover the pre-existing report, not this model adapter.
The new synthetic R checks are provided for user execution and remain unrun by
the assistant. Actual GPU integration, crew routing, resource use and rendered
plot correctness must be confirmed by the user.
