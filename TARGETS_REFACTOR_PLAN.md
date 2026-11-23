# Targets dependency repair and hardware-efficiency refactor

Status: scope narrowed after implementation discussion. Only the normal-target
dependency repair and small fail-closed historical-baseline policy are retained.
The in-progress warm-process supervisor, protocol and scheduler-budget machinery
were removed at the user's request. Hardware changes below remain proposals,
not descriptions of the current implementation. Research date: 2026-11-16.
Supersedes the earlier performance proposal in `TABPFN_PERFORMANCE_PLAN.md`
where recommendations differ.

No R, R-backed language server, targets workflow, model inference, or raw-data
inspection was performed. Research used official online documentation and
release-tagged source. Nix expression evaluation (not package builds) resolved
the existing analysis package names to **targets 1.10.1, crew 1.0.0, tarchetypes
0.12.0**. These are the flake's package declarations, not a verification of the
versions loaded in the user's live session or mutable `.Rlib`.

## 1. What the official guidance establishes

### Dependencies should be ordinary target references

Use an upstream target's symbol in a downstream command. Targets then owns
value loading, dependency edges, invalidation and scheduling. Calling
`tar_read_raw()` on the active pipeline's same store inside an ordinary target
is not a supported substitute. The 1.10.1 implementation explicitly checks this
with `tar_assert_allow_meta()` before accessing storage. Literate report
integration has special handling; it does not justify this snapshot helper.

Sources: [tar_target reference][target], [tar_read reference][read],
[1.10.1 tar_read_raw source][read-source].

### Selection is not a no-refit policy

`tar_make(names = tabpfn_comparison_files)` normally traverses that target's
upstream dependencies and can rebuild outdated upstream targets. `shortcut =
TRUE` assumes upstream targets are current; it is not a sound substitute for an
explicit frozen-baseline policy. Do not recommend it as the routine repair.

Source: [tar_make reference][make].

### Frozen execution requires care

`tar_cue(mode = "never")` suppresses routine invalidation of a previously
successful target. The documentation explicitly allows execution if metadata
is absent or the previous run errored; target-class changes are also listed
among the early invalidation checks. Missing stored output can instead leave
an unloadable target when metadata is present. In short, a cue is not an
absolute no-refit guarantee or protection against deleting files.

Use it together with a command-level fail-closed policy if automatic historical
baseline refitting must be prohibited. Never edit targets' private metadata or
object files to implement that policy.

Sources: [tar_cue reference][cue], [1.10.1 tar_cue source][cue-source].

### Crew parallelism must actually be enabled

`tar_option_set(controller = ...)` supplies the controller used by `tar_make()`.
A `crew_controller_group()` can route individual targets with
`tar_resources_crew(controller = "cpu"/"gpu")`. Named routing is task
placement, not an automatic GPU memory/exclusivity reservation. `workers`
limits concurrent worker processes, not physical GPU count. Two fits can share
one GPU if the workload fits; safety and throughput must be measured.

Sources: [targets crew manual][targets-crew], [crew groups][groups],
[local controller reference][local].

### Persistent crew workers do not make our child processes persistent

Crew 1.0.0 already defaults to `tasks_max = Inf` and `seconds_idle = 300`, so it
can reuse a worker for multiple tasks. This does not reuse the fresh subprocess
launched by `processx::run()` in every fold. Merely increasing `seconds_idle`
or setting `tasks_max = Inf` again will not remove that startup overhead.

The pinned local launcher uses `file.path(R.home("bin"), "Rscript")`.
`r_arguments` supplies arguments, not a different executable. Moving the GPU
runtime directly into a different-R crew worker would need a custom launcher
or environment redesign; it is not a one-line controller option.

Sources: [crew 1.0.0 controller source][local-source],
[crew 1.0.0 launcher source][launcher-source], [processx process][processx].

### Choose task granularity by cost, not by aesthetic preference

The targets manual recommends batching short jobs to amortize scheduling and
transfer overhead while retaining enough independent tasks for concurrency.
Neither one target per fold nor one target per whole study is universally best.
Per-target `packages` overrides avoid loading unrelated libraries. Memory,
retrieval and storage options control R object transport/lifetime, not Python
model/GPU lifetime. The pinned 1.10.1 reference should govern available options;
do not copy newer documented features without checking compatibility.

Sources: [targets performance manual][performance],
[dynamic branching manual][dynamic], [1.10.1 tar_target source][target-source].

## 2. Current-code findings

1. **Confirmed failure:** `tabpfn_cached_sl` calls `tabpfn_read_cached_sl()`,
   which reads the current store with `tar_read_raw()` during target execution.
   Its broad error handler then misleadingly describes the rejected access as
   an unavailable baseline. This is a dependency-design bug, not evidence of
   missing cached results.
2. **Current scheduling disconnect:** `_targets.R` currently contains
   `# controller = crew_controllers`. CPU/GPU controllers are defined, but are
   not attached to the pipeline there. Resource routing alone cannot activate
   them. This may be a deliberate debugging change by the user; leave it intact
   until the scheduling repair is approved. It does not prove how a previous
   invocation was configured.
3. **Confirmed repeated work:** the three contrast targets each loop through
   ten folds, launching a fresh GPU R/Python process for every fold. A complete
   run entails thirty child startups, repeated checkpoint hashing, runtime
   imports and estimator construction. Their relative runtime shares remain
   unmeasured.
4. **Confirmed serialization when enabled:** the declared GPU pool has one
   worker and each contrast's fold loop is sequential. Even restoring crew
   leaves only one active fit unless its configured worker limit changes.
5. **Potential oversubscription:** the normal pool uses almost all available
   cores, while GPU children request three host threads each. Increasing GPU
   concurrency without reserving CPU capacity can slow both pools.
6. **Excess worker setup:** the GPU orchestration targets inherit the full
   global SuperLearner/report package list. Their GPU subprocess is separate
   and has a much smaller R environment. Target-specific package declarations
   can reduce unnecessary orchestration-worker package loading, though exact
   savings require measurement.
7. **Insufficient telemetry:** existing diagnostics have total fold time and
   PyTorch memory peaks, but not stage timings or live child progress. A
   dispatched message does not establish that GPU inference has begun.

The later `qs` messages are not an independently established root cause; no
serializer/dependency change is proposed on that basis. No package upgrade is
required for the dependency repair or the initial scheduling/process refactor.

## 3. Recommended design

### A. Fix the dependency graph first, with no model reruns

Replace the snapshot read helper with an ordinary aggregation target, preserving
its existing output shape if useful:

```r
# Design sketch only: not an implementation patch.
tar_target(
  tabpfn_cached_sl,
  list(results = list(
    AD = roc_result_AD_vs_Control,
    FTD = roc_result_FTD_vs_Control,
    LBD = roc_result_LBD_vs_Control
  ))
)
```

Remove its `mode = "always"` cue and the same-store read helper. Dependency
hashes, not hidden store reads, should decide when this aggregation/plot updates.

For the three historical control-comparison result targets only:

- Keep their exact existing names and existing stored values.
- Make historical-baseline mode explicit, enabled by default for this comparison.
- In that mode, apply `tar_cue(mode = "never")` AND make any attempted execution
  stop with an actionable frozen-baseline error **before calling any fitting
  function**. This protects missing/errored metadata and similar edge cases.
- Retain normal fitting behind an explicit opt-in recomputation mode. Other
  SuperLearner analyses retain their current execution policy.
- Do not globally apply `mode = "never"`, change target classes to dynamic
  patterns, rename the baseline targets, or clear the store.
- Explain that this intentionally compares a historical baseline with the
  current TabPFN analysis, not necessarily the same cohort/folds/software.

The command guard is our project policy, not a built-in targets feature.
It must block fitting unconditionally while frozen, rather than merely checking
whether an arbitrary file exists. If a successful baseline cannot load, fail;
do not silently fit a replacement.

This is the preferred normal-target approach. A separately exported versioned
baseline file is an optional later archival solution, not required for the fix.
If ever used, export outside the running pipeline and track it with a file
input target. File tracking detects changes; it does not make artifacts immutable.

Scope this patch to the baseline aggregation, frozen-baseline declarations,
tests and documentation. Keep the GPU fit commands/helpers unchanged so fixing
the plot's dependencies need not invalidate completed TabPFN predictions.

### B. Restore controlled target-level concurrency

When approved, reconnect the existing controller group. Keep normal targets on
`cpu` and fitting targets on `gpu`. Begin the bounded GPU concurrency comparison
with one versus two workers, not unrestricted workers. Keep folds sequential
inside each target so concurrency has one clear owner.

Reserve host CPU capacity for GPU preprocessing/PyTorch instead of letting
both pools independently assume all cores are available. Make worker/thread
budgets explicit and configurable; compare throughput and memory, not just
instantaneous GPU utilization. Crew does not enforce a VRAM limit itself.

Override packages for new orchestration/aggregation targets with a minimal
explicit list; use namespace-qualified calls. Avoid changing global package
options for existing analyses. Retain current storage format and conservatively
retain transport/memory defaults until measurements justify another setting.

### C. Run one warm GPU subprocess per contrast

Keep the existing analysis R and isolated GPU R/Python installations. A small
runner module should expose one interface, conceptually:

```text
predict_tabpfn_batch(prepared_data, fold_manifest, fitting_config)
  -> fold probabilities + diagnostics
```

Internally:

1. Transfer the prepared comparison data and fold manifest once via private
   temporary files, not thirty separately serialized train/test payloads.
2. Start the GPU R subprocess once per contrast (three startups, not thirty).
3. Verify the checkpoint and initialize Python/PyTorch/GPU once at child startup.
4. For each fold, build a fresh estimator with the existing seed and four-member
   settings; train only on its training rows and predict its held-out rows.
5. Release the estimator/fold resources between folds. Never reuse trained
   context or fitted preprocessing across folds. Checkpoint tensor reuse is not
   promised: validate supported pinned TabPFN behavior separately if needed.
6. Return ordinary R predictions/diagnostics only. Compute ROC metrics in the
   original analysis R, avoiding new GPU-environment package requirements.
7. Keep aggregate-only fold/stage progress visible in dedicated worker/job logs.
   Child stdout does not automatically become the main targets console when
   crew runs the target remotely; document where users monitor it.
8. Supervise the process externally: preserve the per-fold timeout, add a
   bounded startup timeout and ensure interrupt/error cleanup. A single
   600-second timeout for all ten folds would be an incorrect replacement.
   Use process progress events to enforce fold deadlines, without introducing
   a general persistent service/RPC backend.

Retain three contrast-sized fit targets initially. Their fit costs and the
reference probe justify testing this coarse batching before adding complexity.
If real timings show poor load balance, test small fold batches later. More
batches trade better caching/retry/load balancing for more startup/transfer;
choose empirically. Native GPU crew workers or persistent child pools across
separate targets are deferred, not prerequisites for useful gains.

### D. Separate fitting and presentation invalidation

Keep preprocessing/fold manifests, prediction batches, ROC summaries, plot
objects and file exports conceptually separate. Use normal dependencies for
all of them. Preserve the same preprocessing, ordered row positions, explicit
fold assignments and learner seeds when relocating work. Stable seeds must not
be derived from worker number, completion order or newly chosen batch position.

Fit configuration should not contain plot paths/styles or observation-only
monitoring settings unnecessarily. Changing plot code should rebuild plots,
not prediction batches. Export targets remain `format = "file"` so deleting
a PNG/CSV regenerates the export from cached predictions.

## 4. Implementation order and acceptance gates

1. **Dependency repair:** implement A only. Add synthetic tests for normal
   dependency edges and frozen-baseline behavior (successful baseline, missing
   metadata, errored metadata, missing object). The selected plot build must
   never call SuperLearner and must reuse available TabPFN predictions.
2. **Observability:** add timings for process startup, checkpoint verification,
   imports, GPU initialization, fitting, synchronized prediction and transfer;
   record end-to-end time, process RSS and GPU allocation/reservation. Whole-device
   memory/utilization should be measured separately by the user. Logs must not
   include participant identifiers, features, labels or prediction vectors.
3. **Warm-process change:** implement C with a fake prediction adapter in
   synthetic user-run tests, covering exact training/validation membership,
   seeds, class mapping, failure propagation, deadlines and cleanup. Compare
   warm-process and existing process-per-fold predictions on fixed fixtures.
4. **Concurrency tuning:** with crew connected, compare one versus two warm
   contrast jobs. Adopt two only if actual wall time improves without excessive
   memory, failures or prediction drift. Tune host threads only afterward.
5. **Handoff:** retain the single entry target `tabpfn_comparison_files`, the
   two-panel plot and numeric curve export. Document how frozen baselines are
   deliberately refreshed, how to monitor progress, and how to retry failures.

The reference project's parallel probe measured 1.77x pooled throughput for two
GPU jobs, with 3.40 GiB sampled whole-device VRAM, on its small fixture. That is
supporting evidence for a bounded trial, not an ADDF performance/capacity promise.
See `~/collab/tabpfn-test/docs/findings/2026-10-10-pilot-parallelism.md`.

Record current completed prediction outputs before benchmarking using supported
user-run operations outside target execution. Do not delete/refit successful
results just to repair baseline aggregation. Runner code changes normally
invalidate dependent fit targets; benchmark deliberately rather than implying
that a performance refactor cannot trigger new TabPFN work. Never refit
SuperLearner as a side effect.

All actual R/synthetic/model/targets acceptance tests must be run by the user
under AGENTS.md. Assistant-side verification is static only. Do not declare
runtime correctness or speedup from documentation or Air/Nix evaluation.

## 5. Decisions needing confirmation before implementation

- Approve frozen historical mode for exactly the three control baseline result
  targets, with a fail-closed guard and explicit opt-in to recomputation.
- Approve reconnecting the commented controller group and a bounded two-worker
  GPU trial; retain separate configurable CPU/thread budgets.
- Approve contrast-sized warm batches as the first runtime refactor, keeping
  the current Nix package pins and isolated GPU R installation.

[read]: https://docs.ropensci.org/targets/reference/tar_read.html
[read-source]: https://github.com/ropensci/targets/blob/1.10.1/R/tar_read_raw.R
[target]: https://docs.ropensci.org/targets/reference/tar_target.html
[target-source]: https://github.com/ropensci/targets/blob/1.10.1/R/tar_target.R
[cue]: https://docs.ropensci.org/targets/reference/tar_cue.html
[cue-source]: https://github.com/ropensci/targets/blob/1.10.1/R/tar_cue.R
[make]: https://docs.ropensci.org/targets/reference/tar_make.html
[targets-crew]: https://books.ropensci.org/targets/crew.html
[groups]: https://wlandau.github.io/crew/articles/groups.html
[local]: https://wlandau.github.io/crew/reference/crew_controller_local.html
[local-source]: https://github.com/wlandau/crew/blob/1.0.0/R/crew_controller_local.R
[launcher-source]: https://github.com/wlandau/crew/blob/1.0.0/R/crew_launcher_local.R
[performance]: https://books.ropensci.org/targets/performance.html
[dynamic]: https://books.ropensci.org/targets/dynamic.html
[processx]: https://processx.r-lib.org/reference/process.html
