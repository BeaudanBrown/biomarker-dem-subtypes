# TabPFN performance investigation and improvement plan

Superseded by [the source-backed targets refactor plan](TARGETS_REFACTOR_PLAN.md),
which also addresses the failed baseline read and the currently disconnected
crew controller. Retained below as the initial source-inspection proposal.

Status: proposed, not implemented. Source inspection only; no R executed and no
raw data inspected. Low hardware utilisation has been reported, but no timed
reproduction is available. Structural overhead below is confirmed by code;
its contribution to the reported run is not measured. Do not modify the active
run or interpret a dispatched target as proof that inference has started.

## Current execution shape

The graph contains three fit targets (one per diagnostic contrast), routed to a
single GPU worker. Each target loops over ten folds. For every fold the existing
analysis R worker writes temporary RDS inputs, waits for a newly launched GPU R
subprocess, and reads its output. The child verifies the checkpoint, initializes
reticulate/Python/PyTorch, constructs a new estimator, fits, predicts and exits.
There are thirty fresh subprocesses for a complete run. Normal CPU targets may
run concurrently, but the selected comparison branch has little CPU-only work
to fill those gaps.

## Candidate bottlenecks, in priority order

1. **Repeated runtime startup (likely important, not yet timed).**
   `tabpfn_run_fold()` launches Rscript inside the fold loop. Each child imports
   R dependencies/Python/PyTorch and initializes GPU access from scratch.
   Short inference bursts can be separated by substantial host/startup work.
2. **Deliberate serialization (confirmed concurrency ceiling).**
   `_targets.R` configures one GPU worker; `tabpfn_auroc()` runs folds serially.
   A second worker could overlap two contrasts, but would not eliminate startup
   overhead or parallelize folds within one contrast. Only three coarse jobs
   exist, so adding more than three workers cannot help this graph.
3. **Host thread/preprocessing limits (confirmed configuration; impact unknown).**
   PyTorch has three host threads, TabPFN preprocessing has one job, and BLAS
   limits are one. This bounds CPU parallelism but does not establish CPU
   saturation. More threads could increase contention rather than speed.
4. **Repeated checkpoint verification/loading (confirmed repeated work).**
   The checkpoint is SHA-256 hashed in every child and a new estimator is
   constructed every fold. Filesystem caching can reduce disk reads, but does
   not remove hashing or Python/model initialization. Quantify these separately.
5. **Small GPU work units (hypothesis).**
   Each fold predicts only its held-out subset. GPU work can be short relative
   to initialization and host preprocessing. Actual ADDF sizes/costs have not
   been inspected. High GPU occupancy is not itself the optimization objective.
6. **Worker imports and temporary-file serialization (secondary candidates).**
   GPU crew workers inherit the large global analysis package list, although
   their prediction subprocess needs only the isolated GPU R dependencies.
   Every fold also writes/reads compressed RDS payloads. Measure before changing
   transport or package handling; these may be minor on small tables.

The current diagnostic records contain total fold time and GPU allocation peaks,
not stage timings, live progress or whole-device utilization. Captured subprocess
stdout/stderr is not streamed by `processx::run()`. An idle-looking parent worker
can simply be waiting on its child. These gaps prevent a diagnosis of the live
stall from the dispatch message alone.

## Stage 1: observability before tuning

Add aggregate-only progress and timings for worker/task entry, input transfer,
checkpoint verification, Python import, GPU initialization, model fit, prediction,
diagnostics and result transfer. Synchronize GPU work at timed boundaries so
asynchronous work is charged to the correct stage. Include total end-to-end time.
Do not log features, labels, predictions, participant IDs, or sensitive paths.
Record process RSS and whole-device memory/utilization separately from PyTorch
allocator peaks. Make progress visible while the child runs, not just afterward.

User-run bounded probe: one fixed fold with unchanged inputs, four members,
settings and seed, first serially and then two independent jobs concurrently.
Keep outputs local; user supplies only aggregate timing/memory summaries and
prediction-agreement checks. Do not overwrite completed fit targets merely to
profile. Preserve failure records and timeouts.

## Stage 2: low-risk concurrency experiment

Temporarily test two GPU crew workers with three host threads each, retaining
all existing fold/model settings. Compare total wall time/throughput, host RAM,
whole-device VRAM and failures against serial execution. Verify prediction/row
alignment using the same seeds; report any drift rather than assuming equality.
The reference project's two-worker probe improved throughput by 1.77x on its
small fixture, but that is not an ADDF speedup or capacity guarantee.

Adopt two workers only if the actual bounded comparison improves throughput
with acceptable resource use. If a target is waiting before subprocess launch,
investigate scheduling instead; increasing GPU concurrency is not a fix for an
unexplained queue/launch stall.

## Stage 3: amortize process startup if timings support it

Preferred first structural change: one GPU subprocess per contrast, running its
ten folds sequentially. This reduces child startups from thirty to three and
lets Python/PyTorch/GPU initialization be reused. Verify the checkpoint once at
child start. Keep the isolated GPU R installation and pin checks.

Create a fresh fold-specific estimator using only its training rows and retain
fold-specific seeds. Never reuse a trained estimator/context or fitted
preprocessor across folds. Release estimator resources between folds; verify
memory does not accumulate. Retain prediction-only outputs, stage diagnostics,
per-fold timeout enforcement and fail-fast behavior. A single contrast-wide
process timeout is not equivalent to the current per-fold timeout.

Do not promise that this automatically reuses loaded checkpoint tensors:
constructor/loading behavior must be checked against the pinned TabPFN API.
Compare predictions with the current execution path on a fixed fixture before
adopting the change. Combine warm-process execution with two-worker scheduling
only after testing each independently.

## Deferred unless measured costs justify them

- Fold-level targets for independent caching/retry and better load balancing.
  These do not themselves remove fresh-process overhead. A persistent backend
  across targets would require an explicit lifecycle/protocol and memory policy.
- Increasing preprocessing jobs or host threads: tune one parameter at a time.
- Alternate fit modes, checkpoint/model reuse APIs, or GPU batching: require
  pinned-library verification, resource measurements and prediction checks.
- IPC replacing RDS or trimming worker package imports: only if stage timings
  show a meaningful share of total runtime.

## Invariants and acceptance

No SuperLearner refits, changes to existing preprocessing, extra model members,
new tuning, hosted prediction, or weight downloads. Keep all caches and data
local. No Nix dependency upgrades are needed for the initial instrumentation,
worker-count experiment, or contrast-level subprocess restructuring.

Optimize lower end-to-end time with unchanged predictions within an explicitly
checked tolerance, preserved folds/seeds, bounded memory and clear failures—not
GPU utilization percentage alone. Actual R/model probes must be run by the user
under AGENTS.md; assistant checks remain static and data-free.
