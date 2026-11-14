# Numerical experiments

The runtime ABI and `.pmx` byte layout remain unchanged. Numerical choices
belong to `PaintMixPrecompute`; optimization packages never enter the runtime
consumer dependency closure. Existing promoted payloads and legacy sidecars
remain readable. New candidates have content-derived model IDs and schema-2
sidecars. No cleanup experiment promotes a candidate.

## Independent choices

Surrogate fitting separates physical coefficients, their optimization
coordinates (`SoftplusParameters` or `DirectParameters`), push/pull objectives,
derivatives, Optimization.jl algorithms, continuation, and held-out diagnostics.
The default remains softplus plus L-BFGS. Direct coefficients use the same
physical initialization and a positive lower bound equal to `epsilon`.
L-BFGS-B runs through `OptimizationLBFGSB`, not an interior barrier wrapper.

```julia
using PaintMixPrecompute
using Optimization, OptimizationOptimJL, OptimizationLBFGSB
cfg = load_config()
db, spectra, quad = pipeline_inputs(cfg)
fit = fit_surrogates(cfg, spectra, quad; algorithm = BFGS())
# Any compatible Optimization.jl algorithm and AD backend can be supplied.
problem = SurrogateProblem(FitSettings(cfg), spectra, quad)
x = initial_parameters(problem.parameters, spectra)
prob = optimization_problem(problem, x, 1.0; adtype = AutoForwardDiff())
sol = solve(prob, LBFGS(); maxiters = 100)
```

`fit_surrogates` accepts `algorithm`, `adtype`, and `solve_kwargs`. Algorithms
must support the problem's bounds and requested derivatives. Its default
specialized ForwardDiff gradient uses 16-wide chunks; a manual Hessian is
available for second-order methods. `derivative = "finite_difference"` is an
independent check. The stage enforces its wall-clock budget even for backends
that ignore `maxtime`, and records termination and evaluation counts.

Inverse generation independently selects an evaluator, a point solver, and a
seed/traversal strategy. `ActiveFaceLM()` retains the allocation-free default;
`OptimizationInverse(algorithm; kwargs...)` substitutes Optimization.jl on the
same active faces. The analytic evaluator Jacobian is composed with the
softmax derivative, including for `ForwardFloatLUT` evaluators. Backend setup
can allocate on this experimental path; it is not used for production's 16.7
million default point solves. Enumeration explores local face optima and does
not certify global optimality.

```julia
st = TableSettings(cfg; profile = "dev", strategy = :bulk, threads = 4)
inverse = generate_tables_inverse(SpectralModel(spectra, quad), st;
    point_solver = OptimizationInverse(BFGS()))
```

`TableSettings` resolves the actual grid, coarse grid, fine iteration cap,
evaluator (`:spectral` or `:float_lut`), strategy (`:reference`, `:bulk`, or
`:coarse`), and worker limit. Worker count is honored and excluded from
numerical identity. Each slab records iterations, evaluations, residuals, and
unconverged points; resumed slabs are labeled instead of inventing statistics.

## Controlled comparison

```sh
julia --project=precompute precompute/scripts/compare.jl RESULTS fit softplus-lbfgs release
julia --project=precompute precompute/scripts/compare.jl RESULTS fit softplus-lbfgsb release
julia --project=precompute precompute/scripts/compare.jl RESULTS fit direct-lbfgsb release
# After each fit, generate full 256^3 candidate tables and validate.
julia -t 8 --project=precompute precompute/scripts/compare.jl RESULTS all softplus-lbfgs release
julia -t 8 --project=precompute precompute/scripts/compare.jl RESULTS all softplus-lbfgsb release
julia -t 8 --project=precompute precompute/scripts/compare.jl RESULTS all direct-lbfgsb release
julia --project=precompute precompute/scripts/compare.jl RESULTS summary
```

The three rows distinguish solver effects from coordinate effects. They share
physical initial coefficients, input spectra, positivity floor, objective,
quadrature, continuation, iteration/time limits, and validation inputs. The
paper formulation uses direct positive K/S with L-BFGS-B. PY74 instead of PY73,
linear-sRGB conventions, and concentration-space quadrature remain explicit
implementation choices; this experiment does not claim exact reproduction.
Native gradient stopping criteria also change meaning with coefficient
coordinates; use the common physical push/pull and held-out metrics when
judging quality, rather than comparing convergence flags alone.
Per-wavelength common K/S scale is unidentifiable; compare RGB predictions
and quality distributions rather than assuming coefficient distances measure
mixing differences. Gauge normalization is not silently introduced.

`comparison.toml` and `comparison.md` report timing, evaluation counts,
termination, independent gamut and perceptual diagnostics, table quality, and
acceptance gates. Wall-clock and iteration budgets are both recorded; timing
comparisons require the same CPU allocation. A threshold missed within the
budget remains a failed target, not an implied convergence success.

## Identities, caches, and validation

`fit_job_id` hashes resolved numerical settings, input digests, implementation,
and numerical package versions. `table_job_id` includes fitted physical K/S,
actual CLI overrides and grids, model conventions, and table implementation.
Canonical serialization orders dictionary keys. Paths and timestamps do not
identify numerical results. Custom algorithms should be recorded explicitly
by callers managing their own experimental artifacts.

Physical fit files have coefficient digests. Schema-2 checkpoints use
`job.toml` plus atomic, SHA-256-checked slab completion records. Legacy
checkpoints are not reused; legacy release payloads are unaffected.
`--stage validate` and `--stage promote` verify the candidate's artifact
manifest, payload digest, model ID, sidecar, and resolved table-job identity.
They never attach newly invented provenance to an unrelated existing payload.

Validation corpus version 2 uses a specified SplitMix64 stream and hashes the
actual RGB samples, pairs, and concentrations. It covers simplex faces,
edges, vertices, near-pure pigments, and interior points. The same corpus is
required across comparison rows. End-to-end reports distinguish unclipped
linear and Oklab errors from clipped display-sRGB errors. The decomposition
reports surrogate bias, forward interpolation, forward bytes, concentration
quantization, inverse objective excess, and inverse concentration bytes.
These distributions are not additive. Exact residual round trips alone do
not establish mixing quality or a correct inverse.

## SLURM execution

Set `PAINTMIX_WORK` to a separate disk-backed scratch workspace. Login nodes
are for file transfers and submission only. `precompute/slurm/compute.sh` requires
`PAINTMIX_WORK`, sets scratch `TMPDIR` and a private depot, and runs Julia
with the allocated CPU count. It defaults to one CPU, 16 GiB, and 12 hours;
it never requests an exclusive node. Override memory and CPU count per job.
Set `PAINTMIX_JULIA` to an appropriate module-provided Julia executable when
available; otherwise the helper uses the staged Julia 1.13.3 runtime.

Run each single-threaded fit with one CPU. Run `scripts/scaling.jl` with a
small CPU allocation before choosing table workers; it measures forward,
reference-inverse, and coarse-to-fine workloads. A full-node reservation
requires evidence that the workload scales to it. Keep fit and table jobs
separate so a single-threaded fit does not hold a table job's idle CPUs.
Freeze a source snapshot for the comparison: editing implementation files
correctly changes cache identities. Save compact reports in the project;
keep full tables, checkpoints, environments, and logs under scratch.

## Recorded cleanup run

Compact evidence is retained in [`results/cleanup/`](results/cleanup/).
The original and cleaned tiny forward/reference/bulk tables have identical
Float64 content digests. The complete default fine-slab kernel also measured
zero steady-state allocations on 4,096 points. After correcting a task-local
scratch-binding race, the fine-pass probe took 22.4, 11.8, 5.91, and 4.66
seconds at 1, 2, 4, and 8 workers. Four CPUs per table job were selected; no
exclusive node was reserved. All numerical checks and comparison jobs after
the initial local baseline test ran on SLURM compute nodes using a private
scratch workspace and Julia 1.13.3.
