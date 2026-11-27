# PaintMixPrecompute

Offline generation of the `.pmx` table pair the runtime `PaintMix` package
evaluates. Nothing in this project is loadable from the runtime path: it owns
the optimizer, the automatic differentiation, the spectral model, and the
input data.

    precompute/
      Project.toml        runtime-independent deps: Optim, ForwardDiff, DuckDB, DataFrames, PaintMix
      config/default.toml settled conventions and generation parameters
      inputs/             selected sources, checksums, and the CIE observer
      src/
        PaintMixPrecompute.jl  module, config loading and validation
        inputdb.jl             spreadsheet -> DuckDB -> validated inputs
        spectra.jl             wavelength grid, D65/observer quadrature
        kubelka_munk.jl        equations (1)-(7) in generic arithmetic
        experiments.jl         resolved settings and canonical experiment identities
        surrogates.jl          equations (15)-(17), parameterizations and K/S fitting
        solvers.jl             Optimization.jl point adapters and traversal strategies
        unmix.jl               equation (9), the simplex-constrained inverse
        tables.jl              float table generation, quantization, padding
        validate.jl            quality, padding, continuity, round-trip reports
        export.jl              `.pmx` writing plus the provenance sidecar
        pipeline.jl            reusable stages, verified cache loading, orchestration
      scripts/
        generate.jl         CLI for the whole pipeline
        import_inputs.jl    spreadsheet -> DuckDB
      test/runtests.jl      precompute-side tests
      output/               checkpoints and candidate artifacts (ignored)

The dependency direction is one-way: `precompute -> PaintMix`. The runtime
package never depends on anything here.

The forward table and the coarse concentration field interpolate through
`PaintMix.LUTKernels`, the runtime package's documented shared-kernel
interface. These kernels are the shared numerical interface; payload export also calls
the runtime model validation and file-format APIs.

## Naming

Public names are ASCII (`theta_parameters`, `initial_theta`,
`alpha_schedule`). Scratch variables in the math and the solvers use Greek
letters (`θ`, `λ`, `δ`, `α`, `φ`) to match the paper's notation. Nothing
stored, exported, or hashed uses a Greek name.

## Running the pipeline

```sh
# Import the measurement workbooks and the CIE tables into DuckDB.
julia --project=precompute precompute/scripts/import_inputs.jl --force

# Fit the surrogates, generate both tables, validate, and promote.
julia -t auto --project=precompute precompute/scripts/generate.jl --profile dev --no-promote
julia -t auto --project=precompute precompute/scripts/generate.jl --profile release
```

`--stage inputs|fit|tables|validate|promote|all` stops after a stage;
`--solver reference|bulk|coarse` selects the inverse generation strategy
(default `coarse`); `--evaluator spectral|float_lut` independently selects its
forward evaluator. Artifacts live under `precompute/output/<profile>/`, with
physical fits keyed separately from table checkpoints. Reuse verifies resolved
settings and content digests. See [EXPERIMENTS.md](EXPERIMENTS.md) for the
Optimization.jl interface, controlled comparisons, and SLURM execution.

The import step reads the workbook and the observer CSV through DuckDB.jl
and DataFrames.jl, so no external client is needed. DuckDB's `excel`
extension is installed on first use, which needs network access once.

## What each step does

**Step 3 — spectral reference (`inputdb.jl`, `spectra.jl`,
`kubelka_munk.jl`).** `import_inputs.jl` reads the selected ranges of
`data/Final_artist_database.xlsx` through DuckDB's `excel` reader and
writes a normalized database recording, for every sample, its file, sheet,
and cell. `PigmentSpectra(cfg, db)` turns that into checked matrices;
`Quadrature(cfg, db)` builds the trapezoidal weights, the D65 and CIE 1931
2-degree samples, and the pinned XYZ-to-linear-sRGB matrix. `mix_rgb`
evaluates equations (1)-(7) in generic arithmetic, and `mix_rgb_jacobian!` is
the analytic `dRGB/dc` used by the inverse solver; a test checks the analytic
Jacobian against finite differences.

**Step 4 — surrogates and inverse (`surrogates.jl`, `unmix.jl`, `solvers.jl`).**
`fit_surrogates` defaults to equation (15) with `Optim.LBFGS` through
Optimization.jl over
`K = epsilon + softplus(thetaK)`, warm-started as `alpha` halves from
`alpha_initial` until the gamut fits inside the RGB cube (`Epush <=
push_tolerance`). `unmix_reference` is the accuracy reference: it enumerates
all 15 nonempty pigment subsets, largest first, and stops early when the
interior solve is exact. `unmix_bulk!` is the cheaper continuation path used
by `inverse_slab!`.

Both solvers share one Levenberg-Marquardt core, `lm_active!`, over the
subset's softmax coordinates. The damped normal equations are at most 3 x 3,
so they are solved with `StaticArrays` `SMatrix \ SVector` rather than by
hand-rolled Gaussian elimination: the static LU is unrolled, allocates
nothing, and measured about ten times faster per solve. `dRGB/dc` stays
analytic: a `ForwardDiff` Jacobian costs about seven times more per call,
because it pushes four partials through the whole wavelength loop, and the
inverse table needs the Jacobian at every LM iteration. A test cross-checks
`unmix_reference` against `LeastSquaresOptim.LevenbergMarquardt` on the same
objective, driven by `ForwardDiff`, so the hand-written damping schedule,
stationarity test, and Jacobian convention have an independent oracle.

**Step 5 — tables and export (`tables.jl`, `export.jl`).**
`generate_forward` evaluates `mix_rgb` at every concentration-grid vertex,
padding off-simplex vertices by Euclidean projection onto the simplex.
`generate_inverse` fills the RGB table slab by slab, with deterministic
traversal, per-slab checkpoints, and atomic completion markers.
`quantize_simplex` rounds all four concentrations jointly to integers summing
to 255 and stores the first three. `export_model` writes the payload through
`PaintMix.write_model` and a TOML sidecar recording the resolved configuration hash,
input checksums, surrogate parameters and convergence history, environment
versions, seeds, precision, and quantization and padding rules.

**Step 6 — validation and promotion (`validate.jl`).** `quality_report`
compares the byte tables with the spectral reference over a fixed corpus;
`padding_report` measures the padding rule's error on and just inside every
simplex face; `roundtrip_report`, `quantization_report`,
`continuity_report`, and `behavior_report` cover the remaining acceptance
rows. `acceptance_gates` turns the reports into the pass/fail table shown by
`generate.jl`, and `promote_validated` copies into `data/default/` only when
every gate passes.

## Measured results

On this machine (12 physical cores, `julia -t 14`, `surface_divisions = 20`):

| Stage | Cost |
| --- | --- |
| Surrogate fit, 24 alpha steps | ~13 min single-threaded |
| Forward table, 256³ | 17 s |
| Inverse table, coarse n = 64 plus fine polish | a few minutes; the fine pass is ~45 µs per vertex single-threaded |
| Validation | < 1 s |

Legacy slab entry points and the pipeline use one traversal, seed-policy
dispatch, and optional diagnostic accumulator. Coefficient layout conversion
also lives at one boundary; diagnostic push, pull, and boundary extrema share
one spectral pass.

The inverse solver is allocation-free in the hot path: `unmix_bulk!` and
`_solve_adaptive!` allocate nothing per solve, and a test asserts it. Getting
there took a few rounds. An earlier round removed the `ntuple` closures in the
softmax kernel. The active set is now a fixed-width `(count, NTuple{4,Int})`
pair: it travels through the adaptive face reduction and the reference subset
loop without an arity chain, and `_polish_face!` is the single place that
unpacks it into the static-length tuple the softmax kernel needs. `_project_onto`
captures an immutable reciprocal instead of the sum its loop mutates. The old
`Tuple{Int}` / `NTuple{2,Int}` / ... chain cost about 96 bytes per grid vertex on
the bulk path, which is gigabytes of garbage at 256³. `unmix_reference` still
allocates a little per solve; the test bounds it as a regression guard rather
than claiming zero.

A later change replaced the fixed-arity tuple arithmetic with `SVector` and
`SMatrix` throughout both packages, which retired the four-way `_softmax_k`
and three-way `_damped` arity families and the `_sub3`/`_mul3`/`_theta_plus`
helpers. The allocation guarantee is unchanged.

The release validation record is written to
`precompute/output/release/<model-id>.toml`. For the promoted model it shows
that the provisional LUT-versus-spectral quality targets are not met, while
the hard invariants (simplex quantization, exact round trip,
on-simplex padding, qualitative behavior) pass; the cause is the 8-bit
concentration storage and the steep near-vertex Kubelka-Munk curve. See the
root `README.md`. `generate.jl` refuses to promote a failing payload unless
`--force` is passed, and then records the overridden gate names in the
sidecar.

## Recorded deviations from the plan

  * **Surface quadrature weights.** Equations (16) and (17) integrate over
    the RGB gamut boundary with area element `ds`. The fit uses equally
    spaced concentration samples on the simplex faces with fixed weights,
    which is what the paper describes doing. `rgb_surface_weights` computes
    the literal area element and `quadrature_report` records the ratio
    between the two weightings, so the deviation is measured rather than
    assumed. Using the area element as a weight would make the integral a
    moving-boundary one and would have to be differentiated through.
  * **Solver and coordinates.** The default retains L-BFGS on softplus
    coordinates. Direct positive coefficients and true L-BFGS-B are separate
    comparison candidates; any compatible Optimization.jl algorithm can be
    supplied to the fitting and inverse interfaces.
  * **Inverse objective.** The spectral model and the float forward table
    are both available as evaluators (`SpectralModel` and
    `ForwardFloatLUT`); the default profile solves against the spectral
    model using reference coarse vertices and active-face fine polishing.
  * **Observer tables.** The source workbook carries no observer data, so the
    CIE 1931 2-degree functions and D65 are an explicit input file with its
    own checksum. See `inputs/README.md`.
  * **PY74, not PY73.** The selected source provides PY74; see
    `inputs/README.md`.

## Verification

Tests are split into configuration, input, spectral, surrogate, inverse,
table, validation, and experiment suites. Run `julia --project=precompute/test
precompute/test/runtests.jl`; synthetic tests need no measurement files.
`baseline.jl` records tiny forward/reference/bulk table digests and can run
unchanged against two source snapshots. The runtime consumer check continues
to enforce StaticArrays as the only non-stdlib runtime dependency.


## Explorer traces

Tracing is opt-in. Use the same profile, solver choices, configuration and
output directory for the fit and tables stages. The tables stage can reuse
verified slab checkpoints when adding inverse diagnostics to an existing run.

```sh
julia --project=precompute precompute/scripts/generate.jl \
  --profile dev --stage fit --trace-fit --no-promote
julia -t auto --project=precompute precompute/scripts/generate.jl \
  --profile dev --stage tables --trace-inverse --no-promote
julia --project=explorer explorer/run.jl \
  --payload precompute/output/dev/default.pmx
```

Use both flags in one `--stage all` run to capture both traces and validate
the payload. Requesting a fit trace reruns an untraced cached fit to capture
its intermediate coefficients. A traced fit can be reused with either
flag setting. The inverse flag uses the unquantized tables, including any
resumed slabs.

Both artifacts are schema-1 TOML documents. Cache entries and exported
sidecars carry versioned, relative filenames and SHA-256 digests. Exported
traces identify the fit, input checksums, resolved configuration, physical
surrogate and exact payload. Keep the sidecar and its trace files together
when moving a payload. Legacy sidecars can omit the `traces` section.

- **Fit.** Each snapshot stores the four rows of physical K and S in the
  recorded pigment and wavelength order, with alpha, loss terms, cumulative
  solve time and termination. Step zero records the optimizer's physical
  initialization. Integration data support replay without the input database.
  The artifact stores the weighted fitting boundary once and a separate
  deterministic validation sample with fitting nodes excluded. Sampled
  maximum cube violation and Oklab drift from the measured model provide
  evidence at those concentrations. Capture copies the coefficients once per
  solve, then evaluates the validation sample after continuation. It supports
  at most 64 solves plus initialization, 48 boundary divisions, and 256
  wavelengths. It records continuation boundaries exclusively.
- **Inverse.** Up to 12 probes cover near-black, a central mixture, near-pure
  mixtures, the cube corner with the largest reference SSE among the eight
  corners, and the largest runtime reconstruction residuals on a fixed
  64-color corpus. The diagnostics use the generation evaluator and two
  deterministic seeds with the normal active-face LM kernel. They preserve
  seed and face start and final states after the budget of 64 accepted
  iterates per probe is reached. Each face is capped at 64 LM iterations.
  Reference enumeration uses up to four restarts per face and provides a
  numerical SSE comparison. The continuous fit, demonstration rounding,
  actual float and byte vertex data, and interpolation and residual terms
  have distinct fields. Capture runs up to 12 diagnostic and 20 reference
  point solves alongside table generation.

Each trace file is capped at 16 MiB. Diagnostics run outside the table
traversal and are stored separately from the `.pmx` payload.
