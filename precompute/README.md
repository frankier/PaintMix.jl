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
        surrogates.jl          equations (15)-(17), the K/S fit
        unmix.jl               equation (9), the simplex-constrained inverse
        tables.jl              float table generation, quantization, padding
        validate.jl            quality, padding, continuity, round-trip reports
        export.jl              `.pmx` writing plus the provenance sidecar
      scripts/
        generate.jl         CLI for the whole pipeline
        import_inputs.jl    spreadsheet -> DuckDB
      test/runtests.jl      precompute-side tests
      output/               checkpoints and candidate artifacts (ignored)

The dependency direction is one-way: `precompute -> PaintMix`. The runtime
package never depends on anything here.

## Running the pipeline

```sh
# Import the measurement workbooks and the CIE tables into DuckDB.
julia --project=precompute precompute/scripts/import_inputs.jl --force

# Fit the surrogates, generate both tables, validate, and promote.
julia -t auto --project=precompute precompute/scripts/generate.jl --profile dev --no-promote
julia -t auto --project=precompute precompute/scripts/generate.jl --profile release
```

`--stage inputs|fit|tables|validate|promote|all` stops after a stage;
`--solver reference|bulk` selects the inverse solver. Expensive artifacts are
cached under `precompute/output/<profile>/` and reused only when the
configuration hash, input checksums, grid size, and surrogate settings all
match, so changing the configuration invalidates every cache.

The import step reads the workbook and the observer CSV through DuckDB.jl
and DataFrames.jl, so no external client is needed. DuckDB's `excel`
extension is installed on first use, which needs network access once.

## What each step does

**Step 3 — spectral reference (`inputdb.jl`, `spectra.jl`,
`kubelka_munk.jl`).** `import_inputs.jl` reads the selected ranges of
`data/Final_artist_database.xlsx` through DuckDB's `excel` reader and
writes a normalized database recording, for every sample, its file, sheet,
and cell. `load_spectra` turns that into checked matrices;
`build_quadrature` builds the trapezoidal weights, the D65 and CIE 1931
2-degree samples, and the pinned XYZ-to-linear-sRGB matrix. `mix_rgb`
evaluates equations (1)-(7) in generic arithmetic, and `mix_rgb_jacobian!` is
the analytic `dRGB/dc` used by the inverse solver; a test checks the analytic
Jacobian against finite differences.

**Step 4 — surrogates and inverse (`surrogates.jl`, `unmix.jl`).**
`fit_surrogates` solves equation (15) with `Optim.LBFGS` over
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
`PaintMix.write_model` and a TOML sidecar recording the configuration hash,
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

The inverse solver is allocation-free in the hot path: `unmix_bulk!` and
`_solve_adaptive!` allocate nothing per solve, and a test asserts it. Getting
there took two rounds. An earlier round removed the `ntuple` closures in the
softmax kernel. This round removed two remaining boxes: the active set now
travels as a concrete `(count, NTuple{4,Int})` pair instead of a
`Tuple{Int}` / `NTuple{2,Int}` / ... chain, and `_project_onto` captures an
immutable reciprocal instead of the sum its loop mutates. Together those
cost about 96 bytes per grid vertex on the bulk path, which is gigabytes of
garbage at 256³. `unmix_reference` still allocates a little inside its
subset recursion; the test bounds it as a regression guard rather than
claiming zero.

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
  * **Solver substitution.** The paper uses L-BFGS-B; the fit uses
    unconstrained L-BFGS on softplus-transformed parameters. This is
    intentional and recorded in the sidecar and in the fit history.
  * **Inverse objective.** The spectral model and the float forward table
    are both available as evaluators (`SpectralModel` and
    `ForwardFloatLUT`); the production profile solves against the spectral
    model with `reference_slab!`. `inverse_table_slab!` is the cheaper,
    table-based alternative, kept because it is the same code path with a
    different evaluator and is covered by tests.
  * **Observer tables.** The source workbook carries no observer data, so the
    CIE 1931 2-degree functions and D65 are an explicit input file with its
    own checksum. See `inputs/README.md`.
  * **PY74, not PY73.** The selected source provides PY74; see
    `inputs/README.md`.
