# PaintMix.jl

[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

Practical pigment mixing in RGB, after Sochorová and Jamriška, *Practical
Pigment Mixing for Digital Painting* (2021), for Julia 1.13+.

The runtime package evaluates a pair of precomputed lookup tables:

```julia
using PaintMix
using StaticArrays: SVector

blue = SVector(0.02, 0.09, 0.42)      # linear-light sRGB
yellow = SVector(0.71, 0.62, 0.02)

mix(blue, yellow, 0.5)         # paint-like mixture, not an RGB average
mix(blue, yellow, 0.0) === blue
```

Colors are **linear-light** sRGB with the sRGB primaries and the D65 white
point. For image data, convert explicitly:

```julia
linear_from_srgb8(SVector(0x40, 0x80, 0xc0))   # bytes -> linear Float32
srgb8_from_linear(SVector(0.1, 0.4, 0.9))      # linear -> bytes, clips and rounds
```

## What the package does and does not do

It encodes a color as four pigment concentrations plus a signed RGB residual,
mixes latents linearly, and decodes:

```
encode(x)          = Latent(U(x), x - M(U(x)))
decode(z)          = M(z.c) + z.r
mix(a, b, t)       = decode((1 - t) * encode(a) + t * encode(b))
weighted_mix(xs,ws)= decode(sum(ws[i] * encode(xs[i])) / sum(ws))
```

Encoding costs two trilinear lookups and decoding one, so a binary mix costs
five unless latents are cached. There is no optimizer, spectral model, or
file-format dependency in the runtime path; its only dependency is
StaticArrays, for the fixed-size channel vectors and the `SMatrix` solves.

Endpoints are exact (`mix(a, b, 0) === a`), fractions outside `[0, 1]` clamp,
and `decode` never clips: residuals can push a channel outside `[0, 1]`, which
is how the model reproduces colors that no mixture of the four pigments can
produce. Clipping happens only in `srgb8_from_linear`.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/` | the runtime package: tables, interpolation, latent algebra |
| `data/default/` | the promoted release payload and the format contract |
| `precompute/` | `PaintMixPrecompute`: spectra, surrogate fit, inverse solver, export |
| `build/` | the C ABI, `juliac` compilation, and the generated C/Python/R bindings |
| `benchmarks/` | runtime benchmarks |
| `explorer/` | standalone interactive viewer and experiment comparisons |
| `test/` | fast runtime tests on synthetic tables |
| `PLAN.md` | the implementation plan this repository follows |

The root project is a Julia workspace; `precompute`, `precompute/test`, `test`, `build`, and
`benchmarks` are members with their own dependencies. The dependency direction
is one-way: `precompute -> PaintMix` and `build -> PaintMix`.

## Status

Plan steps 1-6 and local release verification are implemented:

  * the workspace, with `precompute`, `build`, `test`, and `benchmarks` as
    members and a one-way dependency direction into `PaintMix`;
  * the settled conventions and the selected input data
    (`precompute/inputs/README.md`, `precompute/config/default.toml`);
  * the runtime: trilinear lookup on byte tables, latents with residuals,
    exact binary endpoints, bulk and weighted mixing, the `.pmx` format with
    simplex validation, and the sRGB adapters;
  * the compiled boundary: a 12-entrypoint C ABI, generated `paintmix.h`,
    `paintmix_py`, and `paintmix` (R) bindings, and C, Python, and R clients
    checked against Julia on
    tiny, synthetic full-size, and promoted release payloads;
  * the precompute pipeline: spreadsheet import into DuckDB, the equations
    (1)-(7) spectral reference with an analytic Jacobian, the surrogate fit of
    equations (15)-(17), the active-face inverse solver of equation (9),
    resumable coarse-to-fine table generation, joint simplex quantization,
    `.pmx` export with a full provenance sidecar, and the validation reports;
  * a promoted 256³ release payload in `data/default/default.pmx`
    (model id `f78bfa4fe793ca748d7647605cb37147`).

The promoted release payload passes trimmed compilation and relocated C,
Python and R checks, including execution without Julia on the client PATH.
The [follow-up record](precompute/results/followups/README.md) retains the
measured toolchain and hashes. Bundled distribution and the wider platform
matrix remain release work.

### Release model

Historical release validation (12 physical cores, `julia -t 14`):

| Stage | Result |
| --- | --- |
| Surrogate fit (20 surface divisions) | 13 min single-threaded; final `Epush` 7.9e-8, max dense cube violation 0.0037, max Oklab deviation 0.043 |
| Forward table (256³) | 17 s |
| Inverse table (coarse n = 64 + fine polish) | ~2-3 min on top of a resumed run |
| Validation | < 1 s |
| Payload | 96 MiB, model id `f78bfa4fe793ca748d7647605cb37147` |

**The provisional quality targets are not met**, and the override is recorded
in `precompute/output/release/<model-id>.toml`:

| Report | Measured | Plan target |
| --- | --- | --- |
| LUT-vs-spectral linear channel error | mean 0.026, p99 0.18, max 0.39 | p99 ≤ 2/255, max ≤ 8/255 |
| LUT-vs-spectral encoded channel error | mean 0.023, p99 0.25, max 0.59 | p99 ≤ 2/255, max ≤ 8/255 |
| LUT-vs-spectral Oklab distance | mean 0.0084, p99 0.067, max 1.09 | p99 ≤ 0.01 |
| Padding on the simplex face | p99 0.0026, max 0.0081 | p99 ≤ 2/255 |
| Round trip | float64 0, float32 3e-8 | ≤ 1e-12 / ≤ 2e-6 |

Every hard invariant passes: the simplex quantization,
the exact round trip, the on-simplex padding, and the qualitative behavior
(blue + yellow is a green, magenta + yellow an orange, white tints stay
monotone).

One-byte concentration storage contributes strongly to the mixing tail:
near a pure pigment, a concentration step of 1/255 can move the Kubelka-Munk
RGB curve by about 0.3. The residual preserves an exact round trip, while
linear interpolation of residuals cannot follow that nonlinear curve.
A separate, controlled direct-coefficient experiment reduces mixing Oklab
p99 from 0.0356 to 0.0217 with 16-bit concentrations. The packed simplex
forward lattice saves storage with little change in this metric. All variants
still miss the mixing quality targets; forward and inverse field resolution
near pure pigments remain work. See the
[full-size comparison](precompute/results/followups/README.md).

The LUT-versus-spectral Oklab error and the paper's surrogate-fitting bias
measure different approximation layers and should not be equated. New
comparison runs report surrogate bias, inverse quality, interpolation,
quantization, and complete mixing separately; see
[`precompute/EXPERIMENTS.md`](precompute/EXPERIMENTS.md).

Evidence collected so far:

| Claim | Where |
| --- | --- |
| Release payload promoted with the quality override recorded | `precompute/output/release/` (ignored), `data/default/README.md` |
| Full-size payload survives `juliac --trim`, relocates, and matches the file | `build/EMBEDDING.md` |
| C, Python, and R clients match Julia across randomized inputs and error cases | `build/smoke/`, all checks pass |
| Zero steady-state allocations in every kernel | `test/test_allocations.jl`, `benchmarks/README.md` |
| Structs have concretely-typed fields, and Aqua checks pass | `test/test_quality.jl` |
| A fresh consumer environment installs `PaintMix` and its StaticArrays closure | `test/consumer_env.jl` |
| The compiled library links only against the Julia runtime | `ldd build/out/paintmix.so` |
| Spectral Jacobian matches finite differences | `precompute/test/runtests.jl` |
| Inverse solver recovers concentrations and respects the simplex | `precompute/test/runtests.jl` |

## Development

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # workspace
julia --project=. -e 'using Pkg; Pkg.test()'          # runtime tests (fast)
julia --project=precompute/test precompute/test/runtests.jl
julia --project=test test/consumer_env.jl             # dependency isolation

julia --project=build build/compile.jl --profile tiny # shared library + bindings
bash build/smoke/run.sh                               # C + Python + R vs Julia

julia --project=benchmarks benchmarks/run.jl --n 256

# Precompute pipeline: import the workbooks, then generate a payload.
julia --project=precompute precompute/scripts/generate.jl --check-inputs
julia --project=precompute precompute/scripts/import_inputs.jl --force
julia -t auto --project=precompute precompute/scripts/generate.jl --profile dev --no-promote
julia -t auto --project=precompute precompute/scripts/generate.jl --profile release --force
```

## Numerical experiments and cleanup

The offline pipeline separates coefficient coordinates, objectives,
derivatives, Optimization.jl algorithms, continuation, and table traversal.
Defaults and the promoted payload remain unchanged. Candidate comparisons
cover softplus/L-BFGS, softplus/L-BFGS-B, direct-coefficient/L-BFGS-B and an
opt-in gauge normalization. A common held-out physical stopping criterion
records target attainment separately from native optimizer termination.
New physical fits, checkpoints, and sidecars have versioned identities and
content digests; old release payloads remain readable. See
[the experiment guide](precompute/EXPERIMENTS.md) for APIs and commands, error
layers, and resource-sized SLURM jobs. Explorer accepts `--comparison PATH`
to show candidate metrics side by side on `/fit`.

Experimental precision and simplex interpolation candidates use a separate
[`.pmxe` format](data/EXPERIMENTAL.md). `read_model` loads them for the Julia
mixing API; the default model and compiled ABI retain production `.pmx`.
The [experiment guide](precompute/EXPERIMENTS.md) describes regeneration and
controlled comparisons.
