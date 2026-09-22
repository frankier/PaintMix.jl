# PaintMix.jl

Practical pigment mixing in RGB, after Sochorová and Jamriška, *Practical
Pigment Mixing for Digital Painting* (2021), for Julia 1.13+.

The runtime package evaluates a pair of precomputed lookup tables:

```julia
using PaintMix

blue = (0.02, 0.09, 0.42)      # linear-light sRGB
yellow = (0.71, 0.62, 0.02)

mix(blue, yellow, 0.5)         # paint-like mixture, not an RGB average
mix(blue, yellow, 0.0) === blue
```

Colors are **linear-light** sRGB with the sRGB primaries and the D65 white
point. For image data, convert explicitly:

```julia
linear_from_srgb8((0x40, 0x80, 0xc0))   # bytes -> linear Float32
srgb8_from_linear((0.1, 0.4, 0.9))      # linear -> bytes, clips and rounds
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
file-format dependency in the runtime path.

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
| `build/` | the C ABI, `juliac` compilation, and the generated C/Python bindings |
| `benchmarks/` | runtime benchmarks |
| `test/` | fast runtime tests on synthetic tables |
| `PLAN.md` | the implementation plan this repository follows |

The root project is a Julia workspace; `precompute`, `test`, `build`, and
`benchmarks` are members with their own dependencies. The dependency direction
is one-way: `precompute -> PaintMix` and `build -> PaintMix`.

## Status

Plan steps 1-6 are implemented:

  * the workspace, with `precompute`, `build`, `test`, and `benchmarks` as
    members and a one-way dependency direction into `PaintMix`;
  * the settled conventions and the selected input data
    (`precompute/inputs/README.md`, `precompute/config/default.toml`);
  * the runtime: trilinear lookup on byte tables, latents with residuals,
    exact binary endpoints, bulk and weighted mixing, the `.pmx` format with
    simplex validation, and the sRGB adapters;
  * the compiled boundary: a 12-entrypoint C ABI, generated `paintmix.h` and
    `paintmix_py` bindings, and C and Python clients checked against Julia on
    both a tiny and a full-size (96 MiB) payload;
  * the precompute pipeline: spreadsheet import into DuckDB, the equations
    (1)-(7) spectral reference with an analytic Jacobian, the surrogate fit of
    equations (15)-(17), the active-face inverse solver of equation (9),
    resumable coarse-to-fine table generation, joint simplex quantization,
    `.pmx` export with a full provenance sidecar, and the validation reports;
  * a promoted 256³ release payload in `data/default/default.pmx`
    (model id `f78bfa4fe793ca748d7647605cb37147`).

The compiled library under `build/` still embeds the step-2 *synthetic*
full-size payload, which is what validated the trimming and relocation
behavior. Recompiling it with the release payload is plan step 7 and is not
part of the work recorded here.

### Release model

Measured on this machine (12 physical cores, `julia -t 14`):

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

The mixing tail is a format limitation, not a solver failure. The inverse
table stores each concentration in one byte, so near a pure pigment the
concentration is quantized to 1/255 while the Kubelka-Munk RGB curve changes
by ~0.3 over that step. The residual makes the round trip exact, but the
*linear* residual interpolation between two colors does not follow the
nonlinear curve, which is what the p99 and max channel errors measure. The
coarse-to-fine inverse itself matches the reference solver exactly on a
held-out sample; the forward table's own error on random simplex points is
mean 0.0014, p99 0.0034. Removing the tail needs finer concentration storage
(a new `.pmx` variant) or a simplex-aware interpolation, both of which change
the runtime format.

The mean Oklab distance of 0.0084 is at the just-noticeable-difference scale
that the paper reports for its surrogate fit, so the typical mixture is about
as close to the spectral reference as the published model claims to be; the
plan's p99 target of 0.01 is stricter than the paper's own distribution,
which has a tail past 1 JND.

Evidence collected so far:

| Claim | Where |
| --- | --- |
| Release payload promoted with the quality override recorded | `precompute/output/release/` (ignored), `data/default/README.md` |
| Full-size payload survives `juliac --trim`, relocates, and matches the file | `build/EMBEDDING.md` |
| C and Python clients match Julia across randomized inputs and error cases | `build/smoke/`, all checks pass |
| Zero steady-state allocations in every kernel | `test/test_allocations.jl`, `benchmarks/README.md` |
| Structs have concretely-typed fields, and Aqua checks pass | `test/test_quality.jl` |
| A fresh consumer environment installs only `PaintMix` | `test/consumer_env.jl` |
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
bash build/smoke/run.sh                               # C + Python vs Julia

julia --project=benchmarks benchmarks/run.jl --n 256

# Precompute pipeline: import the workbooks, then generate a payload.
julia --project=precompute precompute/scripts/generate.jl --check-inputs
julia --project=precompute precompute/scripts/import_inputs.jl --force
julia -t auto --project=precompute precompute/scripts/generate.jl --profile dev --no-promote
julia -t auto --project=precompute precompute/scripts/generate.jl --profile release --force
```
