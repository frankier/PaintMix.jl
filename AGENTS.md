# AGENTS.md

## What this is

`PaintMix.jl` implements Sochorová and Jamriška, *Practical Pigment Mixing for
Digital Painting* (2021). The runtime package evaluates two precomputed lookup
tables: an inverse table (linear sRGB -> four pigment concentrations) and a
forward table (concentrations -> linear sRGB). Encoding adds a signed RGB
residual. The module is a pure evaluator with no non-stdlib dependencies, no
optimizer, and no spectral model.

The root project is a Julia workspace. `precompute`, `build`, `test`,
and `benchmarks` are members with their own dependencies. The dependency
direction is one-way: `precompute -> PaintMix` and `build -> PaintMix`.

## File layout

| Path | Contents |
| --- | --- |
| `src/PaintMix.jl` | module, exports, includes |
| `src/types.jl` | models, headers, latents, payload errors |
| `src/lookup.jl` | indexing, byte decoding, trilinear interpolation |
| `src/srgb.jl` | linear-light sRGB <-> 8-bit adapters |
| `src/mixing.jl` | `encode`, `decode`, `mix`, bulk and weighted mixing |
| `src/tables.jl` | `.pmx` reading, writing, validation, default model |
| `data/default/` | the promoted release payload (untracked build artifact) |
| `precompute/` | `PaintMixPrecompute`: spectra, surrogate fit, inverse solver, export |
| `build/` | the C ABI, `juliac` compilation, generated C/Python bindings |
| `benchmarks/` | runtime benchmarks, which also assert allocations |
| `test/` | fast runtime tests on in-memory synthetic tables |

## Docs

* `README.md` — scope, usage, status, measured results, all commands.
* `PLAN.md` — the implementation plan and the conventions this repository follows.
* `precompute/README.md` — pipeline stages, caching, and the input database.
* `build/README.md` — the C ABI contract and compilation.
* `build/EMBEDDING.md` — how the trimmed library embeds the payload.
* `data/default/README.md` — the `.pmx` byte layout and what it does not say.
* `benchmarks/README.md`, `precompute/inputs/README.md` — data and timings.

Read the relevant file before changing that area. When the format or a
convention changes, update the README that documents it in the same change.

## Style

* Be concise. Prefer short functions, short comments, and short names.
  Comment the convention or the reason, never the mechanics.
* Use the standard library when possible. The runtime package must keep an
  empty `deps` section in `Project.toml`; `test/runtests.jl` enforces this.
  Add a dependency only in the member project that needs it.
* Follow the SciML style guide, <https://docs.sciml.ai/SciMLStyle/dev/>:
  `snake_case` functions, `CamelCase` types, `!` for mutating methods, no
  type piracy, and docstrings on exported names.
* Make reasonable use of multiple dispatch. Separate the design into
  orthogonal parts and dispatch on one argument at a time, instead of writing
  branchy functions or `if x isa T` checks. Keep kernels type-generic and free
  of steady-state allocation.

## Commands

```sh
julia --project=. -e 'using Pkg; Pkg.test()'                 # runtime tests
julia --project=precompute/test -e 'include("runtests.jl")'  # precompute tests
julia --project=benchmarks benchmarks/run.jl --n 64          # benchmarks
julia --project=build build/compile.jl --profile tiny        # shared library
bash build/smoke/run.sh                                      # C + Python vs Julia
```

Test with synthetic tables. The measurement workbooks and the release payload
are untracked; tests that need them must skip when they are absent.

## Formatting

Run Runic.jl like so:

```sh
runic --check --diff .  # check
runic --inplace .      # apply
```

`prek run --all-files` runs the same check through `prek.toml`. CI runs
`fredrikekre/runic-action`. Both install Runic in their own environment, so
neither touches the runtime project's empty `[deps]` section.
