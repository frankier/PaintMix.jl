# PaintMixLib — the C ABI surface

This project is the entry project for `juliac`. It is **not** a Julia-facing
package: it exists to expose a small, fixed set of `Base.@ccallable` functions
that evaluate the same tables the Julia `PaintMix` package uses.

The compiled library embeds one validated `.pmx` payload at
`build/data/payload.pmx`, read while the image is built, so there is no file
parsing, download, or table path at run time.

## Contract

  * **Scalars.** `f32` entrypoints use `float` and `f64` entrypoints use
    `double` throughout. No mixed-precision calls exist; a caller that wants
    `Float32` semantics calls the `f32` variants.
  * **Color order.** RGB, in that order. Arrays are flat and channel-fastest:
    `r0, g0, b0, r1, g1, b1, …`. Latents are seven scalars per color:
    `c1, c2, c3, c4, r, g, b`, in that order, so an array of latents is
    `7 * count` elements long.
  * **Coordinates.** Linear-light sRGB, D65, nominally `[0, 1]`. Values
    outside that range are accepted and never clipped by the mixing
    entrypoints.
  * **Buffers.** Every output buffer is caller-owned and passed as a
    `CVector{:borrowed,T}`, which is `(dims::NTuple{1,Int64}, data::Ptr{T})`.
    The library never allocates and never retains a pointer.
  * **Aliasing.** Inputs may alias each other. An output buffer may not
    overlap an input buffer read by the same call.
  * **Lengths.** `dims[1]` is the element count. A buffer shorter than the
    entrypoint requires is `PM_ERR_LENGTH`; a null `data` pointer with a
    positive count is `PM_ERR_NULL`.
  * **Zero length.** `count == 0` is a success for bulk mixing and writes
    nothing; for weighted mixing it is `PM_ERR_TOTAL`, because zero weights
    sum to zero. Null pointers are allowed when the corresponding count is
    zero.
  * **Errors.** Every entrypoint returns `JLWStatus`
    (`{ code::Int32, message::UInt8[256] }`) or a struct with a `JLWStatus`
    field. `code == 0` is success. No Julia exception crosses the boundary:
    inputs are validated before any buffer is touched, so a failed call
    leaves outputs unchanged.
  * **Threads.** The tables are read-only after load. Concurrent calls from
    any number of threads are safe, and there is no initialization step.

## Building

```sh
julia --project=build/build-env -e 'using Pkg; Pkg.instantiate()'
julia --project=build -e 'using Pkg; Pkg.instantiate()'
julia --project=build build/compile.jl --profile tiny      # fast development build
julia --project=build build/compile.jl --payload data/default/default.pmx
bash build/smoke/run.sh                                    # C + Python vs Julia
```

`compile.jl` copies the chosen payload to `build/data/payload.pmx`, which is
the path `library.jl` reads while the image is built. Without `--bundle` the
library still needs Julia's `lib` directory on the loader path; `--bundle`
produces the self-contained `juliac --bundle` tree instead (hundreds of
megabytes).

## Where the payload is validated

The payload is checked once, in `compile.jl`, by `PaintMix.read_model`: the
inverse table's simplex invariant. `library.jl` then embeds it with
`validate = false`, because the image build and the ABI-metadata probe run
unoptimized and the O(n^3) loop costs minutes there; see
[EMBEDDING.md](EMBEDDING.md) for the measurement. Every O(1) structural check
still runs while the image is built.

## Smoke tests

`build/smoke/run.sh` generates reference values with the Julia package from
the same payload, then checks the compiled library from C (`client.c`, linked
against `paintmix.h`) and from Python (`client.py`, importing the generated
`paintmix_py` package). Both clients cover mixing, encoding, decoding,
weighted mixing, the error contract, and the model-info identity.
