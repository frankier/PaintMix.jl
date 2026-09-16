# Embedding a full-size payload in a trimmed library

Plan step 2 requires proving, before any expensive spectral work, that the
release-sized table pair survives `juliac --trim` and relocation. This is the
record of that experiment: what was built, what was measured, and what is
still not proven.

## Method

```sh
# 100 663 424 bytes: n = 256, 3 * 256^3 bytes per table, plus the 128-byte header
julia --project=build build/scripts/make_dummy_payload.jl full build/data/payload-full.pmx

julia --project=build build/compile.jl \
    --payload build/data/payload-full.pmx --out build/out-full

build/smoke/run.sh build/out-full build/data/payload-full.pmx
```

The dummy payload is deterministic and synthetic. It is a payload-shaped
object, not a model: its purpose is to exercise the size, the serialization,
and the loader.

The payload is embedded by `build/library.jl` reading it while the image is
built and constructing the model from it. There is no sidecar file, no
download, and no table path at run time.

## Results

Measured on x86-64 Linux, Julia 1.13.0, `juliac` 0.3.8, `trim = :safe`.

| Quantity | Value |
| --- | --- |
| Payload | 100 663 424 bytes (96 MiB of tables) |
| Payload generation | ~1.3 s |
| `juliac` image | 98.977 MiB |
| Linked `paintmix.so` | 98.006 MiB |
| Whole build (`compile.jl`, warm depot) | ~21 s |
| Model-info `grid_n` / id | match the payload file |
| Client checks (C / Python) | all pass, no failures |

The client checks are the load-bearing evidence. Trimming rewrites the image
and the loader relocates it, so the table bytes in memory are not the bytes
on disk in any naive sense. The clients do not compare raw bytes; they
recompute mixing, encoding, and decoding results and compare them against
values the Julia package computed from the payload file. Tables that the
loader mangled would change those results.

## One toolchain-cost finding, and how it was handled

`JuliaLibWrapping.build_library` always includes the entry file once, in an
unoptimized subprocess (`--compile=min -O0`), to collect `@api` metadata. Our
library has no `@api` functions, so that subprocess produces nothing — but it
still executed the top-level payload load. With full simplex validation that
cost **12 minutes 22 seconds** per build at n = 256, because the byte loop
runs interpreted.

The fix keeps validation where it belongs and keeps the embedded path cheap:

  * `compile.jl` validates the payload in an optimized process
    (`PaintMix.read_model`, which checks the simplex invariant) and prints
    the identity it validated.
  * `library.jl` then calls
    `PaintMix.model_from_bytes(bytes; validate = false)`, which still
    performs every O(1) structural check — magic, version, dimensions,
    offsets in range, storage conventions — but skips the O(n^3) loop.
  * Table copies use `copyto!`, a `memmove` even when interpreted.

Build time went from 12m22s to 21s with no loss of verification: the check
moved from an interpreted subprocess into the driver.

## What this does and does not prove

Proven: a payload of the exact release size round-trips through `juliac
--trim`, links into a ~98 MiB shared library, loads without any external
file, and yields mixing, encoding, and decoding results that match the Julia
package in two foreign-language clients.

Not proven here: that the toolchain behaves the same on other platforms, with
a bundled runtime (`--bundle`), or on a machine without the depot that built
it. Step 7 covers the release matrix. The dummy payload is also not a
substitute for the production model: nothing here says anything about table
quality, only about transport.
