# `data/default/` — the packaged release payload

This directory holds the validated table pair that ships with `PaintMix`.
`default.pmx` is the file `PaintMix.default_model()` loads; it is produced by
`PaintMixPrecompute` and promoted here only after the acceptance gates pass.

It is a build artifact rather than source: it is not committed, and a
checkout without it fails loudly and tells you the path it expected rather
than falling back to a substitute model. Rebuild it with

```sh
julia -t auto --project=precompute precompute/scripts/generate.jl --profile release
```

which fits the surrogates, generates both `256^3` tables, writes
`precompute/output/release/<model-id>.toml` with the full provenance and
validation record, and copies the payload here when every gate passes. The
model id in the header is derived from that provenance, so a payload can
always be traced back to the configuration and inputs that produced it.

## The shipped model

The promoted payload has model id `f78bfa4fe793ca748d7647605cb37147`. It was
built from the four pigments and the CIE/D65 inputs recorded in
`precompute/inputs/`, using the surrogate fit and inverse settings in
`precompute/config/default.toml`.

Its validation record reports that the provisional LUT-versus-spectral
quality targets in `PLAN.md` are *not* met: the mean encoded-sRGB channel
error is 0.023 and the p99 is 0.25, against a p99 target of `2/255`. The hard
invariants do pass — joint simplex quantization, the exact encoder/decoder
round trip, the on-simplex padding error (p99 0.0026),
and the qualitative mixing behavior. The tail comes from storing each
concentration in one byte: near a pure pigment the concentration step is
`1/255` while the Kubelka-Munk RGB curve moves by roughly 0.3, so the linear
residual interpolation between two colors cannot follow the curve. This is a
property of the 8-bit `.pmx` format, not of the solver; see the root
`README.md` for the measured numbers. `generate.jl` refuses to promote a
payload that fails a gate unless `--force` is passed, and records the
overridden gate names in the sidecar.

## Payload layout

A `.pmx` payload is a fixed 128-byte little-endian header followed by two raw
byte tables. `PaintMix.model_from_bytes` is the reference implementation;
`precompute/src/export.jl` writes payloads through `PaintMix.write_model` so
there is exactly one definition of the format.

| Offset | Size | Field | Meaning |
| --- | --- | --- | --- |
| 0 | 8 | magic | `PAINTMIX` |
| 8 | 2 | format_version | `1` |
| 10 | 2 | header_bytes | `128` |
| 12 | 1 | storage | `0` = `UInt8`. `1`/`2` are reserved for `f32`/`f64` |
| 13 | 1 | channels | `3` |
| 14 | 1 | table_count | `2` |
| 15 | 1 | color_space | `0` = linear-light sRGB, sRGB primaries, D65 |
| 16 | 4 | grid_n | edge length `n` of both tables |
| 20 | 4 | flags | bit 0: forward table padded by simplex projection; bit 1: inverse table quantized by largest remainder |
| 24 | 16 | model_id | identifies the pigment set, conventions, and generator settings |
| 40 | 1 | byte_scale | `0` = `b/255` |
| 41 | 1 | interpolation | `0` = trilinear |
| 42 | 1 | index_order | `0` = channel-fastest, `ch + 3 * (i + n * (j + n * k))` |
| 43 | 13 | reserved | zero |
| 56 | 8 | inverse_offset | byte offset of the inverse payload |
| 64 | 8 | inverse_bytes | `3 * n^3` |
| 72 | 8 | forward_offset | |
| 80 | 8 | forward_bytes | `3 * n^3` |
| 88 | 40 | reserved | zero |
| 128 | | inverse table | RGB `->` concentrations |
| | | forward table | concentrations `->` RGB |

Both tables are dense `n x n x n x 3` byte arrays. Values are `b / 255`, and
interpolation is trilinear with exact treatment of the last grid plane (a
coordinate of exactly `1` uses the final cell and weight one, never reading
past the buffer).

## What the bytes do not say

The header records every convention the runtime needs to evaluate the tables.
These are the parts a consumer has to be told out of band:

  * **Pigment order.** The four concentrations are, in order, the entries of
    `[[pigments]]` in `precompute/config/default.toml`: PB15:4, PR122, PY74,
    PW6 for the default model. All four sum to one; the inverse table stores
    the first three and the fourth is `1 - c1 - c2 - c3`. A payload with a
    different pigment set has a different `model_id` and is a different model,
    not a variant of this one.
  * **Storage order in the inverse table.** The three stored bytes at a vertex
    are `c1, c2, c3`; `c4` is implied.
  * **Cube versus simplex.** The inverse table's domain is the whole RGB cube,
    so every vertex is a legitimate query and stores a valid simplex point.
    The forward table's domain of interest is `c1 + c2 + c3 <= 1`, but
    trilinear lookup can read cube corners outside that region. Those vertices
    are padded: the four-component concentration vector at an off-simplex
    vertex is projected onto the simplex before the offline evaluation.
    **This rule is an implementation choice**: the paper does not specify it.
    `precompute/src/tables.jl` must quantify the resulting error near every
    simplex face, edge, and vertex before a release payload is promoted.
  * **Residuals.** Residuals are computed against the same quantized,
    interpolated forward table the runtime decodes with, and are stored as
    signed floating-point values. They are never clamped or quantized to
    bytes, so a latent can decode outside `[0, 1]`.
  * **Provenance.** `precompute/output/<model-id>.toml` records the
    configuration hash, input checksums, surrogate parameters, Julia and
    package versions, seeds, objective settings, precision,
    quantization/padding rules, and validation results for the payload. The
    header carries only `model_id`; the sidecar carries the story.
