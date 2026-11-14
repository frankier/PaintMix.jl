# Experimental `.pmxe` payloads

This opt-in format is for evaluating concentration precision and forward
interpolation. It does not replace `.pmx`, the promoted default model, or the
compiled C ABI. The runtime still depends only on StaticArrays.

```julia
using PaintMix
candidate = read_model("candidate.pmxe")  # ExperimentalModel
mixed = mix(candidate, a, b, 0.5)
write_model("copy.pmxe", candidate)
```

The Julia scalar, bulk and weighted APIs accept both model types. Latents
belong to the particular model that encoded them. The compiled library
continues to embed production `.pmx` payloads and rejects `.pmxe` input.

## Choices

The inverse remains a dense RGB cube with trilinear interpolation. Its
three stored concentrations use either UInt8 divided by 255 or UInt16
divided by 65535. All four concentrations are quantized jointly, including
the implied fourth, to the corresponding integer sum.

Forward RGB remains three bytes per vertex, scaled by 255. Mode 0 uses the
existing dense concentration cube, Euclidean simplex projection for padding,
and trilinear interpolation. Mode 1 stores a packed concentration simplex
with tetrahedral interpolation. These are independent of inverse precision.

The packed coordinates are cumulative concentrations:
`(u,v,w) = (c1,c1+c2,c1+c2+c3)`, so `0 <= u <= v <= w <= 1`.
A lattice vertex `(i,j,k)` with `0 <= i <= j <= k < n` represents
`(i,j-i,k-j,n-1-k)/(n-1)`. Vertices are channel-fastest, with index
`k*(k+1)*(k+2)/6 + j*(j+1)/2 + i`. There are `n*(n+1)*(n+2)/6`
vertices. Sorting the three cell fractions selects a tetrahedron; ties
increment higher axes first. All four accessed vertices remain in the ordered
domain, including on faces and at pure pigments. No padded vertices are read.

At n=256, the packed forward table is 8,487,168 bytes; the dense forward
is 50,331,648 bytes. A 16-bit inverse is 100,663,296 bytes. The packed layout
saves forward storage but remains subject to linear interpolation error and
8-bit RGB quantization. Wider inverse storage cannot remove inverse-field
interpolation or solver errors.

## Experimental version 1

All integer fields and UInt16 samples are little-endian. The header is 128
bytes; both tables follow without alignment padding. Grid edges from 2 to
1024 are supported. Unknown versions, conventions, reserved fields, incorrect
lengths, corrupt checksums and invalid inverse simplex sums are rejected.
This experimental schema can evolve independently of production format v1.

| Offset | Size | Value |
| --- | --- | --- |
| 0 | 8 | `PAINTMXE` |
| 8 | 2 | experimental version 1 |
| 10 | 2 | header bytes 128 |
| 12 | 1 | inverse bits: 8 or 16 |
| 13 | 1 | forward mode: 0 dense cube, 1 packed simplex |
| 14 | 1 | channel count 3 |
| 15 | 1 | color space 0: linear sRGB, D65 |
| 16 | 4 | grid edge n |
| 20 | 4 | reserved, zero |
| 24 | 16 | content/convention model ID |
| 40 | 4 | inverse scale: 255 or 65535 |
| 44 | 4 | forward scale: 255 |
| 48 | 1 | inverse interpolation 0: trilinear |
| 49 | 1 | index order 0: channel-fastest |
| 50 | 6 | reserved, zero |
| 56 | 8 | inverse offset: 128 |
| 64 | 8 | inverse bytes: `3*n^3*bits/8` |
| 72 | 8 | forward offset: `128+inverse_bytes` |
| 80 | 8 | forward bytes: `3*n^3` or `n*(n+1)*(n+2)/2` |
| 88 | 4 | CRC32C of the entire file with this field zero |
| 92 | 36 | reserved, zero |

CRC32C detects accidental corruption. Experiment reports additionally carry
SHA-256 payload digests and source identities. Pigment names, physical K/S,
input digests and fitting settings remain sidecar metadata.
