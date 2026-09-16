# Runtime benchmarks

```sh
julia --project=benchmarks benchmarks/run.jl --n 64 --samples 20000
julia --project=benchmarks benchmarks/run.jl --payload ../build/data/payload-full.pmx
```

`--n` builds a synthetic table of that size; `--payload` measures a real
payload including its load and validation time. Every measurement asserts the
allocation count, because "the mixing path does not allocate" is a contract,
not an aspiration.

## Measured

x86-64 Linux, Julia 1.13.0, single thread, warm cache. Real release geometry
(n = 256, 96 MiB of tables) with a synthetic payload:

| Operation | Per call | Per color | Allocations |
| --- | --- | --- | --- |
| `encode` | 62 ns | 62 ns | 0 |
| `decode` | 28 ns | 28 ns | 0 |
| `mix` (cold latents) | 148 ns | 148 ns | 0 |
| `mix!` (preallocated) | 147 ns | 147 ns | 0 |
| `bulk_mix!`, 50 000 random pairs | 58 ms | 1160 ns | 0 |
| `weighted_mix!`, 50 000 colors | 28 ms | 557 ns | 0 |
| payload load + full validation | 0.38 s | — | — |

The scalar figures use one fixed color, so its eight table corners stay in L1.
The batch figures use uniformly random colors across the whole table, so
almost every sample misses cache: that is the ~8x difference per color between
`mix` and `bulk_mix!`, and it is the number that matters for painting, where
neighboring pixels share colors far more than random ones do.

`encode` costs two lookups (inverse, then forward for the residual) and
`decode` one, which is why `mix` is about 2.4x `encode`.

The release payload must be regenerated before the load-time figure means
anything for a real model: this run used a synthetic table with the same
geometry, and validation time depends only on size.
