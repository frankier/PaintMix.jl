# Cleanup verification and release comparison

All three candidates were fitted and generated at 256³ on SLURM
compute nodes, using Julia 1.13.3 and the same measured inputs, positivity
floor, continuation schedule, numerical budgets, and version-2 validation
corpus. Defaults and the promoted `data/default/default.pmx` remain unchanged.

## Results

| Coefficients / optimizer | Fit seconds | Mixing Oklab p99 | Encoded channel p99 | Maximum gamut violation |
| --- | ---: | ---: | ---: | ---: |
| Softplus / L-BFGS (current default) | 891.7 | 0.0604 | 0.1970 | 0.00194 |
| Softplus / L-BFGS-B | 137.8 | 0.0615 | 0.2098 | 0.00307 |
| Direct positive K/S / L-BFGS-B | 239.8 | 0.0356 | 0.0982 | 0.01637 |

Changing the optimizer greatly reduced fitting time in this run. Direct
coefficients improved the measured mixing errors but left more of the fitted
gamut outside the RGB cube. Native convergence flags are not a common
physical stopping criterion across coordinate systems. These are local fits
with finite iteration budgets, not global-optimum or exact-paper-reproduction
claims.

Every candidate passes simplex quantization, Float32/Float64 round trips,
on-simplex padding, and the three qualitative behavior gates. Every candidate
fails the three provisional LUT-versus-spectral quality targets (channel p99,
channel max, and Oklab p99). No candidate was promoted.

The [complete comparison](comparison.md) and [machine-readable report](comparison.toml)
include independent surrogate diagnostics, approximation-layer errors,
evaluation counts, timing, acceptance gates, and artifact identities. The
float-table generation timings exclude payload quantization/export and
validation; full job wall times are in the SLURM accounting record.

Physical coefficients, resolved settings, fit history, validation, payload
SHA-256, and per-slab solver statistics are retained in the complete sidecars:

- [Current default](candidates/softplus-lbfgs.toml)
- [Softplus / L-BFGS-B](candidates/softplus-lbfgsb.toml)
- [Direct coefficients / L-BFGS-B](candidates/direct-lbfgsb.toml)

The large working payloads and checkpoints remain in the execution
workspace's scratch storage. Scratch is temporary; these compact repository
reports preserve the physical results and provenance without committing
hundreds of megabytes of tables. Configuration paths in the retained reports
are repository-relative, and compute-node labels are anonymized. Numerical
results, job IDs, and content digests are preserved.

## Verification

| Check | Result |
| --- | --- |
| Runtime tests | 5,658 passed |
| Fresh runtime consumer | 4 checks passed; 17 resolved packages, no optimizer dependencies |
| Exact isolated-consumer CI command | Passed the runtime suite |
| Precompute tests, two Julia workers | 3,577 passed, including scratch ownership, serial/parallel equality, generic optimizer coverage, cache integrity, and stage reuse |
| Explorer with the real legacy sidecar | 533 passed |
| Explorer with the real direct-coefficient payload and comparison report | Passed; physical coefficients reconstruct exactly and all three rows render |
| Runic, whole project | Passed |
| Original versus final tiny tables | Forward, reference inverse, and bulk inverse Float64 digests match exactly |
| Complete default fine-slab kernel | Zero allocations over 4,096 points |
| Full release candidates | All three 256³ payloads generated and validated successfully |

The optional browser E2E lane was skipped. Runtime and ABI layouts did not
change; the compiled-library smoke lane is available explicitly in CI.

The original and final numerical digests are in [baseline.toml](baseline.toml)
and [final-baseline.toml](final-baseline.toml). The real-data allocation
measurement is in [allocations.toml](allocations.toml).

## Resource sizing

The corrected [scaling measurement](corrected-scaling.toml) took 22.4, 11.8,
5.91, and 4.66 seconds for the representative fine pass with 1, 2, 4, and 8
workers. Fits used one CPU each; release table jobs used four CPUs each. No
exclusive node was reserved. A task-local scratch binding is enforced by
regression tests, and superseded checkpoints from the earlier scheduler
revision were excluded by their implementation-dependent job identities.

[SLURM accounting](slurm-accounting.txt) records completed job IDs, allocated
CPUs, compute nodes, wall time, total CPU time, and peak resident memory.
See the [experiment guide](../../EXPERIMENTS.md) for APIs and reproduction.
