# Cleanup follow-ups

The five follow-ups were implemented in order on `main`: resolved solver
objects/options and dependency identities; shared layout/diagnostic/slab
implementations; gauge coordinates and a common physical stopping criterion;
experimental precision/interpolation payloads with Julia runtime loading;
and verification of the real release payload and both browser lanes.
The original dirty Explorer work is separately saved on `explorer-fixups`
(`cf6ef9d`). No candidate was promoted and the default payload is unchanged.

The fitting/table snapshot is `de0dfea`; the precision/runtime snapshot is
`18017b4`. The experiment guide documents commands and conventions. Full
candidate payloads and sidecars are retained under the ignored
`precompute/output/followups/` tree. The verified compiled output is under
`build/out-release/`; it is non-bundled and needs the tested Julia runtime's
shared libraries. These generated artifacts are not committed.

## Physical stopping experiment

All four release-profile fits use the same physical initialization, 20-division
fit quadrature, 200-iteration step cap, 7,200-second budget, and continuation
to alpha >= 1e-6 with a 40-halving cap. Success requires maximum distance
outside the linear-RGB cube <= 1e-3 on 3,792 disjoint held-out concentrations.
Their common digest is recorded in every row. Each candidate has full 256³
tables and validation. Authoritative physical K/S, coordinates and complete
continuation traces are retained in [fits/](fits/).

**No fit reached the held-out target.** Every fit exhausted the continuation
schedule. Native optimizer success and small dense-grid error did not imply
held-out success, particularly near pure pigments. These are completed-budget
timings, not times to equal achieved accuracy. The softplus fits ran on a
Skylake CPU cohort; direct and gauge fits ran on a Cascade Lake cohort.
Consequently the wall times do not establish a strict solver speed ranking.

| Variant | Fit seconds | Final held-out violation | Dense violation | Mixing Oklab p99 |
| --- | ---: | ---: | ---: | ---: |
| direct-lbfgsb | 218.7 | 0.00601325 | 0.00601325 | 0.0347657 |
| gauge-lbfgsb | 108.0 | 0.00910577 | 0.00217818 | 0.0617079 |
| softplus-lbfgs | 838.4 | 0.0123582 | 0.00122532 | 0.0615893 |
| softplus-lbfgsb | 132.4 | 0.00401187 | 0.000346489 | 0.0650923 |

All hard table gates pass; the provisional mixing quality gates fail for all
four candidates. [comparison.toml](comparison.toml) retains the full metrics,
job/model IDs, native termination, actual held-out violations and target flags.
The quantization fraction was recounted from each unchanged payload using the
corrected vertex counter from `18017b4`; the previous diagnostic and an explicit
correction marker are retained for audit. Acceptance gates are unchanged.
The gauge removes one unidentifiable common K/S scale per wavelength, but it
does not resolve the near-vertex gamut problem within this budget.

## Precision and interpolation

The four precision rows share the earlier direct-coefficient candidate's
physical model and verified floating inverse checkpoints. No fitting or
inverse solving changes between these rows. Forward RGB remains 8-bit and
the inverse grid remains 256³. The common corpus has 512 colors, 20,000 mixes
and boundary-rich concentration probes. Full `.pmxe` files were written,
read through the Julia runtime, checked and measured in one single-CPU job.

| Inverse / forward | Payload MiB | Mixing Oklab p99 | Encoded channel p99 | Forward channel p99 | Concentration-only channel p99 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 8-cube | 96.000 | 0.0356442 | 0.0981935 | 0.0266552 | 0.0649834 |
| 8-simplex | 56.094 | 0.0356485 | 0.0979922 | 0.0247742 | 0.0649834 |
| 16-cube | 144.000 | 0.0217084 | 0.0633512 | 0.0266552 | 0.000671807 |
| 16-simplex | 104.094 | 0.0217286 | 0.0628512 | 0.0247742 | 0.000671807 |

Wider inverse storage greatly reduces the isolated concentration error and
reduces complete mixing error, but all four variants still miss the
provisional mixing quality targets. The packed simplex forward lattice
saves about five sixths of forward storage and changes mixing accuracy only
slightly on this model. With a 16-bit inverse it produces a 104.094 MiB
candidate versus 96 MiB for the original byte cube pair.

All four variants allocate zero bytes over the warm 20,000-mix probe and
have Float32 round-trip maximum error 1.49e-8. The minimum of three warm
probes measured 346-365 ns/mix; this is a fixed-corpus measurement, not a
cross-machine throughput guarantee. [precision.toml](precision.toml) retains
source/payload digests, complete error distributions, timings and versions.
The approximation layers are not additive.

## Verification

| Check | Result |
| --- | --- |
| Runtime tests, including experimental format and allocation checks | 6,274 passed |
| Precompute, two workers, including CLI precision workflow and promotion guard | 4,719 passed |
| Fresh consumer in a clean Julia environment | 4 passed; 17 packages, no optimizer/build dependencies |
| Default Float64 forward/reference/bulk table digests | Exactly match the original baseline |
| Whole-project Runic | Passed |
| Real release payload in trimmed library | n=256; model f78bfa4fe793ca748d7647605cb37147 |
| Real-payload C / Python / R smoke | 1,781 / 1,783 / 1,784 passed; zero failures |
| Relocated C / Python / R, Julia absent from PATH | Same checks passed |
| Explorer fast lane with the unchanged legacy release sidecar | 533 passed |
| Main Chromium E2E | 30 passed |
| Saved Explorer fixups Chromium E2E | 35 passed, including stronger stroke assertions |

The compiled image was 98.874 MiB and the linked library 98.057 MiB.
Compilation took 15.691 seconds and linking 0.570 seconds, excluding setup.
The tested toolchain is Julia 1.13.3, JuliaC 0.3.10, JLWInterop 0.2.0,
StaticArrays 1.9.24 and wrapping commit
`4f21b6de6a61590bf8e59c9d6c88b4a62722e177`. Python 3.13.13 / NumPy 2.5.4
and R 4.6.1 / rdyncall 0.10.2 exercised the generated bindings. The browser
lane used its pinned Julia Playwright revision, driver 1.61.1 and Chromium
revision 1228. [release-verification.toml](release-verification.toml) records
payload/library hashes and exact tool metadata.

Numerical work and checks ran on SLURM compute allocations. Fits used one
CPU, tables four, and threaded tests/browser checks two. No exclusive node
was reserved. Scheduler accounting retains actual job states with anonymous
compute labels. Some setup attempts failed before the final checks: missing
modules, NumPy and browser libraries were supplied in private environments.
The smoke and relocation clients passed in job 3090837; an appended consumer
check in the R loader environment crashed. The clean-environment repeat
3090853 passed. Browser libraries were then isolated to browser child
processes; the final browser jobs passed. The failed setup/job statuses are
retained rather than represented as wholly successful allocations.

## Build dependency pins

The build projects now pin JuliaC 0.3.10, JLWInterop 0.2.0, library StaticArrays
1.9.24 and the immutable wrapping revision. The separate tools manifest is
tracked to lock its transitive dependencies. All 28 verified tool source trees
remain unchanged. The shared workspace lock advances StaticArrays from
1.9.22 to 1.9.24; runtime consumer compatibility remains `1`.

The pinned configuration passed the real-payload C/Python/R smoke and
relocation checks again. Runtime tests (6,274), two-worker precompute tests
(4,719) and fresh-consumer checks (4; 17 packages) passed. Original Float64
baseline table digests remain identical. [build-pins.toml](build-pins.toml)
retains the resolved versions, source hashes and repeated verification.

## Remaining work

1. Improve boundary sampling and stopping behavior until the common physical
   target is reached; use fresh independent probes when changing fit samples.
2. Reduce forward/inverse field interpolation error near pure pigments. The
   16-bit and packed-lattice experiments show that concentration precision
   alone is insufficient for the quality gates.
3. Decide whether to stabilize the experimental format, then extend Explorer
   and compiled bindings to it and run a bundled, multi-platform release matrix.
4. Review/merge `explorer-fixups` when desired. It remains separate from main.
