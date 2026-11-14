# Expose sampled inverse-solver traces for explorer RGB-to-concentration explanations

The explorer needs to explain the bridge from fitted pigments to the inverse
lookup table (paper Eq. 9 and §3.2). `/tables` shows finished concentration,
basin and residual slices and aggregate validation, but cannot explain the
optimisation behind a selected RGB point.

`precompute/src/solvers.jl::SlabSummary` aggregates iteration/evaluation
counts, unconverged count and maximum SSE. `UnmixResult` retains the winning
concentrations and summary statistics, but not seed comparisons, accepted
iterates, active-face transitions or the reference-solver gap. These cannot
be recovered from quantized payload bytes.

## Requested work

- Provide a separate opt-in diagnostic solve for a bounded, deterministic
  probe set: near-black, pure-pigment neighbourhoods, unattainable RGB and
  cases selected from quality tails.
- Record target RGB, evaluator/strategy/settings identity, initial seeds,
  candidate active faces/restarts, accepted concentration iterates,
  reconstructed linear RGB, SSE, stationarity/termination and selected winner.
  Include a reference-solver comparison where feasible, labelled as numerical
  evidence rather than a global-optimality guarantee.
- Version the artifact and tie it to config/input/surrogate/payload identity.
  Limit probes, iterations and bytes. Keep the runtime evaluator and default
  allocation-free inverse hot path unchanged; never trace all 256³ vertices.
- Let `/tables` show a recorded probe through target → constrained solve →
  floating concentrations → joint byte rounding → runtime interpolation and
  residual. Distinguish continuous solver, interpolation and quantization
  errors. Label a fresh diagnostic solve separately from generation history
  when the selected point was not recorded.

## Acceptance

Synthetic interior/boundary/unattainable probes preserve the simplex and their
trace winners agree with the normal solver. Quantization is a separate step.
Existing payloads without traces work. Check default hot-path allocations and
throughput. Document capture and regeneration commands.
