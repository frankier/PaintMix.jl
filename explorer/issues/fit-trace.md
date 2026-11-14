# Record bounded surrogate-fit snapshots for explorer continuation replay

The explorer has endpoint spectra and scalar continuation curves, but cannot
show how surrogate pigments and their mixture gamut move during optimisation
(paper §3.1, Fig. 8, Eqs. 15–17).

`precompute/src/surrogates.jl::fit_surrogates` records alpha, Epush, Epull,
objective, cumulative seconds and solve statistics after each continuation
solve. `SurrogateFit`, the cache in `pipeline.jl`, and the export in
`export.jl` retain only the final physical K/S coefficients. The explorer
cannot reconstruct intermediate spectra or gamuts from those scalar summaries.

## Requested work

- Add an opt-in, versioned trace artifact with initial pigments and physical
  K/S snapshots after each continuation solve, joined to alpha, step, terms,
  elapsed time, termination, and model/config/input identity. Use physical
  coefficients so every parameterization can be replayed.
- Record wavelength/pigment order, fixed boundary concentration sample
  identities/weights, and an independent validation corpus or deterministic
  recipe. Include worst cube violation and perceptual drift per snapshot;
  distinguish sampled containment from a proof.
- Optionally sample accepted optimiser iterations through compatible
  callbacks, with explicit stride/size limits. Avoid retaining every objective
  evaluation or duplicating large RGB arrays.
- Carry trace references through cache, resume, export and experiment reports;
  preserve loading of existing sidecars without a trace.
- Use the artifact in `/fit` for a step scrubber synchronizing spectra,
  RGB/CIE gamut movement, and loss terms. Label missing traces explicitly.

## Acceptance

A small synthetic fit can replay initial, intermediate and final physical
coefficients. Its final snapshot agrees with the exported fit. Traces respect
size limits and determinism apart from timings. Old/no-trace sidecars still
render. Document the schema, capture cost, and regeneration commands.
