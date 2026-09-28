"""
    Explorer

Interactive, server-rendered viewer for the PaintMix input database, the
precompute intermediates, and the runtime payload.

This is a development tool. It is a standalone project, not a member of the
root workspace, and it is the only place in the repository that depends on a
web server, a browser toolkit, or a template engine. The dependency direction
is one way: `Explorer -> PaintMixPrecompute -> PaintMix`.

Phase 0 exposes the data layer and the server skeleton (`/`, `/provenance`,
`/healthz`); phase 1 adds the `/spectra` figure, phase 2 the `/cie`
chromaticity diagram, and phase 3 the `/palette` cube and mixer. Later phases
add one figure page at a time.
"""
module Explorer

using Bonnie
using Bonito
using HTTP
using OteraEngine
using Oxygen
using PaintMix
using PaintMixPrecompute
using Printf: @sprintf
using StaticArrays: SVector
using TOML
using WGLMakie
using WGLMakie: Makie

# Module-local Oxygen state. `@oxidize` binds the route table to this module
# instead of Oxygen's global context, so several viewers (and the test suite)
# can run in one process without sharing routes.
@oxidize

include("colorimetry.jl")
include("data.jl")
include("palette.jl")
include("templates.jl")
include("apps/spectra.jl")
include("apps/cie.jl")
include("apps/palette.jl")
include("routes.jl")

export ExplorerData,
    DerivedCurves,
    SpectraReadout,
    SpectraReadoutRow,
    PigmentChromaticity,
    CieProbe,
    encoded,
    linear,
    xy_of_linear,
    hex_from_linear,
    spectral_locus,
    srgb_primaries,
    d65_xy,
    wavelength_linear,
    pigment_chromaticities,
    nominal_xy,
    mixture_gamut,
    convex_hull,
    clip_convex,
    probe_xy,
    NOMINAL_COLORS,
    MixCurve,
    MixerResult,
    NominalRow,
    pigment_rgb,
    pigment_rgbs,
    pigment_label,
    nominal_rgb,
    mix_curve,
    ramp_hex,
    ramp_gradient,
    pairwise_ramps,
    mixer_result,
    mixer_ramp,
    nominal_rows,
    palette_cube_app,
    palette_mixer_app,
    palette_ramp_matrix,
    palette_caption,
    SPECTRA_QUANTITY_KEYS,
    quantity_label,
    quantity_log,
    quantity_matrix,
    fitted_matrix,
    wavelength_index,
    spectra_readout,
    spectra_app,
    spectra_caption,
    cie_app,
    cie_caption,
    load_explorer_data,
    serve_explorer,
    close_explorer,
    register_routes!,
    EXPLORER_DATA

end # module Explorer
