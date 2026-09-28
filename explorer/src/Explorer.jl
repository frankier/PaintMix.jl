"""
    Explorer

Interactive, server-rendered viewer for the PaintMix input database, the
precompute intermediates, and the runtime payload.

This is a development tool. It is a standalone project, not a member of the
root workspace, and it is the only place in the repository that depends on a
web server, a browser toolkit, or a template engine. The dependency direction
is one way: `Explorer -> PaintMixPrecompute -> PaintMix`.

Phase 0 exposes the data layer and the server skeleton: the `/`, `/provenance`,
and `/healthz` routes. Later phases add one figure page at a time.
"""
module Explorer

using Bonnie
using Bonito
using HTTP
using OteraEngine
using Oxygen
using PaintMix
using PaintMixPrecompute
using StaticArrays: SVector
using TOML

# Module-local Oxygen state. `@oxidize` binds the route table to this module
# instead of Oxygen's global context, so several viewers (and the test suite)
# can run in one process without sharing routes.
@oxidize

include("colorimetry.jl")
include("data.jl")
include("templates.jl")
include("routes.jl")

export ExplorerData,
    DerivedCurves,
    encoded,
    linear,
    xy_of_linear,
    load_explorer_data,
    serve_explorer,
    close_explorer,
    register_routes!,
    EXPLORER_DATA

end # module Explorer
