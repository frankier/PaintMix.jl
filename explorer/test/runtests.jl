# Explorer test entry point. The fast lane is browser-free. The browser lane is
# opt-in and heavy; see `test_e2e.jl`.
#
#   julia --project=explorer -e 'using Pkg; Pkg.test()'
#   EXPLORER_E2E=1 julia --project=explorer -e 'using Pkg; Pkg.test()'

using Test

@testset "Explorer" begin
    include("test_colorimetry.jl")
    include("test_spectra.jl")
    include("test_cie.jl")
    include("test_palette.jl")
    include("test_paint.jl")
    include("test_data.jl")
    include("test_server.jl")
end

include("test_e2e.jl")
