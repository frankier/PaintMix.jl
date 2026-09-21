# Fast runtime tests. They use only in-memory synthetic tables, so they run
# in CI without spectral data, an optimizer, or the precompute pipeline.
#
# Precompute-side numerical tests belong to `precompute/test/`; the compiled
# ABI smoke tests belong to `build/smoke/`.

using Pkg: Pkg
using Test
using PaintMix

include("fixtures/synthetic.jl")
using .SyntheticFixtures

@testset "PaintMix" begin
    @testset "runtime dependency isolation" begin
        # The runtime package must stay evaluable with no non-stdlib
        # dependencies; test/consumer_env.jl checks a fresh consumer
        # environment, and this catches the change that would break it.
        project = Pkg.TOML.parsefile(joinpath(pkgdir(PaintMix), "Project.toml"))
        @test isempty(get(project, "deps", Dict{String, Any}()))
    end
    include("test_lookup.jl")
    include("test_mixing.jl")
    include("test_colorspace.jl")
    include("test_format.jl")
    include("test_allocations.jl")
end
