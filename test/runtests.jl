# Fast runtime tests. They use only in-memory synthetic tables, so they run
# in CI without spectral data, an optimizer, or the precompute pipeline.
#
# Precompute-side numerical tests belong to `precompute/test/`; the compiled
# ABI smoke tests belong to `build/smoke/`.

using Pkg: Pkg
using Test
using PaintMix
using StaticArrays: SVector

include("fixtures/synthetic.jl")
using .SyntheticFixtures

@testset "PaintMix" begin
    @testset "runtime dependency isolation" begin
        # The runtime depends on StaticArrays and nothing else; precompute,
        # build, and benchmark dependencies must not leak in.
        # test/consumer_env.jl checks the resolved closure of a fresh
        # consumer environment.
        project = Pkg.TOML.parsefile(joinpath(pkgdir(PaintMix), "Project.toml"))
        @test sort(collect(keys(get(project, "deps", Dict{String, Any}())))) == ["StaticArrays"]
    end
    include("test_lookup.jl")
    include("test_mixing.jl")
    include("test_colorspace.jl")
    include("test_format.jl")
    include("test_allocations.jl")
    include("test_quality.jl")
end
