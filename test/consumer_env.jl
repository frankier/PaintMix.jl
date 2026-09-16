#!/usr/bin/env julia
# Check that installing the runtime package alone pulls in nothing from the
# precompute or build toolchains.
#
#   julia --project=test consumer_env.jl
#
# The workspace shares dependency resolution, not dependency visibility, so a
# workspace test proves nothing about what a consumer gets. This script builds
# a fresh environment in a temporary directory, `develop`s only PaintMix into
# it, resolves, and reports the resulting dependency names. Offline: the
# runtime package has no non-stdlib dependencies to resolve.

using Pkg
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))

const FORBIDDEN = ["Optim", "ForwardDiff", "JuliaLibWrapping", "JuliaC", "JLWInterop", "DuckDB", "XLSX"]

@testset "fresh consumer environment" begin
    dir = mktempdir()
    Pkg.activate(dir; io = devnull)
    Pkg.develop(PackageSpec(; path = ROOT); io = devnull)
    Pkg.instantiate(; io = devnull)
    deps = Pkg.dependencies()
    names = sort([d.name for d in values(deps) if d.name !== nothing])
    @test "PaintMix" in names
    found = intersect(names, FORBIDDEN)
    @test isempty(found)
    println("consumer environment resolved $(length(names)) packages: ",
        join(names, ", "))
    @test length(names) <= 6  # PaintMix and a handful of stdlibs
end
