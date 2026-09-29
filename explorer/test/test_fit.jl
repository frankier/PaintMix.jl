# Fit-page tests. The history and diagnostics accessors run on a hand-built
# sidecar dictionary, so they never touch disk. The real sidecar is checked
# when it is present.

using Test
using PaintMix
using PaintMixPrecompute
using Explorer
using TOML
using WGLMakie: Makie

function fit_history_entries()
    return [
        Dict{String, Any}(
            "step" => 1, "alpha" => 1.0e5, "objective" => 0.5, "Epush" => 0.4,
            "Epull" => 1.0e-3, "seconds" => 1.0, "iterations" => 10,
            "converged" => false,
        ),
        Dict{String, Any}(
            "step" => 2, "alpha" => 5.0e4, "objective" => 0.1, "Epush" => 1.0e-8,
            "Epull" => 2.0e-3, "seconds" => 2.0, "iterations" => 12,
            "converged" => true,
        ),
        Dict{String, Any}(
            "step" => 3, "alpha" => 2.5e4, "objective" => 0.05, "Epush" => 1.0e-9,
            "Epull" => 3.0e-3, "seconds" => 3.0, "iterations" => 8,
            "converged" => true,
        ),
    ]
end

function fit_sidecar()
    return Dict{String, Any}(
        "surrogate" => Dict{String, Any}(
            "history" => fit_history_entries(),
            "diagnostics" => Dict{String, Any}(
                "Epush_fit" => 1.0e-9, "Epull_dense" => 1.0e-3,
                "max_cube_violation" => 0.0, "max_oklab_deviation" => 0.02,
                "dense_divisions" => 40,
            ),
        ),
    )
end

@testset "fit" begin
    sidecar = fit_sidecar()
    cfg = Dict{String, Any}(
        "surrogate" => Dict{String, Any}("push_tolerance" => 1.0e-8)
    )
    data = Explorer.ExplorerData(;
        model = Explorer.synthetic_model(), sidecar = sidecar, cfg = cfg
    )

    @testset "history accessors" begin
        h = Explorer.fit_history(data)
        @test h isa Explorer.FitHistory
        @test h.steps == [1, 2, 3]
        @test h.epush[2] == 1.0e-8
        @test h.converged == [false, true, true]
        @test Explorer.push_tolerance(data) == 1.0e-8
        @test Explorer.push_met_step(h, 1.0e-8) == 2
        @test Explorer.push_met_step(h, 1.0e-12) === nothing
        @test Explorer.fit_diagnostics(data)["Epush_fit"] == 1.0e-9
    end

    @testset "missing data degrades" begin
        @test Explorer.fit_history(nothing) === nothing
        @test Explorer.fit_history(Dict{String, Any}()) === nothing
        @test Explorer.fit_diagnostics(nothing) === nothing
        @test Explorer.fit_diagnostics(Dict{String, Any}()) === nothing
        @test Explorer.push_tolerance(Dict{String, Any}()) == 1.0e-8
        # A history entry missing a field becomes NaN rather than an error.
        partial = Dict{String, Any}(
            "surrogate" => Dict{String, Any}("history" => [Dict{String, Any}()])
        )
        hp = Explorer.fit_history(partial)
        @test hp !== nothing
        @test isnan(hp.epush[1])
    end

    @testset "server-rendered fit cards" begin
        caption = Explorer.fit_caption(data)
        @test occursin("Continuation history", caption)
        @test occursin("first step meeting tolerance", caption)
        @test occursin("Fit diagnostics", caption)
        @test occursin("Chromaticity displacement", caption)
        @test occursin("Not available: no fitted surrogate.", caption)

        empty = Explorer.fit_caption(
            Explorer.ExplorerData(; model = Explorer.synthetic_model())
        )
        @test occursin("Not available: no sidecar.", empty)
    end

    @testset "figures build" begin
        h = Explorer.fit_history(data)
        @test Explorer.fit_history_figure(h, 1.0e-8) isa Makie.Figure
        c = synthetic_curves()
        @test Explorer.fit_spectra_figure(c) isa Makie.Figure
    end

    # The real sidecar, when present, must satisfy the same accessors.
    release_dir = joinpath(Explorer.OUTPUT, "release")
    sidecars = isdir(release_dir) ?
        filter(f -> endswith(f, ".toml"), readdir(release_dir)) : String[]
    if isempty(sidecars)
        @test_skip "no release sidecar present; real history not checked"
    else
        @testset "real release sidecar" begin
            real = TOML.parsefile(joinpath(release_dir, first(sort(sidecars))))
            d = Explorer.ExplorerData(;
                model = Explorer.synthetic_model(), sidecar = real
            )
            h = Explorer.fit_history(d)
            @test h !== nothing
            @test length(h.steps) == length(h.epush)
            @test all(>=(0), h.epush)
            @test Explorer.fit_diagnostics(d) !== nothing
            @test !isempty(Explorer.acceptance_gate_rows(d))
            caption = Explorer.fit_caption(d)
            @test occursin("Continuation history", caption)
        end
    end
end
