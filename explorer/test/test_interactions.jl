using Test
using Bonito
using Explorer
using WGLMakie: Makie

isdefined(@__MODULE__, :PlotFixtures) || include("fixtures.jl")

@testset "plot interaction modes" begin
    fig = Makie.Figure()
    ax = Makie.Axis(fig[1, 1]; limits = (0.0, 0.8, 0.0, 0.9))
    session = Bonito.Session(Bonito.NoConnection(); asset_server = Bonito.NoServer())
    try
        # The paint canvas uses this policy for every built-in gesture.
        Explorer.navigation!(ax, false)
        @test all(!first(Makie.interactions(ax)[name]) for name in Explorer.AXIS_NAVIGATION)

        mode = Bonito.Observable("pick")
        Explorer.cie_tools(session, ax, mode)
        data = PlotFixtures.explorer_data()
        Explorer._cie_probe!(session, ax, data, data.quad, data.derived, mode)
        handler = last(Makie.interactions(ax)[:probe])
        event = Makie.MouseEvent(
            Makie.MouseEventTypes.leftclick, 1.0, Makie.Point2d(0.3, 0.3), Makie.Point2f(100, 100),
            0.0, Makie.Point2d(0.3, 0.3), Makie.Point2f(100, 100),
        )
        @test Makie.process_interaction(handler, event, ax) == Makie.Consume(true)
        mode[] = "navigate"
        @test all(first(Makie.interactions(ax)[name]) for name in Explorer.AXIS_NAVIGATION)
        @test Makie.process_interaction(handler, event, ax) == Makie.Consume(false)
        mode[] = "pick"
        @test all(!first(Makie.interactions(ax)[name]) for name in Explorer.AXIS_NAVIGATION)
    finally
        close(session)
    end
end
