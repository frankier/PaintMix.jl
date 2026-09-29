# Fast server tests: real HTTP and a real websocket, no browser. The websocket
# canary mirrors Bonnie's `test/test_canary.jl` so a broken Bonito embedding
# surface is caught here.

using Test
using HTTP
using HTTP.WebSockets: WebSockets
using Bonito
using Bonnie
using Explorer
using PaintMix

include("helpers.jl")

const SYNTHETIC_ID = "00000000000000000000000000000000"

@testset "server" begin
    @testset "routes render through Bonnie" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        base = "http://127.0.0.1:$port"
        try
            resp = HTTP.get("$base/"; retry = false)
            @test resp.status == 200
            body = String(resp.body)
            @test startswith(body, "<!doctype html>")
            @test occursin("Bonito.init_session", body)
            @test occursin("/bonito/assets/", body)
            @test occursin(SYNTHETIC_ID, body)
            # A forgotten `|> safe` would render the bootstrap as text.
            @test !occursin("&lt;script", body)

            css = HTTP.get("$base/static/explorer.css"; retry = false)
            @test css.status == 200
            @test occursin("--bg", String(css.body))

            resp = HTTP.get("$base/nope"; status_exception = false, retry = false)
            @test resp.status == 404
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/healthz JSON contract" begin
        data = Explorer.ExplorerData(;
            model = Explorer.synthetic_model(), source = "synthetic",
            notes = ["a note"],
        )
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            resp = HTTP.get("http://127.0.0.1:$port/healthz"; retry = false)
            @test resp.status == 200
            @test occursin("application/json", string(resp.headers))
            body = String(resp.body)
            @test occursin("\"status\":\"ok\"", body)
            @test occursin("\"model_id\":\"$SYNTHETIC_ID\"", body)
            @test occursin("\"sidecar\":false", body)
            @test occursin("a note", body)
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "degraded page shows the not-available state" begin
        data = Explorer.ExplorerData(;
            model = Explorer.synthetic_model(), source = "synthetic",
            notes = ["sidecar disabled (--no-sidecar)"],
        )
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/provenance"; retry = false).body)
            @test occursin("Not available: no sidecar.", body)
            @test occursin("sidecar disabled (--no-sidecar)", body)
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/spectra degrades without the spectral reference" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model(), source = "synthetic")
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/spectra"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("Spectral reference not available", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/spectra"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/spectra renders the figure against real spectra" begin
        ref = spectral_reference()
        if ref === nothing
            @test_skip "input database absent; spectra figure not served"
        else
            curves = Explorer.derive_curves(
                ref.spectra, ref.quad, ref.spectral, nothing,
                Explorer.synthetic_model(), ref.db,
            )
            data = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
                derived = curves,
            )
            port = free_port()
            viewer = Explorer.serve_explorer(data; port = port, async = true)
            try
                body = String(HTTP.get("http://127.0.0.1:$port/spectra"; retry = false).body)
                @test occursin("Bonito.init_session", body)
                @test occursin("canvas", body)
                @test occursin("Spectra caption", body)
                @test occursin("kins", body)
            finally
                Explorer.close_explorer(viewer)
            end
        end
    end

    @testset "/cie degrades without the spectral reference" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model(), source = "synthetic")
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/cie"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("Spectral reference not available", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/cie"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/cie renders the diagram against real spectra" begin
        ref = spectral_reference()
        if ref === nothing
            @test_skip "input database absent; CIE figure not served"
        else
            curves = Explorer.derive_curves(
                ref.spectra, ref.quad, ref.spectral, nothing,
                Explorer.synthetic_model(), ref.db,
            )
            data = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
                derived = curves,
            )
            port = free_port()
            viewer = Explorer.serve_explorer(data; port = port, async = true)
            try
                body = String(HTTP.get("http://127.0.0.1:$port/cie"; retry = false).body)
                @test occursin("Bonito.init_session", body)
                @test occursin("canvas", body)
                @test occursin("CIE 1931", body)
                @test occursin("Nominal paper colors", body)
                @test occursin("Click to probe", body)
            finally
                Explorer.close_explorer(viewer)
            end
        end
    end

    @testset "/palette renders the cube and mixer from the model alone" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/palette"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("canvas", body)
            @test occursin("Pairwise mixture ramps", body)
            @test occursin("Nominal paper colors", body)
            @test occursin("linear-gradient", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/palette"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/paint renders the canvas and mixer from the model alone" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/paint"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("canvas", body)
            @test occursin("How the dab works", body)
            @test occursin("Mixer ramp", body)
            @test occursin("linear-gradient", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/paint"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/tables renders slices from the model alone" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/tables"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("canvas", body)
            @test occursin("Forward plane", body)
            @test occursin("fourth concentration", body)
            @test occursin("Acceptance gates", body)
            @test occursin("Not available: no sidecar.", body)
            @test occursin("not available: the input database is absent", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/tables"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/fit degrades without spectra and sidecar" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/fit"; retry = false).body)
            @test occursin("Bonito.init_session", body)
            @test occursin("Fit data not available", body)
            @test occursin("Not available: no sidecar.", body)
            @test !occursin("&lt;script", body)

            fig = HTTP.get("http://127.0.0.1:$port/fig/fit"; retry = false)
            @test fig.status == 200
            @test occursin("Bonito.init_session", String(fig.body))
        finally
            Explorer.close_explorer(viewer)
        end
    end

    @testset "/tables renders the spectral-error map against real spectra" begin
        ref = spectral_reference()
        if ref === nothing
            @test_skip "input database absent; tables spectral map not served"
        else
            data = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
            )
            port = free_port()
            viewer = Explorer.serve_explorer(data; port = port, async = true)
            try
                body = String(HTTP.get("http://127.0.0.1:$port/tables"; retry = false).body)
                @test occursin("Bonito.init_session", body)
                @test occursin("available at 96² resolution", body)
                @test occursin("Acceptance gates", body)
            finally
                Explorer.close_explorer(viewer)
            end
        end
    end

    @testset "websocket canary: session handshake" begin
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        try
            body = String(HTTP.get("http://127.0.0.1:$port/"; retry = false).body)
            id = root_session_id(body)
            sessions = viewer.handle.context.sessions
            @test length(sessions) == 1
            session = Bonnie.lookup(sessions, id)
            @test session !== nothing

            WebSockets.open("ws://127.0.0.1:$port/bonito/ws/$id") do ws
                WebSockets.send(
                    ws, client_message(
                        session, Dict{String, Any}(
                            "msg_type" => Bonito.JSDoneLoading, "exception" => "nothing",
                            "session" => id,
                        )
                    )
                )
                @test wait_for(() -> Bonito.isready(session; throw = false))
                @test isopen(session)
            end
        finally
            Explorer.close_explorer(viewer)
        end
    end
end
