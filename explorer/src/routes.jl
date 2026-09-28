# Oxygen routes and JSON endpoints.
#
# The route table lives in this module's `CONTEXT` (see `@oxidize` in
# `Explorer.jl`), so the viewer is isolated from other Oxygen users in the
# same process. Registration is explicit and idempotent, which keeps it out of
# precompilation and lets tests start several viewers.

"""
    EXPLORER_DATA

The process-wide [`ExplorerData`](@ref) the route handlers read. Set it with
[`serve_explorer`](@ref); tests may set it directly.
"""
const EXPLORER_DATA = Ref{Union{Nothing, ExplorerData}}(nothing)

const ROUTES_REGISTERED = Ref(false)
const BONNIE = Ref{Any}(nothing)

"""
    explorer_data() -> ExplorerData

The loaded data, or an error naming the missing setup call.
"""
function explorer_data()
    d = EXPLORER_DATA[]
    d === nothing && error(
        "Explorer data is not loaded; call `serve_explorer` or set `EXPLORER_DATA[]`"
    )
    return d
end

# The Bonnie integration handle, built once per process. Kept lazy so its
# session registry and timers are never captured by precompilation.
function bonnie()
    handle = BONNIE[]
    handle === nothing || return handle
    handle = Bonnie.setup!(Val(:oxygen); app = @__MODULE__)
    BONNIE[] = handle
    return handle
end

"""
    status_app(d) -> Bonito.App

A tiny live Bonito app. Phase 0 has no figures yet; this exercises the
bootstrap and websocket path that every later figure page depends on.
"""
function status_app(d::ExplorerData)
    return App() do
        slider = Bonito.Slider(1:10)
        return Bonito.DOM.div(
            Bonito.DOM.p(
                "Live Bonito session probe (model ", PaintMix.model_id(d.model), ")"
            ),
            Bonito.DOM.div(slider, Bonito.DOM.div(slider.value)),
        )
    end
end

"""
    register_routes!()

Register the viewer's routes on this module's Oxygen context. Idempotent.
"""
function register_routes!()
    ROUTES_REGISTERED[] && return
    ROUTES_REGISTERED[] = true

    # Re-read static files on each request: the CSS is edited during development.
    dynamicfiles(STATIC, "/static")

    @get "/" function (req::HTTP.Request)
        d = explorer_data()
        head = head_content()
        body = app_html(status_app(d))
        return Oxygen.html(
            render_page(
                "index"; title = "PaintMix explorer", head = head,
                body = index_body(d) * body
            )
        )
    end

    @get "/provenance" function (req::HTTP.Request)
        return Oxygen.html(
            render_page(
                "provenance"; title = "PaintMix explorer — provenance",
                body = provenance_body(explorer_data())
            )
        )
    end

    @get "/spectra" function (req::HTTP.Request)
        d = explorer_data()
        head = head_content()
        body = app_html(spectra_app(d)) * spectra_caption(d)
        return Oxygen.html(
            render_page(
                "spectra"; title = "PaintMix explorer — spectra", head = head,
                body = body
            )
        )
    end

    # Standalone figure page for debugging and iframes, one per app.
    @get "/fig/spectra" function (req::HTTP.Request)
        return Bonnie.app_page(
            spectra_app(explorer_data()); title = "PaintMix explorer — spectra figure"
        )
    end

    @get "/healthz" function (req::HTTP.Request)
        d = explorer_data()
        return json(
            Dict(
                "status" => "ok",
                "model_id" => PaintMix.model_id(d.model),
                "grid_n" => PaintMix.grid_n(d.model),
                "source" => d.source,
                "sidecar" => d.sidecar !== nothing,
                "notes" => d.notes,
            )
        )
    end

    return
end

"""
    serve_explorer(data; host = "127.0.0.1", port = 8080, async = true)

Install `data`, register the routes, and start the server. Returns
`(; server, handle, data)`. Stop it with [`close_explorer`](@ref).
"""
function serve_explorer(
        data::ExplorerData;
        host::AbstractString = "127.0.0.1", port::Integer = 8080, async::Bool = true
    )
    EXPLORER_DATA[] = data
    register_routes!()
    handle = bonnie()
    server = serve(;
        host = host, port = port, middleware = [handle.middleware], async = async,
        show_banner = false, docs = false, metrics = false,
    )
    return (; server, handle, data)
end

"""
    close_explorer(viewer)

Stop the server and close every remaining Bonito session.
"""
function close_explorer(viewer)
    close(viewer.handle.context.sessions)
    isopen(viewer.server) && close(viewer.server)
    return
end
