# `/paint`: the brush canvas and the mixer ramp beside it.
#
# The canvas is a coarse bitmap of latents streamed as a display matrix. A
# drag dabs the disc through `PaintMix.mix`'s own blend, so overlapping
# strokes compose the way the runtime mixes. The naive-RGB toggle lerps the
# display instead and re-encodes, which is the comparison the page is for.
# Everything needs only the runtime payload.

const PAINT_DEFAULT_SIZE = 128
const PAINT_MAX_RADIUS = 32

# Frames are throttled while a drag is in flight; the release flushes.
const PAINT_FRAME_INTERVAL = 0.05

"""
    paint_app(d::ExplorerData) -> Bonito.App

The interactive paint canvas: a drag-painted bitmap with pigment-weight,
radius, flow, brush-source, and paint-vs-RGB controls, plus the mixer ramp.
"""
function paint_app(d::ExplorerData)
    return App() do session
        paint_dom(session, d)
    end
end

function paint_dom(session, d::ExplorerData)
    model = d.model
    canvas = paint_canvas(model; width = PAINT_DEFAULT_SIZE, height = PAINT_DEFAULT_SIZE)
    frame = Bonito.Observable(canvas.display)

    fig = Makie.Figure(; size = (560, 560))
    ax = Makie.Axis(fig[1, 1]; title = "Paint canvas")
    Makie.image!(ax, frame; interpolate = false)
    Makie.hidedecorations!(ax)
    Makie.hidespines!(ax)

    labels = [pigment_label(d, i) for i in 1:4]
    colors = pigment_rgbs(model)

    weight_sliders = [Bonito.Slider(0.0:0.01:1.0; value = i == 1 ? 1.0 : 0.0) for i in 1:4]
    radius_slider = Bonito.Slider(1:PAINT_MAX_RADIUS; value = 8)
    flow_slider = Bonito.Slider(0.05:0.05:1.0; value = 0.5)
    rgb_box = Bonito.Checkbox(false)
    source_dropdown = Bonito.Dropdown(["pigment mixture", "picked color"]; index = 1)
    color_obs = Bonito.Observable("#c42761")
    clear_button = Bonito.Button("Clear")

    weight_brush = Makie.lift(
        (w1, w2, w3, w4) -> brush_from_weights(model, (w1, w2, w3, w4)),
        weight_sliders[1].value, weight_sliders[2].value, weight_sliders[3].value,
        weight_sliders[4].value,
    )
    picked_brush = Makie.lift(hex -> brush_from_hex(model, hex), color_obs)
    brush = Makie.lift(
        (src, wb, pb) -> src == "picked color" ? pb : wb, source_dropdown.value,
        weight_brush, picked_brush,
    )

    # --- interaction -------------------------------------------------------
    dragging = Ref(false)
    last_frame = Ref(0.0)
    function flush_frame!(force::Bool)
        now = time()
        (force || now - last_frame[] >= PAINT_FRAME_INTERVAL) || return
        last_frame[] = now
        frame[] = canvas_image(canvas)
        return
    end
    function paint_at(pos)
        h, _ = size(canvas.latent)
        col, row = paint_pixel(pos, h)
        mode = rgb_box.value[] ? Val(:rgb) : Val(:paint)
        dab!(
            canvas, model, brush[], col, row, radius_slider.value[],
            flow_slider.value[], mode,
        )
        flush_frame!(false)
        return
    end

    mouse = Makie.events(ax.scene)
    Bonito.on(session, mouse.mousebutton) do event
        if event.button == Makie.Mouse.left
            if event.action == Makie.Mouse.press
                dragging[] = true
                paint_at(Makie.mouseposition(ax.scene))
            elseif event.action == Makie.Mouse.release
                dragging[] = false
                flush_frame!(true)
            end
        end
    end
    Bonito.on(session, mouse.mouseposition) do pos
        dragging[] && paint_at(pos)
    end
    Bonito.on(session, clear_button.value) do _
        clear!(canvas, model)
        flush_frame!(true)
    end

    # --- readouts ----------------------------------------------------------
    picked = Bonito.Observable{Any}(_pick_table(picked_result(model, color_obs[])))
    Bonito.on(session, color_obs) do hex
        picked[] = _pick_table(picked_result(model, hex))
    end

    weight_controls = [
        Bonito.DOM.div(
            Bonito.DOM.label(
                Bonito.DOM.span(
                    ; class = "swatch", style = "background:$(hex_from_linear(colors[i]))"
                ),
                " ", labels[i],
            ),
            weight_sliders[i],
        ) for i in 1:4
    ]
    ramp_style = Makie.lift(
        (w1, w2, w3, w4) -> ramp_gradient(mixer_ramp(model, (w1, w2, w3, w4))),
        weight_sliders[1].value, weight_sliders[2].value, weight_sliders[3].value,
        weight_sliders[4].value,
    )
    controls = Bonito.DOM.div(
        Bonito.DOM.h3("Brush"),
        Bonito.DOM.label(source_dropdown, " source"),
        Bonito.DOM.h3("Pigment weights"),
        weight_controls...,
        Bonito.DOM.h3("Brush radius (px)"),
        radius_slider,
        Bonito.DOM.h3("Flow"),
        flow_slider,
        Bonito.DOM.h3("Blend"),
        Bonito.DOM.label(rgb_box, " linear RGB (naive)"),
        Bonito.DOM.h3("Picked color"),
        Bonito.DOM.input(
            ; type = "color", value = color_obs,
            onchange = js"event => $(color_obs).notify(event.srcElement.value);",
        ),
        Bonito.DOM.div(picked; class = "readout-wrap"),
        Bonito.DOM.h3("Actions"),
        clear_button,
        Bonito.DOM.h3("Mixer ramp"),
        Bonito.DOM.div(; class = "ramp mixer-ramp", style = ramp_style),
        Bonito.DOM.p(
            "The ramp runs from paper white to the weighted pigment mixture.";
            class = "muted",
        ),
    )
    return Bonito.DOM.div(
        Bonito.DOM.div(fig; class = "figure-panel"), controls; class = "spectra-grid"
    )
end

function _pick_table(::Nothing)
    return Bonito.DOM.p("Enter a #rrggbb color."; class = "muted")
end

function _pick_table(r::PickResult)
    swatch = Bonito.DOM.span(; class = "swatch", style = "background:$(r.hex)")
    row(label, value) = Bonito.DOM.tr(Bonito.DOM.th(label), Bonito.DOM.td(value))
    return Bonito.DOM.table(
        Bonito.DOM.caption("Picked color"),
        Bonito.DOM.tbody(
            row("hex", Bonito.DOM.span(swatch, " ", r.hex)),
            row(
                "concentrations",
                @sprintf(
                    "%.3f, %.3f, %.3f, %.3f", r.concentrations[1],
                    r.concentrations[2], r.concentrations[3], r.concentrations[4]
                ),
            ),
            row(
                "residual",
                @sprintf("%.4f, %.4f, %.4f", r.residual[1], r.residual[2], r.residual[3]),
            ),
        );
        class = "readout",
    )
end

"""
    paint_caption(d::ExplorerData) -> String

The static caption under the canvas: the dab formula, the default size, and
the brush sources. Needs only the runtime payload.
"""
function paint_caption(d::ExplorerData)
    io = IOBuffer()
    print(io, "<section class=\"cards\"><article class=\"card\">")
    print(io, "<h2>How the dab works</h2>")
    print(
        io,
        "<p>Each covered pixel lerps its latent toward the brush latent, ",
        "then decodes:</p>",
    )
    print(
        io,
        "<p><code>z &larr; Latent((1&minus;&alpha;)z.c + &alpha;z_brush.c, ",
        "(1&minus;&alpha;)z.r + &alpha;z_brush.r)</code>, then <code>decode</code>.</p>",
    )
    print(
        io,
        "<p class=\"muted\">This is exactly <code>mix</code>, so overlapping ",
        "strokes compose without re-encoding and the signed residual survives ",
        "the stroke.</p>",
    )
    print(io, "</article>")

    print(io, "<article class=\"card\"><h2>Canvas</h2>")
    print(
        io, html_table(
            [
                ("default size", "$(PAINT_DEFAULT_SIZE) x $(PAINT_DEFAULT_SIZE) px"),
                ("brush sources", "pigment mixture or picked color"),
                ("paint mode", "lerp the latents, decode"),
                ("RGB mode", "lerp the linear display, re-encode"),
                ("frame throttle", "$(round(Int, 1 / PAINT_FRAME_INTERVAL)) Hz while dragging"),
            ]
        )
    )
    print(
        io,
        "<p class=\"muted\">The display matrix is the only large buffer that ",
        "crosses the websocket, which is why the default edge is ",
        "$(PAINT_DEFAULT_SIZE).</p>",
    )
    print(io, "</article></section>")
    return String(take!(io))
end
