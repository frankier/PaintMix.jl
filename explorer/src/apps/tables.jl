# `/tables`: 2D slices through the inverse and forward lookup tables, and the
# validation reports from the sidecar.
#
# The slice figure works from the runtime payload alone. The spectral-error
# map and the validation cards additionally need the spectral reference and the
# sidecar; when either is absent the page shows a note instead of failing.

"""
    tables_slice_data(d, axis, level) -> NamedTuple

Every slice map for one fixed channel and level: the four concentration
fields, the basin map, the residual magnitude, and the spectral-error map
(`nothing` without the spectral reference).
"""
function tables_slice_data(d::ExplorerData, axis::Symbol, level::Real)
    fields = concentration_fields(d.model, axis, level)
    return (;
        fields,
        basin = basin_map(fields),
        residual = residual_field(d.model, axis, level),
        spectral = d.spectral === nothing ? nothing :
            spectral_error_field(d.model, d.spectral, axis, level),
    )
end

"""
    forward_image(plane) -> Matrix{RGBAf}

The forward-plane colors as a display image; the region outside the simplex is
transparent.
"""
function forward_image(plane::ForwardPlane)
    img = Matrix{Makie.RGBAf}(undef, size(plane.rgb)...)
    for idx in eachindex(plane.rgb)
        if plane.valid[idx]
            c = linear_color(plane.rgb[idx])
            img[idx] = Makie.RGBAf(c.r, c.g, c.b, 1.0)
        else
            img[idx] = Makie.RGBAf(0.0, 0.0, 0.0, 0.0)
        end
    end
    return img
end

"""
    tables_app(d::ExplorerData) -> Bonito.App

The slice figure: concentration fields, basin map, residual field, forward
plane, and the reduced-resolution spectral-error map.
"""
function tables_app(d::ExplorerData)
    return App() do session
        tables_dom(session, d)
    end
end

function tables_dom(session, d::ExplorerData)
    model = d.model
    axis_dropdown = Bonito.Dropdown(["B", "G", "R"]; index = 1)
    level_slider = Bonito.Slider(0.0:0.01:1.0; value = 0.5)
    c4_slider = Bonito.Slider(0.0:0.05:0.95; value = 0.0)

    axis_obs = Makie.lift(
        label -> Dict("R" => :r, "G" => :g, "B" => :b)[label], axis_dropdown.value
    )
    slices = Makie.lift(
        (axis, level) -> tables_slice_data(d, axis, level), axis_obs, level_slider.value
    )
    plane = Makie.lift(c4 -> forward_plane(model, c4), c4_slider.value)

    labels = [pigment_label(d, i) for i in 1:4]
    fig = Makie.Figure(; size = (1040, 560))
    xlabel = Makie.lift(a -> slice_axis_labels(a)[1], axis_obs)
    ylabel = Makie.lift(a -> slice_axis_labels(a)[2], axis_obs)

    for slot in 1:4
        ax = Makie.Axis(
            fig[1, 2slot - 1]; title = "c$(slot) — $(labels[slot])",
            xlabel = xlabel, ylabel = ylabel,
        )
        p = Makie.heatmap!(
            ax, Makie.lift(s -> s.fields[slot], slices);
            colorrange = (0.0, 1.0), colormap = :viridis,
        )
        Makie.Colorbar(fig[1, 2slot], p)
    end

    basin_cmap = Makie.cgrad(
        [linear_color(pigment_rgb(model, i)) for i in 1:4]; categorical = true
    )
    basin_ax = Makie.Axis(fig[2, 1]; title = "basin", xlabel = xlabel, ylabel = ylabel)
    basin_plot = Makie.heatmap!(
        basin_ax, Makie.lift(s -> s.basin, slices);
        colormap = basin_cmap, colorrange = (1, 4),
    )
    Makie.Colorbar(fig[2, 2], basin_plot; ticks = (1:4, labels))

    residual_ax = Makie.Axis(
        fig[2, 3]; title = "residual ‖x − M(U(x))‖", xlabel = xlabel, ylabel = ylabel
    )
    residual_plot = Makie.heatmap!(
        residual_ax, Makie.lift(s -> s.residual, slices); colormap = :magma
    )
    Makie.Colorbar(fig[2, 4], residual_plot)

    fwd_ax = Makie.Axis(
        fig[2, 5]; title = "forward M(c), c₄ fixed", xlabel = "c₁", ylabel = "c₂"
    )
    Makie.image!(fwd_ax, Makie.lift(forward_image, plane); interpolate = false)

    spectral_ax = Makie.Axis(
        fig[2, 7]; title = "forward vs spectral", xlabel = xlabel, ylabel = ylabel
    )
    if d.spectral === nothing
        Makie.text!(
            spectral_ax, 0.5, 0.5; text = "spectral reference not available",
            align = (:center, :center), space = :relative,
        )
    else
        spectral_plot = Makie.heatmap!(
            spectral_ax, Makie.lift(s -> s.spectral, slices); colormap = :inferno
        )
        Makie.Colorbar(fig[2, 8], spectral_plot)
    end

    controls = Bonito.DOM.div(
        Bonito.DOM.h3("Slice"),
        Bonito.DOM.label("fixed channel ", axis_dropdown),
        Bonito.DOM.div(Bonito.DOM.span(level_slider.value), class = "readout-lambda"),
        level_slider,
        Bonito.DOM.h3("Forward plane"),
        Bonito.DOM.label("fourth concentration c₄ ", c4_slider),
        Bonito.DOM.p(
            "Rows and columns of every map are the two free channels; the " *
                "fixed channel is the one selected above.";
            class = "muted",
        ),
    )
    return Bonito.DOM.div(
        Bonito.DOM.div(fig; class = "figure-panel"),
        Bonito.DOM.div(controls; class = "spectra-controls"),
    )
end

"""
    tables_caption(d::ExplorerData) -> String

The server-rendered validation cards: acceptance gates, the overridden gates,
and the per-report summaries. Works without the sidecar.
"""
function tables_caption(d::ExplorerData)
    io = IOBuffer()
    print(io, "<section class=\"cards\">")
    _slice_card(io, d)
    _gates_card(io, d)
    _quality_card(io, d)
    print(io, "</section><section class=\"cards\">")
    _padding_card(io, d)
    _continuity_card(io, d)
    print(io, "</section><section class=\"cards\">")
    _roundtrip_card(io, d)
    _quantization_card(io, d)
    print(io, "</section>")
    return String(take!(io))
end

function _slice_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Slice maps</h2>")
    spectral = d.spectral === nothing ? "not available: the input database is absent" :
        "available at $(TABLE_SPECTRAL_N)² resolution"
    print(
        io, html_table(
            [
                ("grid", "$(PaintMix.grid_n(d.model))³"),
                ("inverse source", "stored vertices, not re-encoded"),
                ("forward plane", "c₄ fixed, c₁ and c₂ free"),
                ("spectral error", spectral),
            ]
        )
    )
    print(
        io,
        "<p class=\"muted\">The concentration fields are the inverse table's ",
        "stored bytes; the residual is ‖x − M(U(x))‖ at the same vertices. The ",
        "spectral-error map compares the forward table with the spectral model ",
        "at the concentrations the runtime inverse recovers.</p>",
    )
    print(io, "</article>")
    return
end

function _gates_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Acceptance gates</h2>")
    rows = acceptance_gate_rows(d)
    if isempty(rows)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        overridden = Set(promoted_gate_names(d))
        print(io, "<table class=\"data\"><tbody>")
        for (name, ok) in rows
            status = ok ? "pass" : "FAIL"
            cls = ok ? "" : (name in overridden ? " class=\"warn\"" : " class=\"fail\"")
            print(
                io, "<tr", cls, "><td>", _esc(name), "</td><td>", status, "</td><td>",
                name in overridden ? "overridden on promotion" : "", "</td></tr>",
            )
        end
        print(io, "</tbody></table>")
    end
    overridden = promoted_gate_names(d)
    if !isempty(overridden)
        print(
            io, "<p class=\"warn\">Promoted with overridden gates: ",
            _esc(join(overridden, ", ")), "</p>",
        )
    end
    print(io, "</article>")
    return
end

function _quality_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>LUT vs spectral quality</h2>")
    q = validation_report(d, "quality")
    if !(q isa AbstractDict)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        ch = Base.get(q, "encoded_channel_error", Base.get(q, "channel_error", Dict()))
        ok = Base.get(q, "oklab_error", Dict())
        print(
            io, html_table(
                [
                    ("pairs", Base.get(q, "pairs", "?")),
                    ("corpus colors", Base.get(q, "corpus_colors", "?")),
                    ("encoded channel mean", number(Base.get(ch, "mean", nothing))),
                    ("encoded channel p99", number(Base.get(ch, "p99", nothing))),
                    ("encoded channel max", number(Base.get(ch, "max", nothing))),
                    ("Oklab mean", number(Base.get(ok, "mean", nothing))),
                    ("Oklab p99", number(Base.get(ok, "p99", nothing))),
                    ("Oklab max", number(Base.get(ok, "max", nothing))),
                ]
            )
        )
    end
    print(io, "</article>")
    return
end

function _padding_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Forward padding error by depth</h2>")
    p = validation_report(d, "padding")
    if !(p isa AbstractDict)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        by_depth = Base.get(p, "by_depth", Dict())
        rows = [
            (
                key, number(Base.get(by_depth[key], "mean", nothing)),
                number(Base.get(by_depth[key], "p99", nothing)),
                number(Base.get(by_depth[key], "max", nothing)),
            ) for key in sort(collect(keys(by_depth)))
        ]
        print(io, html_table(rows; headers = ["depth", "mean", "p99", "max"]))
        print(
            io, "<p class=\"muted\">Rule: ", _esc(string(Base.get(p, "padding_rule", "?"))),
            ". Depth 0 is the simplex boundary itself.</p>",
        )
    end
    print(io, "</article>")
    return
end

function _continuity_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Inverse continuity</h2>")
    c = validation_report(d, "continuity")
    if !(c isa AbstractDict)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        jump = Base.get(c, "concentration_jump", Dict())
        decode = Base.get(c, "decode_error", Dict())
        print(
            io, html_table(
                [
                    ("lines", Base.get(c, "lines", "?")),
                    ("concentration jump mean", number(Base.get(jump, "mean", nothing))),
                    ("concentration jump p99", number(Base.get(jump, "p99", nothing))),
                    ("concentration jump max", number(Base.get(jump, "max", nothing))),
                    ("decode error p99", number(Base.get(decode, "p99", nothing))),
                    ("decode error max", number(Base.get(decode, "max", nothing))),
                ]
            )
        )
        print(io, "<p class=\"muted\">", _esc(string(Base.get(c, "note", ""))), "</p>")
    end
    print(io, "</article>")
    return
end

function _roundtrip_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Round trip</h2>")
    r = validation_report(d, "roundtrip")
    if !(r isa AbstractDict)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        f32 = Base.get(r, "float32", Dict())
        f64 = Base.get(r, "float64", Dict())
        print(
            io, html_table(
                [
                    ("samples", Base.get(r, "samples", "?")),
                    ("Float32 mean", number(Base.get(f32, "mean", nothing))),
                    ("Float32 max", number(Base.get(f32, "max", nothing))),
                    ("Float64 mean", number(Base.get(f64, "mean", nothing))),
                    ("Float64 max", number(Base.get(f64, "max", nothing))),
                ]
            )
        )
    end
    print(io, "</article>")
    return
end

function _quantization_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Quantization</h2>")
    q = validation_report(d, "quantization")
    if !(q isa AbstractDict)
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        print(
            io, html_table(
                [
                    ("vertices", Base.get(q, "vertices", "?")),
                    ("invalid vertices", Base.get(q, "invalid_vertices", "?")),
                    ("sum 255 fraction", number(Base.get(q, "sum_255_fraction", nothing))),
                    ("rule", Base.get(q, "rule", "?")),
                ]
            )
        )
        print(io, "<p class=\"muted\">", _esc(string(Base.get(q, "note", ""))), "</p>")
    end
    print(io, "</article>")
    return
end
