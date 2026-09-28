# `/palette`: the sRGB cube, mixing curves, the pairwise ramp matrix, and the
# mixer ramp.
#
# The cube and the mixer need only the runtime payload, so the page works when
# the input database is absent. The mixing curve for the selected pigment pair
# is drawn inside the cube as a 3D path, with the naive linear-RGB line dashed
# for contrast. The ramp matrix and the nominal table are server-rendered HTML.

# The eight cube corners, indexed by bit 0 = R, bit 1 = G, bit 2 = B.
const _CUBE_CORNERS = ntuple(
    k -> Makie.Point3f(
        Float32((k - 1) & 1), Float32(((k - 1) >> 1) & 1), Float32(((k - 1) >> 2) & 1)
    ), 8
)

# The twelve cube edges: corner pairs that differ in exactly one bit.
function _cube_edges!(ax)
    segments = Makie.Point3f[]
    for a in 1:8, b in (a + 1):8
        count_ones((a - 1) ⊻ (b - 1)) == 1 || continue
        push!(segments, _CUBE_CORNERS[a], _CUBE_CORNERS[b])
    end
    Makie.linesegments!(ax, segments; color = (:black, 0.45), linewidth = 1.2)
    return
end

function _point3f(c)
    return Makie.Point3f(Float32(c[1]), Float32(c[2]), Float32(c[3]))
end

function _path_points(matrix::Matrix{Float64})
    return [
        Makie.Point3f(
            Float32(matrix[1, j]), Float32(matrix[2, j]), Float32(matrix[3, j])
        ) for j in axes(matrix, 2)
    ]
end

"""
    palette_cube_app(d::ExplorerData) -> Bonito.App

The sRGB cube: edges, pure-pigment runtime vertices, nominal paper colors,
the D65 white point, and the mixing path of a selected pigment pair.
"""
function palette_cube_app(d::ExplorerData)
    return App() do session
        palette_cube_dom(session, d)
    end
end

function palette_cube_dom(session, d::ExplorerData)
    model = d.model
    colors = pigment_rgbs(model)
    nominal = ntuple(nominal_rgb, 4)
    labels = [pigment_label(d, i) for i in 1:4]

    fig = Makie.Figure(; size = (620, 620))
    ax = Makie.Axis3(
        fig[1, 1]; xlabel = "R", ylabel = "G", zlabel = "B", title = "sRGB cube",
        limits = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
    )
    _cube_edges!(ax)
    for i in 1:4
        Makie.scatter!(
            ax, [_point3f(colors[i])];
            color = linear_color(colors[i]), marker = :circle, markersize = 15,
            strokecolor = :black, strokewidth = 0.8,
        )
        Makie.scatter!(
            ax, [_point3f(nominal[i])];
            color = linear_color(nominal[i]), marker = :rect, markersize = 12,
            strokecolor = :black, strokewidth = 1.0,
        )
    end
    Makie.scatter!(
        ax, [Makie.Point3f(1, 1, 1)]; marker = :xcross, color = :black, markersize = 11
    )

    pairs = [(i, j) for i in 1:4 for j in (i + 1):4]
    pair_labels = ["$(labels[i]) → $(labels[j])" for (i, j) in pairs]
    dropdown = Bonito.Dropdown(pair_labels; index = 1)
    path = Makie.lift(dropdown.value) do label
        idx = findfirst(==(label), pair_labels)
        i, j = pairs[idx]
        curve = mix_curve(model, colors[i], colors[j]; n = 64)
        return (; paint = _path_points(curve.paint), naive = _path_points(curve.naive))
    end
    Makie.lines!(
        ax, Makie.lift(p -> p.naive, path);
        color = (:gray, 0.9), linewidth = 2, linestyle = :dash,
    )
    Makie.lines!(ax, Makie.lift(p -> p.paint, path); color = :orange, linewidth = 3)

    controls = Bonito.DOM.div(
        Bonito.DOM.h3("Mixing path"),
        dropdown,
        Bonito.DOM.p(
            "Solid: paint mixing through the model. Dashed: naive linear-RGB interpolation.";
            class = "muted",
        ),
        Bonito.DOM.h3("Legend"),
        Bonito.DOM.p("Circles: pure-pigment runtime vertices.", class = "muted"),
        Bonito.DOM.p("Squares: nominal paper colors.", class = "muted"),
        Bonito.DOM.p("Cross: D65 white.", class = "muted"),
    )
    return Bonito.DOM.div(
        Bonito.DOM.div(fig; class = "figure-panel"), controls; class = "spectra-grid"
    )
end

"""
    palette_mixer_app(d::ExplorerData) -> Bonito.App

The mixer ramp: one weight slider per pigment, a paper-to-mixture fraction,
a large swatch, the encoded hex, and the recovered concentrations.
"""
function palette_mixer_app(d::ExplorerData)
    return App() do session
        palette_mixer_dom(session, d)
    end
end

function palette_mixer_dom(session, d::ExplorerData)
    model = d.model
    labels = [pigment_label(d, i) for i in 1:4]
    colors = pigment_rgbs(model)

    weight_sliders = [Bonito.Slider(0.0:0.01:1.0; value = 1.0) for _ in 1:4]
    t_slider = Bonito.Slider(0.0:0.01:1.0; value = 1.0)
    result = Makie.lift(
        (w1, w2, w3, w4, t) -> mixer_result(model, (w1, w2, w3, w4); t = t),
        weight_sliders[1].value, weight_sliders[2].value, weight_sliders[3].value,
        weight_sliders[4].value, t_slider.value,
    )
    ramp = Makie.lift(
        (w1, w2, w3, w4) -> ramp_gradient(mixer_ramp(model, (w1, w2, w3, w4))),
        weight_sliders[1].value, weight_sliders[2].value, weight_sliders[3].value,
        weight_sliders[4].value,
    )

    swatch = Bonito.DOM.div(
        ; class = "mixer-swatch", style = Makie.lift(r -> "background:$(r.hex)", result)
    )
    readout = Bonito.DOM.table(
        Bonito.DOM.tbody(
            Bonito.DOM.tr(
                Bonito.DOM.th("encoded"),
                Bonito.DOM.td(Bonito.DOM.code(Makie.lift(r -> r.hex, result))),
            ),
            Bonito.DOM.tr(
                Bonito.DOM.th("concentrations"),
                Bonito.DOM.td(
                    Makie.lift(
                        r -> @sprintf(
                            "%.3f, %.3f, %.3f, %.3f", r.concentrations[1],
                            r.concentrations[2], r.concentrations[3], r.concentrations[4]
                        ),
                        result,
                    ),
                ),
            ),
        );
        class = "readout",
    )

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
    controls = Bonito.DOM.div(
        Bonito.DOM.h3("Pigment weights"), weight_controls...,
        Bonito.DOM.h3("Paper to mixture"), t_slider,
        Bonito.DOM.div(swatch, class = "mixer-row"),
        Bonito.DOM.div(; class = "ramp mixer-ramp", style = ramp),
        Bonito.DOM.div(readout; class = "readout-wrap"),
    )
    return Bonito.DOM.div(Bonito.DOM.div(controls; class = "spectra-controls"))
end

"""
    palette_ramp_matrix(d::ExplorerData) -> String

The 4×4 pairwise mixture ramp matrix, server-rendered. Each strip is the
sampled paint mixture, so its shape is the model's, not a straight blend.
"""
function palette_ramp_matrix(d::ExplorerData)
    ramps = pairwise_ramps(d.model)
    labels = [pigment_label(d, i) for i in 1:4]
    io = IOBuffer()
    print(io, "<section class=\"card\"><h2>Pairwise mixture ramps</h2>")
    print(io, "<table class=\"ramp-matrix\"><thead><tr><th></th>")
    for label in labels
        print(io, "<th>", _esc(label), "</th>")
    end
    print(io, "</tr></thead><tbody>")
    for i in 1:4
        print(io, "<tr><th>", _esc(labels[i]), "</th>")
        for j in 1:4
            print(
                io, "<td><span class=\"ramp\" style=\"background:",
                ramp_gradient(ramps[i, j]), "\"></span></td>"
            )
        end
        print(io, "</tr>")
    end
    print(io, "</tbody></table>")
    print(
        io,
        "<p class=\"muted\">Each strip runs from the row pigment to the column ",
        "pigment through <code>mix</code>, not a linear RGB blend.</p></section>",
    )
    return String(take!(io))
end

"""
    palette_caption(d::ExplorerData) -> String

The nominal table: every slot's design color against its runtime vertex, with
the per-channel residual and the runtime chromaticity. Works without the
spectral reference; the chromaticity column then shows a dash.
"""
function palette_caption(d::ExplorerData)
    io = IOBuffer()
    print(io, "<section class=\"cards\"><article class=\"card\">")
    print(io, "<h2>Nominal paper colors vs runtime vertices</h2>")
    rows = []
    for r in nominal_rows(d)
        xy = r.xy === nothing ? "—" : @sprintf("%.4f, %.4f", r.xy[1], r.xy[2])
        push!(
            rows, (
                string(r.slot), r.name, r.code,
                @sprintf("%.3f, %.3f, %.3f", r.nominal[1], r.nominal[2], r.nominal[3]),
                r.nominal_hex,
                @sprintf("%.3f, %.3f, %.3f", r.runtime[1], r.runtime[2], r.runtime[3]),
                r.runtime_hex,
                @sprintf(
                    "%+.4f, %+.4f, %+.4f", r.residual[1], r.residual[2], r.residual[3]
                ),
                xy,
            )
        )
    end
    print(
        io, html_table(
            rows; headers = [
                "slot", "nominal name", "code", "nominal linear", "nominal hex",
                "runtime linear", "runtime hex", "residual", "runtime xy",
            ]
        )
    )
    print(
        io,
        "<p class=\"muted\">Residual is the runtime vertex minus the nominal ",
        "paper color. The forward table is 8-bit, so the runtime hex clips ",
        "PY74's negative blue channel.</p>",
    )
    print(
        io,
        "<p class=\"muted\">The workbook's other twenty pigments are not shown: ",
        "the C.I.-name-to-column mapping is unverified (PLAN.md, problem 6).</p>",
    )
    print(io, "</article>")

    if d.quad !== nothing
        primaries = srgb_primaries(d.quad)
        d65 = d65_xy(d.quad)
        print(io, "<article class=\"card\"><h2>sRGB reference</h2>")
        print(
            io, html_table(
                [
                    (
                        "R primary",
                        @sprintf("%.4f, %.4f", primaries[1][1], primaries[1][2]),
                    ),
                    (
                        "G primary",
                        @sprintf("%.4f, %.4f", primaries[2][1], primaries[2][2]),
                    ),
                    (
                        "B primary",
                        @sprintf("%.4f, %.4f", primaries[3][1], primaries[3][2]),
                    ),
                    ("D65 white", @sprintf("%.4f, %.4f", d65[1], d65[2])),
                ]
            )
        )
        print(io, "</article>")
    end
    print(io, "</section>")
    return String(take!(io))
end
