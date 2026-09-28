# `/cie`: the CIE 1931 chromaticity diagram.
#
# One WGLMakie axis carries the spectral locus with a wavelength crown, the
# sRGB chromaticity triangle, the three chromaticity layers of each pigment,
# the nominal paper colors, and the spectral/runtime mixture gamut. Clicking
# the diagram runs the inverse lookup and shows the recovered concentrations
# and residual. Everything degrades to a note when the spectral reference is
# absent.

const CIE_LAYER_LABELS = ("K/S (measured)", "surrogate (fitted)", "forward table (runtime)")

# Mixture-gamut sampling: divisions per simplex axis. 16 keeps the page fast
# while giving a smooth hull.
const CIE_MIXTURE_DIVISIONS = 16

# The displayable fill color of a grid wavelength; the locus point itself is
# the exact, unclamped chromaticity.
wavelength_color(quad::Quadrature, j::Integer) = linear_color(wavelength_linear(quad, j))

"""
    cie_app(d::ExplorerData) -> Bonito.App

The interactive CIE 1931 figure. Renders a note when the spectral reference
is unavailable.
"""
function cie_app(d::ExplorerData)
    return App() do session
        cie_dom(session, d)
    end
end

function cie_dom(session, d::ExplorerData)
    c = d.derived
    quad = d.quad
    spectral = d.spectral
    if quad === nothing || spectral === nothing || isempty(c.wavelength)
        return Bonito.DOM.div(
            Bonito.DOM.p(
                "Spectral reference not available: the input database is absent.";
                class = "muted",
            ),
        )
    end

    locus = spectral_locus(quad)
    primaries = srgb_primaries(quad)
    d65 = d65_xy(quad)
    nominal = nominal_xy(quad)
    colors = ntuple(i -> pigment_color(c, i), 4)

    # State exists before the plots so every plot can attach its visibility.
    pigment_on = [Bonito.Observable(true) for _ in 1:4]
    fitted_available = c.fitted_xy !== nothing
    layer_on = [
        Bonito.Observable(true), Bonito.Observable(fitted_available),
        Bonito.Observable(true),
    ]

    fig = Makie.Figure(; size = (880, 860))
    ax = Makie.Axis(
        fig[1, 1];
        xlabel = "x", ylabel = "y", aspect = Makie.DataAspect(),
        limits = (0.0, 0.8, 0.0, 0.9), title = "CIE 1931 chromaticity",
    )

    _cie_triangle!(ax, primaries, d65)
    _cie_locus!(ax, quad, locus)
    gamut_controls = _cie_gamut!(session, ax, spectral, d.model, quad, primaries)
    _cie_pigments!(ax, c, colors, pigment_on, layer_on)
    _cie_nominal!(ax, nominal)
    probe_controls = _cie_probe!(session, ax, d, quad, c)

    # --- controls ----------------------------------------------------------
    pigment_boxes = [Bonito.Checkbox(true) for _ in 1:4]
    for i in 1:4
        Bonito.on(session, pigment_boxes[i].value) do v
            pigment_on[i][] = v
        end
    end
    pigment_controls = Bonito.DOM.div(
        Bonito.DOM.h3("Pigments"),
        [
            Bonito.DOM.label(
                pigment_boxes[i],
                Bonito.DOM.span(; class = "swatch", style = "background:$(hex_from_linear(c.measured_rgb[i]))"),
                " ", c.codes[i],
            ) for i in 1:4
        ]...,
    )

    layer_boxes = [Bonito.Checkbox(layer_on[k][]) for k in 1:3]
    for k in 1:3
        Bonito.on(session, layer_boxes[k].value) do v
            layer_on[k][] = v
        end
    end
    layer_controls = Bonito.DOM.div(
        Bonito.DOM.h3("Layers"),
        [Bonito.DOM.label(layer_boxes[k], " ", CIE_LAYER_LABELS[k]) for k in 1:3]...,
    )

    controls = Bonito.DOM.div(
        pigment_controls, layer_controls, gamut_controls, probe_controls;
        class = "spectra-controls",
    )
    return Bonito.DOM.div(
        Bonito.DOM.div(fig; class = "figure-panel"), controls; class = "spectra-grid"
    )
end

# sRGB triangle, primary labels, and the D65 marker.
function _cie_triangle!(ax, primaries, d65)
    tri = [Makie.Point2f(primaries[i][1], primaries[i][2]) for i in 1:3]
    Makie.poly!(
        ax, tri;
        color = (:white, 0.0), strokecolor = (:black, 0.5), strokewidth = 1.4,
        linestyle = :dash,
    )
    for (label, p) in (("R", primaries[1]), ("G", primaries[2]), ("B", primaries[3]))
        Makie.text!(
            ax, [Makie.Point2f(p[1] + 0.012, p[2])]; text = [label], fontsize = 12,
            color = (:black, 0.7),
        )
    end
    Makie.scatter!(
        ax, [Makie.Point2f(d65[1], d65[2])];
        marker = :xcross, color = :black, markersize = 11,
    )
    Makie.text!(
        ax, [Makie.Point2f(d65[1] + 0.008, d65[2] + 0.012)]; text = ["D65"],
        fontsize = 11, color = (:black, 0.7),
    )
    return
end

# The locus closed by the line of purples, with a radial crown of wavelength
# ticks and a label every other sample.
function _cie_locus!(ax, quad, locus)
    lx = [p[1] for p in locus]
    ly = [p[2] for p in locus]
    push!(lx, locus[1][1])
    push!(ly, locus[1][2])
    Makie.lines!(ax, lx, ly; color = (:black, 0.8), linewidth = 1.6)

    n = length(locus)
    cx = sum(p -> p[1], locus) / n
    cy = sum(p -> p[2], locus) / n
    segments = Makie.Point2f[]
    segment_colors = Makie.RGBf[]
    label_points = Makie.Point2f[]
    labels = String[]
    for (j, p) in enumerate(locus)
        dx, dy = p[1] - cx, p[2] - cy
        nrm = hypot(dx, dy)
        nrm == 0 && continue
        ux, uy = dx / nrm, dy / nrm
        push!(segments, Makie.Point2f(p[1], p[2]))
        push!(segments, Makie.Point2f(p[1] + 0.035 * ux, p[2] + 0.035 * uy))
        col = wavelength_color(quad, j)
        push!(segment_colors, col, col)
        if isodd(j)
            push!(label_points, Makie.Point2f(p[1] + 0.048 * ux, p[2] + 0.048 * uy))
            push!(labels, string(round(Int, quad.wavelength[j])))
        end
    end
    Makie.linesegments!(ax, segments; color = segment_colors, linewidth = 6)
    Makie.text!(
        ax, label_points; text = labels, fontsize = 8, color = (:black, 0.7),
        align = (:center, :center),
    )
    return
end

# The spectral and runtime mixture gamuts, toggled independently. Returns the
# control widget.
function _cie_gamut!(session, ax, spectral, model, quad, primaries)
    tri = collect(primaries)
    spectral_hull = clip_convex(
        convex_hull(mixture_gamut(spectral, quad; d = CIE_MIXTURE_DIVISIONS)), tri
    )
    runtime_hull = clip_convex(
        convex_hull(mixture_gamut(model, quad; d = CIE_MIXTURE_DIVISIONS)), tri
    )
    spectral_vis = Bonito.Observable(true)
    runtime_vis = Bonito.Observable(true)
    spectral_box = Bonito.Checkbox(true)
    runtime_box = Bonito.Checkbox(true)
    Bonito.on(session, spectral_box.value) do v
        spectral_vis[] = v
    end
    Bonito.on(session, runtime_box.value) do v
        runtime_vis[] = v
    end
    if length(spectral_hull) >= 3
        Makie.poly!(
            ax, [Makie.Point2f(p[1], p[2]) for p in spectral_hull];
            color = (:steelblue, 0.12), strokecolor = (:steelblue, 0.8),
            strokewidth = 1.2, visible = spectral_vis,
        )
    end
    if length(runtime_hull) >= 3
        Makie.poly!(
            ax, [Makie.Point2f(p[1], p[2]) for p in runtime_hull];
            color = (:orange, 0.1), strokecolor = (:orange, 0.9),
            strokewidth = 1.2, visible = runtime_vis,
        )
    end
    return Bonito.DOM.div(
        Bonito.DOM.h3("Mixture gamut"),
        Bonito.DOM.label(spectral_box, " spectral (K/S)"),
        Bonito.DOM.label(runtime_box, " runtime forward table"),
        Bonito.DOM.p(
            "Hull of the simplex surface, clipped to the sRGB triangle."; class = "muted"
        ),
    )
end

# The three chromaticity layers of every pigment, connected so the fit and
# quantization displacement is visible. `layer 1 = measured`, `2 = fitted`,
# `3 = runtime`.
function _cie_pigments!(ax, c, colors, pigment_on, layer_on)
    for i in 1:4
        points = Tuple{Float64, Float64}[c.measured_xy[i]]
        if c.fitted_xy !== nothing
            push!(points, c.fitted_xy[i])
        end
        push!(points, c.runtime_xy[i])
        if length(points) >= 2
            Makie.lines!(
                ax, [p[1] for p in points], [p[2] for p in points];
                color = (colors[i], 0.35), linewidth = 1, visible = pigment_on[i],
            )
        end
        layer_visible(layer) = Makie.lift(
            (p, l) -> p && l, pigment_on[i], layer_on[layer]
        )
        Makie.scatter!(
            ax, [Makie.Point2f(c.measured_xy[i][1], c.measured_xy[i][2])];
            color = colors[i], marker = :circle, markersize = 12,
            strokecolor = :black, strokewidth = 0.8, visible = layer_visible(1),
        )
        if c.fitted_xy !== nothing
            Makie.scatter!(
                ax, [Makie.Point2f(c.fitted_xy[i][1], c.fitted_xy[i][2])];
                color = :transparent, marker = :circle, markersize = 12,
                strokecolor = colors[i], strokewidth = 2.0, visible = layer_visible(2),
            )
        end
        Makie.scatter!(
            ax, [Makie.Point2f(c.runtime_xy[i][1], c.runtime_xy[i][2])];
            marker = :cross, color = colors[i], markersize = 12,
            strokecolor = :black, strokewidth = 0.6, visible = layer_visible(3),
        )
    end
    return
end

# Nominal paper colors as outlined squares.
function _cie_nominal!(ax, nominal)
    for i in 1:4
        p = nominal[i]
        lin = NOMINAL_COLORS[i][2]
        Makie.scatter!(
            ax, [Makie.Point2f(p[1], p[2])];
            marker = :rect, markersize = 13, color = linear_color(lin),
            strokecolor = :black, strokewidth = 1.0,
        )
    end
    return
end

# The click handler, the probe marker, and the readout. Returns the readout
# widget.
function _cie_probe!(session, ax, d, quad, c)
    probe_obs = Bonito.Observable{Any}(nothing)
    probe_pt = Bonito.Observable(Makie.Point2f[])
    Makie.scatter!(
        ax, probe_pt;
        color = :black, marker = :circle, markersize = 11,
        strokecolor = :white, strokewidth = 1.5,
    )
    Bonito.on(session, Makie.events(ax.scene).mousebutton) do event
        if event.button == Makie.Mouse.left && event.action == Makie.Mouse.press
            p = Makie.mouseposition(ax.scene)
            result = probe_xy(d.model, quad, p[1], p[2]; derived = c)
            probe_obs[] = result
            probe_pt[] = result === nothing ? Makie.Point2f[] :
                [Makie.Point2f(result.xy[1], result.xy[2])]
        end
    end
    readout = Bonito.Observable{Any}(_probe_table(nothing))
    Bonito.on(session, probe_obs) do p
        readout[] = _probe_table(p)
    end
    return Bonito.DOM.div(
        Bonito.DOM.h3("Click to probe"),
        Bonito.DOM.p(
            "Click the diagram to run the inverse lookup at that chromaticity.";
            class = "muted",
        ),
        Bonito.DOM.div(readout; class = "readout-wrap"),
    )
end

function _probe_table(::Nothing)
    return Bonito.DOM.p("No probe yet."; class = "muted")
end

function _probe_table(p::CieProbe)
    swatch(hex) = Bonito.DOM.span(; class = "swatch", style = "background:$hex")
    row(label, value) = Bonito.DOM.tr(Bonito.DOM.th(label), Bonito.DOM.td(value))
    return Bonito.DOM.table(
        Bonito.DOM.caption("Click probe"),
        Bonito.DOM.tbody(
            row("target xy", @sprintf("%.4f, %.4f", p.xy[1], p.xy[2])),
            row("in sRGB gamut", p.in_gamut ? "yes" : "no (clamped)"),
            row("target color", Bonito.DOM.span(swatch(p.target_hex), " ", p.target_hex)),
            row(
                "concentrations",
                @sprintf(
                    "%.3f, %.3f, %.3f, %.3f", p.concentrations[1],
                    p.concentrations[2], p.concentrations[3], p.concentrations[4]
                ),
            ),
            row(
                "reconstructed",
                Bonito.DOM.span(swatch(p.reconstructed_hex), " ", p.reconstructed_hex),
            ),
            row(
                "residual",
                @sprintf(
                    "%.4f, %.4f, %.4f", p.residual[1], p.residual[2], p.residual[3]
                ),
            ),
            row(
                "nearest pigment",
                p.nearest_pigment === nothing ? "—" :
                    @sprintf("%s (%.4f)", p.nearest_pigment, p.nearest_distance),
            ),
        );
        class = "readout",
    )
end

"""
    cie_caption(d::ExplorerData) -> String

The static caption under the CIE figure: the sRGB primaries and D65, the
three chromaticity layers of each pigment, and the nominal paper colors.
Empty when the spectral reference is unavailable.
"""
function cie_caption(d::ExplorerData)
    quad = d.quad
    c = d.derived
    (quad === nothing || isempty(c.wavelength)) && return ""
    primaries = srgb_primaries(quad)
    d65 = d65_xy(quad)
    io = IOBuffer()
    print(io, "<section class=\"cards\">")

    print(io, "<article class=\"card\"><h2>Reference points</h2>")
    print(
        io, html_table(
            [
                ("R primary", @sprintf("%.4f, %.4f", primaries[1][1], primaries[1][2])),
                ("G primary", @sprintf("%.4f, %.4f", primaries[2][1], primaries[2][2])),
                ("B primary", @sprintf("%.4f, %.4f", primaries[3][1], primaries[3][2])),
                ("D65 white", @sprintf("%.4f, %.4f", d65[1], d65[2])),
            ]
        )
    )
    print(io, "</article>")

    print(io, "<article class=\"card\"><h2>Pigment chromaticity layers</h2>")
    rows = []
    for p in pigment_chromaticities(c)
        fitted = p.fitted === nothing ? "—" :
            @sprintf("%.4f, %.4f", p.fitted[1], p.fitted[2])
        push!(
            rows, (
                p.code,
                @sprintf("%.4f, %.4f", p.measured[1], p.measured[2]), fitted,
                @sprintf("%.4f, %.4f", p.runtime[1], p.runtime[2]),
            )
        )
    end
    print(io, html_table(rows; headers = ["Pigment", "K/S measured", "fitted", "runtime"]))
    print(io, "</article></section>")

    print(io, "<section class=\"cards\"><article class=\"card\"><h2>Nominal paper colors</h2>")
    nominal = nominal_xy(quad)
    rows = []
    for i in 1:4
        lin = NOMINAL_COLORS[i][2]
        push!(
            rows, (
                NOMINAL_COLORS[i][1],
                @sprintf("%.2f, %.2f, %.2f", lin[1], lin[2], lin[3]),
                hex_from_linear(lin),
                @sprintf("%.4f, %.4f", nominal[i][1], nominal[i][2]),
            )
        )
    end
    print(
        io, html_table(
            rows; headers = ["Name", "nominal linear RGB", "encoded", "chromaticity"]
        )
    )
    print(io, "</article></section>")
    return String(take!(io))
end
