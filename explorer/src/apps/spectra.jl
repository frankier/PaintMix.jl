# `/spectra`: measured, fitted, and derived curves over the wavelength grid.
#
# One WGLMakie figure with one axis per quantity group (K, S, K/S, and R∞ with
# R′). Pigment toggles, quantity toggles, a measured/fitted toggle, a log-scale
# toggle, and a wavelength crosshair drive the figure; a server-rendered readout
# table follows the crosshair and carries the source spreadsheet cell of every
# sample. All observables are created inside the `App` closure, so each page
# load gets its own state.

const SPECTRA_DEFAULT_NM = 550.0

# Display color of a pigment: its measured pure color, clamped (PY74's blue
# channel is negative in linear light) and gamma encoded for an sRGB canvas.
function pigment_color(d::DerivedCurves, i::Integer)
    c = d.measured_rgb[i]
    lin = SVector{3, Float32}(
        clamp(Float32(c[1]), 0.0f0, 1.0f0),
        clamp(Float32(c[2]), 0.0f0, 1.0f0),
        clamp(Float32(c[3]), 0.0f0, 1.0f0),
    )
    s = PaintMix.srgb_from_linear(lin)
    return Makie.RGBf(clamp(s[1], 0.0f0, 1.0f0), clamp(s[2], 0.0f0, 1.0f0), clamp(s[3], 0.0f0, 1.0f0))
end

function pigment_hex(d::DerivedCurves, i::Integer)
    b = encoded(d.measured_rgb[i])
    return @sprintf("#%02x%02x%02x", b[1], b[2], b[3])
end

"""
    spectra_app(d::ExplorerData) -> Bonito.App

The interactive spectra figure. When the spectral reference is unavailable the
app renders a note instead of a figure, so the route never fails.
"""
function spectra_app(d::ExplorerData)
    return App() do session
        spectra_dom(session, d)
    end
end

function spectra_dom(session, d::ExplorerData)
    c = d.derived
    if isempty(c.wavelength)
        return Bonito.DOM.div(
            Bonito.DOM.p(
                "Spectral reference not available: the input database is absent.";
                class = "muted",
            ),
        )
    end

    wavelength = c.wavelength
    default = SPECTRA_DEFAULT_NM in wavelength ? SPECTRA_DEFAULT_NM :
        wavelength[cld(length(wavelength), 2)]
    λ = Bonito.Observable(Float64(default))
    pigment_on = [Bonito.Observable(true) for _ in 1:4]
    quantity_on = [Bonito.Observable(true) for _ in SPECTRA_QUANTITY_KEYS]
    fitted_available = c.fitted_K !== nothing
    fitted_on = Bonito.Observable(fitted_available)
    log_on = Bonito.Observable(true)
    colors = ntuple(i -> pigment_color(c, i), 4)
    keys = (:K, :S, :KS, :Rinf)

    fig = Makie.Figure(; size = (900, 880))
    axes = (
        _spectra_axis(fig, 1, "K", :K; last = false),
        _spectra_axis(fig, 2, "S", :S; last = false),
        _spectra_axis(fig, 3, "K/S", :KS; last = false),
        _spectra_axis(fig, 4, "R∞, R′", :Rinf; last = true),
    )
    for (ax, q) in zip(axes[1:3], keys[1:3])
        for i in 1:4
            qi = _quantity_position(q)
            _plot_curve!(
                ax, c, Val(q), i, colors[i], alpha = 1.0,
                pigment_on = pigment_on[i], quantity_on = quantity_on[qi],
                fitted_on = fitted_on,
            )
        end
    end
    for (q, alpha) in ((:Rinf, 1.0), (:Rprime, 0.45))
        qi = _quantity_position(q)
        for i in 1:4
            _plot_curve!(
                axes[4], c, Val(q), i, colors[i], alpha = alpha,
                pigment_on = pigment_on[i], quantity_on = quantity_on[qi],
                fitted_on = fitted_on,
            )
        end
    end
    for ax in axes
        Makie.vlines!(ax, λ; color = (:black, 0.35), linestyle = :dashdot)
    end

    _apply_spectra_scales!(axes, keys, log_on)
    Bonito.on(session, log_on) do _
        _apply_spectra_scales!(axes, keys, log_on)
    end

    # --- controls ----------------------------------------------------------
    slider = Bonito.Slider(wavelength; value = default)
    Bonito.on(session, slider.value) do v
        λ[] = Float64(v)
    end
    wavelength_control = Bonito.DOM.div(
        Bonito.DOM.h3("Wavelength"),
        slider,
        Bonito.DOM.div(Bonito.DOM.span(slider.value), " nm"; class = "readout-lambda"),
    )

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
                Bonito.DOM.span(; class = "swatch", style = "background:$(pigment_hex(c, i))"),
                " ",
                c.codes[i],
            ) for i in 1:4
        ]...,
    )

    quantity_boxes = [Bonito.Checkbox(true) for _ in SPECTRA_QUANTITY_KEYS]
    for (qi, _) in enumerate(SPECTRA_QUANTITY_KEYS)
        Bonito.on(session, quantity_boxes[qi].value) do v
            quantity_on[qi][] = v
        end
    end
    quantity_controls = Bonito.DOM.div(
        Bonito.DOM.h3("Quantities"),
        [
            Bonito.DOM.label(quantity_boxes[qi], " ", quantity_label(q)) for
                (qi, q) in enumerate(SPECTRA_QUANTITY_KEYS)
        ]...,
    )

    log_box = Bonito.Checkbox(true)
    Bonito.on(session, log_box.value) do v
        log_on[] = v
    end
    option_children = Any[
        Bonito.DOM.h3("Options"), Bonito.DOM.label(log_box, " log K, S, K/S"),
    ]
    if fitted_available
        fitted_box = Bonito.Checkbox(true)
        Bonito.on(session, fitted_box.value) do v
            fitted_on[] = v
        end
        push!(option_children, Bonito.DOM.label(fitted_box, " measured + fitted"))
    end
    options = Bonito.DOM.div(option_children...)

    readout = Bonito.Observable{Any}(nothing)
    function refresh_readout()
        pigments = ntuple(i -> pigment_on[i][], 4)
        quantities = ntuple(i -> quantity_on[i][], length(quantity_on))
        return readout[] = _readout_table(
            spectra_readout(
                c, λ[]; pigments = pigments, quantities = quantities,
                fitted = fitted_on[]
            ),
        )
    end
    Bonito.onany(
        (_...) -> refresh_readout(),
        session, λ, pigment_on..., quantity_on..., fitted_on,
    )
    refresh_readout()

    controls = Bonito.DOM.div(
        wavelength_control, pigment_controls, quantity_controls, options,
        Bonito.DOM.div(readout; class = "readout-wrap");
        class = "spectra-controls",
    )
    return Bonito.DOM.div(
        Bonito.DOM.div(fig; class = "figure-panel"), controls;
        class = "spectra-grid",
    )
end

function _spectra_axis(fig, row::Integer, ylabel::AbstractString, q::Symbol; last::Bool)
    scale = quantity_log(q) ? log10 : identity
    ax = Makie.Axis(
        fig[row, 1]; ylabel = ylabel, xlabel = last ? "λ (nm)" : "", yscale = scale
    )
    last || Makie.hidexdecorations!(ax)
    return ax
end

function _apply_spectra_scales!(axes, keys, log_on)
    for (ax, q) in zip(axes, keys)
        ax.yscale = (log_on[] && quantity_log(q)) ? log10 : identity
    end
    return
end

_quantity_position(q::Symbol) = findfirst(==(q), SPECTRA_QUANTITY_KEYS)

# One measured line, and a dashed fitted line when the surrogate is present.
# Visibility combines the pigment toggle, the quantity toggle, and the
# measured/fitted toggle.
function _plot_curve!(
        ax, c::DerivedCurves, q::Val, i::Integer, color; alpha::Float64,
        pigment_on, quantity_on, fitted_on,
    )
    width = 1.6
    rgba = Makie.RGBAf(color.r, color.g, color.b, alpha)
    visible = Makie.lift((p, qo) -> p && qo, pigment_on, quantity_on)
    Makie.lines!(
        ax, c.wavelength, quantity_matrix(c, q)[i, :];
        color = rgba, linewidth = width, visible = visible,
    )
    fm = fitted_matrix(c, q)
    if fm !== nothing
        fvisible = Makie.lift(
            (p, qo, f) -> p && qo && f, pigment_on, quantity_on, fitted_on
        )
        Makie.lines!(
            ax, c.wavelength, fm[i, :];
            color = rgba, linewidth = width, linestyle = :dash, visible = fvisible,
        )
    end
    return
end

function _format_value(x::Missing)
    return "—"
end

function _format_value(x::Real)
    return @sprintf("%.4g", x)
end

function _readout_table(r::SpectraReadout)
    header = Bonito.DOM.tr(
        Bonito.DOM.th("Pigment"),
        Bonito.DOM.th("Layer"),
        (Bonito.DOM.th(quantity_label(q)) for q in r.quantities)...,
        Bonito.DOM.th("source cells"),
    )
    body = [
        Bonito.DOM.tr(
            Bonito.DOM.td(row.pigment),
            Bonito.DOM.td(string(row.layer)),
            (Bonito.DOM.td(_format_value(v)) for v in row.values)...,
            Bonito.DOM.td(row.cells; class = "cells"),
        ) for row in r.rows
    ]
    caption = Bonito.DOM.caption("λ = " * @sprintf("%.0f", r.wavelength) * " nm")
    return Bonito.DOM.table(
        caption, Bonito.DOM.thead(header), Bonito.DOM.tbody(body...); class = "readout"
    )
end

"""
    spectra_caption(d::ExplorerData) -> String

The static caption under the figure: the wavelength grid, the Saunderson
constants, and the `kins` contradiction recorded in `inputs/README.md`. Empty
when the spectral reference is unavailable.
"""
function spectra_caption(d::ExplorerData)
    c = d.derived
    isempty(c.wavelength) && return ""
    w = c.wavelength
    step = length(w) > 1 ? w[2] - w[1] : 0.0
    rows = [
        ("Grid", @sprintf("%.0f–%.0f nm, step %.0f nm, %d samples", first(w), last(w), step, length(w))),
        ("Integral", "trapezoidal on the configuration grid"),
        ("Lines", "solid = measured, dashed = fitted; on the R∞/R′ axis the lighter line is R′"),
    ]
    note = ""
    if d.spectra !== nothing
        push!(rows, ("Saunderson k1", string(d.spectra.k1)))
        push!(rows, ("Saunderson k2", string(d.spectra.k2)))
        note = string(Base.get(d.spectra.provenance, "kins_note", ""))
    end
    io = IOBuffer()
    print(io, "<section class=\"card caption\"><h2>Spectra caption</h2>")
    print(io, html_table(rows))
    if !isempty(note)
        print(
            io, "<p class=\"warn\" title=\"", _esc(note), "\">kins contradiction: ",
            _esc(note), "</p>",
        )
    end
    print(io, "</section>")
    return String(take!(io))
end
