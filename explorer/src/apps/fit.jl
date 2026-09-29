# `/fit`: the surrogate fit. Measured vs fitted spectra and chromaticity
# displacement, plus the continuation history and diagnostics from the
# sidecar.
#
# The spectra comparison needs the measured and fitted curves; the history
# needs the sidecar. The page shows whichever is available and a note for the
# rest, so it never fails on a partial checkout.

"""
    fit_spectra_figure(c) -> Makie.Figure

Measured (solid) against fitted (dashed) `K`, `S`, `R∞`, and `R′` for every
pigment. Call only when `c.fitted_K` is present.
"""
function fit_spectra_figure(c::DerivedCurves)
    fig = Makie.Figure(; size = (900, 700))
    specs = ((:K, "K"), (:S, "S"), (:Rinf, "R∞"), (:Rprime, "R′"))
    for (idx, (q, label)) in enumerate(specs)
        row = cld(idx, 2)
        col = mod1(idx, 2)
        ax = Makie.Axis(
            fig[row, col]; ylabel = label, xlabel = row == 2 ? "λ (nm)" : "",
            yscale = quantity_log(q) ? log10 : identity,
        )
        for i in 1:4
            color = linear_color(c.measured_rgb[i])
            Makie.lines!(
                ax, c.wavelength, quantity_matrix(c, q)[i, :];
                color = color, linewidth = 1.4,
            )
            fm = fitted_matrix(c, q)
            fm === nothing && continue
            Makie.lines!(
                ax, c.wavelength, fm[i, :];
                color = color, linewidth = 1.4, linestyle = :dash,
            )
        end
    end
    return fig
end

"""
    fit_history_figure(h, tolerance) -> Makie.Figure

`Epush`, `Epull`, the objective, and elapsed seconds against the continuation
step. The `Epush` axis marks the tolerance and the first step that meets it.
"""
function fit_history_figure(h::FitHistory, tolerance::Real)
    fig = Makie.Figure(; size = (1000, 660))
    ax_push = Makie.Axis(
        fig[1, 1]; xlabel = "step", ylabel = "Epush", yscale = log10,
        title = "Cube push",
    )
    Makie.lines!(ax_push, h.steps, h.epush; color = :firebrick)
    Makie.scatter!(ax_push, h.steps, h.epush; color = :firebrick, markersize = 7)
    Makie.hlines!(ax_push, [Float64(tolerance)]; color = :black, linestyle = :dash)
    met = push_met_step(h, tolerance)
    if met !== nothing
        Makie.scatter!(
            ax_push, [h.steps[met]], [h.epush[met]];
            color = :green, markersize = 13, strokecolor = :black, strokewidth = 0.8,
        )
    end
    Makie.vlines!(
        ax_push, [last(h.steps)]; color = (:gray, 0.7), linestyle = :dot
    )

    ax_pull = Makie.Axis(
        fig[1, 2]; xlabel = "step", ylabel = "Epull", yscale = log10,
        title = "Perceptual pull",
    )
    Makie.lines!(ax_pull, h.steps, h.epull; color = :steelblue)
    Makie.scatter!(ax_pull, h.steps, h.epull; color = :steelblue, markersize = 7)

    ax_obj = Makie.Axis(
        fig[2, 1]; xlabel = "step", ylabel = "objective", yscale = log10,
        title = "Objective",
    )
    Makie.lines!(ax_obj, h.steps, h.objective; color = :black)
    Makie.scatter!(ax_obj, h.steps, h.objective; color = :black, markersize = 7)

    ax_time = Makie.Axis(
        fig[2, 2]; xlabel = "step", ylabel = "seconds", title = "Elapsed"
    )
    Makie.lines!(ax_time, h.steps, h.seconds; color = :darkorange)
    Makie.scatter!(ax_time, h.steps, h.seconds; color = :darkorange, markersize = 7)
    return fig
end

"""
    fit_app(d::ExplorerData) -> Bonito.App

The fit page's figures: measured vs fitted spectra, and the continuation
history.
"""
function fit_app(d::ExplorerData)
    return App() do session
        fit_dom(session, d)
    end
end

function fit_dom(session, d::ExplorerData)
    c = d.derived
    h = fit_history(d)
    has_spectra = !isempty(c.wavelength) && c.fitted_K !== nothing
    if !has_spectra && h === nothing
        return Bonito.DOM.div(
            Bonito.DOM.p(
                "Fit data not available: the spectral reference and the " *
                    "provenance sidecar are both absent.";
                class = "muted",
            ),
        )
    end

    panels = Any[]
    if has_spectra
        push!(
            panels,
            Bonito.DOM.div(fit_spectra_figure(c); class = "figure-panel"),
        )
    end
    if h !== nothing
        push!(
            panels,
            Bonito.DOM.div(fit_history_figure(h, push_tolerance(d)); class = "figure-panel"),
        )
    end
    return Bonito.DOM.div(panels...; class = "fit-stack")
end

"""
    fit_caption(d::ExplorerData) -> String

The server-rendered fit cards: history summary, diagnostics, and chromaticity
displacement.
"""
function fit_caption(d::ExplorerData)
    io = IOBuffer()
    print(io, "<section class=\"cards\">")
    _fit_history_card(io, d)
    _fit_diagnostics_card(io, d)
    print(io, "</section>")
    _fit_chromaticity_card(io, d)
    return String(take!(io))
end

function _fit_history_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Continuation history</h2>")
    h = fit_history(d)
    if h === nothing
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        tol = push_tolerance(d)
        met = push_met_step(h, tol)
        print(
            io, html_table(
                [
                    ("steps", length(h.steps)),
                    ("push tolerance", number(tol)),
                    ("first step meeting tolerance", met === nothing ? "never" : met),
                    ("stopped on step", last(h.steps)),
                    ("converged steps", count(h.converged)),
                ]
            )
        )
        rows = [
            (
                h.steps[i], number(h.alpha[i]), number(h.epush[i]),
                number(h.epull[i]), number(h.objective[i]), number(h.seconds[i]),
                h.iterations[i], h.converged[i] ? "yes" : "no",
            ) for i in eachindex(h.steps)
        ]
        print(
            io, html_table(
                rows; headers = [
                    "step", "alpha", "Epush", "Epull", "objective", "seconds",
                    "iterations", "converged",
                ]
            )
        )
    end
    print(io, "</article>")
    return
end

function _fit_diagnostics_card(io, d::ExplorerData)
    print(io, "<article class=\"card\"><h2>Fit diagnostics</h2>")
    diagnostics = fit_diagnostics(d)
    if diagnostics === nothing
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    else
        keys = (
            "Epush_fit", "Epush_dense", "Epull_fit", "Epull_dense",
            "max_cube_violation", "max_oklab_deviation", "dense_divisions",
        )
        rows = [(k, number(Base.get(diagnostics, k, nothing))) for k in keys]
        print(io, html_table(rows))
    end
    print(io, "</article>")
    return
end

function _fit_chromaticity_card(io, d::ExplorerData)
    print(io, "<section class=\"card\"><h2>Chromaticity displacement</h2>")
    c = d.derived
    if c.fitted_xy === nothing
        print(io, "<p class=\"muted\">Not available: no fitted surrogate.</p>")
    else
        rows = []
        for i in 1:4
            push!(
                rows, (
                    c.codes[i],
                    @sprintf("%.4f, %.4f", c.measured_xy[i][1], c.measured_xy[i][2]),
                    @sprintf("%.4f, %.4f", c.fitted_xy[i][1], c.fitted_xy[i][2]),
                    @sprintf("%.4f, %.4f", c.runtime_xy[i][1], c.runtime_xy[i][2]),
                    number(
                        sqrt(
                            (c.fitted_xy[i][1] - c.measured_xy[i][1])^2 +
                                (c.fitted_xy[i][2] - c.measured_xy[i][2])^2
                        ),
                    ),
                )
            )
        end
        print(
            io, html_table(
                rows; headers = [
                    "pigment", "measured xy", "fitted xy", "runtime xy",
                    "fitted − measured",
                ]
            )
        )
        print(
            io,
            "<p class=\"muted\">The fitted surrogate is pulled toward the ",
            "measured pigments but pushed inside the sRGB cube, so its ",
            "chromaticity differs from the raw K/S point. The runtime column ",
            "adds the 8-bit forward-table quantization.</p>",
        )
    end
    print(io, "</section>")
    return
end
