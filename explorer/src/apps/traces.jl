# The same snapshot drives the physical spectra, both gamut views and the
# loss readout. Fixed axes and a ghost of the initial gamut anchor movement.
function fit_trace_dom(session, trace)
    dom = Bonito.DOM
    trace === nothing && return dom.p("Continuation trace unavailable for this payload."; class = "muted")
    rows = trace["snapshots"]
    wavelength = trace["wavelength_nm"]
    samples = trace["boundary"]["concentrations"]
    points = SVector{4, Float64}.(samples[1:max(1, cld(length(samples), 640)):end])
    models = [PaintMixPrecompute.replay_fit_snapshot(trace, i) for i in eachindex(rows)]
    gamuts = [[mix_rgb(m, c) for c in points] for m in models]
    xy = [Makie.Point2f[Makie.Point2f(xy_of_linear(m.quad, p)...) for p in rgb] for (m, rgb) in zip(models, gamuts)]
    limits = ntuple(3) do i
        lo, hi = extrema(p[i] for rgb in gamuts for p in rgb)
        (min(lo, 0.0) - 0.03, max(hi, 1.0) + 0.03)
    end

    slider = Bonito.Slider(0:(length(rows) - 1); value = 0)
    i = Makie.lift(step -> step + 1, slider.value)
    fig = Makie.Figure(; size = (1000, 700))
    pigment_colors = (:steelblue, :firebrick, :goldenrod, :gray38)
    for (col, field, label) in ((1, :K, "Absorption K"), (2, :S, "Scattering S"))
        ax = Makie.Axis(fig[1, col]; title = label, xlabel = "Wavelength (nm)", yscale = log10)
        lo, hi = extrema(v for m in models for v in getfield(m, field))
        Makie.ylims!(ax, lo / 1.2, hi * 1.2)
        for p in 1:4
            y = Makie.lift(j -> getfield(models[j], field)[:, p], i)
            Makie.lines!(ax, wavelength, y; color = pigment_colors[p], label = trace["pigments"][p])
        end
        col == 1 && Makie.axislegend(ax; position = :rt)
    end

    cube = Makie.Axis3(
        fig[2, 1]; title = "Sampled linear RGB gamut",
        xlabel = "R", ylabel = "G", zlabel = "B", limits = (limits[1]..., limits[2]..., limits[3]...)
    )
    _cube_edges!(cube)
    Makie.scatter!(cube, _point3f.(first(gamuts)); color = (:gray, 0.35), markersize = 3)
    Makie.scatter!(cube, Makie.lift(j -> _point3f.(gamuts[j]), i); color = :steelblue, markersize = 4)

    cie = Makie.Axis(fig[2, 2]; title = "Sampled CIE 1931 gamut", xlabel = "x", ylabel = "y", aspect = Makie.DataAspect())
    quad = first(models).quad
    _cie_triangle!(cie, srgb_primaries(quad), d65_xy(quad))
    Makie.scatter!(cie, first(xy); color = (:gray, 0.35), markersize = 3)
    Makie.scatter!(cie, Makie.lift(j -> xy[j], i); color = :steelblue, markersize = 4)
    xy_bounds = filter(p -> all(isfinite, p), reduce(vcat, xy))
    for (dim, set_limits, nominal_max) in ((1, Makie.xlims!, 0.75), (2, Makie.ylims!, 0.85))
        lo, hi = extrema(p[dim] for p in xy_bounds)
        set_limits(cie, min(lo, 0.0) - 0.02, max(hi, nominal_max) + 0.02)
    end

    summary = Makie.lift(j -> fit_step_summary(rows[j], length(rows) - 1), i)
    return dom.section(
        dom.h2("Continuation replay"),
        dom.label("Continuation step ", slider),
        dom.p(summary; role = "status"),
        dom.p("Gray points show the initial gamut. The cube violation and Oklab drift describe the held-out samples."),
        dom.div(fig; class = "figure-panel"); id = "fit-trace"
    )
end

function fit_step_summary(row, last_step)
    step = row["step"] == 0 ? "Initial (0/$last_step)" : "Step $(row["step"])/$last_step"
    return join(
        [
            step, "α $(number(row["alpha"]))", "Epush $(number(row["Epush"]))",
            "Epull $(number(row["Epull"]))", "Objective $(number(row["objective"]))",
            "$(number(row["seconds"])) s", row["termination"],
            "Sampled cube violation $(number(row["sampled_max_cube_violation"]))",
            "Sampled Oklab drift $(number(row["sampled_max_oklab_drift"]))",
        ], " · "
    )
end

_probe_values(v::AbstractVector) = join(number.(v), ", ")
_probe_values(v) = number(v)

function inverse_trace_panel(d, params)
    trace = Base.get(d.traces, "inverse", nothing)
    probes = trace === nothing ? [] : trace["probes"]
    index = tryparse(Int, Base.get(params, "probe", "1"))
    row = index !== nothing && checkbounds(Bool, probes, index) ? probes[index] : nothing
    io = IOBuffer()
    print(io, "<section class=\"card\" id=\"inverse-trace\"><h2>Inverse solve walkthrough</h2>")
    if trace === nothing
        print(io, "<p class=\"muted\">Recorded inverse probes unavailable for this payload.</p>")
    else
        print(io, "<form action=\"/tables#inverse-trace\"><label>Recorded probe <select name=\"probe\">")
        for (j, probe) in enumerate(probes)
            selected = j == index ? " selected" : ""
            print(io, "<option value=\"$j\"$selected>", _esc(probe["label"]), "</option>")
        end
        print(io, "</select></label> <button>Show probe</button></form>")
    end
    print(io, "<form action=\"/tables#inverse-trace\">")
    for (key, channel) in zip(("r", "g", "b"), ("R", "G", "B"))
        print(io, "<label>Linear $channel <input type=\"number\" name=\"$key\" min=\"0\" max=\"1\" step=\"any\" required></label> ")
    end
    print(io, "<button>Inspect RGB</button></form>")

    if any(key -> haskey(params, key), ("r", "g", "b"))
        values = ntuple(i -> tryparse(Float64, Base.get(params, ("r", "g", "b")[i], "")), 3)
        if !all(x -> x !== nothing && 0 <= x <= 1, values)
            print(io, "<p>Enter three finite linear RGB channels from 0 to 1.</p>")
            return String(take!(io)) * "</section>"
        end
        rgb = SVector{3, Float64}(values)
        found = findfirst(p -> p["target"] == collect(rgb), probes)
        if found === nothing
            fresh_inverse_probe(io, d, rgb)
            return String(take!(io)) * "</section>"
        end
        row = probes[found]
    end
    if row !== nothing
        print(io, "<h3>", _esc(row["label"]), "</h3>")
        print(io, "<p>Recorded diagnostic probe. The solve used independent seeds alongside table generation.</p>")
        inverse_probe_result(io, row, trace["evaluator"])
    end
    return String(take!(io)) * "</section>"
end

function fresh_inverse_probe(io, d, rgb)
    spectral = d.fitted
    if spectral === nothing && haskey(d.traces, "fit")
        trace = d.traces["fit"]
        spectral = PaintMixPrecompute.replay_fit_snapshot(trace, length(trace["snapshots"]))
    end
    if spectral === nothing
        latent = PaintMix.encode(d.model, rgb)
        print(io, "<p>Fresh diagnostic solve unavailable. This payload has no fitted spectral coefficients available for replay.</p>")
        print(
            io, html_table(
                [
                    ("Target linear RGB", _probe_values(rgb)), ("Runtime concentrations", _probe_values(PaintMix.concentrations(latent))),
                    ("Runtime residual", _probe_values(PaintMix.residual(latent))),
                ]
            )
        )
        return nothing
    end
    cfg = d.sidecar === nothing ? d.cfg : Base.get(d.sidecar, "resolved_config", d.cfg)
    settings = UnmixSettings(haskey(cfg, "unmix") ? cfg : load_config())
    row = PaintMixPrecompute.diagnose_inverse(spectral, rgb, settings)
    row["runtime"] = PaintMixPrecompute.probe_runtime(d.model, spectral, rgb, SVector{4}(row["concentrations"]))
    print(io, "<h3>Fresh diagnostic solve</h3><p>This request uses the fitted spectral evaluator. Generation history and float table interpolation are unavailable for this RGB.</p>")
    inverse_probe_result(io, row, "spectral")
    return nothing
end

function inverse_probe_result(io, row, evaluator)
    ref = row["reference"]
    print(
        io, html_table(
            [
                ("Target linear RGB", _probe_values(row["target"])), ("Evaluator", evaluator), ("Strategy", row["strategy"]),
                ("Floating concentrations", _probe_values(row["concentrations"])), ("Reconstructed linear RGB", _probe_values(row["rgb"])),
                ("Continuous solver SSE", number(row["sse"])), ("Termination", row["termination"]),
                ("Reference SSE and signed gap", _probe_values([ref["sse"], ref["sse_gap"]])),
            ]
        )
    )
    print(io, "<p>The reference is a numerical comparison across active faces with bounded local restarts.</p>")
    runtime = row["runtime"]
    stages = (
        ("rounded_bytes_of_diagnostic_winner", "Joint byte rounding of diagnostic concentrations"),
        ("rounded_winner_rgb", "Spectral RGB after joint rounding"),
        ("rounded_winner_max_channel_delta", "Joint rounding RGB delta (max channel)"),
        ("float_inverse_concentrations", "Interpolated float inverse concentrations"),
        ("float_inverse_sse", "Float inverse SSE"),
        ("inverse_interpolation_and_vertex_solve_delta", "Inverse interpolation and vertex solve delta (max RGB channel)"),
        ("inverse_byte_quantization_concentration_delta", "Inverse byte quantization delta (max concentration)"),
        ("runtime_concentrations", "Interpolated runtime concentrations"),
        ("forward_interpolation_max_channel_delta", "Forward interpolation delta (max RGB channel)"),
        ("forward_byte_quantization_max_channel_delta", "Forward byte quantization delta (max RGB channel)"),
        ("runtime_rgb", "Runtime forward RGB"),
        ("runtime_residual", "Signed runtime residual"),
    )
    print(io, html_table([(label, _probe_values(runtime[key])) for (key, label) in stages if haskey(runtime, key)]))
    if haskey(runtime, "inverse_cell_vertices")
        vertices = runtime["inverse_cell_vertices"]
        print(io, "<details><summary>Inverse cell vertices and stored byte rounding</summary>")
        print(
            io, html_table(
                [(number(v["weight"]), _probe_values(v["float_concentrations"]), _probe_values(v["stored_bytes"])) for v in vertices];
                headers = ["Interpolation weight", "Float vertex concentrations", "Stored concentration bytes"]
            ), "</details>"
        )
    end
    for (i, seed) in enumerate(row["seeds"])
        label = i == row["winning_seed"] ? "Seed $i (winner)" : "Seed $i"
        print(io, "<details><summary>$label · ", _esc(_probe_values(seed["concentrations"])), "</summary>")
        for face in seed["faces"]
            active = join(face["active"], ", ")
            print(io, "<h4>Active pigment slots ", _esc(active), "</h4>")
            entries = [("initial", _probe_values(face["initial_concentrations"]), _probe_values(face["initial_rgb"]), number(face["initial_sse"]))]
            for point in face["accepted"]
                push!(entries, (string(point["iteration"]), _probe_values(point["concentrations"]), _probe_values(point["rgb"]), number(point["sse"])))
            end
            push!(entries, ("final", _probe_values(face["concentrations"]), _probe_values(face["rgb"]), number(face["sse"])))
            print(io, html_table(entries; headers = ["Iteration", "Concentrations", "Linear RGB", "SSE"]))
            print(io, "<p>", _esc(face["termination"]), " · simplex projected gradient ", number(face["simplex_projected_gradient_inf"]), "</p>")
        end
        print(io, "</details>")
    end
    if row["omitted_accepted_iterates"] > 0
        print(io, "<p>", row["omitted_accepted_iterates"], " accepted iterates omitted at the capture limit.</p>")
    end
    return nothing
end
