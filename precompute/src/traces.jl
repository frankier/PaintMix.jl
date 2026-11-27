# Trace documents are separate from the numerical caches and payload format.
const TRACE_MAX_BYTES = 16 * 1024^2
const FIT_TRACE_MAX_STEPS = 64
const INVERSE_TRACE_ITERATES = 64
const INVERSE_TRACE_PROBES = 12

function write_trace(dir, doc)
    io = IOBuffer(; maxsize = TRACE_MAX_BYTES + 1)
    TOML.print(io, doc; sorted = true)
    bytes = take!(io)
    length(bytes) <= TRACE_MAX_BYTES || throw(InputError("trace exceeds the size limit"))
    digest = bytes2hex(sha256(bytes))
    file = "$(doc["kind"])-trace-$digest.toml"
    mkpath(dir)
    mktemp(dir) do path, out
        write(out, bytes)
        close(out)
        mv(path, joinpath(dir, file); force = true)
    end
    return Dict("file" => file, "sha256" => digest, "schema_version" => 1)
end

function read_trace(dir, reference, kind)
    digest = reference["sha256"]
    occursin(r"^[a-f0-9]{64}$", digest) || throw(InputError("invalid trace digest"))
    file = "$kind-trace-$digest.toml"
    reference["file"] == file || throw(InputError("invalid trace filename"))
    path = joinpath(dir, file)
    bytes = open(path) do io
        read(io, TRACE_MAX_BYTES + 1)
    end
    length(bytes) <= TRACE_MAX_BYTES || throw(InputError("trace exceeds the size limit"))
    bytes2hex(sha256(bytes)) == digest || throw(InputError("trace digest mismatch"))
    doc = TOML.parse(String(bytes))
    doc["kind"] == kind && doc["schema_version"] == reference["schema_version"] == 1 ||
        throw(InputError("unsupported trace schema"))
    return doc
end

coefficient_digest(K, S) = experiment_id(Dict("K" => K, "S" => S))

function check_fit_trace_limits(st, W)
    st.alpha_halvings < FIT_TRACE_MAX_STEPS && st.surface_divisions <= 48 && W <= 256 ||
        throw(ConfigError("fit tracing supports up to 64 solves, 48 surface divisions and 256 wavelengths"))
    return nothing
end

function build_fit_trace(cfg, spectra, problem, snapshots, history)
    sq = problem.quadrature
    # The corpus includes interior and near-boundary points. Exclude the
    # actual fitting nodes including shared edges and pure vertices.
    corpus = validation_corpus(; colors = 128, pairs = 1)
    heldout = filter(c -> all(t -> maximum(abs, c - t) > 1.0e-14, sq.points), corpus.concentrations)
    targets = [linear_srgb_to_oklab(mix_rgb(problem.base, c)) for c in heldout]
    x = initial_parameters(problem.parameters, spectra)
    push, pull = objective_terms(problem, x)
    alpha = first(alpha_schedule(FitSettings(cfg)))
    initial = Dict(
        "step" => 0, "alpha" => alpha, "Epush" => push, "Epull" => pull,
        "objective" => push + alpha * pull, "seconds" => 0.0,
        "termination" => "Initial"
    )
    rows = map(snapshots, [initial; history]) do (K, S), entry
        model = SpectralModel(problem.base; K, S)
        worst_push, worst_drift = 0.0, 0.0
        for (c, target) in zip(heldout, targets)
            rgb = mix_rgb(model, c)
            worst_push = max(worst_push, cube_outside_penalty(rgb))
            worst_drift = max(worst_drift, oklab_distance_squared(linear_srgb_to_oklab(rgb), target))
        end
        merge(
            entry, Dict(
                "absorption" => _matrix_rows(K), "scattering" => _matrix_rows(S),
                "sampled_max_cube_violation" => sqrt(worst_push),
                "sampled_max_oklab_drift" => sqrt(worst_drift)
            )
        )
    end
    quad = problem.base.quad
    return Dict{String, Any}(
        "kind" => "fit", "schema_version" => 1,
        "config_sha256" => experiment_id(cfg), "pigments" => collect(spectra.codes),
        "input_coefficients_sha256" => coefficient_digest(spectra.K, spectra.S),
        "wavelength_nm" => quad.wavelength,
        "integration" => Dict(
            "weights" => quad.weights, "x_bar" => quad.x_bar, "y_bar" => quad.y_bar,
            "z_bar" => quad.z_bar, "d65" => quad.d65, "xyz_to_rgb" => collect(vec(quad.xyz_to_rgb)),
            "norm" => quad.norm, "k1" => problem.base.k1, "k2" => problem.base.k2
        ),
        "boundary" => Dict("concentrations" => collect.(sq.points), "weights" => sq.weights, "faces" => sq.face),
        "validation" => Dict(
            "corpus_sha256" => corpus.digest, "concentrations" => collect.(heldout),
            "interpretation" => "finite samples of cube violation and Oklab drift from the measured model"
        ),
        "surrogate_sha256" => coefficient_digest(last(snapshots)...), "snapshots" => rows
    )
end

"""Replay one fit step with the integration data recorded in the trace."""
function replay_fit_snapshot(trace, index)
    row = trace["snapshots"][index]
    q = trace["integration"]
    quad = Quadrature{Float64}(;
        wavelength = trace["wavelength_nm"], weights = q["weights"], x_bar = q["x_bar"],
        y_bar = q["y_bar"], z_bar = q["z_bar"], d65 = q["d65"],
        xyz_to_rgb = SMatrix{3, 3, Float64}(q["xyz_to_rgb"]), norm = q["norm"]
    )
    K = hcat(row["absorption"]...)
    S = hcat(row["scattering"]...)
    return SpectralModel(K, S, quad, q["k1"], q["k2"])
end

mutable struct InverseTraceScratch
    core::SolverScratch
    seeds::Vector{Dict{String, Any}}
    winner::Int
    remaining::Int
    omitted::Int
end

Base.getproperty(s::InverseTraceScratch, name::Symbol) =
    name in (:function_evaluations, :gradient_evaluations) ? getproperty(getfield(s, :core), name) : getfield(s, name)
Base.setproperty!(s::InverseTraceScratch, name::Symbol, value) =
    name in (:function_evaluations, :gradient_evaluations) ? setproperty!(getfield(s, :core), name, value) : setfield!(s, name, value)

function _trace_seed!(trace::InverseTraceScratch, seed)
    push!(trace.seeds, Dict("concentrations" => collect(seed), "faces" => Dict{String, Any}[]))
    return nothing
end

_trace_winner!(trace::InverseTraceScratch) = (trace.winner = length(trace.seeds); nothing)

struct FaceIterates
    trace::InverseTraceScratch
    rows::Vector{Dict{String, Any}}
end

function _trace_iterate!(capture::FaceIterates, iteration, c, rgb, sse)
    if capture.trace.remaining == 0
        capture.trace.omitted += 1
        return nothing
    end
    push!(capture.rows, Dict("iteration" => iteration, "concentrations" => collect(c), "rgb" => collect(rgb), "sse" => sse))
    capture.trace.remaining -= 1
    return nothing
end

function _polish!(trace::InverseTraceScratch, model, rgb, active, seed, settings)
    theta = _theta_from_c(active, seed)
    capture = FaceIterates(trace, Dict{String, Any}[])
    c, sse, iterations, converged = lm_active!(trace.core, model, rgb, active, theta, settings, capture)
    initial = _c_from_active(active, _softmax_k(SVector(theta..., zero(eltype(seed)))))
    initial_rgb = _mix(model, initial)
    fitted_rgb = _mix(model, c)
    _jacobian4!(trace.core.J, model, c)
    gradient = SVector{4}(transpose(trace.core.J) * (2 * (fitted_rgb - rgb)))
    stationarity = maximum(abs, c - project_to_simplex(c - gradient))
    push!(
        last(trace.seeds)["faces"], Dict(
            "active" => collect(active), "initial_concentrations" => collect(initial),
            "initial_rgb" => collect(initial_rgb), "initial_sse" => _sse(rgb, initial_rgb),
            "accepted" => capture.rows, "concentrations" => collect(c),
            "rgb" => collect(fitted_rgb), "sse" => sse, "iterations" => iterations,
            "termination" => length(active) == 1 ? "vertex" : (converged ? "stationary" : "stalled_or_limit"),
            "simplex_projected_gradient_inf" => stationarity
        )
    )
    return c, sse, iterations, converged
end

_trace_pure(i) = SVector{4, Float64}(ntuple(j -> Float64(i == j), 4))

function diagnostic_seeds(model, rgb)
    uniform = SVector(0.25, 0.25, 0.25, 0.25)
    nearest = argmin(i -> _sse(rgb, _mix(model, _trace_pure(i))), 1:4)
    return uniform, 0.1 * uniform + 0.9 * _trace_pure(nearest)
end

diagnostic_settings(settings) =
    UnmixSettings{Float64}(min(settings.max_iterations, 64), settings.tolerance, settings.face_threshold, min(settings.restarts, 4))

"""Replay the normal bulk kernel on two deterministic seeds with bounded LM tracing."""
function diagnose_inverse(model, rgb, settings; seeds = diagnostic_seeds(model, rgb))
    length(seeds) in 1:2 || throw(ArgumentError("diagnostic solve supports one or two seeds"))
    settings = diagnostic_settings(settings)
    trace = InverseTraceScratch(SolverScratch(), Dict{String, Any}[], 0, INVERSE_TRACE_ITERATES, 0)
    winner = unmix_bulk!(trace, model, rgb, seeds, settings)
    ref = unmix_reference(model, rgb; settings)
    return Dict{String, Any}(
        "target" => collect(rgb), "evaluations" => winner.function_evaluations,
        "gradient_evaluations" => winner.gradient_evaluations, "iterations" => winner.iterations, "settings" => settings_dict(settings),
        "strategy" => "active-face LM with deterministic independent seeds",
        "seeds" => trace.seeds, "winning_seed" => trace.winner,
        "concentrations" => collect(winner.c), "rgb" => collect(_mix(model, winner.c)),
        "sse" => winner.sse, "termination" => string(winner.termination),
        "accepted_iterate_limit" => INVERSE_TRACE_ITERATES,
        "omitted_accepted_iterates" => trace.omitted,
        "reference" => Dict(
            "concentrations" => collect(ref.c), "rgb" => collect(_mix(model, ref.c)), "sse" => ref.sse,
            "termination" => string(ref.termination), "sse_gap" => winner.sse - ref.sse,
            "interpretation" => "numerical active-face enumeration with bounded local restarts"
        )
    )
end

function _table_concentrations(xyz)
    c = PaintMix._repair_simplex(xyz...)
    return SVector(c..., max(0.0, 1.0 - sum(c)))
end

function probe_runtime(model, spectral, rgb, continuous, ft = nothing)
    latent = PaintMix.encode(model, rgb)
    runtime_c = PaintMix.concentrations(latent)
    runtime_rgb = PaintMix.forward_rgb(model, runtime_c)
    bytes = _inverse_bytes(continuous[1], continuous[2], continuous[3])
    rounded = SVector{4, Int}(bytes..., 255 - sum(Int, bytes))
    c_rounded = rounded / 255
    spectral_rgb = _mix(spectral, runtime_c)
    continuous_rgb = _mix(spectral, continuous)
    rounded_rgb = _mix(spectral, c_rounded)
    result = Dict{String, Any}(
        "rounded_bytes_of_diagnostic_winner" => collect(rounded),
        "rounded_winner_rgb" => collect(rounded_rgb),
        "rounded_winner_max_channel_delta" => maximum(abs, rounded_rgb - continuous_rgb),
        "runtime_concentrations" => collect(runtime_c),
        "runtime_rgb" => collect(runtime_rgb), "runtime_residual" => collect(rgb - runtime_rgb),
        "spectral_rgb_at_runtime_concentrations" => collect(spectral_rgb)
    )
    ft === nothing && return result
    corners, fx, fy, fz = cell(ft.inverse, ft.n, rgb..., 1.0)
    byte_corners, _, _, _ = cell(model.inverse.data, ft.n, rgb..., 1.0)
    xyz = trilinear3(corners, fx, fy, fz)
    c = _table_concentrations(xyz)
    vertices = map(1:8) do index
        i = index - 1
        floating, bytes = corners[index], byte_corners[index]
        weight = (isodd(i) ? fx : 1 - fx) *
            (isodd(i >> 1) ? fy : 1 - fy) * (isodd(i >> 2) ? fz : 1 - fz)
        return Dict(
            "weight" => weight, "float_concentrations" => collect(_table_concentrations(floating)),
            "stored_bytes" => [Int.(bytes)..., 255 - sum(Int, bytes)]
        )
    end
    float_forward_rgb = _mix(forward_float_lut(ft.n, ft.forward), runtime_c)
    float_inverse_rgb = _mix(spectral, c)
    merge!(
        result, Dict(
            "inverse_cell_vertices" => vertices,
            "float_inverse_concentrations" => collect(c),
            "float_inverse_rgb" => collect(float_inverse_rgb),
            "float_inverse_sse" => _sse(rgb, float_inverse_rgb),
            "inverse_interpolation_and_vertex_solve_delta" => maximum(abs, float_inverse_rgb - continuous_rgb),
            "inverse_byte_quantization_concentration_delta" => maximum(abs, runtime_c - c),
            "forward_interpolation_rgb" => collect(float_forward_rgb),
            "forward_interpolation_max_channel_delta" => maximum(abs, float_forward_rgb - spectral_rgb),
            "forward_byte_quantization_max_channel_delta" => maximum(abs, runtime_rgb - float_forward_rgb)
        )
    )
    return result
end

function build_inverse_trace(model, spectral, evaluator, ft, provenance, settings)
    corpus = validation_corpus(; colors = 64, pairs = 1)
    reference_sse(rgb) = unmix_reference(evaluator, rgb; settings = diagnostic_settings(settings)).sse
    corner = argmax(reference_sse, corpus.colors[1:8])
    uniform = SVector(0.25, 0.25, 0.25, 0.25)
    seeds = [
        ("Near black", SVector(0.001, 0.001, 0.001)),
        ("Interior mixture", clamp.(_mix(evaluator, uniform), 0, 1)),
        ("Cube corner with largest sampled reference SSE", corner),
        [("Near pigment $i", clamp.(_mix(evaluator, 0.02 * uniform + 0.98 * _trace_pure(i)), 0, 1)) for i in 1:4]...,
    ]
    # Use stable corpus order to break ties and retain distinct target RGBs.
    ordered = sort(corpus.colors; by = rgb -> -sum(abs2, rgb - PaintMix.forward_rgb(model, PaintMix.concentrations(PaintMix.encode(model, rgb)))))
    seen = Set(last.(seeds))
    for rgb in ordered
        length(seeds) == INVERSE_TRACE_PROBES && break
        rgb in seen && continue
        push!(seeds, ("Residual quality tail", rgb))
        push!(seen, rgb)
    end
    probes = map(seeds) do (label, rgb)
        row = diagnose_inverse(evaluator, rgb, settings)
        row["label"] = label
        row["runtime"] = probe_runtime(model, spectral, rgb, SVector{4}(row["concentrations"]), ft)
        row
    end
    return Dict{String, Any}(
        "kind" => "inverse", "schema_version" => 1, "probe_limit" => INVERSE_TRACE_PROBES,
        "corpus_sha256" => corpus.digest, "evaluator" => provenance["table_settings"]["evaluator"],
        "capture" => "deterministic diagnostic solves alongside the generated table",
        "generation_table_settings" => provenance["table_settings"], "generation_point_solver" => provenance["point_solver"],
        "config_sha256" => provenance["resolved_config_sha256"], "inputs" => provenance["inputs"],
        "fit_job_id" => provenance["fit_job_id"], "table_job_id" => provenance["table_job_id"],
        "payload_sha256" => provenance["payload_sha256"], "model_id" => PaintMix.model_id(model),
        "surrogate_sha256" => provenance["surrogate_sha256"], "probes" => probes
    )
end
