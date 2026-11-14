# Resolved numerical settings and deterministic experiment identities.

"""Validated settings for one surrogate fit, independent of the CLI."""
Base.@kwdef struct FitSettings
    parameterization::Symbol = :softplus
    optimizer::Symbol = :lbfgs
    derivative::Symbol = :forwarddiff
    epsilon::Float64 = 1.0e-6
    max_iterations::Int = 200
    time_limit_seconds::Float64 = 7200.0
    surface_divisions::Int = 20
    alpha_initial::Float64 = 1.0e5
    alpha_final::Float64 = 1.0e-2
    alpha_halvings::Int = 40
    push_tolerance::Float64 = 1.0e-8
    gradient_tolerance::Float64 = 1.0e-8
end

function FitSettings(cfg::AbstractDict)
    s = cfg["surrogate"]
    st = FitSettings(;
        parameterization = Symbol(get(s, "parameterization", "softplus")),
        optimizer = Symbol(get(s, "optimizer", "lbfgs")),
        derivative = Symbol(get(s, "derivative", "forwarddiff")),
        epsilon = Float64(s["epsilon"]), max_iterations = Int(s["max_iterations"]),
        time_limit_seconds = Float64(get(s, "time_limit_seconds", 7200.0)),
        surface_divisions = Int(s["surface_divisions"]),
        alpha_initial = Float64(s["alpha_initial"]), alpha_final = Float64(s["alpha_final"]),
        alpha_halvings = Int(s["alpha_halvings"]),
        push_tolerance = Float64(get(s, "push_tolerance", 1.0e-8)),
        gradient_tolerance = Float64(get(s, "gradient_tolerance", 1.0e-8)),
    )
    st.parameterization in (:softplus, :direct) || throw(ConfigError("unknown parameterization"))
    haskey(OPTIMIZER_FACTORIES, st.optimizer) || throw(ConfigError("unknown surrogate optimizer"))
    st.derivative in (:forwarddiff, :finite_difference) || throw(ConfigError("unknown derivative method"))
    for x in (st.epsilon, st.time_limit_seconds, st.alpha_initial, st.alpha_final, st.gradient_tolerance)
        isfinite(x) && x > 0 || throw(ConfigError("fit scales and tolerances must be positive and finite"))
    end
    isfinite(st.push_tolerance) && st.push_tolerance >= 0 || throw(ConfigError("invalid push tolerance"))
    st.alpha_initial >= st.alpha_final || throw(ConfigError("alpha_initial must be >= alpha_final"))
    st.max_iterations >= 1 && st.surface_divisions >= 1 && st.alpha_halvings >= 0 ||
        throw(ConfigError("invalid fit iteration or sampling limits"))
    return st
end

"""Resolved table settings; execution threads do not affect table identity."""
struct TableSettings
    n::Int
    strategy::Symbol
    evaluator::Symbol
    coarse_n::Int
    fine_max_iterations::Int
    threads::Int
end

function TableSettings(
        cfg::AbstractDict; profile = "dev", strategy = :coarse,
        evaluator = :spectral, coarse_n = nothing, threads = Threads.nthreads()
    )
    profile in ("dev", "release") || throw(ConfigError("profile must be dev or release"))
    n = Int(cfg["grid"][profile == "dev" ? "dev_n" : "release_n"])
    g = get(cfg, "generation", Dict{String, Any}())
    cn = min(something(coarse_n, Int(get(g, "coarse_n", 64))), max(2, n ÷ 2))
    strategy in (:reference, :bulk, :coarse) || throw(ConfigError("unknown generation strategy"))
    evaluator in (:spectral, :float_lut) || throw(ConfigError("unknown inverse evaluator"))
    threads >= 1 || throw(ConfigError("threads must be >= 1"))
    strategy == :coarse && !(2 <= cn < n) && throw(ConfigError("coarse generation needs 2 <= coarse_n < grid_n"))
    fine = Int(get(g, "fine_max_iterations", 15))
    fine >= 1 || throw(ConfigError("fine_max_iterations must be positive"))
    return TableSettings(n, strategy, evaluator, cn, fine, Int(threads))
end

"""Convert immutable settings to their portable, explicit representation."""
settings_dict(s) = Dict{String, Any}(string(k) => _setting_value(getfield(s, k)) for k in fieldnames(typeof(s)))
_setting_value(x::Symbol) = string(x)
_setting_value(x) = x

"""Canonical serialization for identities, independent of dictionary insertion order."""
function canonical_text(x)
    io = IOBuffer()
    _canonical(io, x)
    return String(take!(io))
end
function _canonical(io, x::AbstractDict)
    print(io, '{')
    for k in sort!(collect(keys(x)); by = string)
        print(io, repr(string(k)), ':')
        _canonical(io, x[k])
        print(io, ',')
    end
    return print(io, '}')
end
function _canonical(io, x::AbstractArray)
    print(io, "array", size(x), '[')
    for v in x
        _canonical(io, v)
        print(io, ',')
    end
    return print(io, ']')
end
_canonical(io, x::Symbol) = _canonical(io, string(x))
_canonical(io, x::Union{AbstractString, Real, Nothing}) = print(io, repr(x))

"""SHA-256 of a canonical experiment specification."""
experiment_id(x) = bytes2hex(sha256(canonical_text(x)))

function implementation_hash(files)
    return experiment_id(Dict(f => bytes2hex(sha256(read(joinpath(@__DIR__, f)))) for f in files))
end

"""Derive a development profile before hashing or starting any stage."""
function profile_config(cfg::AbstractDict, profile::AbstractString)
    profile in ("dev", "release") || throw(ConfigError("profile must be dev or release"))
    out = deepcopy(cfg)
    if profile == "dev"
        g = get(cfg, "generation", Dict{String, Any}())
        for (key, cap) in (("surface_divisions", 8), ("max_iterations", 40), ("alpha_halvings", 12))
            out["surrogate"][key] = min(Int(out["surrogate"][key]), Int(get(g, "dev_" * key, cap)))
        end
    end
    return validate_config(out)
end

"""Identity of the physical fit problem, solver, inputs, and implementation."""
function fit_job_id(cfg, db; solver = fit_solver_spec(cfg), source_id = "")
    spec = Dict{String, Any}(k => cfg[k] for k in ("spectra", "color", "saunderson", "pigments"))
    spec["fit"] = settings_dict(FitSettings(cfg))
    spec["solver"] = solver_description(solver)
    spec["solver_source"] = source_id
    spec["inputs"] = Dict(f.role => f.sha256 for f in db.source_files)
    spec["observer"] = db.observer_source.sha256
    spectra, quad = PigmentSpectra(cfg, db), Quadrature(cfg, db)
    spec["physical_inputs"] = Dict(
        "K" => spectra.K, "S" => spectra.S,
        "weights" => quad.weights, "observer" => hcat(quad.x_bar, quad.y_bar, quad.z_bar), "illuminant" => quad.d65
    )
    spec["implementation"] = implementation_hash(["solver_specs.jl", "experiments.jl", "surrogates.jl", "kubelka_munk.jl", "spectra.jl"])
    spec["environment"] = numerical_environment(solver)
    return experiment_id(spec)
end

"""Identity of a table job, including physical coefficients and resolved overrides."""
function table_job_id(cfg, fit, st::TableSettings; point_solver = inverse_solver_spec(cfg), source_id = "")
    settings = settings_dict(st)
    delete!(settings, "threads")
    return experiment_id(
        Dict(
            "settings" => settings, "K" => fit.K, "S" => fit.S,
            "environment" => numerical_environment(point_solver),
            "point_solver" => solver_description(point_solver), "solver_source" => source_id,
            "runtime_kernels" => bytes2hex(sha256(read(joinpath(pkgdir(PaintMix), "src", "lutkernels.jl")))),
            "unmix" => cfg["unmix"], "quantization" => cfg["quantization"],
            "grid" => Dict(k => v for (k, v) in cfg["grid"] if !(k in ("dev_n", "release_n"))),
            "spectra" => cfg["spectra"], "color" => cfg["color"], "saunderson" => cfg["saunderson"],
            "implementation" => implementation_hash(["solver_specs.jl", "tables.jl", "solvers.jl", "unmix.jl", "export.jl", "kubelka_munk.jl"]),
        )
    )
end

"""Validated quality sampling and acceptance policy, separate from model choices."""
struct ValidationSettings
    channel_p99::Float64
    channel_max::Float64
    oklab_p99::Float64
    roundtrip_f32::Float64
    roundtrip_f64::Float64
    colors::Int
    pairs::Int
    seed::Int
end
function ValidationSettings(cfg::AbstractDict)
    v = get(cfg, "validation", Dict{String, Any}())
    st = ValidationSettings(
        Float64(get(v, "lut_channel_p99", 2 / 255)),
        Float64(get(v, "lut_channel_max", 8 / 255)), Float64(get(v, "lut_oklab_p99", 0.01)),
        Float64(get(v, "roundtrip_f32", 2.0e-6)), Float64(get(v, "roundtrip_f64", 1.0e-12)),
        Int(get(v, "colors", 512)), Int(get(v, "pairs", 20000)), Int(get(v, "seed", 20210913))
    )
    all(x -> isfinite(x) && x >= 0, (st.channel_p99, st.channel_max, st.oklab_p99, st.roundtrip_f32, st.roundtrip_f64)) ||
        throw(ConfigError("validation tolerances must be finite and nonnegative"))
    st.colors >= 12 && st.pairs >= 1 && st.seed >= 0 || throw(ConfigError("invalid validation sample counts or seed"))
    return st
end
