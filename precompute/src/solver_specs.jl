# Configuration and identity live outside optimizer kernels.

const OPTIMIZER_FACTORIES = Dict{Symbol, Any}(
    :lbfgs => Optim.LBFGS, :bfgs => Optim.BFGS, :newton => Optim.Newton,
    :nelder_mead => Optim.NelderMead, :lbfgsb => OptimizationLBFGSB.LBFGSB,
)

"""Register a keyword constructor for a named Optimization.jl algorithm."""
function register_optimizer!(name::Union{Symbol, AbstractString}, constructor)
    OPTIMIZER_FACTORIES[Symbol(name)] = constructor
    return nothing
end

"""Construct a registered algorithm with explicit constructor options."""
function optimizer(name::Union{Symbol, AbstractString}; kwargs...)
    factory = get(OPTIMIZER_FACTORIES, Symbol(name), nothing)
    factory === nothing && throw(ConfigError("unknown optimizer $name; register its constructor or supply a SolverSpec"))
    return factory(; kwargs...)
end

"""An algorithm, derivative backend, and solve options with a portable cache identity.

Opaque callbacks or custom state require an explicit `identity` that changes
whenever their behavior changes. Ordinary algorithm fields and options are
recorded automatically; no object addresses or file paths enter the identity.
"""
struct SolverSpec{A, D, K}
    algorithm::A
    adtype::D
    solve_kwargs::K
    identity::Union{Nothing, String}
end
SolverSpec(algorithm; adtype = nothing, solve_kwargs = (;), identity = nothing) =
    SolverSpec(algorithm, adtype, solve_kwargs, identity === nothing ? nothing : String(identity))

_spec_value(x::Union{Real, AbstractString, Symbol}, opaque) = _setting_value(x)
_spec_value(::Nothing, opaque) = Dict("type" => "Nothing")
_spec_value(x::Base.RefValue, opaque) = Dict("ref" => _spec_value(x[], opaque))
_spec_value(x::Type, opaque) = Dict("type" => string(x))
_spec_value(x::Module, opaque) = Dict("module" => string(x))
_spec_value(x::NamedTuple, opaque) = Dict(string(k) => _spec_value(v, opaque) for (k, v) in pairs(x))
_spec_value(x::AbstractDict, opaque) = Dict(string(k) => _spec_value(v, opaque) for (k, v) in pairs(x))
_spec_value(x::Tuple, opaque) = [_spec_value(v, opaque) for v in x]
_spec_value(x::AbstractArray, opaque) = Dict("shape" => collect(size(x)), "values" => [_spec_value(v, opaque) for v in vec(x)])
function _spec_value(x, opaque)
    T = typeof(x)
    owner = parentmodule(T)
    custom_code = !(owner in (Base, Core, @__MODULE__)) && Base.PkgId(owner).uuid === nothing
    if x isa Ptr || x isa IO || ismutabletype(T) || custom_code
        opaque || throw(ArgumentError("$(T) has opaque state or custom code; supply an explicit SolverSpec identity"))
        return Dict("opaque_type" => string(T))
    end
    return Dict("type" => string(T), "fields" => Dict(string(k) => _spec_value(getfield(x, k), opaque) for k in fieldnames(T)))
end

"""Portable description of the actual algorithm, AD backend, and solve options."""
function solver_description(s::SolverSpec)
    opaque = s.identity !== nothing
    return Dict{String, Any}(
        "algorithm" => _spec_value(s.algorithm, opaque), "adtype" => _spec_value(s.adtype, opaque),
        "solve_options" => _spec_value(s.solve_kwargs, opaque),
        "identity" => something(s.identity, "automatic"),
    )
end
solver_description(s) = _spec_value(s, false)

_options(d::AbstractDict) = (; (Symbol(k) => v for (k, v) in pairs(d))...)
_options(d::NamedTuple) = d

"""Resolve the fit solver from configuration without constraining programmatic choices."""
function fit_solver_spec(cfg)
    section = get(get(cfg, "solvers", Dict()), "fit", Dict())
    name = get(section, "method", cfg["surrogate"]["optimizer"])
    algorithm_options = _options(get(section, "algorithm_options", Dict()))
    solve_options = _options(get(section, "solve_options", Dict()))
    return SolverSpec(optimizer(name; algorithm_options...); solve_kwargs = solve_options)
end

"""Resolve a point solver independently of the inverse traversal strategy."""
function inverse_solver_spec(cfg)
    section = get(get(cfg, "solvers", Dict()), "inverse", Dict())
    name = get(section, "method", "active_face_lm")
    name == "active_face_lm" && return ActiveFaceLM()
    algorithm_options = _options(get(section, "algorithm_options", Dict()))
    options = _options(get(section, "solve_options", Dict()))
    return OptimizationInverse(optimizer(name; algorithm_options...); options...)
end

function _development_source_hash(source)
    files = String[]
    for (dir, _, names) in walkdir(joinpath(source, "src"))
        append!(files, [joinpath(dir, n) for n in names if endswith(n, ".jl")])
    end
    return experiment_id(Dict(relpath(f, source) => bytes2hex(sha256(read(f))) for f in files))
end

_algorithm_component(c) = c
_algorithm_component(c::SolverSpec) = c.algorithm

"""Exact versions and source identities for the selected numerical dependency closure."""
function numerical_environment(components...)
    deps = Pkg.dependencies()
    roots = [Base.PkgId(m).uuid for m in (Optim, Optimization, OptimizationOptimJL, OptimizationLBFGSB, ForwardDiff, StaticArrays, PaintMix)]
    for c in components
        values = c isa SolverSpec ? (_algorithm_component(c), c.adtype) : (_algorithm_component(c),)
        for value in values
            id = Base.PkgId(parentmodule(typeof(value))).uuid
            id === nothing || id == Base.PkgId(@__MODULE__).uuid || push!(roots, id)
        end
    end
    seen = Set{Base.UUID}()
    result = Dict{String, Any}()
    while !isempty(roots)
        id = pop!(roots)
        id in seen && continue
        push!(seen, id)
        info = get(deps, id, nothing)
        info === nothing && continue
        record = Dict{String, Any}("uuid" => string(id), "version" => string(info.version))
        if info.tree_hash !== nothing
            record["tree_hash"] = string(info.tree_hash)
        elseif info.source !== nothing && isdir(joinpath(info.source, "src"))
            record["source_hash"] = _development_source_hash(info.source)
        end
        result[info.name] = record
        append!(roots, values(info.dependencies))
    end
    return Dict("julia" => string(VERSION), "packages" => result)
end

"""Load trusted Julia solver configuration and fingerprint its source.

The file must define `configure_solvers(cfg)` returning a named tuple with
`fit` and `inverse` fields. It can import any installed Optimization backend.
"""
function load_solver_file(path::AbstractString, cfg)
    scope = Module(:PaintMixSolverConfig)
    Base.include(scope, abspath(path))
    isdefined(scope, :configure_solvers) || throw(ConfigError("solver file must define configure_solvers(cfg)"))
    choices = Base.invokelatest() do
        getfield(scope, :configure_solvers)(cfg)
    end
    return choices, bytes2hex(sha256(read(path)))
end
