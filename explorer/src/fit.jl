# `/fit` data: the surrogate continuation history and diagnostics recorded in
# the provenance sidecar. Pure data; the app file plots it.

"""
    FitHistory

One vector per column of `sidecar["surrogate"]["history"]`, one entry per
continuation step. `epush` is the cube-push penalty, `epull` the perceptual
pull toward the measured pigments, and `objective` their `alpha`-weighted sum.
"""
struct FitHistory
    steps::Vector{Int}
    alpha::Vector{Float64}
    objective::Vector{Float64}
    epush::Vector{Float64}
    epull::Vector{Float64}
    seconds::Vector{Float64}
    iterations::Vector{Int}
    converged::Vector{Bool}
end

function _history_float(entry::AbstractDict, key::AbstractString)
    v = Base.get(entry, key, nothing)
    return v isa Real ? Float64(v) : NaN
end

"""
    fit_history(sidecar) -> Union{Nothing, FitHistory}

The continuation history, or `nothing` when the sidecar is absent or records
none. Missing fields in an entry become `NaN` rather than an error, so a
partial sidecar still plots.
"""
function fit_history(sidecar::Union{Nothing, AbstractDict})
    sidecar === nothing && return nothing
    surrogate = Base.get(sidecar, "surrogate", nothing)
    surrogate isa AbstractDict || return nothing
    entries = Base.get(surrogate, "history", nothing)
    entries isa AbstractVector && !isempty(entries) || return nothing
    count = length(entries)
    steps = Vector{Int}(undef, count)
    alpha = Vector{Float64}(undef, count)
    objective = Vector{Float64}(undef, count)
    epush = Vector{Float64}(undef, count)
    epull = Vector{Float64}(undef, count)
    seconds = Vector{Float64}(undef, count)
    iterations = Vector{Int}(undef, count)
    converged = Vector{Bool}(undef, count)
    for (i, entry) in enumerate(entries)
        entry isa AbstractDict || continue
        s = Base.get(entry, "step", i)
        steps[i] = s isa Real ? Int(s) : i
        alpha[i] = _history_float(entry, "alpha")
        objective[i] = _history_float(entry, "objective")
        epush[i] = _history_float(entry, "Epush")
        epull[i] = _history_float(entry, "Epull")
        seconds[i] = _history_float(entry, "seconds")
        iterations[i] = Int(something(Base.get(entry, "iterations", 0), 0))
        converged[i] = Base.get(entry, "converged", false) === true
    end
    return FitHistory(
        steps, alpha, objective, epush, epull, seconds, iterations, converged
    )
end

fit_history(d::ExplorerData) = fit_history(d.sidecar)

"""
    push_tolerance(cfg) -> Float64

The `Epush` value at which the continuation stops, from the configuration, or
the `1e-8` default when it is absent.
"""
function push_tolerance(cfg::AbstractDict)
    surrogate = Base.get(cfg, "surrogate", nothing)
    surrogate isa AbstractDict || return 1.0e-8
    return Float64(Base.get(surrogate, "push_tolerance", 1.0e-8))
end

push_tolerance(d::ExplorerData) = push_tolerance(d.cfg)

"""
    push_met_step(h, tolerance) -> Union{Nothing, Int}

The index of the first history step whose `Epush` reaches `tolerance`, or
`nothing` when it never does.
"""
function push_met_step(h::FitHistory, tolerance::Real)
    for i in eachindex(h.epush)
        isfinite(h.epush[i]) && h.epush[i] <= tolerance && return i
    end
    return nothing
end

"""
    fit_diagnostics(sidecar) -> Union{Nothing, Dict}

The `surrogate.diagnostics` section, or `nothing` when the sidecar is absent.
"""
function fit_diagnostics(sidecar::Union{Nothing, AbstractDict})
    sidecar === nothing && return nothing
    surrogate = Base.get(sidecar, "surrogate", nothing)
    surrogate isa AbstractDict || return nothing
    diagnostics = Base.get(surrogate, "diagnostics", nothing)
    return diagnostics isa AbstractDict ? diagnostics : nothing
end

fit_diagnostics(d::ExplorerData) = fit_diagnostics(d.sidecar)
