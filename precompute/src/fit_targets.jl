struct PushTarget
    tolerance::Float64
end
struct GamutTarget
    tolerance::Float64
    points::Vector{SVector{4, Float64}}
    digest::String
end

_fit_stopping(st, problem) = _fit_stopping(Val(st.stop_rule), st, problem)
_fit_stopping(::Val{:push}, st, problem) = PushTarget(st.push_tolerance)
function _fit_stopping(::Val{:heldout_gamut}, st, problem)
    points = unique(vcat(SurfaceQuadrature(41).points, validation_corpus(; colors = 256, pairs = 1).concentrations))
    # Exclude fitting nodes; retain near-vertex and interior probes.
    filter!(c -> all(t -> maximum(abs, c - t) > 1.0e-14, problem.quadrature.points), points)
    return GamutTarget(st.gamut_tolerance, points, experiment_id([collect(c) for c in points]))
end
_fit_progress(rule::PushTarget, problem, x, epush) = (reached = epush <= rule.tolerance, values = Dict{String, Any}())
function _fit_progress(rule::GamutTarget, problem, x, epush)
    base = problem.base
    K, S = _parameter_matrices_wl(problem.parameters, x, size(base.K, 1))
    worst = maximum(c -> cube_outside_penalty(mix_rgb_params(K, S, base.k1, base.k2, base.quad, c)), rule.points)
    violation = sqrt(worst)
    reached = isfinite(violation) && violation <= rule.tolerance
    return (reached, values = Dict{String, Any}("heldout_gamut_violation" => violation, "heldout_target_reached" => reached))
end
_target_description(rule::PushTarget) = Dict("rule" => "push", "tolerance" => rule.tolerance)
_target_description(rule::GamutTarget) = Dict("rule" => "heldout_gamut", "tolerance" => rule.tolerance, "samples" => length(rule.points), "corpus_sha256" => rule.digest)
