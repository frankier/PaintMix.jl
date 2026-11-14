# Plan step 4: surrogate pigment fitting, equations (15)-(17).
#
# The original pigments `P*` mix to colors outside the sRGB cube. The
# surrogate pigments `Q*` are fitted so that their gamut lies inside the cube
# while staying as close as possible in Oklab:
#
#     arg min_Q  E_push(Q) + alpha * E_pull(Q, P*)   s.t.  K, S > 0     (15)
#     E_push(Q) = integral_dOmega  max(0, phi(mix_Q(c)))^2 ds           (16)
#     E_pull(Q) = integral_dOmega  ||psi(mix_Q(c)) - psi(mix_P*(c))||^2 ds (17)
#
# Positivity is handled by reparameterizing `K = epsilon + softplus(thetaK)`
# and the same for `S`, so unconstrained `Optim.LBFGS` can be used. The paper
# uses L-BFGS-B; this is an intentional, recorded solver substitution.

"""
    SurfaceQuadrature{T}

A quadrature over the four faces of the concentration simplex, at equally
spaced concentrations, as the paper describes for equations (16) and (17).

  * `points`: concentration 4-vectors, one coordinate exactly zero.
  * `weights`: equal per-sample weights; every face carries the same total.
  * `face`: index `1:4` of the zero coordinate, for reporting.

The weights are the *concentration*-space measure, not the RGB surface area
element. See [`rgb_surface_weights`](@ref) for that diagnostic and
`precompute/README.md` for why the fit uses this simpler, smooth measure.
"""
struct SurfaceQuadrature{T <: AbstractFloat}
    points::Vector{SVector{4, T}}
    weights::Vector{T}
    face::Vector{Int}
end

"""
    SurfaceQuadrature(d) -> SurfaceQuadrature
    SurfaceQuadrature(cfg) -> SurfaceQuadrature

The equally spaced surface quadrature for `d` divisions per axis, or read
`surrogate.surface_divisions` from `cfg`.
"""
SurfaceQuadrature(cfg::AbstractDict) =
    SurfaceQuadrature(Int(get(cfg["surrogate"], "surface_divisions", 20)))

function SurfaceQuadrature(d::Integer)
    d >= 1 || throw(ArgumentError("surface quadrature needs at least one division"))
    points = SVector{4, Float64}[]
    faces = Int[]
    step = 1 / d
    for f in 1:4
        free = [i for i in 1:4 if i != f]
        for a in 0:d, b in 0:(d - a)
            c = zeros(4)
            c[free[1]] = a * step
            c[free[2]] = b * step
            c[free[3]] = (d - a - b) * step
            push!(points, SVector(c[1], c[2], c[3], c[4]))
            push!(faces, f)
        end
    end
    n = length(points)
    weights = fill(1 / n, n)
    return SurfaceQuadrature(points, weights, faces)
end

"""
    surface_targets(model, sq) -> Vector{SVector{3,Float64}}

Oklab coordinates of `mix_P*(c)` at every quadrature sample. The pull targets
do not depend on the fitted parameters, so they are computed once.
"""
function surface_targets(model::SpectralModel, sq::SurfaceQuadrature)
    return [linear_srgb_to_oklab(mix_rgb(model, c)) for c in sq.points]
end

"""
    Epush(model, sq) -> Float64

Equation (16): the weighted sum of squared cube violations over the
quadrature samples of the model's gamut boundary.

The paper writes this as an integral over `dOmega` of `max(0, phi)^2`; here
`dOmega` is approximated by the equally spaced concentration samples, so the
weights are uniform. `phi^2` is evaluated with
[`cube_outside_penalty`](@ref), which is the same value without the square
root.
"""
function Epush(model::SpectralModel, sq::SurfaceQuadrature)
    total = 0.0
    @inbounds for (i, c) in enumerate(sq.points)
        total += sq.weights[i] * cube_outside_penalty(mix_rgb(model, c))
    end
    return total
end

"""
    Epull(model, sq, targets) -> Float64

Equation (17): the weighted sum of squared Oklab differences between the
fitted model and the original pigments at corresponding concentrations.
`targets` is the output of [`surface_targets`](@ref).
"""
function Epull(model::SpectralModel, sq::SurfaceQuadrature, targets)
    total = 0.0
    @inbounds for (i, c) in enumerate(sq.points)
        total += sq.weights[i] * oklab_distance_squared(
            linear_srgb_to_oklab(mix_rgb(model, c)), targets[i]
        )
    end
    return total
end

"""
    rgb_surface_weights(model, sq) -> Vector{Float64}

The RGB-surface area element `|t1 x t2|` of `dOmega` at every quadrature
sample, where `t1, t2` are the images of two tangent directions of the
concentration face.

This is the *literal* `ds` of equations (16) and (17). It depends on the
fitted parameters, so using it as a weight makes the integral a moving-boundary
one; the fit therefore uses the fixed concentration measure and this function
exists to quantify the difference. Reported by
[`quadrature_report`](@ref).
"""
function rgb_surface_weights(model::SpectralModel{T}, sq::SurfaceQuadrature) where {T}
    W = size(model.K, 1)
    J = Matrix{T}(undef, 3, 4)
    out = Vector{Float64}(undef, length(sq.points))
    @inbounds for (i, c) in enumerate(sq.points)
        mix_rgb_jacobian!(J, model, c)
        f = sq.face[i]
        free = Int[j for j in 1:4 if j != f]
        # Tangent directions of the face in concentration space.
        d = free[3]
        u1 = J[1, free[1]] - J[1, d]
        u2 = J[2, free[1]] - J[2, d]
        u3 = J[3, free[1]] - J[3, d]
        v1 = J[1, free[2]] - J[1, d]
        v2 = J[2, free[2]] - J[2, d]
        v3 = J[3, free[2]] - J[3, d]
        out[i] = sqrt(
            (u2 * v3 - u3 * v2)^2 + (u3 * v1 - u1 * v3)^2 + (u1 * v2 - u2 * v1)^2
        )
    end
    return out
end

"""
    quadrature_report(model, sq) -> Dict{String,Any}

Compare the concentration-space quadrature weights with the RGB-surface area
weights, so the deviation the fit accepts is measured rather than assumed.
"""
function quadrature_report(model::SpectralModel, sq::SurfaceQuadrature)
    jac = rgb_surface_weights(model, sq)
    uni = sq.weights
    total_u = sum(uni)
    total_j = sum(jac)
    ratios = map(eachindex(uni)) do i
        total_j > 0 ? (uni[i] * total_j) / (jac[i] * total_u) : NaN
    end
    nu = count(!isnan, ratios)
    return Dict{String, Any}(
        "samples" => length(uni),
        "weighting" => "concentration-simplex-area",
        "rgb_surface_jacobian" => Dict{String, Any}(
            "total" => total_j,
            "mean" => total_j / length(jac),
            "max" => maximum(jac),
            "min" => minimum(jac),
        ),
        "normalized_ratio" => Dict{String, Any}(
            "mean" => nu == 0 ? 0.0 : sum(r for r in ratios if !isnan(r)) / nu,
            "min" => nu == 0 ? 0.0 : minimum(r for r in ratios if !isnan(r)),
            "max" => nu == 0 ? 0.0 : maximum(r for r in ratios if !isnan(r)),
        ),
        "note" =>
            "equations (16)-(17) integrate over the RGB gamut boundary; the fit " *
            "uses equally spaced concentrations, as the paper describes, and these " *
            "ratios bound the resulting weighting deviation",
    )
end

# --- positivity reparameterization ----------------------------------------

"""
    softplus(x)

`log(1 + exp(x))`, written in the stable form
`max(x, 0) + log1p(exp(-abs(x)))` so large magnitudes do not overflow.
"""
@inline function softplus(x::T) where {T <: Real}
    return max(x, zero(T)) + log1p(exp(-abs(x)))
end

"""
    inv_softplus(y)

The inverse of [`softplus`](@ref): `log(expm1(y))`, with `y` floored at a
small positive value so the initial parameters are finite even when a
measured coefficient is below the configured floor.
"""
@inline function inv_softplus(y::T) where {T <: Real}
    v = max(y, T(1.0e-12))
    return v + log(-expm1(-v))
end

@inline _unpack(θ::Real, epsilon::Real) = epsilon + softplus(θ)

"""
    SurrogateFit{T}

The fitted surrogates plus the record of how they were obtained.

  * `K`, `S`: `4 x W` pigment-major fitted coefficients.
  * `theta`: the optimization coordinates (the parameterization is recorded in diagnostics).
  * `history`: one entry per continuation step, with `alpha`, objective
    terms, iteration count, and convergence flag.
  * `diagnostics`: worst violations and perceptual deviations.
  * `quadrature`: the [`SurfaceQuadrature`](@ref) used.
"""
struct SurrogateFit{T <: AbstractFloat}
    K::Matrix{T}
    S::Matrix{T}
    theta::Vector{T}
    history::Vector{Dict{String, Any}}
    diagnostics::Dict{String, Any}
    quadrature::SurfaceQuadrature{T}
end

"""
    initial_theta(spectra, epsilon) -> Vector

Parameters that reproduce the original pigments exactly:
`theta = inv_softplus(K - epsilon)` and likewise for `S`.
"""
function initial_theta(spectra::PigmentSpectra{T}, epsilon::Real) where {T}
    W = length(spectra.wavelength)
    theta = Vector{T}(undef, 8W)
    @inbounds for i in 1:4, j in 1:W
        theta[(i - 1) * W + j] = inv_softplus(spectra.K[i, j] - epsilon)
        theta[4W + (i - 1) * W + j] = inv_softplus(spectra.S[i, j] - epsilon)
    end
    return theta
end

"""
    theta_parameters(theta, W, epsilon) -> (K, S)

Split a parameter vector into `4 x W` pigment-major `K` and `S` matrices.
"""
function theta_parameters(theta::AbstractVector{T}, W::Integer, epsilon::Real) where {T}
    K = Matrix{T}(undef, 4, W)
    S = Matrix{T}(undef, 4, W)
    @inbounds for i in 1:4, j in 1:W
        K[i, j] = _unpack(theta[(i - 1) * W + j], epsilon)
        S[i, j] = _unpack(theta[4W + (i - 1) * W + j], epsilon)
    end
    return K, S
end

"""
    theta_parameters_wl(theta, W, epsilon) -> (K, S)

Like [`theta_parameters`](@ref) but in the wavelength-major `W x 4` layout
that [`mix_rgb_params`](@ref) indexes. The pigment-major layout is what the
provenance sidecar records; this one is what the objective evaluates.
"""
function theta_parameters_wl(theta::AbstractVector{T}, W::Integer, epsilon::Real) where {T}
    K = Matrix{T}(undef, W, 4)
    S = Matrix{T}(undef, W, 4)
    @inbounds for i in 1:4, j in 1:W
        K[j, i] = _unpack(theta[(i - 1) * W + j], epsilon)
        S[j, i] = _unpack(theta[4W + (i - 1) * W + j], epsilon)
    end
    return K, S
end

"""
    surrogate_objective(theta, base, sq, targets, epsilon, alpha) -> Real

`E_push(Q(theta)) + alpha * E_pull(Q(theta), P*)`, the objective of equation
(15) for one continuation step. Accepts `ForwardDiff.Dual` elements in
`theta`.
"""
function surrogate_objective(theta, base, sq, targets, epsilon, alpha)
    problem = SurrogateProblem(SoftplusParameters(Float64(epsilon)), base, sq, targets)
    return _objective_value(problem, theta, alpha)
end

# Compatibility with the original derivative-checking entry point.
surrogate_objective(theta, base, sq, targets, epsilon, alpha, scratch::AbstractMatrix) =
    surrogate_objective(theta, base, sq, targets, epsilon, alpha)

"""
    alpha_schedule(cfg) -> Vector{Float64}

The continuation schedule: `alpha_initial`, halved `alpha_halvings` times,
never below `alpha_final`.
"""
function alpha_schedule(cfg::AbstractDict)
    s = cfg["surrogate"]
    α0 = Float64(s["alpha_initial"])
    αmin = Float64(s["alpha_final"])
    n = Int(s["alpha_halvings"])
    out = Float64[]
    α = α0
    for _ in 0:n
        α < αmin && break
        push!(out, α)
        α /= 2
    end
    isempty(out) && push!(out, αmin)
    return out
end

# --- problem definition, parameterization, and derivative adapters ----------

"""Positive coefficients represented as `epsilon + softplus(theta)`."""
struct SoftplusParameters
    epsilon::Float64
end

"""Physical coefficients with explicit positive lower bounds."""
struct DirectParameters
    epsilon::Float64
end

parameterization(st::FitSettings) = parameterization(Val(st.parameterization), st.epsilon)
parameterization(::Val{:softplus}, epsilon) = SoftplusParameters(epsilon)
parameterization(::Val{:direct}, epsilon) = DirectParameters(epsilon)

"""Initialize either parameterization from the same floored physical coefficients."""
initial_parameters(p::SoftplusParameters, spectra) = initial_theta(spectra, p.epsilon)
function initial_parameters(p::DirectParameters, spectra)
    # Use the same initialization as softplus, including its numerical floor.
    K, S = theta_parameters(initial_theta(spectra, p.epsilon), length(spectra.wavelength), p.epsilon)
    return vcat(vec(permutedims(K)), vec(permutedims(S)))
end

parameter_matrices(p::SoftplusParameters, x, W) = theta_parameters(x, W, p.epsilon)
function parameter_matrices(::DirectParameters, x, W)
    return permutedims(reshape(x[1:(4W)], W, 4)), permutedims(reshape(x[(4W + 1):end], W, 4))
end
_parameter_matrices_wl(p::SoftplusParameters, x, W) = theta_parameters_wl(x, W, p.epsilon)
_parameter_matrices_wl(::DirectParameters, x, W) =
    (reshape(view(x, 1:(4W)), W, 4), reshape(view(x, (4W + 1):(8W)), W, 4))
_parameter_bounds(::SoftplusParameters, x) = (;)
_parameter_bounds(p::DirectParameters, x) = (lb = fill(p.epsilon, length(x)), ub = fill(Inf, length(x)))

"""The physical surrogate fitting problem, independent of an optimization algorithm."""
struct SurrogateProblem{P, B, Q, V}
    parameters::P
    base::B
    quadrature::Q
    targets::V
end

function SurrogateProblem(st::FitSettings, spectra, quad)
    base = SpectralModel(spectra, quad)
    sq = SurfaceQuadrature(st.surface_divisions)
    return SurrogateProblem(parameterization(st), base, sq, surface_targets(base, sq))
end

"""Physical push and pull terms, evaluated before choosing their relative weight."""
function objective_terms(problem::SurrogateProblem, x)
    base = problem.base
    K, S = _parameter_matrices_wl(problem.parameters, x, size(base.K, 1))
    push = zero(eltype(x))
    pull = zero(eltype(x))
    for (i, c) in enumerate(problem.quadrature.points)
        rgb = mix_rgb_params(K, S, base.k1, base.k2, base.quad, c)
        w = problem.quadrature.weights[i]
        push += w * cube_outside_penalty(rgb)
        pull += w * oklab_distance_squared(linear_srgb_to_oklab(rgb), problem.targets[i])
    end
    return push, pull
end

function _objective_value(problem, x, alpha)
    push, pull = objective_terms(problem, x)
    return push + alpha * pull
end

"""Build an Optimization.jl problem; `adtype` accepts any compatible AD backend."""
function optimization_problem(problem::SurrogateProblem, x, alpha; adtype = Optimization.AutoForwardDiff())
    f = Optimization.OptimizationFunction((u, p) -> _objective_value(problem, u, p), adtype)
    return Optimization.OptimizationProblem(f, copy(x), alpha; _parameter_bounds(problem.parameters, x)...)
end

# This signal is handled only at the stage boundary. Checking the objective
# also bounds backends that do not implement Optimization's maxtime keyword.
struct FitTimeLimit <: Exception end

_default_algorithm(::Val{:lbfgs}) = Optim.LBFGS()
_default_algorithm(::Val{:lbfgsb}) = OptimizationLBFGSB.LBFGSB()

function _finite_gradient!(G, f, x, lower)
    u = copy(x)
    for i in eachindex(x)
        h = cbrt(eps(Float64)) * max(1.0, abs(x[i]))
        lo = max(x[i] - h, lower)
        u[i] = x[i] + h
        hi_value = f(u)
        u[i] = lo
        G[i] = (hi_value - f(u)) / (x[i] + h - lo)
        u[i] = x[i]
    end
    return G
end
_gradient_floor(::SoftplusParameters) = -Inf
_gradient_floor(p::DirectParameters) = p.epsilon

function _solve_keywords(algorithm, tolerance)
    return (;)
end
_solve_keywords(::Optim.LBFGS, tolerance) = (g_tol = tolerance,)
_solve_keywords(::OptimizationLBFGSB.LBFGSB, tolerance) = (reltol = tolerance,)

"""Solve one continuation step with an arbitrary compatible Optimization.jl algorithm."""
function solve_fit_step(
        problem::SurrogateProblem, x, alpha, algorithm, st::FitSettings;
        deadline = Inf, adtype = nothing, solve_kwargs = (;)
    )
    nf, ng = Ref(0), Ref(0)
    best_x, best_f = copy(x), Ref(Inf)
    raw = u -> _objective_value(problem, u, alpha)
    f = function (u, p)
        time() <= deadline || throw(FitTimeLimit())
        nf[] += 1
        value = raw(u)
        if value isa AbstractFloat && isfinite(value) && value < best_f[]
            best_f[] = value
            copyto!(best_x, u)
        end
        return value
    end
    if adtype === nothing
        gc = ForwardDiff.GradientConfig(raw, x, ForwardDiff.Chunk{min(length(x), 16)}())
        grad! = function (G, u, p)
            time() <= deadline || throw(FitTimeLimit())
            ng[] += 1
            if st.derivative == :finite_difference
                _finite_gradient!(G, raw, u, _gradient_floor(problem.parameters))
            else
                ForwardDiff.gradient!(G, raw, u, gc)
            end
            return nothing
        end
        hess! = (H, u, p) -> ForwardDiff.hessian!(H, raw, u)
        optf = Optimization.OptimizationFunction(f; grad = grad!, hess = hess!)
    else
        optf = Optimization.OptimizationFunction(f, adtype)
    end
    prob = Optimization.OptimizationProblem(optf, copy(x); _parameter_bounds(problem.parameters, x)...)
    options = merge((maxiters = st.max_iterations,), _solve_keywords(algorithm, st.gradient_tolerance), solve_kwargs)
    try
        sol = Optimization.solve(prob, algorithm; options...)
        return Vector{Float64}(sol.u), Dict{String, Any}(
                "iterations" => sol.stats.iterations, "converged" => Optimization.SciMLBase.successful_retcode(sol),
                "termination" => string(sol.retcode), "function_evaluations" => adtype === nothing ? nf[] : sol.stats.fevals,
                "gradient_evaluations" => adtype === nothing ? ng[] : sol.stats.gevals, "objective" => Float64(sol.objective),
            )
    catch err
        err isa FitTimeLimit || rethrow()
        return best_x, Dict{String, Any}(
                "iterations" => 0, "converged" => false,
                "termination" => "TimeLimit", "function_evaluations" => nf[],
                "gradient_evaluations" => ng[], "objective" => best_f[]
            )
    end
end

"""Continuation schedule independent of optimizer and coefficient coordinates."""
function alpha_schedule(st::FitSettings)
    return [st.alpha_initial / 2.0^i for i in 0:st.alpha_halvings if st.alpha_initial / 2.0^i >= st.alpha_final]
end

"""
    fit_surrogates(cfg, spectra, quad; algorithm, adtype, solve_kwargs, log)

Fit physical surrogate coefficients through Optimization.jl. Pass any compatible
algorithm and AD backend; bounds follow the parameterization. The default
remains softplus with L-BFGS. Continuation and diagnostics are solver-independent.
"""
function fit_surrogates(
        cfg::AbstractDict, spectra::PigmentSpectra, quad::Quadrature;
        algorithm = nothing, adtype = nothing, solve_kwargs = (;), log = nothing
    )
    st = FitSettings(cfg)
    if algorithm === nothing
        st.parameterization == :direct && st.optimizer == :lbfgs && throw(ConfigError("direct coefficients require optimizer=lbfgsb or an explicit bounded algorithm"))
        algorithm = _default_algorithm(Val(st.optimizer))
    end
    problem = SurrogateProblem(st, spectra, quad)
    x = initial_parameters(problem.parameters, spectra)
    history = Dict{String, Any}[]
    started = time()
    deadline = started + st.time_limit_seconds
    for (step, alpha) in enumerate(alpha_schedule(st))
        x, entry = solve_fit_step(problem, x, alpha, algorithm, st; deadline, adtype, solve_kwargs)
        epush, epull = objective_terms(problem, x)
        merge!(
            entry, Dict(
                "step" => step, "alpha" => alpha, "Epush" => epush,
                "Epull" => epull, "seconds" => time() - started
            )
        )
        push!(history, entry)
        if log !== nothing
            println(
                log, "  alpha=", alpha, " Epush=", epush, " Epull=", epull,
                " termination=", entry["termination"]
            )
            flush(log)
        end
        epush <= st.push_tolerance && break
        entry["termination"] == "TimeLimit" && break
    end
    K, S = parameter_matrices(problem.parameters, x, length(spectra.wavelength))
    model = SpectralModel(problem.base; K, S)
    diagnostics = surrogate_diagnostics(model, problem.base, problem.quadrature, problem.targets)
    diagnostics["parameterization"] = string(st.parameterization)
    diagnostics["algorithm"] = string(algorithm)
    diagnostics["derivative"] = adtype === nothing ? string(st.derivative) : string(adtype)
    diagnostics["termination"] = isempty(history) ? "NotStarted" : last(history)["termination"]
    return SurrogateFit(K, S, x, history, diagnostics, problem.quadrature)
end

"""
    surrogate_diagnostics(model, base, sq, targets) -> Dict{String,Any}

Independent-ish checks of a fitted surrogate: the worst cube violation and
worst Oklab deviation over the fit quadrature, plus a denser surface sample
that the fit did not use, plus the quadrature weighting report.

Finite sampling is evidence, not a proof of continuous gamut containment.
"""
function surrogate_diagnostics(
        model::SpectralModel, base::SpectralModel, sq::SurfaceQuadrature, targets
    )
    denser = SurfaceQuadrature(40)
    push_fit = Epush(model, sq)
    pull_fit = Epull(model, sq, targets)
    denser_targets = surface_targets(base, denser)
    push_dense = Epush(model, denser)
    pull_dense = Epull(model, denser, denser_targets)
    worst_push = 0.0
    worst_pull = 0.0
    deviations = Float64[]
    @inbounds for (i, c) in enumerate(denser.points)
        rgb = mix_rgb(model, c)
        worst_push = max(worst_push, cube_outside_penalty(rgb))
        d = oklab_distance_squared(linear_srgb_to_oklab(rgb), denser_targets[i])
        worst_pull = max(worst_pull, d)
        push!(deviations, sqrt(d))
    end
    # Radial interior samples complement the boundary quadrature. Neither is
    # a proof of continuous containment.
    for cface in denser.points[1:8:end], depth in (0.1, 0.5, 0.9)
        c = (1 - depth) * cface + depth * SVector(0.25, 0.25, 0.25, 0.25)
        rgb = mix_rgb(model, c)
        worst_push = max(worst_push, cube_outside_penalty(rgb))
        d = oklab_distance_squared(linear_srgb_to_oklab(rgb), linear_srgb_to_oklab(mix_rgb(base, c)))
        worst_pull = max(worst_pull, d)
        push!(deviations, sqrt(d))
    end
    return Dict{String, Any}(
        "Epush_fit" => push_fit,
        "Epull_fit" => pull_fit,
        "Epush_dense" => push_dense,
        "Epull_dense" => pull_dense,
        "max_cube_violation_squared" => worst_push,
        "max_cube_violation" => sqrt(max(worst_push, 0.0)),
        "max_oklab_deviation" => sqrt(max(worst_pull, 0.0)),
        "oklab_deviation" => Dict(
            "mean" => sum(deviations) / length(deviations),
            "p99" => sort(deviations)[ceil(Int, 0.99 * length(deviations))],
            "max" => maximum(deviations), "samples" => length(deviations)
        ),
        "dense_divisions" => 40,
        "quadrature" => quadrature_report(model, sq),
    )
end
