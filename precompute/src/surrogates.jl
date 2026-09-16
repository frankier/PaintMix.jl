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
struct SurfaceQuadrature{T<:AbstractFloat}
    points::Vector{NTuple{4,T}}
    weights::Vector{T}
    face::Vector{Int}
end

"""
    surface_quadrature(cfg_or_divisions) -> SurfaceQuadrature

Build the equally spaced surface quadrature. Accepts either the configuration
dictionary (reading `surrogate.surface_divisions`) or the division count.
"""
surface_quadrature(divisions::Integer) = _surface_quadrature(Int(divisions))
surface_quadrature(cfg::AbstractDict) =
    _surface_quadrature(Int(get(cfg["surrogate"], "surface_divisions", 20)))

function _surface_quadrature(d::Int)
    d >= 1 || throw(ArgumentError("surface quadrature needs at least one division"))
    points = NTuple{4,Float64}[]
    faces = Int[]
    step = 1 / d
    for f in 1:4
        free = [i for i in 1:4 if i != f]
        for a in 0:d, b in 0:(d - a)
            c = zeros(4)
            c[free[1]] = a * step
            c[free[2]] = b * step
            c[free[3]] = (d - a - b) * step
            push!(points, (c[1], c[2], c[3], c[4]))
            push!(faces, f)
        end
    end
    n = length(points)
    weights = fill(1 / n, n)
    return SurfaceQuadrature(points, weights, faces)
end

"""
    surface_colors(model, sq) -> Vector{NTuple{3,T}}

`mix_Q(c)` at every quadrature sample. Used for diagnostics and for the
`E_push`/`E_pull` reference implementations.
"""
function surface_colors(model::SpectralModel{T}, sq::SurfaceQuadrature) where {T}
    return [mix_rgb(model, c) for c in sq.points]
end

"""
    surface_targets(model, sq) -> Vector{NTuple{3,Float64}}

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
        u = (J[:, free[1]] .- J[:, free[3]])
        v = (J[:, free[2]] .- J[:, free[3]])
        cr = (
            u[2] * v[3] - u[3] * v[2],
            u[3] * v[1] - u[1] * v[3],
            u[1] * v[2] - u[2] * v[1],
        )
        out[i] = sqrt(cr[1]^2 + cr[2]^2 + cr[3]^2)
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
    ratios = [total_j > 0 ? (uni[i] * total_j) / (jac[i] * total_u) : NaN for i in eachindex(uni)]
    nu = count(!isnan, ratios)
    return Dict{String,Any}(
        "samples" => length(uni),
        "weighting" => "concentration-simplex-area",
        "rgb_surface_jacobian" => Dict{String,Any}(
            "total" => total_j,
            "mean" => total_j / length(jac),
            "max" => maximum(jac),
            "min" => minimum(jac),
        ),
        "normalized_ratio" => Dict{String,Any}(
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
@inline function softplus(x::T) where {T<:Real}
    return max(x, zero(T)) + log1p(exp(-abs(x)))
end

"""
    inv_softplus(y)

The inverse of [`softplus`](@ref): `log(expm1(y))`, with `y` floored at a
small positive value so the initial parameters are finite even when a
measured coefficient is below the configured floor.
"""
@inline function inv_softplus(y::T) where {T<:Real}
    return log(expm1(max(y, T(1.0e-12))))
end

@inline _unpack_absorb(θ::Real, epsilon::Real) = epsilon + softplus(θ)
@inline _unpack_scatter(θ::Real, epsilon::Real) = epsilon + softplus(θ)

"""
    SurrogateFit{T}

The fitted surrogates plus the record of how they were obtained.

  * `K`, `S`: `4 x W` pigment-major fitted coefficients.
  * `theta`: the unconstrained parameter vector that produced them.
  * `history`: one entry per continuation step, with `alpha`, objective
    terms, iteration count, and convergence flag.
  * `diagnostics`: worst violations and perceptual deviations.
  * `quadrature`: the [`SurfaceQuadrature`](@ref) used.
"""
struct SurrogateFit{T<:AbstractFloat}
    K::Matrix{T}
    S::Matrix{T}
    theta::Vector{T}
    history::Vector{Dict{String,Any}}
    diagnostics::Dict{String,Any}
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
        K[i, j] = _unpack_absorb(theta[(i - 1) * W + j], epsilon)
        S[i, j] = _unpack_scatter(theta[4W + (i - 1) * W + j], epsilon)
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
        K[j, i] = _unpack_absorb(theta[(i - 1) * W + j], epsilon)
        S[j, i] = _unpack_scatter(theta[4W + (i - 1) * W + j], epsilon)
    end
    return K, S
end

"""
    surrogate_objective(theta, base, sq, targets, epsilon, alpha) -> Real

`E_push(Q(theta)) + alpha * E_pull(Q(theta), P*)`, the objective of equation
(15) for one continuation step. Accepts `ForwardDiff.Dual` elements in
`theta`.
"""
function surrogate_objective(
        theta::AbstractVector, base::SpectralModel, sq::SurfaceQuadrature,
        targets, epsilon::Real, alpha::Real, scratch::AbstractMatrix,
    )
    W = size(base.K, 1)
    K, S = theta_parameters_wl(theta, W, epsilon)
    push = 0.0
    pull = 0.0
    @inbounds for (i, c) in enumerate(sq.points)
        rgb = mix_rgb_params(K, S, base.k1, base.k2, base.quad, c)
        push += sq.weights[i] * cube_outside_penalty(rgb)
        pull += sq.weights[i] * oklab_distance_squared(linear_srgb_to_oklab(rgb), targets[i])
    end
    return push + alpha * pull
end

# Without a scratch matrix the objective allocates the two `4 x W` parameter
# matrices per call. ForwardDiff's gradient calls it once per chunk, so that
# is acceptable; the extra argument exists for the finite-difference checks.

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

"""
    fit_surrogates(cfg, spectra, quad; log = nothing) -> SurrogateFit

Solve equation (15) with a warm-started `Optim.LBFGS` over the transformed
parameters, halving `alpha` from `alpha_initial` towards `alpha_final`.

Returns the fitted surrogates together with the convergence history and
diagnostics. The caller is responsible for caching and for deciding whether
the diagnostics pass the acceptance gates; this function reports rather than
promotes.
"""
function fit_surrogates(
        cfg::AbstractDict, spectra::PigmentSpectra{T}, quad::Quadrature{T};
        log = nothing,
    ) where {T}
    validate_config(cfg)
    s = cfg["surrogate"]
    epsilon = Float64(s["epsilon"])
    max_iter = Int(s["max_iterations"])
    time_limit = Float64(get(s, "time_limit_seconds", 3600.0))
    sq = surface_quadrature(cfg)
    base = spectral_model(spectra, quad)
    targets = surface_targets(base, sq)

    θ = initial_theta(spectra, epsilon)
    history = Vector{Dict{String,Any}}()
    W = size(base.K, 1)
    push_tol = Float64(get(s, "push_tolerance", 1.0e-8))

    started = time()
    for (step, α) in enumerate(alpha_schedule(cfg))
        f = x -> surrogate_objective(x, base, sq, targets, epsilon, α, Matrix{T}(undef, 4, W))
        # The ForwardDiff tag is derived from the function object, so the
        # gradient configuration must be rebuilt whenever the closure changes
        # (here, when alpha does).
        gcfg = ForwardDiff.GradientConfig(f, θ, ForwardDiff.Chunk{min(length(θ), 16)}())
        g! = (G, x) -> ForwardDiff.gradient!(G, f, x, gcfg)
        remaining = time_limit - (time() - started)
        remaining <= 0 && break
        opts = Optim.Options(
            iterations = max_iter,
            time_limit = remaining,
            store_trace = false,
            show_trace = false,
        )
        result = Optim.optimize(f, g!, θ, Optim.LBFGS(), opts)
        θ = Vector{Float64}(Optim.minimizer(result))
        K, S = theta_parameters(θ, W, epsilon)
        model = with_parameters(base, K, S)
        entry = Dict{String,Any}(
            "step" => step,
            "alpha" => α,
            "iterations" => Optim.iterations(result),
            "converged" => Optim.converged(result),
            "objective" => Optim.minimum(result),
            "Epush" => Epush(model, sq),
            "Epull" => Epull(model, sq, targets),
            "seconds" => time() - started,
        )
        push!(history, entry)
        if log !== nothing
            println(
                log,
                @sprintf(
                    "  alpha=%-10.4g Epush=%-12.6e Epull=%-12.6e iters=%-4d converged=%s",
                    α, entry["Epush"], entry["Epull"], entry["iterations"], entry["converged"]
                ),
            )
            flush(log)
        end
        entry["Epush"] <= push_tol && break
    end

    K, S = theta_parameters(θ, W, epsilon)
    model = with_parameters(base, K, S)
    diagnostics = surrogate_diagnostics(model, base, sq, targets)
    return SurrogateFit(K, S, θ, history, diagnostics, sq)
end

"""
    surrogate_diagnostics(model, base, sq, targets) -> Dict{String,Any}

Independent-ish checks of a fitted surrogate: the worst cube violation and
worst Oklab deviation over the fit quadrature, plus a denser surface sample
that the fit did not use, plus the quadrature weighting report.

Finite sampling is evidence, not a proof of continuous gamut containment.
"""
function surrogate_diagnostics(model::SpectralModel, base::SpectralModel, sq::SurfaceQuadrature, targets)
    denser = _surface_quadrature(40)
    push_fit = Epush(model, sq)
    pull_fit = Epull(model, sq, targets)
    denser_targets = surface_targets(base, denser)
    push_dense = Epush(model, denser)
    pull_dense = Epull(model, denser, denser_targets)
    worst_push = 0.0
    worst_pull = 0.0
    @inbounds for (i, c) in enumerate(denser.points)
        rgb = mix_rgb(model, c)
        worst_push = max(worst_push, cube_outside_penalty(rgb))
        d = oklab_distance_squared(linear_srgb_to_oklab(rgb), denser_targets[i])
        worst_pull = max(worst_pull, d)
    end
    return Dict{String,Any}(
        "Epush_fit" => push_fit,
        "Epull_fit" => pull_fit,
        "Epush_dense" => push_dense,
        "Epull_dense" => pull_dense,
        "max_cube_violation_squared" => worst_push,
        "max_cube_violation" => sqrt(max(worst_push, 0.0)),
        "max_oklab_deviation" => sqrt(max(worst_pull, 0.0)),
        "dense_divisions" => 40,
        "quadrature" => quadrature_report(model, sq),
    )
end
