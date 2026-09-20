# Plan step 4: the RGB-to-concentration inverse, equation (9).
#
#     unmix(RGB) = arg min_c || mix(c) - RGB ||^2   s.t.  c >= 0, sum(c) = 1
#
# Two solvers share one Levenberg-Marquardt core:
#
#   * `unmix_reference` enumerates all 15 nonempty pigment subsets, so a true
#     boundary solution is found on the lower-dimensional face that contains
#     it. It is the accuracy reference and the definition of correct.
#   * `unmix_bulk` is the throughput path used to fill 256^3 table vertices.
#     It starts from neighbouring solutions and reduces the active face
#     adaptively. Its output is checked against the reference solver on
#     held-out colors.
#
# Both stay in the paper's RGB least-squares objective. Switching to Oklab
# would be a declared model variant.

"""
    UnmixSettings

Numerical settings for the inverse solver: iteration cap, stationarity
tolerance on the step, the concentration below which a pigment is considered
inactive, and how many deterministic restarts the reference solver makes.
"""
struct UnmixSettings{T<:AbstractFloat}
    max_iterations::Int
    tolerance::T
    face_threshold::T
    restarts::Int
end

function unmix_settings(cfg::AbstractDict)
    u = cfg["unmix"]
    T = Float64
    return UnmixSettings{T}(
        Int(u["max_iterations"]),
        T(get(u, "tolerance", 1.0e-10)),
        T(get(u, "face_threshold", 1.0e-6)),
        Int(u["restarts"]),
    )
end

"""
    UnmixResult{T}

One inverse solve: the concentrations, the squared RGB residual, the active
pigment indices, and how the solver finished.
"""
struct UnmixResult{T<:AbstractFloat}
    c::NTuple{4,T}
    sse::T
    active::NTuple{4,Int}   # active indices, zero-padded
    nactive::Int
    iterations::Int
    converged::Bool
    restarts::Int
end

"""
    SolverScratch

Reusable buffers for one solver instance. One per thread. The kernels are
allocation-free once a scratch exists; the public solvers are too, and the
precompute tests check that. That is what makes 256^3 solves of 16.7 million
points practical.

The linear algebra is deferred to StaticArrays (`SMatrix`/`SVector`) rather
than being hand-rolled: the systems are at most 3 x 3, so `A \\ g` is
unrolled and allocation-free. The buffers stay `Matrix`/`Vector` because they
are written once per iteration by the generic `_jacobian4!` kernels.
"""
mutable struct SolverScratch
    J::Matrix{Float64}     # 3 x 4 analytic dRGB/dc
    Jt::Matrix{Float64}    # 3 x 3 analytic dRGB/dtheta for the active face
    A::Matrix{Float64}     # 3 x 3 normal matrix
    r::Vector{Float64}     # 3 residual
    g::Vector{Float64}     # 3 gradient
    d::Vector{Float64}     # 3 step
    c::Vector{Float64}     # 4 concentrations
    ct::Vector{Float64}    # 4 trial concentrations
    rt::Vector{Float64}    # 3 trial residual
end

SolverScratch() = SolverScratch(
    zeros(3, 4), zeros(3, 3), zeros(3, 3), zeros(3), zeros(3), zeros(3),
    zeros(4), zeros(4), zeros(3),
)

# --- small helpers ---------------------------------------------------------

@inline function _softmax_k(logits::NTuple{1,T}) where {T}
    return (one(T),)
end

@inline function _softmax_k(logits::NTuple{2,T}) where {T}
    m = max(logits[1], logits[2])
    a = exp(logits[1] - m)
    b = exp(logits[2] - m)
    s = a + b
    return (a / s, b / s)
end

@inline function _softmax_k(logits::NTuple{3,T}) where {T}
    m = max(logits[1], max(logits[2], logits[3]))
    a = exp(logits[1] - m)
    b = exp(logits[2] - m)
    c = exp(logits[3] - m)
    s = a + b + c
    return (a / s, b / s, c / s)
end

@inline function _softmax_k(logits::NTuple{4,T}) where {T}
    m = max(max(logits[1], logits[2]), max(logits[3], logits[4]))
    a = exp(logits[1] - m)
    b = exp(logits[2] - m)
    c = exp(logits[3] - m)
    d = exp(logits[4] - m)
    s = a + b + c + d
    return (a / s, b / s, c / s, d / s)
end

@inline _zeros_tuple(::Val{1}, ::Type{T}) where {T} = (zero(T),)
@inline _zeros_tuple(::Val{2}, ::Type{T}) where {T} = (zero(T), zero(T))
@inline _zeros_tuple(::Val{3}, ::Type{T}) where {T} = (zero(T), zero(T), zero(T))

@inline _theta_plus(a::NTuple{1,T}, b::NTuple{1,T}) where {T} = (a[1] + b[1],)
@inline _theta_plus(a::NTuple{2,T}, b::NTuple{2,T}) where {T} =
    (a[1] + b[1], a[2] + b[2])
@inline _theta_plus(a::NTuple{3,T}, b::NTuple{3,T}) where {T} =
    (a[1] + b[1], a[2] + b[2], a[3] + b[3])

# Damped normal matrix `A + lambda * diag(max(A_ii, floor))`, returned as an
# `SMatrix` so the caller can solve it with `\`. `lambda` scales the
# diagonal only; the floor keeps a zero diagonal entry from making the
# system singular.
@inline _damp_diag(a::T, λ::T) where {T} = a + λ * max(a, T(1.0e-12))

@inline function _damped(sc::SolverScratch, λ::T, ::Val{1}) where {T}
    a = sc.A
    return SMatrix{1,1,T}(_damp_diag(a[1, 1], λ))
end

@inline function _damped(sc::SolverScratch, λ::T, ::Val{2}) where {T}
    a = sc.A
    return SMatrix{2,2,T}(
        _damp_diag(a[1, 1], λ), a[2, 1],
        a[1, 2], _damp_diag(a[2, 2], λ),
    )
end

@inline function _damped(sc::SolverScratch, λ::T, ::Val{3}) where {T}
    a = sc.A
    return SMatrix{3,3,T}(
        _damp_diag(a[1, 1], λ), a[2, 1], a[3, 1],
        a[1, 2], _damp_diag(a[2, 2], λ), a[3, 2],
        a[1, 3], a[2, 3], _damp_diag(a[3, 3], λ),
    )
end

# Subset probabilities (length K, last logit fixed to zero) -> 4-vector.
@inline function _c_from_active(active::NTuple{K,Int}, p::NTuple{K,T}) where {K,T}
    pos = ntuple(i -> _findpos(active, i), Val(4))
    return ntuple(Val(4)) do i
        pos[i] == 0 ? zero(T) : p[pos[i]]
    end
end

@inline function _findpos(active::NTuple{K,Int}, i::Int) where {K}
    @inbounds for j in 1:K
        active[j] == i && return j
    end
    return 0
end

# Concentrations on a subset -> free logits (last fixed to zero).
function _theta_from_c(active::NTuple{K,Int}, c::NTuple{4,T}) where {K,T}
    last = c[active[K]]
    return ntuple(Val(K - 1)) do i
        ci = c[active[i]]
        ci <= zero(T) && return T(-30)
        return log(max(ci, T(1.0e-300)) / max(last, T(1.0e-300)))
    end
end

# Allocation-free active-set packing: returns `(n, i1, i2, i3, i4)` with the
# active indices in `i1..in` and zeros elsewhere. Returning a fixed 5-tuple
# keeps this off the heap, which matters because it runs once per seed per
# table vertex.
@inline function _active_pack(c::NTuple{4,T}, threshold::T) where {T}
    n = 0
    i1 = 0
    i2 = 0
    i3 = 0
    i4 = 0
    @inbounds for i in 1:4
        if c[i] > threshold
            n += 1
            if n == 1
                i1 = i
            elseif n == 2
                i2 = i
            elseif n == 3
                i3 = i
            else
                i4 = i
            end
        end
    end
    if n == 0
        m = 1
        @inbounds for i in 2:4
            c[i] > c[m] && (m = i)
        end
        return 1, m, 0, 0, 0
    end
    return n, i1, i2, i3, i4
end

# Pack an active set into a concrete `NTuple{4,Int}`, zero padded. Keeping
# the arity out of the type is deliberate: a `Tuple{Int}` / `NTuple{2,Int}`
# / ... chain infers as `Tuple{Vararg{Int}}`, which made every active set a
# boxed value and cost `unmix_bulk!` about 96 bytes per grid vertex. A
# fixed-width tuple plus a separate count allocates nothing.
@inline function _pad4(active::NTuple{K,Int}) where {K}
    return ntuple(i -> i <= K ? active[i] : 0, Val(4))
end

@inline _drop_active(a::NTuple{2,Int}, j::Int) = j == 1 ? (a[2],) : (a[1],)
@inline _drop_active(a::NTuple{3,Int}, j::Int) =
    j == 1 ? (a[2], a[3]) : j == 2 ? (a[1], a[3]) : (a[1], a[2])
@inline _drop_active(a::NTuple{4,Int}, j::Int) =
    j == 1 ? (a[2], a[3], a[4]) :
    j == 2 ? (a[1], a[3], a[4]) :
    j == 3 ? (a[1], a[2], a[4]) : (a[1], a[2], a[3])

@inline function _sse(rgb::NTuple{3,T}, m::NTuple{3,S}) where {T,S}
    d1 = m[1] - rgb[1]
    d2 = m[2] - rgb[2]
    d3 = m[3] - rgb[3]
    return d1 * d1 + d2 * d2 + d3 * d3
end

# Solve the damped normal equations. The dimension is at most three (the
# concentration simplex is 3-dimensional), so StaticArrays' unrolled,
# partially pivoted LU is allocation-free and measured about ten times
# faster per solve than the hand-rolled elimination it replaced. A singular
# system comes back as `nothing`: `\` yields `NaN`/`Inf` rather than
# throwing, so the finiteness check is what detects it, and the caller
# treats that as a rejected step.
@inline function _solve_damped(A::SMatrix{M,M,T}, g::SVector{M,T}) where {M,T}
    if M == 1 && abs(A[1, 1]) < eps(T)
        return nothing
    end
    x = A \ g
    @inbounds for i in 1:M
        isfinite(x[i]) || return nothing
    end
    return x
end

# --- Levenberg-Marquardt on one active face --------------------------------

"""
    lm_active!(sc, model, rgb, active, theta0, settings) -> (c, sse, iters, converged)

Levenberg-Marquardt on the softmax coordinates of one pigment subset. `theta0`
has length `K - 1` where `K = length(active)`; the last logit is fixed to
zero, which removes softmax's additive gauge.

Returns the best concentrations found, the squared residual, the iteration
count, and whether the step stationarity test passed.
"""
function lm_active!(
        sc::SolverScratch, model, rgb::NTuple{3,T},
        active::NTuple{K,Int}, theta0::NTuple{M,T}, settings::UnmixSettings,
    ) where {T,K,M}
    @assert M == K - 1 "softmax needs K - 1 free logits"
    θ = theta0
    p = _softmax_k((θ..., zero(T)))
    c = _c_from_active(active, p)
    m = _mix(model, c)
    sse = _sse(rgb, m)
    K == 1 && return c, sse, 0, true

    λ = T(1.0e-3)
    iterations = 0
    converged = false

    while iterations < settings.max_iterations
        iterations += 1
        # Analytic dRGB/dc at the current point; the softmax chain rule is
        # applied on top of it below.
        _jacobian4!(sc.J, model, c)
        @inbounds for i in 1:3
            sc.r[i] = m[i] - rgb[i]
        end
        # Jtheta = Jc[:, active] * dc/dtheta, a 3 x M matrix.
        @inbounds for i in 1:M, row in 1:3
            sc.Jt[row, i] = zero(T)
        end
        @inbounds for i in 1:M
            for j in 1:K
                d = j <= M ? p[j] * ((i == j ? one(T) : zero(T)) - p[i]) :
                    -p[K] * p[i]
                d == zero(T) && continue
                col = active[j]
                for row in 1:3
                    sc.Jt[row, i] += sc.J[row, col] * d
                end
            end
        end
        # Normal equations: A = Jt'Jt, g = Jt'r.
        @inbounds for i in 1:M, j in 1:M
            s = zero(T)
            for row in 1:3
                s += sc.Jt[row, i] * sc.Jt[row, j]
            end
            sc.A[i, j] = s
        end
        @inbounds for i in 1:M
            s = zero(T)
            for row in 1:3
                s += sc.Jt[row, i] * sc.r[row]
            end
            sc.g[i] = s
        end

        accepted = false
        δ = _zeros_tuple(Val(M), T)
        for _ in 1:12
            A = _damped(sc, λ, Val(M))
            g = SVector{M,T}(ntuple(i -> -sc.g[i], Val(M)))
            x = _solve_damped(A, g)
            x === nothing && break
            δ = ntuple(i -> x[i], Val(M))
            θt = _theta_plus(θ, δ)
            pt = _softmax_k((θt..., zero(T)))
            ct = _c_from_active(active, pt)
            mt = _mix(model, ct)
            st = _sse(rgb, mt)
            if st < sse
                θ = θt
                p = pt
                c = ct
                m = mt
                sse = st
                λ = max(λ / 3, T(1.0e-12))
                accepted = true
                break
            end
            λ *= 10
        end
        if !accepted
            converged = true
            break
        end
        normδ = zero(T)
        @inbounds for i in 1:M
            normδ += δ[i] * δ[i]
        end
        if sqrt(normδ) < settings.tolerance
            converged = true
            break
        end
    end

    return c, sse, iterations, converged
end

"""
    _polish!(sc, model, rgb, active, c, settings) -> (c, sse, iters, converged)

Run LM from concentrations rather than logits, converting to the subset's
free logits first.
"""
function _polish!(
        sc::SolverScratch, model, rgb::NTuple{3,T},
        active::NTuple{K,Int}, c::NTuple{4,T}, settings::UnmixSettings,
    ) where {T,K}
    θ = _theta_from_c(active, c)
    return lm_active!(sc, model, rgb, active, θ, settings)
end

"""
    _solve_adaptive!(sc, model, rgb, seed, settings) -> (c, sse, n, active, iters, converged)

Solve from `seed`, then repeatedly drop the smallest concentration while it
is below `settings.face_threshold` and re-solve on the reduced face. This is
what keeps boundary optima from being approached through saturating softmax
logits. `active` is a zero-padded `NTuple{4,Int}` whose first `n` entries are
the active pigment indices.
"""
function _solve_adaptive!(
        sc::SolverScratch, model, rgb::NTuple{3,T},
        seed::NTuple{4,T}, settings::UnmixSettings,
    ) where {T}
    n, i1, i2, i3, i4 = _active_pack(seed, settings.face_threshold)
    if n == 4
        return _adaptive_face!(sc, model, rgb, (i1, i2, i3, i4), seed, settings)
    elseif n == 3
        return _adaptive_face!(sc, model, rgb, (i1, i2, i3), seed, settings)
    elseif n == 2
        return _adaptive_face!(sc, model, rgb, (i1, i2), seed, settings)
    end
    return _adaptive_face!(sc, model, rgb, (i1,), seed, settings)
end

# Recursion on the active-set size keeps every frame specialized on `K`, so
# the reduced face is a concrete tuple type rather than a union. The *return*
# value carries the active set as a fixed-width `(n, NTuple{4,Int})` pair,
# which is what keeps the boundary of the solver allocation-free.
function _adaptive_face!(
        sc::SolverScratch, model, rgb::NTuple{3,T},
        active::NTuple{K,Int}, c::NTuple{4,T}, settings::UnmixSettings,
    ) where {T,K}
    c, sse, iters, converged = _polish!(sc, model, rgb, active, c, settings)
    K == 1 && return c, sse, K, _pad4(active), iters, converged
    smallest = 0
    smallest_val = T(Inf)
    for (j, idx) in enumerate(active)
        if c[idx] < smallest_val
            smallest_val = c[idx]
            smallest = j
        end
    end
    smallest_val > settings.face_threshold &&
        return c, sse, K, _pad4(active), iters, converged
    reduced = _drop_active(active, smallest)
    c2, sse2, n2, a2, iters2, converged2 =
        _adaptive_face!(sc, model, rgb, reduced, c, settings)
    return c2, sse2, n2, a2, iters + iters2, converged | converged2
end

# --- public solvers --------------------------------------------------------

"""
    unmix_reference(model, rgb; settings) -> UnmixResult

Equation (9) by active-face enumeration: every one of the 15 nonempty pigment
subsets is optimized (size-one subsets are evaluated directly), and the best
objective wins. Seeds are deterministic: the uniform distribution on the
subset, the global uniform distribution projected onto it, and the best
solution found so far, for `settings.restarts` rounds.

This is the accuracy reference for [`unmix_bulk`](@ref). It is a local method
over a non-convex objective, so `converged` and the objective value are
reported rather than a global-optimality claim.
"""
function unmix_reference(
        model, rgb::NTuple{3,T};
        settings::UnmixSettings{T} = UnmixSettings{T}(50, T(1.0e-10), T(1.0e-6), 4),
        scratch::SolverScratch = SolverScratch(),
    ) where {T}
    state = _ReferenceState{T}(
        (T(0.25), T(0.25), T(0.25), T(0.25)), T(Inf), 0, false, 0,
    )
    _reference_pass!(state, scratch, model, rgb, settings, _all_subsets())
    nactive, a1, a2, a3, a4 = _active_pack(state.c, zero(T))
    return UnmixResult(state.c, state.sse, (a1, a2, a3, a4), nactive, state.iters,
        state.converged, state.restarts)
end

# Mutable accumulator for one reference solve. Keeping it in a struct rather
# than a closure lets the subset pass stay allocation-free.
mutable struct _ReferenceState{T<:AbstractFloat}
    c::NTuple{4,T}
    sse::T
    iters::Int
    converged::Bool
    restarts::Int
end

# Walk the fixed subset tuple one element at a time. The tuple is
# heterogeneous, so the element type has to stay in the signature: with an
# abstract `subsets::Tuple`, `first(subsets)` infers as `Any`, every active
# set is boxed, and a reference solve allocates some hundreds of bytes. With
# `S<:Tuple` each level specializes on what is left of the tuple, so the
# active set stays concrete. The early exit matters: once the interior solve
# is exact there is no reason to enumerate the faces.
function _reference_pass!(
        state::_ReferenceState, sc::SolverScratch, model, rgb::NTuple{3,T},
        settings::UnmixSettings, subsets::S,
    ) where {T,S<:Tuple}
    state.sse < T(1.0e-18) && return nothing
    _subset_pass!(state, sc, model, rgb, first(subsets), settings)
    return _reference_pass!(state, sc, model, rgb, settings, Base.tail(subsets))
end

_reference_pass!(
    state::_ReferenceState, sc::SolverScratch, model, rgb::NTuple{3,T},
    settings::UnmixSettings, ::Tuple{},
) where {T} = nothing

function _subset_pass!(
        state::_ReferenceState{T}, sc::SolverScratch, model, rgb::NTuple{3,T},
        active::NTuple{K,Int}, settings::UnmixSettings,
    ) where {T,K}
    if K == 1
        c = ntuple(i -> i == active[1] ? one(T) : zero(T), Val(4))
        s = _sse(rgb, _mix(model, c))
        if s < state.sse
            state.c = c
            state.sse = s
            state.converged = true
        end
        return nothing
    end
    for r in 1:max(settings.restarts, 1)
        seed = if r == 1
            _uniform_on(active, T)
        elseif r == 2 && isfinite(state.sse)
            _project_onto(active, state.c)
        else
            _random_on(active, r)
        end
        state.restarts = max(state.restarts, r)
        c, s, it, conv = _polish!(sc, model, rgb, active, seed, settings)
        state.iters += it
        state.converged |= conv
        if s < state.sse
            state.c = c
            state.sse = s
        end
    end
    return nothing
end

@inline function _uniform_on(active::NTuple{K,Int}, ::Type{T}) where {K,T}
    n = one(T) / K
    return ntuple(i -> _findpos(active, i) == 0 ? zero(T) : n, Val(4))
end

function _project_onto(active::NTuple{K,Int}, c::NTuple{4,T}) where {K,T}
    s = zero(T)
    @inbounds for idx in active
        s += c[idx]
    end
    s <= 0 && return _uniform_on(active, T)
    # `inv_s` is assigned once, so the closure below captures an immutable
    # binding. Capturing `s` itself would box it, because `s` is mutated by
    # the loop above.
    inv_s = inv(s)
    return ntuple(Val(4)) do i
        _findpos(active, i) == 0 ? zero(T) : c[i] * inv_s
    end
end

# Deterministic perturbation, not RNG: restart `r` shifts the uniform seed in
# a fixed pattern so runs are reproducible.
function _random_on(active::NTuple{K,Int}, r::Int) where {K}
    T = Float64
    raw = ntuple(Val(4)) do i
        j = _findpos(active, i)
        j == 0 && return zero(T)
        return T(0.5 + 0.37 * sin(12.9898 * (i + 78.233 * r)))
    end
    s = sum(raw)
    return ntuple(i -> raw[i] / s, Val(4))
end

# The 15 nonempty subsets of {1,2,3,4}, in a fixed order. The full simplex
# comes first: for a point inside the gamut the interior solve is exact and
# the enumeration can stop immediately, which is what makes the accuracy path
# affordable. Faces follow, largest first.
function _all_subsets()
    return (
        (1, 2, 3, 4),
        (1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4),
        (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),
        (1,), (2,), (3,), (4,),
    )
end

"""
    unmix_bulk!(sc, model, rgb, seeds, settings) -> UnmixResult

The throughput path: solve from each seed in `seeds` with the adaptive
active-face reduction and keep the best objective. Seeds normally come from
already-solved neighbouring grid vertices; the uniform seed is the fallback.
"""
function unmix_bulk!(
        sc::SolverScratch, model, rgb::NTuple{3,T},
        seeds::NTuple{N,NTuple{4,T}}, settings::UnmixSettings{T},
    ) where {T,N}
    best_c = seeds[1]
    best_sse = T(Inf)
    n0, a1, a2, a3, a4 = _active_pack(seeds[1], settings.face_threshold)
    best_active = (a1, a2, a3, a4)
    best_n = n0
    total_iters = 0
    converged = false
    for seed in seeds
        c, s, n, a, it, conv = _solve_adaptive!(sc, model, rgb, seed, settings)
        total_iters += it
        converged |= conv
        if s < best_sse
            best_c, best_sse = c, s
            best_active = a
            best_n = n
        end
    end
    return UnmixResult(best_c, best_sse, best_active, best_n, total_iters, converged, N)
end
