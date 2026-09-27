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
struct UnmixSettings{T <: AbstractFloat}
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
struct UnmixResult{T <: AbstractFloat}
    c::SVector{4, T}
    sse::T
    active::NTuple{4, Int}   # active indices, zero-padded
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
end

SolverScratch() =
    SolverScratch(zeros(3, 4), zeros(3, 3), zeros(3, 3), zeros(3), zeros(3))

# --- small helpers ---------------------------------------------------------

# Softmax over a static logit vector. Its length is in the type, so one
# method covers every active-set size the solver uses.
@inline function _softmax_k(logits::SVector{M, T}) where {M, T}
    e = exp.(logits .- maximum(logits))
    return e ./ sum(e)
end

# Damped normal matrix `A + lambda * diag(max(A_ii, floor))`, returned as an
# `SMatrix` so the caller can solve it with `\`. `lambda` scales the
# diagonal only; the floor keeps a zero diagonal entry from making the
# system singular.
@inline _damp_diag(a::T, λ::T) where {T} = a + λ * max(a, T(1.0e-12))

@inline function _damped(sc::SolverScratch, λ::T, ::Val{M}) where {M, T}
    a = sc.A
    return SMatrix{M, M, T}(
        ntuple(Val(M * M)) do k
            # The flat index is column-major, matching `SMatrix` storage.
            j, i = divrem(k - 1, M)
            i += 1
            j += 1
            i == j ? _damp_diag(a[i, j], λ) : a[i, j]
        end
    )
end

# Subset probabilities (length K, last logit fixed to zero) -> 4-vector.
@inline function _c_from_active(active::NTuple{K, Int}, p::SVector{K, T}) where {K, T}
    pos = ntuple(i -> _findpos(active, i), Val(4))
    return SVector{4, T}(
        ntuple(Val(4)) do i
            pos[i] == 0 ? zero(T) : p[pos[i]]
        end
    )
end

@inline function _findpos(active::NTuple{K, Int}, i::Int) where {K}
    @inbounds for j in 1:K
        active[j] == i && return j
    end
    return 0
end

# Concentrations on a subset -> free logits (last fixed to zero).
function _theta_from_c(active::NTuple{K, Int}, c::SVector{4, T}) where {K, T}
    last = c[active[K]]
    return SVector{K - 1, T}(
        ntuple(Val(K - 1)) do i
            ci = c[active[i]]
            ci <= zero(T) && return T(-30)
            return log(max(ci, T(1.0e-300)) / max(last, T(1.0e-300)))
        end
    )
end

const Active4 = NTuple{4, Int}

# Allocation-free active-set packing: returns `(n, active)` with the active
# indices in the first `n` slots and zeros after. The fixed width is
# deliberate: a `Tuple{Int}` / `NTuple{2,Int}` / ... chain infers as
# `Tuple{Vararg{Int}}`, which made every active set a boxed value and cost
# `unmix_bulk!` about 96 bytes per grid vertex. A fixed-width tuple plus a
# separate count allocates nothing.
@inline function _active_pack(c::AbstractVector{T}, threshold::T) where {T}
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
        return 1, (m, 0, 0, 0)
    end
    return n, (i1, i2, i3, i4)
end

@inline function _active_count(a::Active4)
    n = 0
    @inbounds for i in 1:4
        n += a[i] != 0
    end
    return n
end

# Drop the pigment in slot `j` of an `n`-element active set and shift the
# remaining indices left. The zero padding absorbs the shifted hole, so no
# arity-specific methods are needed.
@inline function _drop_active(a::Active4, n::Int, j::Int)
    return ntuple(i -> i < n ? a[i + (i >= j)] : 0, Val(4))
end

@inline function _sse(rgb::SVector{3, T}, m::SVector{3, S}) where {T, S}
    return sum(abs2, m - rgb)
end

# Solve the damped normal equations. The dimension is at most three (the
# concentration simplex is 3-dimensional), so StaticArrays' unrolled,
# partially pivoted LU is allocation-free and measured about ten times
# faster per solve than the hand-rolled elimination it replaced. A singular
# system comes back as `nothing`: `\` yields `NaN`/`Inf` rather than
# throwing, so the finiteness check is what detects it, and the caller
# treats that as a rejected step.
@inline function _solve_damped(A::SMatrix{M, M, T}, g::SVector{M, T}) where {M, T}
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
        sc::SolverScratch, model, rgb::SVector{3, T},
        active::NTuple{K, Int}, theta0::SVector{M, T}, settings::UnmixSettings,
    ) where {T, K, M}
    @assert M == K - 1 "softmax needs K - 1 free logits"
    θ = theta0
    p = _softmax_k(SVector{M + 1, T}(θ..., zero(T)))
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
        δ = zero(SVector{M, T})
        for _ in 1:12
            A = _damped(sc, λ, Val(M))
            g = -SVector{M, T}(ntuple(i -> sc.g[i], Val(M)))
            x = _solve_damped(A, g)
            x === nothing && break
            δ = x
            θt = θ + δ
            pt = _softmax_k(SVector{M + 1, T}(θt..., zero(T)))
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
        if sqrt(sum(abs2, δ)) < settings.tolerance
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
        sc::SolverScratch, model, rgb::SVector{3, T},
        active::NTuple{K, Int}, c::SVector{4, T}, settings::UnmixSettings,
    ) where {T, K}
    θ = _theta_from_c(active, c)
    return lm_active!(sc, model, rgb, active, θ, settings)
end

"""
    _solve_adaptive!(sc, model, rgb, seed, settings) ->
        (c, sse, n, active, iters, converged)

Solve from `seed`, then repeatedly drop the smallest concentration while it
is below `settings.face_threshold` and re-solve on the reduced face. This is
what keeps boundary optima from being approached through saturating softmax
logits. `active` is a zero-padded `NTuple{4,Int}` whose first `n` entries are
the active pigment indices.
"""
function _solve_adaptive!(
        sc::SolverScratch, model, rgb::SVector{3, T},
        seed::SVector{4, T}, settings::UnmixSettings,
    ) where {T}
    n, active = _active_pack(seed, settings.face_threshold)
    return _adaptive_face!(sc, model, rgb, active, n, seed, settings)
end

# The softmax kernel is static-length, so the padded active set is unpacked
# into a concrete `NTuple{n,Int}` here, one branch per face. This is the only
# place the arity is recovered.
@inline function _polish_face!(
        sc::SolverScratch, model, rgb::SVector{3, T},
        active::Active4, n::Int, c::SVector{4, T}, settings::UnmixSettings,
    ) where {T}
    n == 4 && return _polish!(sc, model, rgb, active, c, settings)
    n == 3 &&
        return _polish!(sc, model, rgb, (active[1], active[2], active[3]), c, settings)
    n == 2 && return _polish!(sc, model, rgb, (active[1], active[2]), c, settings)
    return _polish!(sc, model, rgb, (active[1],), c, settings)
end

# Repeatedly drop the smallest concentration while it is below
# `settings.face_threshold` and re-solve on the reduced face. The active set
# stays a fixed-width `(n, NTuple{4,Int})` pair throughout, so the loop is
# allocation-free without an arity recursion.
function _adaptive_face!(
        sc::SolverScratch, model, rgb::SVector{3, T},
        active::Active4, n::Int, c::SVector{4, T}, settings::UnmixSettings,
    ) where {T}
    iters = 0
    converged = false
    while true
        c, sse, k, conv = _polish_face!(sc, model, rgb, active, n, c, settings)
        iters += k
        converged |= conv
        n == 1 && return c, sse, 1, active, iters, converged
        smallest = 0
        smallest_val = T(Inf)
        @inbounds for i in 1:n
            if c[active[i]] < smallest_val
                smallest_val = c[active[i]]
                smallest = i
            end
        end
        smallest_val > settings.face_threshold &&
            return c, sse, n, active, iters, converged
        active = _drop_active(active, n, smallest)
        n -= 1
    end
    return
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
        model, rgb::SVector{3, T};
        settings::UnmixSettings{T} = UnmixSettings{T}(50, T(1.0e-10), T(1.0e-6), 4),
        scratch::SolverScratch = SolverScratch(),
    ) where {T}
    state = _ReferenceState{T}(
        SVector(T(0.25), T(0.25), T(0.25), T(0.25)), T(Inf), 0, false, 0,
    )
    _reference_pass!(state, scratch, model, rgb, settings, _all_subsets())
    nactive, active = _active_pack(state.c, zero(T))
    return UnmixResult(
        state.c, state.sse, active, nactive, state.iters,
        state.converged, state.restarts
    )
end

# Mutable accumulator for one reference solve. Keeping it in a struct rather
# than a closure lets the subset pass stay allocation-free.
mutable struct _ReferenceState{T <: AbstractFloat}
    c::SVector{4, T}
    sse::T
    iters::Int
    converged::Bool
    restarts::Int
end

# The subset list is a homogeneous `NTuple{15, NTuple{4,Int}}`, so a plain
# loop over it keeps the active set concrete. The early exit matters: once
# the interior solve is exact there is no reason to enumerate the faces.
function _reference_pass!(
        state::_ReferenceState, sc::SolverScratch, model, rgb::SVector{3, T},
        settings::UnmixSettings, subsets,
    ) where {T}
    for active in subsets
        state.sse < T(1.0e-18) && return nothing
        _subset_pass!(state, sc, model, rgb, active, settings)
    end
    return nothing
end

function _subset_pass!(
        state::_ReferenceState{T}, sc::SolverScratch, model, rgb::SVector{3, T},
        active::Active4, settings::UnmixSettings,
    ) where {T}
    n = _active_count(active)
    if n == 1
        c = SVector{4, T}(ntuple(i -> i == active[1] ? one(T) : zero(T), Val(4)))
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
            _uniform_on(active, n, T)
        elseif r == 2 && isfinite(state.sse)
            _project_onto(active, n, state.c)
        else
            _random_on(active, r)
        end
        state.restarts = max(state.restarts, r)
        c, s, it, conv = _polish_face!(sc, model, rgb, active, n, seed, settings)
        state.iters += it
        state.converged |= conv
        if s < state.sse
            state.c = c
            state.sse = s
        end
    end
    return nothing
end

@inline function _uniform_on(active::Active4, n::Int, ::Type{T}) where {T}
    w = one(T) / n
    return SVector{4, T}(ntuple(i -> _findpos(active, i) == 0 ? zero(T) : w, Val(4)))
end

function _project_onto(active::Active4, n::Int, c::SVector{4, T}) where {T}
    s = zero(T)
    @inbounds for i in 1:n
        s += c[active[i]]
    end
    s <= 0 && return _uniform_on(active, n, T)
    # `inv_s` is assigned once, so the closure below captures an immutable
    # binding. Capturing `s` itself would box it, because `s` is mutated by
    # the loop above.
    inv_s = inv(s)
    return SVector{4, T}(
        ntuple(Val(4)) do i
            _findpos(active, i) == 0 ? zero(T) : c[i] * inv_s
        end
    )
end

# Deterministic perturbation, not RNG: restart `r` shifts the uniform seed in
# a fixed pattern so runs are reproducible.
function _random_on(active::Active4, r::Int)
    T = Float64
    raw = SVector{4, T}(
        ntuple(Val(4)) do i
            j = _findpos(active, i)
            j == 0 && return zero(T)
            return T(0.5 + 0.37 * sin(12.9898 * (i + 78.233 * r)))
        end
    )
    return raw / sum(raw)
end

# The 15 nonempty subsets of {1,2,3,4}, zero padded so the whole list is one
# concrete type, in a fixed order. The full simplex comes first: for a point
# inside the gamut the interior solve is exact and the enumeration can stop
# immediately, which is what makes the accuracy path affordable. Faces
# follow, largest first.
function _all_subsets()
    return (
        (1, 2, 3, 4),
        (1, 2, 3, 0), (1, 2, 4, 0), (1, 3, 4, 0), (2, 3, 4, 0),
        (1, 2, 0, 0), (1, 3, 0, 0), (1, 4, 0, 0), (2, 3, 0, 0), (2, 4, 0, 0), (3, 4, 0, 0),
        (1, 0, 0, 0), (2, 0, 0, 0), (3, 0, 0, 0), (4, 0, 0, 0),
    )
end

"""
    unmix_bulk!(sc, model, rgb, seeds, settings) -> UnmixResult

The throughput path: solve from each seed in `seeds` with the adaptive
active-face reduction and keep the best objective. Seeds normally come from
already-solved neighbouring grid vertices; the uniform seed is the fallback.
"""
function unmix_bulk!(
        sc::SolverScratch, model, rgb::SVector{3, T},
        seeds::NTuple{N, SVector{4, T}}, settings::UnmixSettings{T},
    ) where {T, N}
    best_c = seeds[1]
    best_sse = T(Inf)
    best_n, best_active = _active_pack(seeds[1], settings.face_threshold)
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
