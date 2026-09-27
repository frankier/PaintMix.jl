# Plan step 3: equations (1)-(7) in generic floating-point arithmetic.
#
# Everything here is written so a `ForwardDiff.Dual` flows through unchanged:
# no `Float64` conversion, no clipping, no value-dependent branching that
# changes the formula. The surrogate fit differentiates these functions; the
# inverse solver calls the analytic Jacobian in `mix_rgb_jacobian!`.

"""
    km_reflectance(q) -> R

Equation (2) in the cancellation-free form

    R = 1 / (1 + q + sqrt(q*q + 2*q)),   q = K / S,

which is algebraically identical to the paper's `1 + q - sqrt(q^2 + 2q)` but
does not lose precision when `q` is small.
"""
@inline function km_reflectance(q::T) where {T <: Real}
    return one(T) / (one(T) + q + sqrt(q * q + 2 * q))
end

"""
    saunderson(R, k1, k2) -> R′

Equation (6): the modified reflectance that accounts for surface reflection.
The paper's convention has no added specular term, so `kins` does not appear.
"""
@inline function saunderson(R::T, k1, k2) where {T <: Real}
    oneT = one(T)
    return ((oneT - k1) * (oneT - k2) * R) / (oneT - k2 * R)
end

"""
    mix_rgb(model, c) -> RGB

`mix(c)` of equation (7): equations (1) and (2), the Saunderson correction
(6), the tristimulus integrals (3)-(5) with the cached trapezoidal weights,
and the pinned XYZ-to-linear-sRGB matrix.

`c` is four concentrations. It must be non-negative and sum to one; callers
on the simplex boundary pass an exact zero. Nothing is clipped: mixtures of
the original pigments may legitimately produce channels outside `[0, 1]`.
"""
@inline function mix_rgb(
        m::SpectralModel{T}, c::AbstractVector{S}
    ) where {T <: AbstractFloat, S <: Real}
    return mix_rgb_params(m.K, m.S, m.k1, m.k2, m.quad, c)
end

"""
    mix_rgb_params(K, S, k1, k2, quad, c) -> RGB

`mix(c)` of equation (7) from raw pigment parameters. This is the form the
surrogate fit calls: during `ForwardDiff` it holds `Dual`-valued `K` and `S`
with a `Float64` quadrature, so it cannot build a `SpectralModel`.
"""
@inline function mix_rgb_params(
        K::AbstractMatrix, S::AbstractMatrix, k1, k2,
        q::Quadrature, c::AbstractVector{S2},
    ) where {S2 <: Real}
    Tacc = promote_type(S2, eltype(K), eltype(S), eltype(q.weights), typeof(k1), typeof(k2))
    W = length(q.weights)
    c1, c2, c3, c4 = c
    X = zero(Tacc)
    Y = zero(Tacc)
    Z = zero(Tacc)
    @inbounds for l in 1:W
        k = c1 * K[l, 1] + c2 * K[l, 2] + c3 * K[l, 3] + c4 * K[l, 4]
        s = c1 * S[l, 1] + c2 * S[l, 2] + c3 * S[l, 3] + c4 * S[l, 4]
        ratio = k / s
        r = km_reflectance(ratio)
        rc = saunderson(r, k1, k2)
        w = q.weights[l] * q.d65[l] * rc
        X += w * q.x_bar[l]
        Y += w * q.y_bar[l]
        Z += w * q.z_bar[l]
    end
    return (q.xyz_to_rgb * SVector(X, Y, Z)) * q.norm
end

# The docstring above documents the model form; the raw-parameter form is
# what the surrogate fit calls.

"""
    mix_rgb_simplex(model, c1, c2, c3) -> RGB

`mix_rgb` with the fourth concentration implied as `1 - c1 - c2 - c3`.
"""
@inline function mix_rgb_simplex(
        m::SpectralModel{T}, c1::S, c2::S, c3::S
    ) where {T <: AbstractFloat, S <: Real}
    return mix_rgb(m, SVector(c1, c2, c3, one(S) - c1 - c2 - c3))
end

"""
    mix_rgb_jacobian!(J, model, c)

The analytic `3 x 4` Jacobian `d mix(c) / d c` into the first `3 x 4` block of
`J` (which must have at least 3 rows and 4 columns). Exact derivatives, used
by the bulk inverse solver and as the reference for finite-difference checks
of the `ForwardDiff` path.
"""
function mix_rgb_jacobian!(
        J::AbstractMatrix{T}, m::SpectralModel{T}, c::AbstractVector{T}
    ) where {T <: AbstractFloat}
    size(J, 1) >= 3 && size(J, 2) >= 4 ||
        throw(ArgumentError("Jacobian buffer must be at least 3 x 4, got $(size(J))"))
    K = m.K
    Sc = m.S
    q = m.quad
    W = length(q.weights)
    c1, c2, c3, c4 = c
    k1 = m.k1
    k2 = m.k2
    @inbounds for i in 1:4
        J[1, i] = zero(T)
        J[2, i] = zero(T)
        J[3, i] = zero(T)
    end
    @inbounds for l in 1:W
        k = c1 * K[l, 1] + c2 * K[l, 2] + c3 * K[l, 3] + c4 * K[l, 4]
        s = c1 * Sc[l, 1] + c2 * Sc[l, 2] + c3 * Sc[l, 3] + c4 * Sc[l, 4]
        ratio = k / s
        root = sqrt(ratio * ratio + 2 * ratio)
        r = one(T) / (one(T) + ratio + root)
        rc = saunderson(r, k1, k2)
        # dR/dq = -R^2 * (1 + (q + 1) / sqrt(q^2 + 2q))
        dR_dq = -(r * r) * (one(T) + (ratio + one(T)) / root)
        dRc_dR = ((one(T) - k1) * (one(T) - k2)) / ((one(T) - k2 * r) * (one(T) - k2 * r))
        fac = q.weights[l] * dRc_dR * dR_dq
        xw = fac * q.x_bar[l] * q.d65[l]
        yw = fac * q.y_bar[l] * q.d65[l]
        zw = fac * q.z_bar[l] * q.d65[l]
        inv_s = one(T) / s
        for i in 1:4
            dq = (K[l, i] - ratio * Sc[l, i]) * inv_s
            J[1, i] += xw * dq
            J[2, i] += yw * dq
            J[3, i] += zw * dq
        end
    end
    @inbounds for i in 1:4
        rgb = (q.xyz_to_rgb * SVector(J[1, i], J[2, i], J[3, i])) * q.norm
        J[1, i] = rgb[1]
        J[2, i] = rgb[2]
        J[3, i] = rgb[3]
    end
    return J
end

# --- perceptual difference and the cube penalty ---------------------------

# Oklab (Ottosson 2020), the color space equation (17) measures distances in.
const _OKLAB_M1 = @SMatrix [
    0.4122214708 0.5363325363 0.0514459929
    0.2119034982 0.6806995451 0.1073969566
    0.0883024619 0.2817188376 0.6299787005
]
const _OKLAB_M2 = @SMatrix [
    0.2104542553 0.793617785 -0.0040720468
    1.9779984951 -2.428592205 0.4505937099
    0.0259040371 0.7827717662 -0.808675766
]

"""
    linear_srgb_to_oklab(rgb) -> RGB

Convert linear-light sRGB to Oklab. Accepts out-of-gamut values; the LMS
cube root is signed so that negatives propagate.
"""
@inline function linear_srgb_to_oklab(rgb::SVector{3, S}) where {S <: Real}
    # `cbrt` accepts negative arguments, which out-of-gamut linear RGB needs.
    # Its derivative is singular at zero; that is inherent, and the tests
    # check gradients away from it.
    return _OKLAB_M2 * cbrt.(_OKLAB_M1 * rgb)
end

"""
    oklab_distance_squared(a, b) -> Real

Squared Euclidean Oklab distance, the integrand of equation (17).
"""
@inline oklab_distance_squared(a::AbstractVector, b::AbstractVector) = sum(abs2, a - b)

"""
    cube_signed_distance(p) -> Real

Signed distance from `p` to the surface of the unit cube: negative inside,
positive outside. This is the paper's `phi` in equation (16), exposed for
diagnostics.

Inside, the distance to the closest face is `-min_i min(p_i, 1 - p_i)`.
Outside it is the Euclidean distance to the cube, computed as a square root.
The fit itself never calls this: [`cube_outside_penalty`](@ref) is the
squared form and avoids the root.
"""
function cube_signed_distance(p::SVector{3, S}) where {S <: Real}
    d = sum(abs2, max.(p .- one(S), zero(S))) + sum(abs2, max.(-p, zero(S)))
    d == zero(S) && return -minimum(min.(p, one(S) .- p))
    return sqrt(d)
end

"""
    cube_outside_penalty(p) -> Real

`max(0, phi(p))^2`: the squared Euclidean distance from `p` to the unit cube,
or zero when `p` is inside. This is the integrand of equation (16).

It is written as a sum of squared per-channel violations, which is exactly
`phi(p)^2` outside the cube but needs no square root and stays smooth for
automatic differentiation.
"""
@inline function cube_outside_penalty(p::SVector{3, S}) where {S <: Real}
    return sum(abs2, max.(p .- one(S), zero(S))) + sum(abs2, max.(-p, zero(S)))
end

# --- evaluator interface ---------------------------------------------------
#
# The inverse solver works against either the spectral model or the float
# forward table. Both expose the same two operations, and the solver is
# generic over the evaluator type so no dispatch happens inside the hot loop.

"""
    _mix(ev, c) -> RGB

Evaluate the forward map at four concentrations.
"""
@inline _mix(m::SpectralModel, c::AbstractVector{S}) where {S <: Real} = mix_rgb(m, c)

"""
    _jacobian4!(J, ev, c)

Fill `J` with the `3 x 4` derivative of the forward map at `c`, in the
*simplex-constrained* sense the softmax chain rule needs: the derivative
along increasing `c[i]` with the fourth concentration taking up the slack is
`J[:, i] - J[:, 4]`. Callers only ever combine columns with directions whose
entries sum to zero, so this convention is exact for both evaluators.
"""
@inline function _jacobian4!(
        J::AbstractMatrix, m::SpectralModel, c::AbstractVector{T}
    ) where {T <: Real}
    return mix_rgb_jacobian!(J, m, c)
end
