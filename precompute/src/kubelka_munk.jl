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
    mix_rgb(model, c) -> NTuple{3}

`mix(c)` of equation (7): equations (1) and (2), the Saunderson correction
(6), the tristimulus integrals (3)-(5) with the cached trapezoidal weights,
and the pinned XYZ-to-linear-sRGB matrix.

`c` is four concentrations. It must be non-negative and sum to one; callers
on the simplex boundary pass an exact zero. Nothing is clipped: mixtures of
the original pigments may legitimately produce channels outside `[0, 1]`.
"""
@inline function mix_rgb(m::SpectralModel{T}, c::NTuple{4, S}) where {T <: AbstractFloat, S <: Real}
    return mix_rgb_params(m.K, m.S, m.k1, m.k2, m.quad, c)
end

"""
    mix_rgb_params(K, S, k1, k2, quad, c) -> NTuple{3}

`mix(c)` of equation (7) from raw pigment parameters. This is the form the
surrogate fit calls: during `ForwardDiff` it holds `Dual`-valued `K` and `S`
with a `Float64` quadrature, so it cannot build a `SpectralModel`.
"""
@inline function mix_rgb_params(
        K::AbstractMatrix, S::AbstractMatrix, k1, k2,
        q::Quadrature, c::NTuple{4, S2},
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
    mt = q.xyz_to_rgb
    n = q.norm
    return (
        (mt[1] * X + mt[2] * Y + mt[3] * Z) * n,
        (mt[4] * X + mt[5] * Y + mt[6] * Z) * n,
        (mt[7] * X + mt[8] * Y + mt[9] * Z) * n,
    )
end

# The docstring above documents the model form; the raw-parameter form is
# what the surrogate fit calls.

"""
    mix_rgb_simplex(model, c1, c2, c3) -> NTuple{3}

`mix_rgb` with the fourth concentration implied as `1 - c1 - c2 - c3`.
"""
@inline function mix_rgb_simplex(
        m::SpectralModel{T}, c1::S, c2::S, c3::S
    ) where {T <: AbstractFloat, S <: Real}
    return mix_rgb(m, (c1, c2, c3, one(S) - c1 - c2 - c3))
end

"""
    mix_rgb_jacobian!(J, model, c)

The analytic `3 x 4` Jacobian `d mix(c) / d c` into the first `3 x 4` block of
`J` (which must have at least 3 rows and 4 columns). Exact derivatives, used
by the bulk inverse solver and as the reference for finite-difference checks
of the `ForwardDiff` path.
"""
function mix_rgb_jacobian!(
        J::AbstractMatrix{T}, m::SpectralModel{T}, c::NTuple{4, T}
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
    mt = q.xyz_to_rgb
    n = q.norm
    @inbounds for i in 1:4
        x, y, z = J[1, i], J[2, i], J[3, i]
        J[1, i] = (mt[1] * x + mt[2] * y + mt[3] * z) * n
        J[2, i] = (mt[4] * x + mt[5] * y + mt[6] * z) * n
        J[3, i] = (mt[7] * x + mt[8] * y + mt[9] * z) * n
    end
    return J
end

# --- perceptual difference and the cube penalty ---------------------------

# Oklab (Ottosson 2020), the color space equation (17) measures distances in.
const _OKLAB_M1 = (
    0.4122214708, 0.5363325363, 0.0514459929,
    0.2119034982, 0.6806995451, 0.1073969566,
    0.0883024619, 0.2817188376, 0.6299787005,
)
const _OKLAB_M2 = (
    0.2104542553, 0.793617785, -0.0040720468,
    1.9779984951, -2.428592205, 0.4505937099,
    0.0259040371, 0.7827717662, -0.808675766,
)

"""
    linear_srgb_to_oklab(rgb) -> NTuple{3}

Convert linear-light sRGB to Oklab. Accepts out-of-gamut values; the LMS
cube root is signed so that negatives propagate.
"""
@inline function linear_srgb_to_oklab(rgb::NTuple{3, S}) where {S <: Real}
    r, g, b = rgb
    m = _OKLAB_M1
    l = m[1] * r + m[2] * g + m[3] * b
    mm = m[4] * r + m[5] * g + m[6] * b
    s = m[7] * r + m[8] * g + m[9] * b
    # `cbrt` accepts negative arguments, which out-of-gamut linear RGB needs.
    # Its derivative is singular at zero; that is inherent, and the tests
    # check gradients away from it.
    l_ = cbrt(l)
    m_ = cbrt(mm)
    s_ = cbrt(s)
    n = _OKLAB_M2
    return (
        n[1] * l_ + n[2] * m_ + n[3] * s_,
        n[4] * l_ + n[5] * m_ + n[6] * s_,
        n[7] * l_ + n[8] * m_ + n[9] * s_,
    )
end

"""
    oklab_distance_squared(a, b) -> Real

Squared Euclidean Oklab distance, the integrand of equation (17).
"""
@inline function oklab_distance_squared(a::NTuple{3, S}, b::NTuple{3, S}) where {S <: Real}
    d1 = a[1] - b[1]
    d2 = a[2] - b[2]
    d3 = a[3] - b[3]
    return d1 * d1 + d2 * d2 + d3 * d3
end

@inline function oklab_distance_squared(a::NTuple{3, S}, b::NTuple{3, S2}) where {S <: Real, S2 <: Real}
    d1 = a[1] - b[1]
    d2 = a[2] - b[2]
    d3 = a[3] - b[3]
    return d1 * d1 + d2 * d2 + d3 * d3
end

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
function cube_signed_distance(p::NTuple{3, S}) where {S <: Real}
    inside = true
    d = zero(S)
    @inbounds for i in 1:3
        v = p[i]
        if v < zero(S)
            inside = false
            d += v * v
        elseif v > one(S)
            inside = false
            t = v - one(S)
            d += t * t
        end
    end
    if inside
        return -min(min(p[1], one(S) - p[1]), min(p[2], one(S) - p[2]), min(p[3], one(S) - p[3]))
    end
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
@inline function cube_outside_penalty(p::NTuple{3, S}) where {S <: Real}
    t1 = p[1] - one(S)
    t2 = p[2] - one(S)
    t3 = p[3] - one(S)
    lo1 = -p[1]
    lo2 = -p[2]
    lo3 = -p[3]
    a = max(t1, zero(S))
    b = max(t2, zero(S))
    c = max(t3, zero(S))
    d = max(lo1, zero(S))
    e = max(lo2, zero(S))
    f = max(lo3, zero(S))
    return a * a + b * b + c * c + d * d + e * e + f * f
end

# --- evaluator interface ---------------------------------------------------
#
# The inverse solver works against either the spectral model or the float
# forward table. Both expose the same two operations, and the solver is
# generic over the evaluator type so no dispatch happens inside the hot loop.

"""
    _mix(ev, c) -> NTuple{3}

Evaluate the forward map at four concentrations.
"""
@inline _mix(m::SpectralModel, c::NTuple{4, S}) where {S <: Real} = mix_rgb(m, c)

"""
    _jacobian4!(J, ev, c)

Fill `J` with the `3 x 4` derivative of the forward map at `c`, in the
*simplex-constrained* sense the softmax chain rule needs: the derivative
along increasing `c[i]` with the fourth concentration taking up the slack is
`J[:, i] - J[:, 4]`. Callers only ever combine columns with directions whose
entries sum to zero, so this convention is exact for both evaluators.
"""
@inline function _jacobian4!(J::AbstractMatrix, m::SpectralModel, c::NTuple{4, T}) where {T <: Real}
    return mix_rgb_jacobian!(J, m, c)
end
