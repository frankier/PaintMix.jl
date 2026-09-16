# Encoding, decoding, and mixing.
#
# Two layers share these kernels:
#
#   * the Julia API (`encode`, `decode`, `mix`, `bulk_mix!`, `weighted_mix!`),
#     which validates and throws, and
#   * the status-code kernels (`bulk_mix_kernel!`, `weighted_mix_kernel!`)
#     used directly by the compiled C ABI, which never throw.
#
# Kernels assume a non-null, correctly sized buffer and finite coordinates;
# the wrappers check those conditions first. Nothing here allocates.

"""
    PM_OK

Status code returned by the bulk kernels on success.
"""
const PM_OK = Int32(0)

"""
    PM_ERR_NULL

Status code: a required pointer argument is null while its element count is
positive.
"""
const PM_ERR_NULL = Int32(1)

"""
    PM_ERR_LENGTH

Status code: a buffer is shorter than its element count requires.
"""
const PM_ERR_LENGTH = Int32(2)

"""
    PM_ERR_NONFINITE

Status code: a color channel, mixing fraction, or weight is not finite.
"""
const PM_ERR_NONFINITE = Int32(3)

"""
    PM_ERR_WEIGHT

Status code: a weight is negative.
"""
const PM_ERR_WEIGHT = Int32(4)

"""
    PM_ERR_TOTAL

Status code: the weights sum to zero, so the weighted average is undefined.
"""
const PM_ERR_TOTAL = Int32(5)

"""
    error_message(code) -> String

A short, fixed description of a status code. Used for foreign-callable error
messages, where the status has to become a bounded string without allocating
on the failure path.
"""
function error_message(code::Integer)
    c = Int32(code)
    c == PM_OK && return "ok"
    c == PM_ERR_NULL && return "null pointer argument"
    c == PM_ERR_LENGTH && return "buffer shorter than element count"
    c == PM_ERR_NONFINITE && return "non-finite input value"
    c == PM_ERR_WEIGHT && return "negative weight"
    c == PM_ERR_TOTAL && return "weights sum to zero"
    return "unknown status $(c)"
end

@inline function _check(rgb::RGB)
    if !(isfinite(rgb[1]) && isfinite(rgb[2]) && isfinite(rgb[3]))
        throw(DomainError(rgb, "color channels must be finite"))
    end
    return nothing
end

@inline function _check_fraction(t::Real)
    isfinite(t) || throw(DomainError(t, "mixing fraction must be finite"))
    return nothing
end

# Project a concentration triple back into the simplex. Interpolating valid
# grid vertices cannot leave the simplex in exact arithmetic, so this only
# repairs floating-point excursions; for a pathological input it is still a
# deterministic projection rather than an error.
@inline function _repair_simplex(c1::T, c2::T, c3::T) where {T<:AbstractFloat}
    a = max(c1, zero(T))
    b = max(c2, zero(T))
    c = max(c3, zero(T))
    s = a + b + c
    if s > one(T)
        f = one(T) / s
        return (a * f, b * f, c * f)
    end
    return (a, b, c)
end

"""
    forward_rgb(model, c) -> RGB{T}

Evaluate the forward table at concentrations `c`, returning unclipped
linear-light sRGB. Coordinates are repaired onto the simplex first. The
fourth concentration is implied; only the first three index the table.
"""
@inline function forward_rgb(
        model::PigmentModel, c1::T, c2::T, c3::T
    ) where {T<:AbstractFloat}
    r1, r2, r3 = _repair_simplex(c1, c2, c3)
    return trilinear(model.forward, r1, r2, r3)
end

@inline forward_rgb(model::PigmentModel, c::Concentrations{T}) where {T<:AbstractFloat} =
    forward_rgb(model, c[1], c[2], c[3])

"""
    encode(model, rgb) -> Latent{T}

Encode a linear-light sRGB color as four concentrations plus a residual.

This performs the inverse lookup, then evaluates the forward table at the
recovered concentrations and stores `rgb - M(c)`. The residual is therefore
measured against the same quantized, interpolated table that
[`decode`](@ref) uses, so the round trip cancels to floating-point accuracy
instead of losing the byte quantization.

`rgb` must be finite; a `DomainError` is thrown otherwise. The returned
latent belongs to `model` and must not be combined with latents from another
model.
"""
@inline function encode(model::PigmentModel, rgb::RGB{T}) where {T<:AbstractFloat}
    _check(rgb)
    x, y, z = rgb
    c1, c2, c3 = trilinear(model.inverse, x, y, z)
    c1, c2, c3 = _repair_simplex(c1, c2, c3)
    c4 = max(zero(T), one(T) - (c1 + c2 + c3))
    m = forward_rgb(model, c1, c2, c3)
    return Latent((c1, c2, c3, c4), (x - m[1], y - m[2], z - m[3]))
end

"""
    decode(model, latent) -> RGB{T}

Decode a latent back to linear-light sRGB: `M(c) + r`.

The result is not clipped. Residuals can push a channel outside `[0, 1]`;
clip explicitly with [`srgb8_from_linear`](@ref) or `clamp` when displaying.

The fourth concentration is reconstructed as `1 - c1 - c2 - c3` rather than
read from the latent, so that interpolation drift cannot move the point off
the simplex.
"""
@inline function decode(model::PigmentModel, z::Latent{T}) where {T<:AbstractFloat}
    c1, c2, c3, _ = z.c
    m = forward_rgb(model, c1, c2, c3)
    r = z.r
    return (m[1] + r[1], m[2] + r[2], m[3] + r[3])
end

@inline function _latent_lerp(a::Latent{T}, b::Latent{T}, t::T) where {T<:AbstractFloat}
    s = one(T) - t
    c = ntuple(i -> s * a.c[i] + t * b.c[i], Val(4))
    r = ntuple(i -> s * a.r[i] + t * b.r[i], Val(3))
    return Latent(c, r)
end

"""
    mix(model, a, b, t) -> RGB{T}

Mix two linear-light sRGB colors as paint: encode both, lerp the latents,
decode.

`t` is the share of `b`. Fractions outside `[0, 1]` clamp to the nearest
endpoint, and `t == 0` or `t == 1` returns `a` or `b` unchanged, so mixing
is exact in the binary endpoints. A `t` of a different real type is
converted to `T`, the color scalar type.
"""
function mix(
        model::PigmentModel, a::RGB{T}, b::RGB{T}, t::T
    ) where {T<:AbstractFloat}
    _check(a)
    _check(b)
    _check_fraction(t)
    t == zero(T) && return a
    t == one(T) && return b
    tc = clamp(t, zero(T), one(T))
    za = encode(model, a)
    zb = encode(model, b)
    return decode(model, _latent_lerp(za, zb, tc))
end

function mix(
        model::PigmentModel, a::RGB{T}, b::RGB{T}, t::Real
    ) where {T<:AbstractFloat}
    return mix(model, a, b, T(t))
end

"""
    mix!(dest, model, a, b, t) -> dest

Write [`mix`](@ref) into the first three elements of `dest`. The caller owns
`dest`, so bulk work can run without allocating a result per sample.
"""
function mix!(
        dest::AbstractVector{T}, model::PigmentModel, a::RGB{T}, b::RGB{T}, t::T
    ) where {T<:AbstractFloat}
    length(dest) >= 3 || throw(ArgumentError("dest must hold at least 3 elements"))
    c = mix(model, a, b, t)
    @inbounds begin
        dest[1] = c[1]
        dest[2] = c[2]
        dest[3] = c[3]
    end
    return dest
end

"""
    bulk_mix_kernel!(dest, model, as, bs, ts, count) -> Int32

Mix `count` color pairs without allocating or throwing. `as`, `bs`, and
`ts` are flat `3 * count`, `3 * count`, and `count` element buffers; `dest`
holds `3 * count` results in the same order.

All inputs are validated before anything is written, so a failed call leaves
`dest` untouched. Returns [`PM_OK`](@ref) or one of the `PM_ERR_*` codes.
"""
function bulk_mix_kernel!(
        dest::AbstractVector{T}, model::PigmentModel,
        as::AbstractVector{T}, bs::AbstractVector{T}, ts::AbstractVector{T},
        count::Integer
    ) where {T<:AbstractFloat}
    n = Int(count)
    n < 0 && return PM_ERR_LENGTH
    length(dest) >= 3n || return PM_ERR_LENGTH
    length(as) >= 3n || return PM_ERR_LENGTH
    length(bs) >= 3n || return PM_ERR_LENGTH
    length(ts) >= n || return PM_ERR_LENGTH
    @inbounds for i in 1:(3n)
        (isfinite(as[i]) && isfinite(bs[i])) || return PM_ERR_NONFINITE
    end
    @inbounds for i in 1:n
        isfinite(ts[i]) || return PM_ERR_NONFINITE
    end
    @inbounds for i in 1:n
        c = mix(model, (as[3i - 2], as[3i - 1], as[3i]),
            (bs[3i - 2], bs[3i - 1], bs[3i]), ts[i])
        dest[3i - 2] = c[1]
        dest[3i - 1] = c[2]
        dest[3i] = c[3]
    end
    return PM_OK
end

"""
    bulk_mix!(dest, model, as, bs, ts) -> dest

Throwing wrapper around [`bulk_mix_kernel!`](@ref) that requires exact
buffer lengths: `as`, `bs`, and `dest` hold `3n` elements and `ts` holds `n`.
"""
function bulk_mix!(
        dest::AbstractVector{T}, model::PigmentModel,
        as::AbstractVector{T}, bs::AbstractVector{T}, ts::AbstractVector{T}
    ) where {T<:AbstractFloat}
    n = length(ts)
    (length(as) == 3n && length(bs) == 3n && length(dest) == 3n) || throw(
        DimensionMismatch(
            "bulk_mix! needs dest, as, bs of length 3n and ts of length n; got " *
                "$(length(dest)), $(length(as)), $(length(bs)), $n"
        )
    )
    _check_status(bulk_mix_kernel!(dest, model, as, bs, ts, n))
    return dest
end

function _check_status(code)
    code == PM_OK || throw(ArgumentError(error_message(code)))
    return nothing
end

"""
    weighted_mix_kernel!(dest, model, colors, weights, count) -> Int32

Encode all `count` colors, accumulate the latents with the given weights,
divide by the weight total, and decode into the first three elements of
`dest`. Zero-allocation and non-throwing, like [`bulk_mix_kernel!`](@ref).

Weights must be finite and non-negative with a strictly positive total. If
exactly one weight is positive, that input color is copied out unchanged.
"""
function weighted_mix_kernel!(
        dest::AbstractVector{T}, model::PigmentModel,
        colors::AbstractVector{T}, weights::AbstractVector{T},
        count::Integer
    ) where {T<:AbstractFloat}
    n = Int(count)
    n < 0 && return PM_ERR_LENGTH
    length(dest) >= 3 || return PM_ERR_LENGTH
    length(colors) >= 3n || return PM_ERR_LENGTH
    length(weights) >= n || return PM_ERR_LENGTH
    total = zero(T)
    single = 0
    @inbounds for i in 1:n
        w = weights[i]
        isfinite(w) || return PM_ERR_NONFINITE
        w < zero(T) && return PM_ERR_WEIGHT
        if w > zero(T)
            total += w
            single = single == 0 ? i : -1
        end
    end
    isfinite(total) || return PM_ERR_NONFINITE
    total > zero(T) || return PM_ERR_TOTAL
    @inbounds for i in 1:(3n)
        isfinite(colors[i]) || return PM_ERR_NONFINITE
    end
    if single > 0
        @inbounds begin
            dest[1] = colors[3single - 2]
            dest[2] = colors[3single - 1]
            dest[3] = colors[3single]
        end
        return PM_OK
    end
    # Seven scalar accumulators: a `Ref` or a mutable container would heap
    # allocate once per call.
    ac1 = zero(T); ac2 = zero(T); ac3 = zero(T); ac4 = zero(T)
    ar1 = zero(T); ar2 = zero(T); ar3 = zero(T)
    @inbounds for i in 1:n
        w = weights[i]
        w == zero(T) && continue
        zi = encode(model, (colors[3i - 2], colors[3i - 1], colors[3i]))
        ac1 += w * zi.c[1]
        ac2 += w * zi.c[2]
        ac3 += w * zi.c[3]
        ac4 += w * zi.c[4]
        ar1 += w * zi.r[1]
        ar2 += w * zi.r[2]
        ar3 += w * zi.r[3]
    end
    inv = one(T) / total
    z = decode(model, Latent(
        (ac1 * inv, ac2 * inv, ac3 * inv, ac4 * inv),
        (ar1 * inv, ar2 * inv, ar3 * inv),
    ))
    @inbounds begin
        dest[1] = z[1]
        dest[2] = z[2]
        dest[3] = z[3]
    end
    return PM_OK
end

"""
    weighted_mix!(dest, model, colors, weights) -> dest

Throwing wrapper around [`weighted_mix_kernel!`](@ref). `colors` holds
`3n` flat channels, `weights` holds `n`, and `dest` holds 3.
"""
function weighted_mix!(
        dest::AbstractVector{T}, model::PigmentModel,
        colors::AbstractVector{T}, weights::AbstractVector{T}
    ) where {T<:AbstractFloat}
    n = length(weights)
    length(colors) == 3n || throw(
        DimensionMismatch("colors must hold 3n channels for n = $n weights")
    )
    length(dest) >= 3 || throw(ArgumentError("dest must hold at least 3 elements"))
    _check_status(weighted_mix_kernel!(dest, model, colors, weights, n))
    return dest
end

"""
    weighted_mix(model, colors, weights) -> RGB{T}

Weighted average of `n` linear-light sRGB colors in latent space: the
intended multi-color operation. Weights must be finite, non-negative, and
total a positive value; they need not sum to one.

Encoding each color, accumulating, and decoding once is not the same as
repeatedly calling [`mix`](@ref): mixing is not associative, and repeated
pairwise mixing loses the residual of the accumulated mixture.
"""
function weighted_mix(
        model::PigmentModel, colors::AbstractVector{<:RGB{T}},
        weights::AbstractVector{T}
    ) where {T<:AbstractFloat}
    n = length(colors)
    length(weights) == n || throw(
        DimensionMismatch("got $(length(colors)) colors but $(length(weights)) weights")
    )
    flat = Vector{T}(undef, 3n)
    @inbounds for i in 1:n
        c = colors[i]
        flat[3i - 2] = c[1]
        flat[3i - 1] = c[2]
        flat[3i] = c[3]
    end
    dest = Vector{T}(undef, 3)
    weighted_mix!(dest, model, flat, weights)
    return (dest[1], dest[2], dest[3])
end

# --- default-model conveniences -------------------------------------------

"""
    encode(rgb) -> Latent

[`encode`](@ref) against [`default_model`](@ref).
"""
encode(rgb::RGB{T}) where {T<:AbstractFloat} = encode(default_model(), rgb)

"""
    decode(latent) -> RGB

[`decode`](@ref) against [`default_model`](@ref).
"""
decode(z::Latent{T}) where {T<:AbstractFloat} = decode(default_model(), z)

"""
    mix(a, b, t) -> RGB

[`mix`](@ref) against [`default_model`](@ref).
"""
mix(a::RGB{T}, b::RGB{T}, t::Real) where {T<:AbstractFloat} =
    mix(default_model(), a, b, T(t))

"""
    weighted_mix(colors, weights) -> RGB

[`weighted_mix`](@ref) against [`default_model`](@ref).
"""
weighted_mix(colors::AbstractVector{<:RGB{T}}, weights::AbstractVector{T}) where {T<:AbstractFloat} =
    weighted_mix(default_model(), colors, weights)
