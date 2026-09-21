# Encoded-sRGB adapters.
#
# The core model works in linear-light sRGB. Image data is almost always
# encoded sRGB, so these helpers do the transfer function and the byte
# rounding explicitly, in one named place. Never feed encoded values into
# `encode`, and never mix encoded concentrations with linear residuals.

"""
    linear_from_srgb(c::RGB{T}) -> RGB{T}

Apply the sRGB electro-optical transfer function (the standard
piecewise `12.92` / `1.055 * x^(1/2.4) - 0.055` curve) channel by channel.

The computation is done in `T`, so a `Float32` input does not silently
round-trip through `Float64`.
"""
@inline function linear_from_srgb(c::RGB{T}) where {T <: AbstractFloat}
    return (_srgb_decode(c[1]), _srgb_decode(c[2]), _srgb_decode(c[3]))
end

"""
    srgb_from_linear(c::RGB{T}) -> RGB{T}

Apply the inverse of [`linear_from_srgb`](@ref). The result is not clipped:
use this for round trips and clip only when converting to bytes.
"""
@inline function srgb_from_linear(c::RGB{T}) where {T <: AbstractFloat}
    return (_srgb_encode(c[1]), _srgb_encode(c[2]), _srgb_encode(c[3]))
end

@inline function _srgb_decode(c::T) where {T <: AbstractFloat}
    return c <= T(0.04045) ? c / T(12.92) : ((c + T(0.055)) / T(1.055))^T(2.4)
end

@inline function _srgb_encode(c::T) where {T <: AbstractFloat}
    return c <= T(0.0031308) ? T(12.92) * c : T(1.055) * c^T(1 / 2.4) - T(0.055)
end

"""
    linear_from_srgb8(rgb::RGB8) -> RGB{Float32}

Convert an encoded-sRGB byte triple to linear light. Bytes scale as
`b / 255`, so `0x00 -> 0.0` and `0xff -> 1.0` exactly, and the transfer
function is then applied in `Float32`.
"""
@inline function linear_from_srgb8(rgb::RGB8)
    s = 1.0f0 / 255.0f0
    return linear_from_srgb((s * Float32(rgb[1]), s * Float32(rgb[2]), s * Float32(rgb[3])))
end

"""
    srgb8_from_linear(c::RGB{T}) -> RGB8

Convert linear light to an encoded-sRGB byte triple.

Rounding rule: clip each channel to `[0, 1]`, apply the sRGB transfer
function, multiply by 255, round half away from zero (`round`), and clamp to
`0x00:0xff`. Clipping is explicit here and nowhere else in the scalar API.
"""
@inline function srgb8_from_linear(c::RGB{T}) where {T <: AbstractFloat}
    return (
        _byte_from_linear(c[1]), _byte_from_linear(c[2]), _byte_from_linear(c[3]),
    )
end

@inline function _byte_from_linear(c::T) where {T <: AbstractFloat}
    e = _srgb_encode(clamp(c, zero(T), one(T)))
    v = round(255 * e)
    return UInt8(clamp(v, zero(v), 255))
end
