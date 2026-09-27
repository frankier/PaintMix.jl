# Generic trilinear lookup kernels.
#
# One kernel serves every table in the repository: the runtime's `UInt8`
# forward and inverse tables, the precompute float forward table, and the
# coarse concentration field. Only the buffer, its element scale, and the
# scalar type vary, so they are arguments rather than separate kernels.
#
# The kernels index a concrete channel-fastest buffer directly and return
# `SVector`s. They allocate nothing and never branch on trailing-singleton
# axes.
#
# Callers must supply finite coordinates. `clamp` propagates NaN, and the
# following `unsafe_trunc` would then be undefined, so non-finite inputs are
# rejected in `mixing.jl` before a lookup happens.

"""
    LUTKernels

Trilinear lookup kernels shared by `PaintMix` and `PaintMixPrecompute`.

The kernels sample a dense channel-fastest grid. A buffer stores three
channels per vertex at the offset `ch + 3 * (i + n * (j + n * k))`, for the
zero-based channel `ch` and grid indices `i`, `j`, `k`.

This is the supported interface for code outside `PaintMix` that needs the
same interpolation as [`trilinear`](@ref PaintMix.trilinear).

# Exports

  * [`cell`](@ref): the eight scaled corners around a position plus weights.
  * [`lerp3`](@ref): linear interpolation of two RGB triples.
  * [`stages`](@ref): the x- and y-lerp stages of the interpolation ladder.
  * [`trilinear3`](@ref): the full trilinear interpolation.
"""
module LUTKernels

using StaticArrays: SVector

export cell, lerp3, stages, trilinear3

# 0-based cell index and in-cell weight for one axis, given `f = x * (n - 1)`
# already clamped to `[0, n - 1]`. `f` is finite, so `unsafe_trunc` is safe.
@inline function _axis(f::T, n::Int) where {T <: AbstractFloat}
    i = unsafe_trunc(Int, f)
    # `x == 1` lands exactly on the last grid plane; use the last valid cell
    # with weight one instead of reading past the end of the buffer.
    if i > n - 2
        i = n - 2
    end
    return i, f - T(i)
end

# The three channels at `base`, scaled by `s`: `one(T) / 255` for the byte
# tables and `one(T)` for the float tables and the coarse field.
@inline function _load3(d::AbstractVector, base::Int, s::T) where {T <: AbstractFloat}
    return @inbounds SVector(s * T(d[base]), s * T(d[base + 1]), s * T(d[base + 2]))
end

# Channel-fastest layout: `ch + 3 * (i + n * (j + n * k))`. Base offsets of
# the eight cell corners for zero-based indices `(i, j, k)`.
@inline function _corner_offsets(n::Int, i::Int, j::Int, k::Int)
    b000 = 3 * (i + n * (j + n * k)) + 1
    b100 = b000 + 3
    b010 = b000 + 3n
    b110 = b010 + 3
    b001 = b000 + 3n * n
    b101 = b001 + 3
    b011 = b001 + 3n
    b111 = b011 + 3
    return b000, b100, b010, b110, b001, b101, b011, b111
end

"""
    cell(d, n, x, y, z, s) -> (corners, fx, fy, fz)

The eight scaled channel triples of the cell containing `(x, y, z)`, plus the
in-cell weights. Each coordinate is clamped to `[0, 1]` and scaled to the
grid. A one-plane table degenerates to the single vertex with zero weights.

`d` is a channel-fastest buffer of at least `3n^3` entries and `s` scales a
stored channel to a color value, e.g. `1/255` for a byte table and `1` for a
float table.
"""
@inline function cell(
        d::AbstractVector, n::Int, x::T, y::T, z::T, s::T
    ) where {T <: AbstractFloat}
    if n == 1
        v = _load3(d, 1, s)
        return (v, v, v, v, v, v, v, v), zero(T), zero(T), zero(T)
    end
    gx = clamp(x, zero(T), one(T)) * T(n - 1)
    gy = clamp(y, zero(T), one(T)) * T(n - 1)
    gz = clamp(z, zero(T), one(T)) * T(n - 1)
    i, fx = _axis(gx, n)
    j, fy = _axis(gy, n)
    k, fz = _axis(gz, n)
    b000, b100, b010, b110, b001, b101, b011, b111 = _corner_offsets(n, i, j, k)
    corners = (
        _load3(d, b000, s), _load3(d, b100, s),
        _load3(d, b010, s), _load3(d, b110, s),
        _load3(d, b001, s), _load3(d, b101, s),
        _load3(d, b011, s), _load3(d, b111, s),
    )
    return corners, fx, fy, fz
end

"""
    lerp3(a, b, w) -> SVector{3, T}

Linear interpolation `a + w * (b - a)` between two RGB triples.
"""
@inline function lerp3(a::SVector{3, T}, b::SVector{3, T}, w::T) where {T <: AbstractFloat}
    return a + w * (b - a)
end

"""
    stages(corners, fx, fy) -> (c00, c10, c01, c11, c0, c1)

The x- and y-lerp stages of the interpolation ladder on the eight `corners`
from [`cell`](@ref). Value-only callers drop the four unused stages; the
Jacobian keeps them.
"""
@inline function stages(corners, fx::T, fy::T) where {T <: AbstractFloat}
    c000, c100, c010, c110, c001, c101, c011, c111 = corners
    c00 = lerp3(c000, c100, fx)
    c10 = lerp3(c010, c110, fx)
    c01 = lerp3(c001, c101, fx)
    c11 = lerp3(c011, c111, fx)
    c0 = lerp3(c00, c10, fy)
    c1 = lerp3(c01, c11, fy)
    return c00, c10, c01, c11, c0, c1
end

"""
    trilinear3(corners, fx, fy, fz) -> SVector{3, T}

The trilinear interpolation of the eight `corners` from [`cell`](@ref) at the
weights `fx`, `fy`, `fz`.
"""
@inline function trilinear3(corners, fx::T, fy::T, fz::T) where {T <: AbstractFloat}
    _, _, _, _, c0, c1 = stages(corners, fx, fy)
    return lerp3(c0, c1, fz)
end

end # module LUTKernels
