# Runtime lookup-table sampling.
#
# The generic trilinear kernels live in the `LUTKernels` submodule, which
# `PaintMixPrecompute` also uses. This file wraps them for the byte tables and
# adds the direct-vertex accessor.

using .LUTKernels: cell, trilinear3

@inline function _channel3(d::Vector{UInt8}, base::Int)::SVector{3, UInt8}
    return @inbounds SVector(d[base], d[base + 1], d[base + 2])
end

"""
    trilinear(lut::ByteLUT, x, y, z) -> RGB{T}

Sample `lut` trilinearly at `(x, y, z)`, each of which is clamped to
`[0, 1]`. Values are returned in `[0, 1]` after scaling bytes by `1/255`.

Coordinates on a grid plane return the stored vertex value exactly. This is
the only interpolation kernel in the runtime; both directions of the model
use it.
"""
@inline function trilinear(
        lut::ByteLUT, x::T, y::T, z::T
    ) where {T <: AbstractFloat}
    corners, fx, fy, fz = cell(lut.data, lut.n, x, y, z, one(T) / T(255))
    return trilinear3(corners, fx, fy, fz)
end

"""
    vertex(lut::ByteLUT, i, j, k) -> RGB8

The stored bytes at zero-based grid indices `(i, j, k)`. Exposed for
validation and tests; not used in the mixing path.
"""
function vertex(lut::ByteLUT, i::Integer, j::Integer, k::Integer)
    n = lut.n
    0 <= i < n && 0 <= j < n && 0 <= k < n || throw(BoundsError(lut, (i, j, k)))
    return _channel3(lut.data, 3 * (Int(i) + n * (Int(j) + n * Int(k))) + 1)
end
