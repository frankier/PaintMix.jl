# Lookup-table indexing and trilinear interpolation.
#
# The kernels take a concrete `ByteLUT` and scalar type, index the raw byte
# buffer directly, and return tuples. They allocate nothing, dispatch on
# nothing, and never branch on trailing-singleton axes.
#
# Callers must supply finite coordinates. `clamp` propagates NaN, and the
# following `unsafe_trunc` would then be undefined, so non-finite inputs are
# rejected in `mixing.jl` before a lookup happens.

# 0-based cell index and in-cell weight for one axis, given `f = x * (n - 1)`
# already clamped to `[0, n - 1]`.
@inline function _axis(f::T, n::Int) where {T <: AbstractFloat}
    i = unsafe_trunc(Int, f)
    # `x == 1` lands exactly on the last grid plane; use the last valid cell
    # with weight one instead of reading past the end of the buffer.
    if i > n - 2
        i = n - 2
    end
    return i, f - T(i)
end

@inline function _channel3(d::Vector{UInt8}, base::Int)::NTuple{3, UInt8}
    return @inbounds (d[base], d[base + 1], d[base + 2])
end

@inline function _lerp3(a::NTuple{3, T}, b::NTuple{3, T}, w::T) where {T <: AbstractFloat}
    s = one(T) - w
    return (s * a[1] + w * b[1], s * a[2] + w * b[2], s * a[3] + w * b[3])
end

@inline function _scale3(v::NTuple{3, UInt8}, s::T) where {T <: AbstractFloat}
    return (s * T(v[1]), s * T(v[2]), s * T(v[3]))
end

"""
    trilinear(lut::ByteLUT, x, y, z) -> NTuple{3,T}

Sample `lut` trilinearly at `(x, y, z)`, each of which is clamped to
`[0, 1]`. Values are returned in `[0, 1]` after scaling bytes by `1/255`.

Coordinates on a grid plane return the stored vertex value exactly. This is
the only interpolation kernel in the runtime; both directions of the model
use it.
"""
@inline function trilinear(
        lut::ByteLUT, x::T, y::T, z::T
    ) where {T <: AbstractFloat}
    n = lut.n
    d = lut.data
    s = one(T) / T(255)
    if n == 1
        return _scale3(_channel3(d, 1), s)
    end
    gx = clamp(x, zero(T), one(T)) * T(n - 1)
    gy = clamp(y, zero(T), one(T)) * T(n - 1)
    gz = clamp(z, zero(T), one(T)) * T(n - 1)
    i, fx = _axis(gx, n)
    j, fy = _axis(gy, n)
    k, fz = _axis(gz, n)
    # Channel-fastest layout: ch + 3 * (i + n * (j + n * k)).
    b000 = 3 * (i + n * (j + n * k)) + 1
    b100 = b000 + 3
    b010 = b000 + 3n
    b110 = b010 + 3
    b001 = b000 + 3n * n
    b101 = b001 + 3
    b011 = b001 + 3n
    b111 = b011 + 3
    c00 = _lerp3(_scale3(_channel3(d, b000), s), _scale3(_channel3(d, b100), s), fx)
    c10 = _lerp3(_scale3(_channel3(d, b010), s), _scale3(_channel3(d, b110), s), fx)
    c01 = _lerp3(_scale3(_channel3(d, b001), s), _scale3(_channel3(d, b101), s), fx)
    c11 = _lerp3(_scale3(_channel3(d, b011), s), _scale3(_channel3(d, b111), s), fx)
    c0 = _lerp3(c00, c10, fy)
    c1 = _lerp3(c01, c11, fy)
    return _lerp3(c0, c1, fz)
end

"""
    vertex(lut::ByteLUT, i, j, k) -> NTuple{3,UInt8}

The stored bytes at zero-based grid indices `(i, j, k)`. Exposed for
validation and tests; not used in the mixing path.
"""
function vertex(lut::ByteLUT, i::Integer, j::Integer, k::Integer)
    n = lut.n
    0 <= i < n && 0 <= j < n && 0 <= k < n || throw(BoundsError(lut, (i, j, k)))
    return _channel3(lut.data, 3 * (Int(i) + n * (Int(j) + n * Int(k))) + 1)
end
