@inline function trilinear(lut::WordLUT, x::T, y::T, z::T) where {T <: AbstractFloat}
    corners, fx, fy, fz = cell(lut.data, lut.n, x, y, z, one(T) / T(65535))
    return trilinear3(corners, fx, fy, fz)
end

@inline _forward_lookup(lut::ByteLUT, x, y, z) = trilinear(lut, x, y, z)
@inline _fraction_before(a, b, f) = f[a] > f[b] || (f[a] == f[b] && a > b)
@inline _axis_step(axis) = SVector(Int(axis == 1), Int(axis == 2), Int(axis == 3))
@inline function _simplex_vertex(lut::SimplexLUT, i, ::Type{T}) where {T}
    o = 3 * _simplex_index(i...) + 1
    scale = one(T) / T(255)
    return @inbounds SVector(T(lut.data[o]), T(lut.data[o + 1]), T(lut.data[o + 2])) .* scale
end

@inline function _forward_lookup(lut::SimplexLUT, x::T, y::T, z::T) where {T <: AbstractFloat}
    d = lut.n - 1
    q = clamp.(SVector(x, x + y, x + y + z), zero(T), one(T)) .* T(d)
    base = min.(floor.(Int, q), d - 1)
    f = q - base
    a, b, c = 1, 2, 3
    _fraction_before(b, a, f) && ((a, b) = (b, a))
    _fraction_before(c, b, f) && ((b, c) = (c, b))
    _fraction_before(b, a, f) && ((a, b) = (b, a))
    # Higher axes win ties, keeping every accessed vertex in the ordered domain.
    v0 = _simplex_vertex(lut, base, T)
    v1 = _simplex_vertex(lut, base + _axis_step(a), T)
    v2 = _simplex_vertex(lut, base + _axis_step(a) + _axis_step(b), T)
    v3 = _simplex_vertex(lut, base + SVector(1, 1, 1), T)
    return (one(T) - f[a]) * v0 + (f[a] - f[b]) * v1 + (f[b] - f[c]) * v2 + f[c] * v3
end
