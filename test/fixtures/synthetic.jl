# Synthetic model fixtures for the runtime tests.
#
# These builders make tiny `PigmentModel`s with known analytic behavior, so
# the lookup kernel, the latent algebra, and the payload format can be tested
# without any spectral data or optimization. They are test scaffolding, not
# part of the package: the release tables come from `PaintMixPrecompute`.

module SyntheticFixtures

using PaintMix
using Random: MersenneTwister

export quantize_simplex,
    project_to_simplex3,
    synthetic_model,
    identity_model,
    permuted_model,
    signed_model,
    offset_model,
    random_model,
    reference_trilinear,
    reference_vertex

"""
    quantize_simplex(c) -> NTuple{4,Int}

Quantize four concentrations to non-negative integers summing to exactly 255
by largest-remainder rounding, with ties broken by index.

Independent rounding would let the first three exceed 255 and imply a
negative fourth; here the fourth is always `255 - b1 - b2 - b3 >= 0`.
"""
function quantize_simplex(c::NTuple{4,<:Real})
    scaled = ntuple(i -> max(Float64(c[i]), 0.0) * 255.0, Val(4))
    floors = ntuple(i -> floor(Int, scaled[i]), Val(4))
    total = sum(floors)
    order = sort(collect(1:4); by = i -> (-(scaled[i] - floors[i]), i))
    out = collect(floors)
    deficit = 255 - total
    k = 1
    while deficit > 0 && k <= 4
        out[order[k]] += 1
        deficit -= 1
        k += 1
    end
    k = 4
    while deficit < 0 && k >= 1
        if out[order[k]] > 0
            out[order[k]] -= 1
            deficit += 1
        end
        k -= 1
    end
    return (out[1], out[2], out[3], out[4])
end

"""
    project_to_simplex3(x, y, z) -> NTuple{3,Float64}

Map a cube point to the concentration simplex. Points inside the simplex are
kept; points outside are scaled toward the origin, which is continuous on
the whole cube and never implies a negative fourth concentration.

The production forward table's padding rule is decided by
`PaintMixPrecompute`; this stand-in exists so fixtures are always valid.
"""
function project_to_simplex3(x::Real, y::Real, z::Real)
    a = max(Float64(x), 0.0)
    b = max(Float64(y), 0.0)
    c = max(Float64(z), 0.0)
    s = a + b + c
    if s <= 1.0
        return (a, b, c)
    end
    return (a / s, b / s, c / s)
end

function _fill_forward!(data::Vector{UInt8}, n::Int, f)
    for k in 0:(n - 1), j in 0:(n - 1), i in 0:(n - 1)
        c = n == 1 ? (0.0, 0.0, 0.0) : (i / (n - 1), j / (n - 1), k / (n - 1))
        rgb = f(c[1], c[2], c[3])
        base = 3 * (i + n * (j + n * k)) + 1
        for ch in 0:2
            data[base + ch] = UInt8(clamp(round(255 * Float64(rgb[ch + 1])), 0, 255))
        end
    end
    return data
end

function _fill_inverse!(data::Vector{UInt8}, n::Int, g)
    for k in 0:(n - 1), j in 0:(n - 1), i in 0:(n - 1)
        rgb = n == 1 ? (0.0, 0.0, 0.0) : (i / (n - 1), j / (n - 1), k / (n - 1))
        c = g(rgb[1], rgb[2], rgb[3])
        q = quantize_simplex((c[1], c[2], c[3], 1 - c[1] - c[2] - c[3]))
        base = 3 * (i + n * (j + n * k)) + 1
        data[base] = UInt8(q[1])
        data[base + 1] = UInt8(q[2])
        data[base + 2] = UInt8(q[3])
    end
    return data
end

"""
    synthetic_model(forward, inverse; n = 2, id = ..., flags = ...) -> PigmentModel

Build a model whose forward table samples `forward(c1, c2, c3)` and whose
inverse table samples `inverse(r, g, b)` on the grid. `inverse` must return
concentrations on the simplex; they are quantized with
[`quantize_simplex`](@ref).
"""
function synthetic_model(
        forward, inverse; n::Integer = 2,
        id::NTuple{16,UInt8} = ntuple(_ -> 0x00, Val(16)),
        flags::UInt32 = FLAG_FORWARD_SIMPLEX_PROJECTED | FLAG_INVERSE_LARGEST_REMAINDER
    )
    m = Int(n)
    payload = 3 * m^3
    fwd = _fill_forward!(Vector{UInt8}(undef, payload), m, forward)
    inv = _fill_inverse!(Vector{UInt8}(undef, payload), m, inverse)
    model = PigmentModel(
        id, ByteLUT(m, inv), ByteLUT(m, fwd), FORMAT_VERSION, flags,
        crc32(fwd), crc32(inv),
    )
    return validate_model(model)
end

"""
    identity_model(n = 2) -> PigmentModel

A synthetic model with a projected-identity forward table and the same
projection as its inverse. Use it for round-trip and algebra tests.
"""
identity_model(n::Integer = 2) = synthetic_model(
    (c1, c2, c3) -> project_to_simplex3(c1, c2, c3),
    project_to_simplex3; n,
)

"""
    permuted_model(n = 2) -> PigmentModel

A forward table sampling the exact affine map `(c1, c2, c3) -> (c3, c1, c2)`.
Every stored byte is `0` or `255`, so trilinear interpolation reproduces the
map exactly for every query, including cube points outside the simplex. This
is the analytic interpolation check.
"""
permuted_model(n::Integer = 2) = synthetic_model(
    (c1, c2, c3) -> (c3, c1, c2),
    project_to_simplex3; n,
)

"""
    signed_model(n = 2) -> PigmentModel

A forward table sampling the exact affine map
`(c1, c2, c3) -> (1 - c1, c2, 1 - c3)`, again byte-exact at the corners.
`_repair_simplex` is the identity on the simplex, so `forward_rgb` can be
compared with the closed form directly.
"""
signed_model(n::Integer = 2) = synthetic_model(
    (c1, c2, c3) -> (1 - c1, c2, 1 - c3),
    project_to_simplex3; n,
)

"""
    offset_model(n = 8) -> PigmentModel

A forward table with a constant offset and a gain below one per channel. A
byte-valued table cannot hold values outside `[0, 1]`, so out-of-gamut
behavior in this model can only come from residuals; the fixture exists so
tests can add a large residual and confirm that `decode` propagates it
without clipping.
"""
offset_model(n::Integer = 8) = synthetic_model(
    (c1, c2, c3) -> (0.9 * c1 + 0.05, 0.9 * c2 + 0.1 * c3, 0.8 * c3 + 0.1 * c1),
    project_to_simplex3; n,
)

"""
    random_model(seed = 1, n = 5) -> PigmentModel

A deterministic pseudo-random model with no analytic meaning, used to fuzz
the kernels and the payload round trip. Its inverse table always satisfies
the simplex constraint.
"""
function random_model(seed::Integer = 1, n::Integer = 5)
    rng = MersenneTwister(seed)
    m = Int(n)
    payload = 3 * m^3
    fwd = rand(rng, UInt8, payload)
    inv = Vector{UInt8}(undef, payload)
    for v in 0:(m^3 - 1)
        base = 3v + 1
        a = rand(rng, 0:255)
        b = rand(rng, 0:(255 - a))
        c = rand(rng, 0:(255 - a - b))
        inv[base] = UInt8(a)
        inv[base + 1] = UInt8(b)
        inv[base + 2] = UInt8(c)
    end
    id = ntuple(i -> rand(rng, UInt8), Val(16))
    model = PigmentModel(
        id, ByteLUT(m, inv), ByteLUT(m, fwd), FORMAT_VERSION,
        FLAG_FORWARD_SIMPLEX_PROJECTED | FLAG_INVERSE_LARGEST_REMAINDER,
        crc32(fwd), crc32(inv),
    )
    return validate_model(model)
end

"""
    reference_vertex(lut, i, j, k, ch) -> UInt8

Read one stored byte by restating the channel-fastest offset formula
independently of the kernel. `ch` is zero-based.
"""
reference_vertex(lut::ByteLUT, i::Int, j::Int, k::Int, ch::Int) =
    lut.data[ch + 1 + 3 * (i + lut.n * (j + lut.n * k))]

"""
    reference_trilinear(lut, x, y, z) -> NTuple{3,Float64}

An independent trilinear implementation: eight weight products rather than
three stages of lerps. Any disagreement with `PaintMix.trilinear` beyond
rounding is a bug in one of the two.
"""
function reference_trilinear(lut::ByteLUT, x::Real, y::Real, z::Real)
    n = lut.n
    if n == 1
        return ntuple(ch -> Float64(reference_vertex(lut, 0, 0, 0, ch - 1)) / 255.0, Val(3))
    end
    gx = clamp(Float64(x), 0.0, 1.0) * (n - 1)
    gy = clamp(Float64(y), 0.0, 1.0) * (n - 1)
    gz = clamp(Float64(z), 0.0, 1.0) * (n - 1)
    i0 = min(floor(Int, gx), n - 2)
    j0 = min(floor(Int, gy), n - 2)
    k0 = min(floor(Int, gz), n - 2)
    fx = gx - i0
    fy = gy - j0
    fz = gz - k0
    return ntuple(Val(3)) do ch
        s = 0.0
        for dk in 0:1, dj in 0:1, di in 0:1
            w = (di == 1 ? fx : 1 - fx) * (dj == 1 ? fy : 1 - fy) *
                (dk == 1 ? fz : 1 - fz)
            v = reference_vertex(lut, i0 + di, j0 + dj, k0 + dk, ch - 1)
            s += w * Float64(v) / 255.0
        end
        s
    end
end

end # module SyntheticFixtures
