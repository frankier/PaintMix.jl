# Core value types: colors, latents, and lookup-table containers.
#
# Everything here is allocation-free and `isbits` except the two contiguous
# byte buffers held by a `PigmentModel`. Scalars are either `Float32` or
# `Float64`; the two are never mixed implicitly.

"""
    RGB{T}

A linear-light sRGB color as a 3-tuple of `T`, nominally in `[0, 1]`.
Values outside that range are legal and meaningful: see [`mix`](@ref) and
[`decode`](@ref), neither of which clips.
"""
const RGB{T<:AbstractFloat} = NTuple{3,T}

"""
    RGB8

An encoded-sRGB color as three bytes. Byte colors carry the sRGB transfer
function; convert with [`linear_from_srgb8`](@ref) before mixing.
"""
const RGB8 = NTuple{3,UInt8}

"""
    Concentrations{T}

Four pigment concentrations as a tuple, non-negative and summing to one.
The index order is fixed by the model's table layout; see
[`PigmentModel`](@ref).
"""
const Concentrations{T<:AbstractFloat} = NTuple{4,T}

"""
    Latent{T}

The paper's latent representation: four pigment concentrations plus the
signed RGB residual that must be added back to the pigment mixture to
recover the encoded color exactly.

`Latent` is an immutable `isbits` value, so it can live in registers, be
stored in arrays, and be passed through the C ABI as seven scalars. The
residual is kept in floating point and is never clamped or quantized.

# Field layout

  * `c`: four concentrations, non-negative and summing to one.
  * `r`: the residual `x - M(c)` in linear-light sRGB, where `x` is the
    encoded color.

A latent belongs to the table model that produced it. Latents from
different models must not be combined; the current API documents the
requirement rather than enforcing it, because a model tag would enlarge the
ABI representation.

`c[4]` is redundant: it always equals `1 - c[1] - c[2] - c[3]`. Decoding
reconstructs it from the first three so that the simplex stays consistent
even after floating-point interpolation.
"""
struct Latent{T<:AbstractFloat}
    c::NTuple{4,T}
    r::NTuple{3,T}
end

"""
    concentrations(z::Latent) -> NTuple{4,T}

The four pigment concentrations of `z`. The fourth is the value stored at
encoding time; decoding always reconstructs it as `1 - c[1] - c[2] - c[3]`.
"""
concentrations(z::Latent) = z.c

"""
    residual(z::Latent) -> NTuple{3,T}

The signed linear-light residual of `z`.
"""
residual(z::Latent) = z.r

function Base.show(io::IO, z::Latent{T}) where {T}
    print(io, "Latent{", T, "}(c = ", z.c, ", r = ", z.r, ")")
    return nothing
end

"""
    ByteLUT

A dense `n x n x n x 3` lookup table stored as raw bytes, one byte per
channel, in channel-fastest order.

The offset of channel `ch` (zero-based, `0:2`) at grid indices `i`, `j`, `k`
(zero-based, `0:n-1`) is

    ch + 3 * (i + n * (j + n * k))

Byte values scale as `b / 255`. The table is a plain `Vector{UInt8}`, so it
is contiguous, relocatable, and cheap to embed in a compiled image. The
struct is immutable but the buffer is not: treat a model's buffers as
read-only.
"""
struct ByteLUT
    n::Int
    data::Vector{UInt8}

    function ByteLUT(n::Integer, data::AbstractVector{UInt8})
        m = Int(n)
        m >= 1 || throw(ArgumentError("lookup table size must be >= 1, got $m"))
        expected = 3 * m^3
        length(data) == expected || throw(
            ArgumentError(
                "lookup table with n = $m needs $expected bytes, got $(length(data))"
            )
        )
        return new(m, convert(Vector{UInt8}, data))
    end
end

Base.:(==)(a::ByteLUT, b::ByteLUT) = a.n == b.n && a.data == b.data
Base.hash(l::ByteLUT, h::UInt) = hash(l.data, hash(l.n, hash(:ByteLUT, h)))

function Base.show(io::IO, l::ByteLUT)
    print(io, "ByteLUT(n = ", l.n, ", ", length(l.data), " bytes)")
    return nothing
end

# Header flag bits. Both rules are generation-time choices that a consumer
# may want to confirm; see `data/default/README.md`.
const FLAG_FORWARD_SIMPLEX_PROJECTED = UInt32(0x00000001)
const FLAG_INVERSE_LARGEST_REMAINDER = UInt32(0x00000002)

"""
    PigmentModel

A validated pair of lookup tables, plus the provenance needed to reject a
mismatched payload.

  * `inverse` maps linear-light sRGB to concentrations, storing the first
    three of the four at every grid vertex. The fourth is reconstructed.
  * `forward` maps concentrations to linear-light sRGB.
  * `id` is the 16-byte model identifier from the payload header.
  * `flags` records generation rules, e.g.
    [`FLAG_FORWARD_SIMPLEX_PROJECTED`](@ref).

Tables are valid on the whole cube `[0, 1]^3`. Vertices outside the
concentration simplex are padding; see `data/default/README.md` for the
documented padding rule.
"""
struct PigmentModel
    id::NTuple{16,UInt8}
    inverse::ByteLUT
    forward::ByteLUT
    format_version::UInt16
    flags::UInt32
end

"""
    grid_n(model::PigmentModel) -> Int

The edge length `n` of both lookup tables.
"""
grid_n(m::PigmentModel) = m.forward.n

"""
    model_id(model::PigmentModel) -> String

The 16-byte model identifier as a lowercase hex string.
"""
function model_id(m::PigmentModel)
    io = IOBuffer()
    for b in m.id
        print(io, string(b, base = 16, pad = 2))
    end
    return String(take!(io))
end

function Base.show(io::IO, m::PigmentModel)
    print(
        io, "PigmentModel(id = ", model_id(m), ", n = ", grid_n(m),
        ", flags = 0x", string(m.flags, base = 16, pad = 8), ")"
    )
    return nothing
end
