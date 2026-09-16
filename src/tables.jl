# The `.pmx` payload format: serialization, validation, and the packaged
# default model.
#
# The format is deliberately dumb: a fixed-width little-endian header
# followed by two raw byte tables. There is no compression, no versioned
# generic container, and no Julia `Serialization`, so a payload can be
# embedded in a compiled image, memory-mapped, or produced by another
# language's tooling.
#
# See `data/default/README.md` for the field semantics that are not encoded
# in the bytes themselves (pigment order, padding rule, provenance).

const FORMAT_VERSION = UInt16(1)

# "PAINTMIX"
const _MAGIC = (0x50, 0x41, 0x49, 0x4e, 0x54, 0x4d, 0x49, 0x58)

const HEADER_BYTES = 128

const STORAGE_U8 = 0x00
const STORAGE_F32 = 0x01
const STORAGE_F64 = 0x02

const COLORSPACE_LINEAR_SRGB_D65 = 0x00

const BYTE_SCALE_255 = 0x00
const INTERP_TRILINEAR = 0x00
const INDEX_CHANNEL_FAST = 0x00
const CHECKSUM_CRC32 = 0x01

const _CRC32_POLY = 0xEDB88320

const _CRC32_TABLE = let t = Vector{UInt32}(undef, 256)
    @inbounds for i in 0:255
        c = UInt32(i)
        for _ in 1:8
            c = (c & one(UInt32)) == one(UInt32) ? (c >> 1) ⊻ _CRC32_POLY : c >> 1
        end
        t[i + 1] = c
    end
    Tuple(t)
end

"""
    crc32(data, len = length(data)) -> UInt32

CRC-32/ISO-HDLC (reflected, polynomial `0xEDB88320`, initial value
`0xffffffff`, final XOR `0xffffffff`), the same function as `zlib.crc32`.

Implemented here so the runtime has no dependencies. Used to checksum
payloads, not to authenticate them.
"""
function crc32(data::AbstractVector{UInt8}, len::Integer = length(data))
    c = typemax(UInt32)
    @inbounds for i in 1:Int(len)
        c = _CRC32_TABLE[((c ⊻ data[i]) & 0xff) + 1] ⊻ (c >> 8)
    end
    return c ⊻ typemax(UInt32)
end

# --- header fields ---------------------------------------------------------

@inline _hdr_u8(b::AbstractVector{UInt8}, o::Int) = b[o + 1]

@inline function _hdr_u16(b::AbstractVector{UInt8}, o::Int)
    return UInt16(b[o + 1]) | (UInt16(b[o + 2]) << 8)
end

@inline function _hdr_u32(b::AbstractVector{UInt8}, o::Int)
    return UInt32(b[o + 1]) | (UInt32(b[o + 2]) << 8) |
        (UInt32(b[o + 3]) << 16) | (UInt32(b[o + 4]) << 24)
end

@inline function _hdr_u64(b::AbstractVector{UInt8}, o::Int)
    lo = UInt64(_hdr_u32(b, o))
    hi = UInt64(_hdr_u32(b, o + 4))
    return lo | (hi << 32)
end

@inline function _put_u8!(b::Vector{UInt8}, o::Int, v::Integer)
    @inbounds b[o + 1] = UInt8(v & 0xff)
    return nothing
end

@inline function _put_u16!(b::Vector{UInt8}, o::Int, v::Integer)
    u = UInt16(v & 0xffff)
    @inbounds begin
        b[o + 1] = UInt8(u & 0xff)
        b[o + 2] = UInt8((u >> 8) & 0xff)
    end
    return nothing
end

@inline function _put_u32!(b::Vector{UInt8}, o::Int, v::Integer)
    u = UInt32(v & 0xffffffff)
    @inbounds for i in 0:3
        b[o + 1 + i] = UInt8((u >> (8i)) & 0xff)
    end
    return nothing
end

@inline function _put_u64!(b::Vector{UInt8}, o::Int, v::Integer)
    u = UInt64(v)
    @inbounds for i in 0:7
        b[o + 1 + i] = UInt8((u >> (8i)) & 0xff)
    end
    return nothing
end

"""
    ModelHeader

The decoded fixed-width header of a `.pmx` payload: format version, grid
size, storage conventions, the model identifier, payload checksums, and
payload offsets. Exposed for tests and tooling; `PigmentModel` keeps only
what the mixing path needs.
"""
struct ModelHeader
    format_version::UInt16
    header_bytes::UInt16
    storage::UInt8
    channels::UInt8
    table_count::UInt8
    color_space::UInt8
    grid_n::Int
    flags::UInt32
    id::NTuple{16,UInt8}
    byte_scale::UInt8
    interpolation::UInt8
    index_order::UInt8
    checksum::UInt8
    inverse_crc32::UInt32
    forward_crc32::UInt32
    inverse_offset::Int
    inverse_bytes::Int
    forward_offset::Int
    forward_bytes::Int
end

struct InvalidPayload <: Exception
    msg::String
end

Base.showerror(io::IO, e::InvalidPayload) = print(io, "InvalidPayload: ", e.msg)

_parse_header(b::AbstractVector{UInt8}) = ModelHeader(
    _hdr_u16(b, 8), _hdr_u16(b, 10), _hdr_u8(b, 12), _hdr_u8(b, 13),
    _hdr_u8(b, 14), _hdr_u8(b, 15), Int(_hdr_u32(b, 16)), _hdr_u32(b, 20),
    ntuple(i -> _hdr_u8(b, 23 + i), Val(16)),
    _hdr_u8(b, 40), _hdr_u8(b, 41), _hdr_u8(b, 42), _hdr_u8(b, 43),
    _hdr_u32(b, 44), _hdr_u32(b, 48),
    Int(_hdr_u64(b, 56)), Int(_hdr_u64(b, 64)),
    Int(_hdr_u64(b, 72)), Int(_hdr_u64(b, 80)),
)

function _check_header(h::ModelHeader)
    h.format_version == FORMAT_VERSION || throw(InvalidPayload(
        "unsupported format version $(h.format_version), runtime understands $(FORMAT_VERSION)"
    ))
    h.header_bytes >= HEADER_BYTES || throw(InvalidPayload(
        "header_bytes $(h.header_bytes) is smaller than the $HEADER_BYTES byte header"
    ))
    h.storage == STORAGE_U8 || throw(InvalidPayload(
        "storage type $(h.storage) is not supported by this runtime (only u8 = 0)"
    ))
    h.channels == 3 || throw(InvalidPayload("channel count must be 3, got $(h.channels)"))
    h.table_count == 2 || throw(InvalidPayload("table count must be 2, got $(h.table_count)"))
    h.color_space == COLORSPACE_LINEAR_SRGB_D65 || throw(InvalidPayload(
        "color space $(h.color_space) is not linear-light sRGB with D65 (0)"
    ))
    h.byte_scale == BYTE_SCALE_255 || throw(InvalidPayload(
        "byte scaling $(h.byte_scale) is not b/255 (0)"
    ))
    h.interpolation == INTERP_TRILINEAR || throw(InvalidPayload(
        "interpolation $(h.interpolation) is not trilinear (0)"
    ))
    h.index_order == INDEX_CHANNEL_FAST || throw(InvalidPayload(
        "index order $(h.index_order) is not channel-fast (0)"
    ))
    h.checksum == CHECKSUM_CRC32 || throw(InvalidPayload(
        "checksum algorithm $(h.checksum) is not CRC-32 (1)"
    ))
    h.grid_n >= 1 || throw(InvalidPayload("grid size must be >= 1, got $(h.grid_n)"))
    expected = 3 * h.grid_n^3
    h.inverse_bytes == expected || throw(InvalidPayload(
        "inverse table has $(h.inverse_bytes) bytes, expected $expected for n = $(h.grid_n)"
    ))
    h.forward_bytes == expected || throw(InvalidPayload(
        "forward table has $(h.forward_bytes) bytes, expected $expected for n = $(h.grid_n)"
    ))
    return nothing
end

function _header_bytes(model::PigmentModel)
    n = grid_n(model)
    payload = 3 * n^3
    b = zeros(UInt8, HEADER_BYTES)
    for (i, m) in enumerate(_MAGIC)
        b[i] = m
    end
    _put_u16!(b, 8, model.format_version)
    _put_u16!(b, 10, HEADER_BYTES)
    _put_u8!(b, 12, STORAGE_U8)
    _put_u8!(b, 13, 3)
    _put_u8!(b, 14, 2)
    _put_u8!(b, 15, COLORSPACE_LINEAR_SRGB_D65)
    _put_u32!(b, 16, n)
    _put_u32!(b, 20, model.flags)
    for (i, v) in enumerate(model.id)
        b[23 + i + 1] = v
    end
    _put_u8!(b, 40, BYTE_SCALE_255)
    _put_u8!(b, 41, INTERP_TRILINEAR)
    _put_u8!(b, 42, INDEX_CHANNEL_FAST)
    _put_u8!(b, 43, CHECKSUM_CRC32)
    _put_u32!(b, 44, crc32(model.inverse.data))
    _put_u32!(b, 48, crc32(model.forward.data))
    # Payload order matches `write_model`: inverse first, then forward.
    _put_u64!(b, 56, HEADER_BYTES)
    _put_u64!(b, 64, payload)
    _put_u64!(b, 72, HEADER_BYTES + payload)
    _put_u64!(b, 80, payload)
    return b
end

"""
    model_from_bytes(bytes; validate = true, checksum = true) -> PigmentModel

Parse a `.pmx` payload held in memory.

`checksum` verifies both table CRC-32 values against the payload's own header;
`validate` additionally checks that every inverse-table vertex satisfies the
simplex constraint (`b1 + b2 + b3 <= 255`). Both default to `true`, which is
what a runtime loader should use.

Making the checks optional exists for one caller: the compiled library embeds
bytes that its build driver has already validated, and re-verifying them
inside the image-building interpreter (which runs without optimization) costs
minutes on a release-sized payload. Turning the checks off is a statement
that the bytes were verified elsewhere, not a relaxation of the format.

Throws `PaintMix.InvalidPayload` on any structural problem. Tables are copied
out of `bytes`, so the caller may release it afterwards.
"""
function model_from_bytes(
        bytes::AbstractVector{UInt8}; validate::Bool = true, checksum::Bool = true
    )
    length(bytes) >= HEADER_BYTES || throw(InvalidPayload(
        "payload is $(length(bytes)) bytes, shorter than the $HEADER_BYTES byte header"
    ))
    for (i, m) in enumerate(_MAGIC)
        bytes[i] == m || throw(InvalidPayload("bad magic: not a PaintMix payload"))
    end
    h = _parse_header(bytes)
    _check_header(h)
    total = max(h.inverse_offset + h.inverse_bytes, h.forward_offset + h.forward_bytes)
    length(bytes) >= total || throw(InvalidPayload(
        "payload is $(length(bytes)) bytes, but the header places tables up to byte $total"
    ))
    inverse = _copy_table(bytes, h.inverse_offset, h.inverse_bytes, h.grid_n)
    forward = _copy_table(bytes, h.forward_offset, h.forward_bytes, h.grid_n)
    icrc = h.inverse_crc32
    fcrc = h.forward_crc32
    if checksum
        icrc = crc32(inverse.data)
        icrc == h.inverse_crc32 || throw(InvalidPayload(
            "inverse table checksum mismatch: payload has $(h.inverse_crc32), computed $icrc"
        ))
        fcrc = crc32(forward.data)
        fcrc == h.forward_crc32 || throw(InvalidPayload(
            "forward table checksum mismatch: payload has $(h.forward_crc32), computed $fcrc"
        ))
    end
    model = PigmentModel(
        h.id, inverse, forward, h.format_version, h.flags, fcrc, icrc,
    )
    validate && validate_model(model)
    return model
end

function _copy_table(bytes::AbstractVector{UInt8}, offset::Int, len::Int, n::Int)
    data = Vector{UInt8}(undef, len)
    # `copyto!` on two `Vector`s is a memmove, so this stays cheap even when
    # the payload is being loaded by an unoptimized image-building process.
    copyto!(data, 1, bytes, offset + 1, len)
    return ByteLUT(n, data)
end

"""
    validate_model(model; simplex = true) -> model

Check the invariants the mixing kernels rely on: both tables have `3n^3`
bytes, every byte pattern is a value (trivially true for `UInt8`, so this
exists for clarity), and — when `simplex` is set — every inverse-table
vertex stores three concentrations whose sum cannot exceed 255, which is
what makes the reconstructed fourth concentration non-negative.

Returns `model`, or throws `ArgumentError` describing the first violation.
"""
function validate_model(model::PigmentModel; simplex::Bool = true)
    n = grid_n(model)
    length(model.forward.data) == 3 * n^3 || throw(ArgumentError(
        "forward table has $(length(model.forward.data)) bytes, expected $(3 * n^3)"
    ))
    length(model.inverse.data) == 3 * n^3 || throw(ArgumentError(
        "inverse table has $(length(model.inverse.data)) bytes, expected $(3 * n^3)"
    ))
    model.forward.n == n || throw(ArgumentError("table sizes disagree"))
    model.inverse.n == n || throw(ArgumentError("table sizes disagree"))
    if simplex
        _check_simplex(model.inverse, n)
    end
    return model
end

function _check_simplex(lut::ByteLUT, n::Int)
    d = lut.data
    vertices = n^3
    @inbounds for v in 0:(vertices - 1)
        base = 3v + 1
        s = Int(d[base]) + Int(d[base + 1]) + Int(d[base + 2])
        s <= 255 || throw(ArgumentError(
            "inverse table vertex $v stores concentrations summing to $s > 255; " *
                "the implied fourth concentration would be negative"
        ))
    end
    return nothing
end

"""
    write_model(io::IO, model)
    write_model(path::AbstractString, model)

Serialize `model` in the `.pmx` format. The header is rewritten from the
model's own bytes, so the stored checksums always match the stored payload.
"""
function write_model(io::IO, model::PigmentModel)
    validate_model(model)
    write(io, _header_bytes(model))
    write(io, model.inverse.data)
    write(io, model.forward.data)
    return nothing
end

write_model(path::AbstractString, model::PigmentModel) =
    open(io -> write_model(io, model), path, "w")

"""
    model_to_bytes(model) -> Vector{UInt8}

Serialize `model` into a new byte vector.
"""
function model_to_bytes(model::PigmentModel)
    return _header_bytes_and_payload(model)
end

function _header_bytes_and_payload(model::PigmentModel)
    out = IOBuffer(; sizehint = HEADER_BYTES + 6 * grid_n(model)^3)
    write_model(out, model)
    return take!(out)
end

"""
    read_model(path) -> PigmentModel

Read, checksum, and validate a `.pmx` payload from disk.
"""
read_model(path::AbstractString) = model_from_bytes(read(path))

"""
    default_payload_path() -> String

Path of the packaged default model: `<package>/data/default/default.pmx`.
The file is produced by `PaintMixPrecompute` and is not present in a
development checkout until a model has been promoted.
"""
function default_payload_path()
    root = pkgdir(@__MODULE__)
    root === nothing && (root = normpath(joinpath(@__DIR__, "..")))
    return joinpath(root, "data", "default", "default.pmx")
end

const _DEFAULT_MODEL = Ref{Union{Nothing,PigmentModel}}(nothing)

"""
    default_model() -> PigmentModel

The packaged four-pigment model, loaded once on first use.

Throws if the payload is missing, which is the case in a checkout where the
precompute pipeline has not produced and promoted a model yet; the error
names the path it expected so the remedy is obvious.
"""
function default_model()
    m = _DEFAULT_MODEL[]
    m === nothing || return m
    path = default_payload_path()
    isfile(path) || throw(ArgumentError(
        "no default pigment model at $path; build one with " *
            "`precompute/scripts/generate.jl` and promote it into data/default/"
    ))
    loaded = read_model(path)
    _DEFAULT_MODEL[] = loaded
    return loaded
end
