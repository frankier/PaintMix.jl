const _EXPERIMENTAL_MAGIC = (0x50, 0x41, 0x49, 0x4e, 0x54, 0x4d, 0x58, 0x45)
const EXPERIMENTAL_FORMAT_VERSION = UInt16(1)
_inverse_bits(::ByteLUT) = 8
_inverse_bits(::WordLUT) = 16
_forward_mode(::ByteLUT) = 0
_forward_mode(::SimplexLUT) = 1
_inverse_scale(lut) = Int(typemax(eltype(lut.data)))

function validate_model(model::ExperimentalModel; simplex::Bool = true)
    n = grid_n(model)
    _experimental_grid(n)
    model.inverse.n == n || throw(ArgumentError("table sizes disagree"))
    simplex && _check_simplex(model.inverse, n, _inverse_scale(model.inverse))
    return model
end

function _experimental_checksum(bytes)
    header = collect(view(bytes, 1:HEADER_BYTES))
    _put_u32!(header, 88, 0)
    crc = Base._crc32c(header)
    io = IOBuffer(bytes; read = true, write = false)
    seek(io, HEADER_BYTES)
    return Base._crc32c(io, length(bytes) - HEADER_BYTES, crc)
end

function model_to_bytes(model::ExperimentalModel)
    validate_model(model)
    n = grid_n(model)
    inverse = collect(reinterpret(UInt8, htol.(model.inverse.data)))
    forward = model.forward.data
    header = zeros(UInt8, HEADER_BYTES)
    header[1:8] .= _EXPERIMENTAL_MAGIC
    _put_u16!(header, 8, EXPERIMENTAL_FORMAT_VERSION)
    _put_u16!(header, 10, HEADER_BYTES)
    _put_u8!(header, 12, _inverse_bits(model.inverse))
    _put_u8!(header, 13, _forward_mode(model.forward))
    _put_u8!(header, 14, 3)
    _put_u32!(header, 16, n)
    header[25:40] .= model.id
    _put_u32!(header, 40, _inverse_scale(model.inverse))
    _put_u32!(header, 44, 255)
    _put_u64!(header, 56, HEADER_BYTES)
    _put_u64!(header, 64, length(inverse))
    _put_u64!(header, 72, HEADER_BYTES + length(inverse))
    _put_u64!(header, 80, length(forward))
    bytes = vcat(header, inverse, forward)
    _put_u32!(bytes, 88, _experimental_checksum(bytes))
    return bytes
end

write_model(io::IO, model::ExperimentalModel) = (write(io, model_to_bytes(model)); nothing)
write_model(path::AbstractString, model::ExperimentalModel) = open(io -> write_model(io, model), path, "w")

"""Read and validate experimental bytes, including lengths, conventions and CRC32C."""
function ExperimentalModel(bytes::AbstractVector{UInt8})
    length(bytes) >= HEADER_BYTES || throw(InvalidPayload("short experimental header"))
    Tuple(bytes[1:8]) == _EXPERIMENTAL_MAGIC || throw(InvalidPayload("bad experimental magic"))
    _hdr_u16(bytes, 8) == EXPERIMENTAL_FORMAT_VERSION || throw(InvalidPayload("unsupported experimental version"))
    _hdr_u16(bytes, 10) == HEADER_BYTES || throw(InvalidPayload("bad experimental header size"))
    bits, mode = _hdr_u8(bytes, 12), _hdr_u8(bytes, 13)
    bits in (8, 16) && mode in (0, 1) || throw(InvalidPayload("unknown experimental storage"))
    _hdr_u8(bytes, 14) == 3 || throw(InvalidPayload("experimental channel count must be three"))
    (all(i -> iszero(bytes[i + 1]), (15, 20, 21, 22, 23)) && all(iszero, view(bytes, 49:56)) && all(iszero, view(bytes, 93:128))) || throw(InvalidPayload("unsupported experimental convention or reserved field"))
    n = Int(_hdr_u32(bytes, 16))
    2 <= n <= 1024 || throw(InvalidPayload("experimental grid out of range"))
    invbytes = 3 * n^3 * (Int(bits) ÷ 8)
    fwdbytes = 3 * (mode == 0 ? n^3 : _simplex_vertices(n))
    _hdr_u32(bytes, 40) == (bits == 8 ? 255 : 65535) && _hdr_u32(bytes, 44) == 255 || throw(InvalidPayload("unsupported experimental scaling"))
    (_hdr_u64(bytes, 56), _hdr_u64(bytes, 64), _hdr_u64(bytes, 72), _hdr_u64(bytes, 80)) == (HEADER_BYTES, invbytes, HEADER_BYTES + invbytes, fwdbytes) || throw(InvalidPayload("invalid experimental offsets or lengths"))
    length(bytes) == HEADER_BYTES + invbytes + fwdbytes || throw(InvalidPayload("experimental file length mismatch"))
    _hdr_u32(bytes, 88) == _experimental_checksum(bytes) || throw(InvalidPayload("experimental checksum mismatch"))
    data = bytes[(HEADER_BYTES + 1):(HEADER_BYTES + invbytes)]
    inverse = bits == 8 ? ByteLUT(n, data) : WordLUT(n, ltoh.(collect(reinterpret(UInt16, data))))
    data = bytes[(HEADER_BYTES + invbytes + 1):end]
    forward = mode == 0 ? ByteLUT(n, data) : SimplexLUT(n, data)
    id = ntuple(i -> bytes[24 + i], Val(16))
    return ExperimentalModel(id, inverse, forward)
end
