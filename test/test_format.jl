# Payload format tests: header conventions, checksums, validation, and the
# packaged default model loader.

# Build a payload for a model that `write_model` would refuse to validate: the
# header is written directly, so the file is well-formed but invalid.
function BytesForTest(model::PigmentModel)
    io = IOBuffer()
    write(io, PaintMix._header_bytes(model))
    write(io, model.inverse.data)
    write(io, model.forward.data)
    return take!(io)
end

@testset "crc32 matches the reference vectors" begin
    @test crc32(UInt8[]) == 0x00000000
    @test crc32(Vector{UInt8}(codeunits("123456789"))) == 0xcbf43926
    @test crc32(Vector{UInt8}(codeunits("The quick brown fox jumps over the lazy dog"))) ==
        0x414fa339
    @test crc32(UInt8[0x00]) == 0xd202ef8d
    # Length argument.
    @test crc32(Vector{UInt8}(codeunits("123456789")), 4) ==
        crc32(Vector{UInt8}(codeunits("1234")))
end

@testset "payload round trip" begin
    model = random_model(17, 4)
    bytes = model_to_bytes(model)
    @test length(bytes) == HEADER_BYTES + 6 * 4^3
    back = model_from_bytes(bytes)
    @test back.inverse.data == model.inverse.data
    @test back.forward.data == model.forward.data
    @test back.id == model.id
    @test back.flags == model.flags
    @test back.format_version == FORMAT_VERSION
    @test model_id(back) == model_id(model)
    @test grid_n(back) == 4
    @test back.inverse_crc32 == crc32(model.inverse.data)
    @test back.forward_crc32 == crc32(model.forward.data)
    # Re-serializing a parsed model is byte-identical.
    @test model_to_bytes(back) == bytes
end

@testset "header conventions are recorded in the payload" begin
    model = identity_model(3)
    h = PaintMix._parse_header(model_to_bytes(model))
    @test h.header_bytes == HEADER_BYTES
    @test h.storage == PaintMix.STORAGE_U8
    @test h.channels == 3
    @test h.table_count == 2
    @test h.color_space == PaintMix.COLORSPACE_LINEAR_SRGB_D65
    @test h.byte_scale == PaintMix.BYTE_SCALE_255
    @test h.interpolation == PaintMix.INTERP_TRILINEAR
    @test h.index_order == PaintMix.INDEX_CHANNEL_FAST
    @test h.checksum == PaintMix.CHECKSUM_CRC32
    @test h.grid_n == 3
    @test h.inverse_offset == HEADER_BYTES
    @test h.forward_offset == HEADER_BYTES + 3 * 3^3
    @test h.inverse_bytes == 3 * 3^3
    @test h.forward_bytes == 3 * 3^3
end

@testset "file round trip" begin
    model = random_model(21, 3)
    mktempdir() do dir
        path = joinpath(dir, "model.pmx")
        write_model(path, model)
        back = read_model(path)
        @test back.forward.data == model.forward.data
        @test back.inverse.data == model.inverse.data
    end
end

@testset "malformed payloads are rejected" begin
    model = random_model(23, 3)
    good = model_to_bytes(model)

    bad_magic = copy(good)
    bad_magic[1] = 0x00
    @test_throws InvalidPayload model_from_bytes(bad_magic)

    version = copy(good)
    PaintMix._put_u16!(version, 8, 99)
    @test_throws InvalidPayload model_from_bytes(version)

    short = good[1:(end - 1)]
    @test_throws InvalidPayload model_from_bytes(short)
    @test_throws InvalidPayload model_from_bytes(UInt8[])

    corrupt = copy(good)
    corrupt[end] ⊻= 0xff
    @test_throws InvalidPayload model_from_bytes(corrupt)

    # A corrupted header checksum field is caught the same way.
    corrupt_crc = copy(good)
    PaintMix._put_u32!(corrupt_crc, 44, 0xdeadbeef)
    @test_throws InvalidPayload model_from_bytes(corrupt_crc)

    wrong_n = copy(good)
    PaintMix._put_u32!(wrong_n, 16, 5)
    @test_throws InvalidPayload model_from_bytes(wrong_n)

    wrong_storage = copy(good)
    PaintMix._put_u8!(wrong_storage, 12, PaintMix.STORAGE_F64)
    @test_throws InvalidPayload model_from_bytes(wrong_storage)

    wrong_scale = copy(good)
    PaintMix._put_u8!(wrong_scale, 40, 0x01)
    @test_throws InvalidPayload model_from_bytes(wrong_scale)

    # `checksum = false` skips the CRC-32 comparison. It is a promise that the
    # caller verified the bytes elsewhere, so a corrupted table loads and the
    # header's checksums are carried through unchanged.
    loaded = model_from_bytes(corrupt; checksum = false)
    @test loaded.forward_crc32 == model.forward_crc32
    @test loaded.forward.data != model.forward.data
    @test loaded.inverse.data == model.inverse.data
    # It is structurally valid, so the simplex check still runs.
    @test validate_model(loaded) === loaded
end

@testset "simplex validation" begin
    good = random_model(29, 3)
    @test validate_model(good) === good

    # Hand-build an inverse table whose first vertex is off the simplex.
    n = 3
    payload = 3 * n^3
    inv = copy(good.inverse.data)
    inv[1] = 0xff
    inv[2] = 0xff
    inv[3] = 0xff
    broken = PigmentModel(
        good.id, ByteLUT(n, inv), good.forward, FORMAT_VERSION, good.flags,
        good.forward_crc32, crc32(inv),
    )
    @test_throws ArgumentError validate_model(broken)
    bytes = BytesForTest(broken)
    @test_throws ArgumentError model_from_bytes(bytes)
    # The checksum is valid, so this really is the simplex check firing.
    @test model_from_bytes(bytes; validate = false) isa PigmentModel
end


@testset "every stored concentration triple is a valid simplex vertex" begin
    for model in (random_model(31, 6), identity_model(4), offset_model(5))
        n = grid_n(model)
        for v in 0:(n^3 - 1)
            base = 3v + 1
            d = model.inverse.data
            @test Int(d[base]) + Int(d[base + 1]) + Int(d[base + 2]) <= 255
        end
    end
end

@testset "default model reports a missing payload clearly" begin
    # The release payload is produced by the precompute pipeline. Until it is
    # promoted, `default_model()` must explain what is missing.
    if !isfile(default_payload_path())
        err = try
            default_model()
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("default.pmx", sprint(showerror, err))
        @test_throws ArgumentError encode((0.1, 0.2, 0.3))
    else
        @test default_model() === default_model()
    end
end
