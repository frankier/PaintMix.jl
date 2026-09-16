# Payload format tests: header conventions, validation, and the packaged
# default model loader.

# Build a payload for a model that `write_model` would refuse to validate: the
# header is written directly, so the file is well-formed but invalid.
function BytesForTest(model::PigmentModel)
    io = IOBuffer()
    write(io, PaintMix._header_bytes(model))
    write(io, model.inverse.data)
    write(io, model.forward.data)
    return take!(io)
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

    wrong_n = copy(good)
    PaintMix._put_u32!(wrong_n, 16, 5)
    @test_throws InvalidPayload model_from_bytes(wrong_n)

    wrong_storage = copy(good)
    PaintMix._put_u8!(wrong_storage, 12, PaintMix.STORAGE_F64)
    @test_throws InvalidPayload model_from_bytes(wrong_storage)

    wrong_scale = copy(good)
    PaintMix._put_u8!(wrong_scale, 40, 0x01)
    @test_throws InvalidPayload model_from_bytes(wrong_scale)

    # Table bytes are not checked. A truncated payload is caught by the
    # header offsets, but a flipped table byte is accepted; only the simplex
    # invariant, checked separately, can reject a table at load time.
    flipped = copy(good)
    flipped[end] ⊻= 0xff
    @test model_from_bytes(flipped) isa PigmentModel
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
        good.id, ByteLUT(n, inv), good.forward, FORMAT_VERSION, good.flags
    )
    @test_throws ArgumentError validate_model(broken)
    bytes = BytesForTest(broken)
    @test_throws ArgumentError model_from_bytes(bytes)
    # Without the simplex check the payload is structurally acceptable.
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
