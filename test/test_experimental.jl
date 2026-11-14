function experimental_fixture(U; packed = true)
    n = 4
    q = U == UInt8 ? U(64) : U(16384)
    inverse = U == UInt8 ? ByteLUT(n, fill(q, 3n^3)) : WordLUT(n, fill(q, 3n^3))
    data = UInt8[]
    for k in 0:(n - 1), j in 0:k, i in 0:j
        append!(data, UInt8.((85i, 85(j - i), 85(k - j))))
    end
    forward = packed ? SimplexLUT(n, data) : ByteLUT(n, fill(UInt8(127), 3n^3))
    return ExperimentalModel(ntuple(i -> UInt8(i), 16), inverse, forward)
end

function experimental_allocations(model, ::Type{T}) where {T}
    a, b = SVector{3, T}(0.2, 0.3, 0.4), SVector{3, T}(0.8, 0.1, 0.2)
    mix(model, a, b, T(0.3))
    return @allocated mix(model, a, b, T(0.3))
end

@testset "experimental simplex interpolation and runtime API" begin
    for U in (UInt8, UInt16), T in (Float32, Float64)
        model = experimental_fixture(U)
        tolerance = T == Float32 ? 2.0e-6 : 1.0e-13
        for a in 0:7, b in 0:(7 - a), c in 0:(7 - a - b)
            concentration = SVector{4, T}(a, b, c, 7 - a - b - c) / T(7)
            got = PaintMix.forward_rgb(model, concentration)
            @test got ≈ concentration[1:3] atol = tolerance
        end
        a, b = SVector{3, T}(0.2, 0.3, 0.4), SVector{3, T}(0.8, 0.1, 0.2)
        @test decode(model, encode(model, a)) ≈ a atol = tolerance
        @test mix(model, a, b, zero(T)) == a
        @test mix(model, a, b, one(T)) == b
        @test mix(model, a, b, T(0.3)) ≈ T(0.7) * a + T(0.3) * b atol = tolerance
        @test weighted_mix(model, [a, b], T[2, 3]) ≈ T(0.4) * a + T(0.6) * b atol = tolerance
        output = zeros(T, 6)
        bulk_mix!(output, model, vcat(a, b), vcat(b, a), fill(T(0.3), 2))
        @test output[1:3] ≈ mix(model, a, b, T(0.3)) atol = tolerance
        @test experimental_allocations(model, T) == 0
    end
end

@testset "experimental portable format and integrity" begin
    for U in (UInt8, UInt16), packed in (false, true)
        model = experimental_fixture(U; packed)
        bytes = model_to_bytes(model)
        loaded = ExperimentalModel(bytes)
        @test loaded.inverse.data == model.inverse.data
        @test loaded.forward.data == model.forward.data
        @test model_to_bytes(loaded) == bytes
        @test model_id(loaded) == model_id(model)
        @test_throws InvalidPayload PigmentModel(bytes)
        for count in (0, 7, 127, 129, length(bytes) - 1)
            @test_throws InvalidPayload ExperimentalModel(bytes[1:count])
        end
        @test_throws InvalidPayload ExperimentalModel(vcat(bytes, UInt8(0)))
        for index in (1, 9, 13, 14, 16, 17, 21, 41, 49, 57, 89, 93, length(bytes))
            bad = copy(bytes)
            bad[index] ⊻= 0x80
            @test_throws InvalidPayload ExperimentalModel(bad)
        end
        mktempdir() do dir
            path = joinpath(dir, "candidate.pmxe")
            write_model(path, model)
            @test read_model(path) isa ExperimentalModel
            @test model_to_bytes(read_model(path)) == bytes
        end
    end
    model = experimental_fixture(UInt16)
    bytes = model_to_bytes(model)
    @test bytes[129:130] == [0x00, 0x40]
    # A valid checksum does not bypass the simplex constraint.
    bytes[129:134] .= 0xff
    PaintMix._put_u32!(bytes, 88, PaintMix._experimental_checksum(bytes))
    @test_throws ArgumentError ExperimentalModel(bytes)
    @test_throws ArgumentError WordLUT(0, UInt16[])
    @test_throws ArgumentError SimplexLUT(4, UInt8[])
end
