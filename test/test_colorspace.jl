# Encoded-sRGB adapter tests: transfer function, byte scaling, and the
# explicit clipping and rounding rule.

@testset "transfer function endpoints and known values" begin
    @test linear_from_srgb((0.0, 0.0, 0.0)) == (0.0, 0.0, 0.0)
    @test linear_from_srgb((1.0, 1.0, 1.0)) == (1.0, 1.0, 1.0)
    @test srgb_from_linear((0.0, 0.0, 0.0)) == (0.0, 0.0, 0.0)
    # The power segment evaluates to 1 - 1 ulp at the white point, which is
    # why the byte adapter, not this function, owns the clamping rule.
    @test maximum(abs.(srgb_from_linear((1.0, 1.0, 1.0)) .- 1.0)) <= 2 * eps(1.0)
    # Mid-encoded-gray is the classic 0.2159 linear value.
    @test linear_from_srgb((0.5, 0.5, 0.5))[1] ≈ 0.21404114048223255 atol = 1e-12
    @test srgb_from_linear((0.21404114048223255, 0.21404114048223255, 0.21404114048223255))[1] ≈
        0.5 atol = 1e-12
    # The linear segment below 0.04045 is continuous with the power segment.
    for t in (0.0, 1e-6, 0.03, 0.04044, 0.04046, 0.5, 0.9, 1.0)
        @test srgb_from_linear(linear_from_srgb((t, t, t)))[1] ≈ t atol = 1e-12
    end
end

@testset "byte scaling is b/255 with exact endpoints" begin
    @test linear_from_srgb8((0x00, 0x00, 0x00)) == (0.0f0, 0.0f0, 0.0f0)
    @test linear_from_srgb8((0xff, 0xff, 0xff)) == (1.0f0, 1.0f0, 1.0f0)
    @test srgb8_from_linear((0.0f0, 0.0f0, 0.0f0)) == (0x00, 0x00, 0x00)
    @test srgb8_from_linear((1.0f0, 1.0f0, 1.0f0)) == (0xff, 0xff, 0xff)
end

@testset "byte round trip is lossless for every byte value" begin
    for b in 0x00:0xff
        @test srgb8_from_linear(linear_from_srgb8((b, b, b))) == (b, b, b)
    end
    # Random triples in Float64 as well.
    for r in 0x00:0x11:0xff, g in 0x00:0x33:0xff, b in 0x00:0x55:0xff
        @test srgb8_from_linear(linear_from_srgb8((r, g, b))) == (r, g, b)
    end
end

@testset "byte conversion is the only clipping step" begin
    @test srgb8_from_linear((-0.5, 0.5, 1.5)) == (0x00, 0xbc, 0xff)
    # The linear tuple is untouched: this is the path a mixer must use when
    # it wants unclipped values.
    lo, hi = -1.0, 2.0
    @test srgb_from_linear((lo, 0.5, hi))[1] < 0
    @test srgb_from_linear((lo, 0.5, hi))[3] > 1
end

@testset "Float32 adapters stay in Float32" begin
    c = linear_from_srgb8((0x80, 0x40, 0x00))
    @test c isa NTuple{3,Float32}
    @test srgb8_from_linear(c) == (0x80, 0x40, 0x00)
    e = srgb_from_linear(c)
    @test e isa NTuple{3,Float32}
    @test maximum(abs.(linear_from_srgb(e) .- c)) <= 1e-6
end
