# Encoded-sRGB adapter tests: transfer function, byte scaling, and the
# explicit clipping and rounding rule.

@testset "transfer function endpoints and known values" begin
    black = SVector(0.0, 0.0, 0.0)
    white = SVector(1.0, 1.0, 1.0)
    @test linear_from_srgb(black) == black
    @test linear_from_srgb(white) == white
    @test srgb_from_linear(black) == black
    # The power segment evaluates to 1 - 1 ulp at the white point, which is
    # why the byte adapter, not this function, owns the clamping rule.
    @test maximum(abs.(srgb_from_linear(white) .- 1.0)) <= 2 * eps(1.0)
    # Mid-encoded-gray is the classic 0.2159 linear value.
    gray = SVector(0.5, 0.5, 0.5)
    @test linear_from_srgb(gray)[1] ≈ 0.21404114048223255 atol = 1.0e-12
    @test srgb_from_linear(SVector(0.21404114048223255, 0.21404114048223255, 0.21404114048223255))[1] ≈
        0.5 atol = 1.0e-12
    # The linear segment below 0.04045 is continuous with the power segment.
    for t in (0.0, 1.0e-6, 0.03, 0.04044, 0.04046, 0.5, 0.9, 1.0)
        @test srgb_from_linear(linear_from_srgb(SVector(t, t, t)))[1] ≈ t atol = 1.0e-12
    end
end

@testset "byte scaling is b/255 with exact endpoints" begin
    @test linear_from_srgb8(SVector(0x00, 0x00, 0x00)) == SVector(0.0f0, 0.0f0, 0.0f0)
    @test linear_from_srgb8(SVector(0xff, 0xff, 0xff)) == SVector(1.0f0, 1.0f0, 1.0f0)
    @test srgb8_from_linear(SVector(0.0f0, 0.0f0, 0.0f0)) == SVector(0x00, 0x00, 0x00)
    @test srgb8_from_linear(SVector(1.0f0, 1.0f0, 1.0f0)) == SVector(0xff, 0xff, 0xff)
end

@testset "byte round trip is lossless for every byte value" begin
    for b in 0x00:0xff
        @test srgb8_from_linear(linear_from_srgb8(SVector(b, b, b))) == SVector(b, b, b)
    end
    # Random triples in Float64 as well.
    for r in 0x00:0x11:0xff, g in 0x00:0x33:0xff, b in 0x00:0x55:0xff
        @test srgb8_from_linear(linear_from_srgb8(SVector(r, g, b))) == SVector(r, g, b)
    end
end

@testset "byte conversion is the only clipping step" begin
    @test srgb8_from_linear(SVector(-0.5, 0.5, 1.5)) == SVector(0x00, 0xbc, 0xff)
    # The linear vector is untouched: this is the path a mixer must use when
    # it wants unclipped values.
    lo, hi = -1.0, 2.0
    @test srgb_from_linear(SVector(lo, 0.5, hi))[1] < 0
    @test srgb_from_linear(SVector(lo, 0.5, hi))[3] > 1
end

@testset "Float32 adapters stay in Float32" begin
    c = linear_from_srgb8(SVector(0x80, 0x40, 0x00))
    @test c isa SVector{3, Float32}
    @test srgb8_from_linear(c) == SVector(0x80, 0x40, 0x00)
    e = srgb_from_linear(c)
    @test e isa SVector{3, Float32}
    @test maximum(abs.(linear_from_srgb(e) .- c)) <= 1.0e-6
end
