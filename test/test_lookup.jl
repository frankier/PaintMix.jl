# Lookup-table kernel tests: index arithmetic, interpolation weights, and
# the exact grid-endpoint rules.

using Random: MersenneTwister

const trilinear = PaintMix.trilinear

@testset "table size validation" begin
    @test_throws ArgumentError ByteLUT(0, UInt8[])
    @test_throws ArgumentError ByteLUT(2, UInt8[0x00])
    @test ByteLUT(2, zeros(UInt8, 24)).n == 2
end

@testset "trilinear matches an independent implementation" begin
    for n in (2, 3, 5, 8), seed in (1, 7)
        lut = random_model(seed, n).forward
        rng = MersenneTwister(seed + 100)
        for _ in 1:200
            x, y, z = rand(rng), rand(rng), rand(rng)
            got = trilinear(lut, x, y, z)
            ref = reference_trilinear(lut, x, y, z)
            @test maximum(abs.(got .- ref)) <= 1e-15
        end
    end
end

@testset "trilinear reproduces stored vertices exactly" begin
    lut = random_model(2, 3).forward
    for k in 0:2, j in 0:2, i in 0:2
        x = i / 2
        y = j / 2
        z = k / 2
        got = trilinear(lut, x, y, z)
        want = ntuple(ch -> Float64(reference_vertex(lut, i, j, k, ch - 1)) / 255, Val(3))
        # A grid plane must reproduce the stored value: no interpolation
        # error. Division and multiplication by 1/255 may differ by one ulp,
        # so the check is exact to within that.
        @test maximum(abs.(got .- want)) <= 1e-15
    end
end

@testset "exact grid endpoints never read past the buffer" begin
    for n in (1, 2, 3, 17), model in (random_model(3, n),)
        lut = model.forward
        for (x, y, z) in (
                (1.0, 1.0, 1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                (1.0, 1.0, 0.0), (0.5, 1.0, 1.0), (1e-300, 1.0, 0.0),
            )
            @test all(isfinite, trilinear(lut, x, y, z))
            @test trilinear(lut, x, y, z) == trilinear(lut, x, y, z)
        end
        # Weight one on the last plane must equal the last vertex exactly.
        last = ntuple(ch -> Float64(reference_vertex(lut, n - 1, n - 1, n - 1, ch - 1)) / 255, Val(3))
        @test maximum(abs.(trilinear(lut, 1.0, 1.0, 1.0) .- last)) <= 1e-15
    end
end

@testset "coordinates outside [0, 1] clamp to the boundary" begin
    lut = random_model(4, 4).forward
    @test trilinear(lut, -0.5, 0.25, 0.75) == trilinear(lut, 0.0, 0.25, 0.75)
    @test trilinear(lut, 0.25, 2.0, 0.75) == trilinear(lut, 0.25, 1.0, 0.75)
    @test trilinear(lut, 0.25, 0.75, -1e9) == trilinear(lut, 0.25, 0.75, 0.0)
end

@testset "affine tables interpolate analytically" begin
    # Both fixtures store only 0x00 and 0xff, so trilinear interpolation of
    # the corner samples is the affine map itself, bit for bit.
    p = permuted_model(2)
    for (c1, c2, c3) in ((0.25, 0.5, 0.125), (1.0, 0.0, 0.0), (0.3, 0.7, 0.9))
        @test trilinear(p.forward, c1, c2, c3) == (c3, c1, c2)
    end
    s = signed_model(2)
    for (c1, c2, c3) in ((0.25, 0.5, 0.125), (1.0, 1.0, 1.0), (0.0, 0.0, 0.0))
        @test trilinear(s.forward, c1, c2, c3) == (1 - c1, c2, 1 - c3)
    end
end

@testset "Float32 tables stay in Float32" begin
    lut = random_model(5, 4).forward
    got = trilinear(lut, 0.3f0, 0.6f0, 0.1f0)
    @test got isa NTuple{3,Float32}
    @test maximum(abs.(got .- trilinear(lut, 0.3, 0.6, 0.1))) <= 1e-6
    @test PaintMix.vertex(lut, 1, 2, 3) isa NTuple{3,UInt8}
    @test_throws BoundsError PaintMix.vertex(lut, 4, 0, 0)
end
