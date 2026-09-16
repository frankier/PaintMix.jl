# Latent algebra tests: residuals, round trips, endpoints, and the bulk and
# weighted mixing kernels.

using Random: MersenneTwister

@testset "latent is an isbits value with tuple fields" begin
    @test isbitstype(Latent{Float32})
    @test isbitstype(Latent{Float64})
    z = Latent((0.25, 0.25, 0.25, 0.25), (0.0, 0.0, 0.0))
    @test concentrations(z) === (0.25, 0.25, 0.25, 0.25)
    @test residual(z) === (0.0, 0.0, 0.0)
    @test occursin("Latent", sprint(show, z))
end

@testset "round trips cancel the table quantization" begin
    rng = MersenneTwister(11)
    for model in (identity_model(2), identity_model(16), random_model(9, 5)), T in (Float32, Float64)
        tol = T === Float32 ? 2.0f-6 : 1e-12
        for _ in 1:200
            x = T(rand(rng))
            y = T(rand(rng))
            z = T(rand(rng))
            got = decode(model, encode(model, (x, y, z)))
            @test eltype(got) === T
            @test maximum(abs.(got .- (x, y, z))) <= tol
        end
    end
end

@testset "encode returns a simplex latent" begin
    rng = MersenneTwister(3)
    model = random_model(3, 6)
    for _ in 1:100
        z = encode(model, (rand(rng), rand(rng), rand(rng)))
        c = concentrations(z)
        @test all(>=(0), c)
        @test sum(c) <= 1 + 8 * eps(Float64)
        # Stored exactly as `1 - (c1 + c2 + c3)`; reassociating the sum can
        # differ by an ulp, and `decode` reconstructs rather than trusting it.
        @test abs(c[4] - (1 - c[1] - c[2] - c[3])) <= 4 * eps(Float64)
    end
end

@testset "same-color mixing and exact binary endpoints" begin
    model = identity_model(16)
    a = (0.31, 0.62, 0.07)
    b = (0.1, 0.2, 0.3)
    for t in (0.0, 0.25, 0.5, 1.0)
        @test mix(model, a, a, t) == a
    end
    @test mix(model, a, b, 0.0) === a
    @test mix(model, a, b, 1.0) === b
    @test mix(model, a, b, -1.0) === a
    @test mix(model, a, b, 2.0) === b
    @test mix(model, a, b, 0) === a
    @test mix(model, a, b, 1) === b
    @test mix(model, a, b, 1 // 2) isa NTuple{3,Float64}
end

@testset "mixing is symmetric under swapping the arguments" begin
    model = identity_model(16)
    a = (0.31, 0.62, 0.07)
    b = (0.11, 0.22, 0.73)
    for t in (0.1, 0.25, 0.5, 0.9)
        @test maximum(abs.(mix(model, a, b, t) .- mix(model, b, a, 1 - t))) <= 1e-12
    end
end

@testset "ties between a color and its own latents are exact" begin
    # `mix(a, b, t)` with t in (0, 1) equals decoding the lerp of the two
    # latents; re-encoding that result and decoding again is stable.
    model = identity_model(8)
    a = (0.9, 0.2, 0.4)
    b = (0.05, 0.8, 0.3)
    m = mix(model, a, b, 0.35)
    @test maximum(abs.(decode(model, encode(model, m)) .- m)) <= 1e-12
end

@testset "decode does not clip" begin
    model = offset_model(8)
    # A byte-valued forward table cannot leave [0, 1]; the excursion has to
    # come from the residual, and decode must keep it.
    z = Latent((0.25, 0.5, 0.25, 0.0), (0.0, 0.9, -0.4))
    got = decode(model, z)
    @test got[2] > 1
    @test got[3] < 0
    @test got == decode(model, z)
    @test srgb8_from_linear(got)[2] == 0xff
    @test srgb8_from_linear(got)[3] == 0x00
    @test all(srgb8_from_linear((-1.0, 2.0, 0.5)) .== (0x00, 0xff, 0xbc))
end

@testset "concentration repair keeps the simplex" begin
    model = identity_model(4)
    z = Latent((0.5, 0.4, 0.3, -0.2), (0.0, 0.0, 0.0))
    got = decode(model, z)
    @test all(isfinite, got)
    neg = Latent((0.5, 0.4, 0.1, 0.0), (0.0, 0.0, 0.0))
    @test all(isfinite, decode(model, neg))
    @test decode(model, neg) == PaintMix.forward_rgb(model, 0.5, 0.4, 0.1)
end

@testset "bulk and scalar kernels agree" begin
    rng = MersenneTwister(5)
    model = random_model(5, 7)
    n = 50
    as = rand(rng, 3n)
    bs = rand(rng, 3n)
    ts = rand(rng, n)
    dest = fill(NaN, 3n)
    @test bulk_mix!(dest, model, as, bs, ts) === dest
    for i in 1:n
        want = mix(model, (as[3i - 2], as[3i - 1], as[3i]),
            (bs[3i - 2], bs[3i - 1], bs[3i]), ts[i])
        @test (dest[3i - 2], dest[3i - 1], dest[3i]) == want
    end
end

@testset "weighted mixing" begin
    model = identity_model(16)
    a = (0.31, 0.62, 0.07)
    b = (0.1, 0.2, 0.3)
    # Two colors reduce to the binary mix.
    @test maximum(
        abs.(weighted_mix(model, [a, b], [1 - 0.3, 0.3]) .- mix(model, a, b, 0.3))
    ) <= 1e-12
    # Scale invariance.
    @test maximum(
        abs.(
            weighted_mix(model, [a, b], [2.0, 5.0]) .-
                weighted_mix(model, [a, b], [20.0, 50.0])
        )
    ) <= 1e-12
    # A single positive weight returns that input color unchanged.
    @test weighted_mix(model, [a, b], [0.0, 3.0]) === b
    colors = [(0.1, 0.2, 0.3), (0.7, 0.1, 0.2), (0.05, 0.05, 0.9)]
    @test weighted_mix(model, colors, [0.0, 0.0, 5.0]) === colors[3]
    # Zero weights are skipped, and accumulation order is the input order.
    @test weighted_mix(model, [a, b, b], [1.0, 0.0, 0.0]) === a
    # n colors flatten correctly.
    flat = reduce(vcat, collect.(colors))
    dest = zeros(3)
    weighted_mix!(dest, model, flat, [1.0, 1.0, 1.0])
    @test maximum(
        abs.((dest[1], dest[2], dest[3]) .- weighted_mix(model, colors, [1.0, 1.0, 1.0]))
    ) <= 1e-15
    @test all(isfinite, weighted_mix(model, colors, [1.0, 0.5, 2.5]))
end

@testset "input validation throws" begin
    model = identity_model(4)
    a = (0.5, 0.5, 0.0)
    b = (0.0, 0.5, 0.5)
    @test_throws DomainError encode(model, (NaN, 0.5, 0.5))
    @test_throws DomainError encode(model, (0.5, Inf, 0.5))
    @test_throws DomainError mix(model, a, b, NaN)
    @test_throws DomainError mix(model, a, b, Inf)
    @test_throws DomainError mix(model, (NaN, 0.0, 0.0), b, 0.5)
    @test_throws ArgumentError weighted_mix(model, [a, b], [-1.0, 2.0])
    @test_throws ArgumentError weighted_mix(model, [a, b], [0.0, 0.0])
    @test_throws ArgumentError weighted_mix(model, [a, b], [NaN, 1.0])
    @test_throws ArgumentError weighted_mix(model, [a, b], [Inf, 1.0])
    @test_throws DimensionMismatch weighted_mix(model, [a, b], [1.0])
    @test_throws DimensionMismatch bulk_mix!(zeros(6), model, zeros(6), zeros(5), zeros(2))
    @test_throws ArgumentError mix!(zeros(2), model, a, b, 0.5)
    @test_throws ArgumentError weighted_mix!(zeros(2), model, zeros(3), [1.0])
end

@testset "status kernels report codes instead of throwing" begin
    model = identity_model(4)
    as = zeros(6)
    bs = zeros(6)
    ts = zeros(2)
    dest = zeros(6)
    @test PaintMix.bulk_mix_kernel!(dest, model, as, bs, ts, 2) == PM_OK
    @test PaintMix.bulk_mix_kernel!(dest, model, as, bs, ts, -1) == PM_ERR_LENGTH
    @test PaintMix.bulk_mix_kernel!(zeros(5), model, as, bs, ts, 2) == PM_ERR_LENGTH
    @test PaintMix.bulk_mix_kernel!(dest, model, zeros(5), bs, ts, 2) == PM_ERR_LENGTH
    @test PaintMix.bulk_mix_kernel!(dest, model, as, bs, zeros(1), 2) == PM_ERR_LENGTH
    bad = copy(as)
    bad[1] = NaN
    @test PaintMix.bulk_mix_kernel!(dest, model, bad, bs, ts, 2) == PM_ERR_NONFINITE
    @test PaintMix.bulk_mix_kernel!(dest, model, as, bs, [Inf, 0.0], 2) == PM_ERR_NONFINITE
    weights = [1.0, -1.0]
    @test PaintMix.weighted_mix_kernel!(zeros(3), model, as, weights, 2) == PM_ERR_WEIGHT
    @test PaintMix.weighted_mix_kernel!(zeros(3), model, as, [0.0, 0.0], 2) == PM_ERR_TOTAL
    @test PaintMix.weighted_mix_kernel!(zeros(3), model, as, [NaN, 1.0], 2) == PM_ERR_NONFINITE
    @test PaintMix.weighted_mix_kernel!(zeros(3), model, zeros(5), [1.0, 1.0], 2) == PM_ERR_LENGTH
    # A failed call leaves the output untouched.
    sentinel = fill(42.0, 6)
    @test PaintMix.bulk_mix_kernel!(sentinel, model, bad, bs, ts, 2) == PM_ERR_NONFINITE
    @test all(==(42.0), sentinel)
    @test error_message(PM_ERR_WEIGHT) == "negative weight"
    @test error_message(PM_OK) == "ok"
    @test occursin("unknown", error_message(9999))
    @test PaintMix.bulk_mix_kernel!(dest, model, as, bs, ts, 0) == PM_OK
end

@testset "Float32 path" begin
    model = identity_model(16)
    a = (0.31f0, 0.62f0, 0.07f0)
    b = (0.1f0, 0.2f0, 0.3f0)
    z = encode(model, a)
    @test z isa Latent{Float32}
    @test decode(model, z) isa NTuple{3,Float32}
    @test mix(model, a, b, 0.5f0) isa NTuple{3,Float32}
    @test maximum(abs.(decode(model, z) .- a)) <= 2.0f-6
    as = Float32[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    dest = zeros(Float32, 6)
    @test bulk_mix!(dest, model, as, as, Float32[0.5, 0.5]) === dest
    @test eltype(dest) === Float32
end
