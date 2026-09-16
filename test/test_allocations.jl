# Steady-state allocation tests. Every kernel used by the compiled ABI must
# allocate nothing once compiled, so bulk work can run inside a render loop.

using Random: Xoshiro

@testset "steady-state allocations are zero" begin
    model = identity_model(16)
    a32 = (0.31f0, 0.62f0, 0.07f0)
    b32 = (0.1f0, 0.2f0, 0.3f0)
    a64 = (0.31, 0.62, 0.07)
    b64 = (0.1, 0.2, 0.3)

    dest64 = zeros(3)
    dest32 = zeros(Float32, 3)
    n = 64
    rng = Xoshiro(2)
    as = rand(rng, 3n)
    bs = rand(rng, 3n)
    ts = rand(rng, n)
    bulk_dest = zeros(3n)
    weights = rand(rng, n)
    wdest = zeros(3)
    as32 = Float32.(as)
    bs32 = Float32.(bs)
    ts32 = Float32.(ts)
    bulk_dest32 = zeros(Float32, 3n)

    # Warm up every specialization before measuring.
    encode(model, a64)
    encode(model, a32)
    decode(model, encode(model, a64))
    mix(model, a64, b64, 0.5)
    mix(model, a32, b32, 0.5f0)
    mix!(dest64, model, a64, b64, 0.5)
    mix!(dest32, model, a32, b32, 0.5f0)
    bulk_mix!(bulk_dest, model, as, bs, ts)
    bulk_mix!(bulk_dest32, model, as32, bs32, ts32)
    weighted_mix!(wdest, model, as, weights)
    PaintMix.trilinear(model.forward, 0.5, 0.5, 0.5)

    @test (@allocated encode(model, a64)) == 0
    @test (@allocated encode(model, a32)) == 0
    @test (@allocated decode(model, encode(model, a64))) == 0
    @test (@allocated mix(model, a64, b64, 0.5)) == 0
    @test (@allocated mix(model, a32, b32, 0.5f0)) == 0
    @test (@allocated mix!(dest64, model, a64, b64, 0.5)) == 0
    @test (@allocated mix!(dest32, model, a32, b32, 0.5f0)) == 0
    @test (@allocated bulk_mix!(bulk_dest, model, as, bs, ts)) == 0
    @test (@allocated bulk_mix!(bulk_dest32, model, as32, bs32, ts32)) == 0
    @test (@allocated weighted_mix!(wdest, model, as, weights)) == 0
    @test (@allocated PaintMix.trilinear(model.forward, 0.5, 0.5, 0.5)) == 0
end
