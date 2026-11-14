@testset "surface quadrature" begin
    sq = SurfaceQuadrature(5)
    @test length(sq.points) == length(sq.weights) == length(sq.face)
    @test sum(sq.weights) ≈ 1.0
    for c in sq.points
        @test sum(c) ≈ 1.0 atol = 1.0e-12
        @test all(>=(0), c)
        @test count(==(0.0), c) >= 1
    end
end


@testset "softplus round trip and alpha schedule" begin
    for x in (-8.0, -1.0, 0.0, 1.0, 6.0)
        @test inv_softplus(softplus(x)) ≈ x atol = 1.0e-9
    end
    cfg = load_config()
    α = alpha_schedule(cfg)
    @test issorted(α; rev = true)
    @test first(α) == cfg["surrogate"]["alpha_initial"]
    @test last(α) >= cfg["surrogate"]["alpha_final"]
end
