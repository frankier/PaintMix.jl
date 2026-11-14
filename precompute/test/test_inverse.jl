@testset "inverse solver recovers concentrations" begin
    cfg, _, _, _, model = synthetic_model()
    settings = UnmixSettings(cfg)
    scratch = SolverScratch()
    for c in (
            SVector(1.0, 0.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0, 0.0),
            SVector(0.0, 0.0, 1.0, 0.0), SVector(0.0, 0.0, 0.0, 1.0),
            SVector(0.25, 0.25, 0.25, 0.25), SVector(0.5, 0.0, 0.5, 0.0),
            SVector(0.6, 0.2, 0.1, 0.1),
        )
        rgb = mix_rgb(model, c)
        r = unmix_reference(model, rgb; settings = settings, scratch = scratch)
        @test r.sse < 1.0e-16
        @test maximum(abs.(collect(r.c) .- collect(c))) < 1.0e-6
    end
    # Bulk solver agrees with the reference on interior and boundary points.
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(0.5, 0.0, 0.5, 0.0),
            SVector(0.8, 0.1, 0.05, 0.05),
        )
        rgb = mix_rgb(model, c)
        ref = unmix_reference(model, rgb; settings = settings, scratch = scratch)
        bulk = unmix_bulk!(scratch, model, rgb, (c,), settings)
        @test bulk.sse <= ref.sse + 1.0e-10
    end
end


@testset "inverse solver allocates nothing per solve" begin
    cfg, _, _, _, model = synthetic_model()
    settings = UnmixSettings(cfg)
    scratch = SolverScratch()
    rgb = mix_rgb(model, SVector(0.4, 0.3, 0.2, 0.1))
    seed = SVector(0.4, 0.3, 0.2, 0.1)
    # Warm up every specialization the measured calls reach. The bulk path
    # runs once per grid vertex, so at 256^3 even a few bytes per solve is
    # gigabytes of garbage; it must stay at exactly zero. A `Union` of
    # tuple arities in the active set used to box here.
    @test unmix_bulk!(scratch, model, rgb, (seed,), settings) isa
        PaintMixPrecompute.UnmixResult
    @test unmix_reference(model, rgb; settings = settings, scratch = scratch) isa
        PaintMixPrecompute.UnmixResult
    @test (@allocated unmix_bulk!(scratch, model, rgb, (seed,), settings)) == 0
    @test (
        @allocated PaintMixPrecompute._solve_adaptive!(
            scratch, model, rgb, seed, settings,
        )
    ) == 0
    # The reference path still allocates a little inside its subset
    # recursion (measured 80 bytes). The bound is a regression guard, not
    # a claim of zero: the old union-typed active set cost 768.
    @test (
        @allocated unmix_reference(model, rgb; settings = settings, scratch = scratch)
    ) <= 256
end


@testset "reference solver agrees with an independent Levenberg-Marquardt" begin
    # `unmix_reference` and `lm_active!` are hand-written: their damping
    # schedule, stationarity test, and Jacobian convention are all local
    # choices that round-tripping cannot falsify. This cross-checks the
    # concentrations themselves against
    # `LeastSquaresOptim.LevenbergMarquardt`, driven by `ForwardDiff`
    # rather than the analytic Jacobian, on the same least-squares
    # objective.
    cfg, _, _, _, model = synthetic_model()
    settings = UnmixSettings(cfg)
    scratch = SolverScratch()
    # Softmax over three free logits spans the interior of the full
    # simplex, so this is the oracle for optima with every concentration
    # above the face threshold. Written out here rather than reusing
    # `lm_active!`'s parameterization, so the oracle stays independent.
    function oracle_softmax(θ1, θ2, θ3)
        m = max(0.0, θ1, θ2, θ3)
        a, b, c, d = exp(θ1 - m), exp(θ2 - m), exp(θ3 - m), exp(-m)
        s = a + b + c + d
        return SVector(a / s, b / s, c / s, d / s)
    end
    oracle_residual(model, rgb) = x -> begin
        m = mix_rgb(model, oracle_softmax(x[1], x[2], x[3]))
        [m[1] - rgb[1], m[2] - rgb[2], m[3] - rgb[3]]
    end
    function oracle_unmix(model, rgb)
        # The oracle's default tolerances (1e-8) stop well short of the
        # reference solver's stationarity tolerance, so tighten them.
        res = LeastSquaresOptim.optimize(
            oracle_residual(model, rgb), zeros(3),
            LeastSquaresOptim.LevenbergMarquardt();
            autodiff = :forward, x_tol = 1.0e-12, f_tol = 1.0e-12, g_tol = 1.0e-12,
        )
        x = res.minimizer
        return oracle_softmax(x[1], x[2], x[3]), res.ssr
    end
    # Colors the model produces: a zero-residual solution exists, so both
    # solvers should report the same concentrations.
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(0.4, 0.3, 0.2, 0.1),
            SVector(0.7, 0.1, 0.1, 0.1), SVector(0.55, 0.25, 0.15, 0.05),
        )
        rgb = mix_rgb(model, c)
        ref = unmix_reference(model, rgb; settings = settings, scratch = scratch)
        oc, ossr = oracle_unmix(model, rgb)
        @test ref.sse < 1.0e-16
        @test ossr < 1.0e-18
        @test maximum(abs.(collect(ref.c) .- collect(oc))) < 1.0e-8
        @test maximum(abs.(collect(ref.c) .- collect(c))) < 1.0e-6
    end
    # Colors the model cannot produce. The optimum is a compromise, so
    # only the objective is comparable: the oracle cannot reach a face
    # exactly, and the reference enumerates faces the oracle does not
    # search. The reference must not do worse.
    for rgb in (
            SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 1.0),
            SVector(0.1, 0.05, 0.9), SVector(0.9, 0.9, 0.05),
        )
        ref = unmix_reference(model, rgb; settings = settings, scratch = scratch)
        _, ossr = oracle_unmix(model, rgb)
        @test ref.sse <= ossr + 1.0e-10
    end
end


@testset "coarse-to-fine inverse matches the reference" begin
    cfg, _, _, _, model = synthetic_model()
    n = 9
    coarse_n = 6
    settings = UnmixSettings(cfg)
    fine = generate_inverse_coarse_to_fine(
        model, n; coarse_n = coarse_n, threads = 1,
        settings = UnmixSettings{Float64}(15, 1.0e-10, 1.0e-6, 1),
        coarse_settings = settings,
    )
    # Every vertex is a valid simplex point and decodes to the requested
    # color at least as well as the coarse seeding allows.
    sc = SolverScratch()
    d = n - 1
    for k in 0:d, j in 0:d, i in 0:d
        o = 3 * (i + n * (j + n * k))
        c1, c2, c3 = fine[o + 1], fine[o + 2], fine[o + 3]
        @test c1 >= -1.0e-9 && c2 >= -1.0e-9 && c3 >= -1.0e-9
        @test c1 + c2 + c3 <= 1 + 1.0e-9
        rgb = SVector(i / d, j / d, k / d)
        res = unmix_reference(model, rgb; settings = settings, scratch = sc)
        # The coarse seed is a starting point, not the answer; the fine
        # solve should land in the reference basin even when it does not
        # reach the same local minimum.
        @test sum(abs2, mix_rgb(model, SVector(c1, c2, c3, 1 - c1 - c2 - c3)) .- rgb) <=
            res.sse + 1.0e-3
    end
    # Slab threading does not change the result.
    fine_threaded = generate_inverse_coarse_to_fine(
        model, n; coarse_n = coarse_n, threads = 2,
        settings = UnmixSettings{Float64}(15, 1.0e-10, 1.0e-6, 1),
        coarse_settings = settings,
    )
    @test fine_threaded == fine
end
