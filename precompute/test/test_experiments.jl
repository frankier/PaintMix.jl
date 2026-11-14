@testset "experiment identity and physical parameterizations" begin
    cfg, db, spectra, quad, model = synthetic_model()
    @test experiment_id(Dict("b" => 2, "a" => 1)) == experiment_id(Dict("a" => 1, "b" => 2))
    @test experiment_id(Dict("a" => [1, 2])) != experiment_id(Dict("a" => [2, 1]))
    @test_throws ConfigError FitSettings(merge(cfg, Dict("surrogate" => merge(cfg["surrogate"], Dict("time_limit_seconds" => Inf)))))
    soft = SoftplusParameters(1.0e-6)
    direct = DirectParameters(1.0e-6)
    xs = initial_parameters(soft, spectra)
    xd = initial_parameters(direct, spectra)
    Ks, Ss = parameter_matrices(soft, xs, length(spectra.wavelength))
    Kd, Sd = parameter_matrices(direct, xd, length(spectra.wavelength))
    @test Ks ≈ Kd
    @test Ss ≈ Sd
    st = FitSettings(cfg)
    problem = SurrogateProblem(st, spectra, quad)
    direct_cfg = deepcopy(cfg)
    direct_cfg["surrogate"]["parameterization"] = "direct"
    direct_cfg["surrogate"]["optimizer"] = "lbfgsb"
    dp = SurrogateProblem(FitSettings(direct_cfg), spectra, quad)
    @test collect(objective_terms(problem, xs)) ≈ collect(objective_terms(dp, xd))
    f = u -> sum(objective_terms(problem, u))
    ad = ForwardDiff.gradient(f, xs)
    fd = similar(ad)
    PaintMixPrecompute._finite_gradient!(fd, f, xs, -Inf)
    @test ad ≈ fd rtol = 1.0e-4 atol = 1.0e-7
    opt = optimization_problem(dp, xd, 1.0)
    @test all(==(direct.epsilon), opt.lb)
    @test all(isinf, opt.ub)
    id = fit_job_id(cfg, db)
    @test fit_job_id(direct_cfg, db) != id
    @test table_job_id(
        cfg, SurrogateFit(Ks, Ss, xs, Dict{String, Any}[], Dict{String, Any}(), SurfaceQuadrature(2)),
        TableSettings(cfg; strategy = :bulk, threads = 1)
    ) ==
        table_job_id(
        cfg, SurrogateFit(Ks, Ss, xs, Dict{String, Any}[], Dict{String, Any}(), SurfaceQuadrature(2)),
        TableSettings(cfg; strategy = :bulk, threads = 2)
    )
end

@testset "Optimization backends and continuation" begin
    cfg, db, spectra, quad, model = synthetic_model()
    cfg["surrogate"]["surface_divisions"] = 2
    cfg["surrogate"]["max_iterations"] = 3
    cfg["surrogate"]["alpha_halvings"] = 0
    cfg["surrogate"]["alpha_initial"] = 1.0
    cfg["surrogate"]["alpha_final"] = 1.0
    for (param, alg) in (("softplus", "lbfgs"), ("softplus", "lbfgsb"), ("direct", "lbfgsb"))
        cfg["surrogate"]["parameterization"] = param
        cfg["surrogate"]["optimizer"] = alg
        fit = fit_surrogates(cfg, spectra, quad)
        @test all(isfinite, fit.K)
        @test all(>(0), fit.S)
        @test length(fit.history) == 1
        @test haskey(fit.history[1], "termination")
        mktempdir() do dir
            path = joinpath(dir, "fit.toml")
            save_surrogate(path, fit, cfg, "job", length(spectra.wavelength))
            loaded = load_surrogate(path, cfg, "job", length(spectra.wavelength))
            @test loaded.K == fit.K
            @test load_surrogate(path, cfg, "other", length(spectra.wavelength)) === nothing
        end
    end
    cfg["surrogate"]["parameterization"] = "softplus"
    cfg["surrogate"]["optimizer"] = "lbfgs"
    fit = fit_surrogates(cfg, spectra, quad; algorithm = PaintMixPrecompute.Optim.BFGS())
    @test all(isfinite, fit.K)
    for algorithm in (PaintMixPrecompute.Optim.NelderMead(), PaintMixPrecompute.Optim.Newton())
        candidate = fit_surrogates(cfg, spectra, quad; algorithm)
        @test all(isfinite, candidate.K)
    end
    automatic = fit_surrogates(cfg, spectra, quad; adtype = PaintMixPrecompute.Optimization.AutoForwardDiff())
    @test all(isfinite, automatic.K)
    cfg["surrogate"]["derivative"] = "finite_difference"
    finite = fit_surrogates(cfg, spectra, quad)
    @test all(isfinite, finite.K)
    st = FitSettings(cfg)
    problem = SurrogateProblem(st, spectra, quad)
    x = initial_parameters(problem.parameters, spectra)
    _, status = solve_fit_step(problem, x, 1.0, PaintMixPrecompute.Optim.LBFGS(), st; deadline = time() - 1)
    @test status["termination"] == "TimeLimit"
    @test !status["converged"]
end

@testset "versioned corpus, decomposition, and integrity" begin
    samples = validation_corpus(; colors = 16, pairs = 512)
    @test samples.digest == validation_corpus(; colors = 16, pairs = 512).digest
    @test any(p -> p[1] == 16 || p[2] == 16, samples.pairs)
    @test all(c -> minimum(c) >= 0 && sum(c) ≈ 1, samples.concentrations)
    cfg, db, spectra, quad, model = synthetic_model()
    n = 4
    fwd = generate_forward(model, n; threads = 1)
    st = TableSettings(n, :reference, :spectral, 2, 3, 1)
    summaries = Dict{Int, Any}()
    inv = generate_tables_inverse(model, st; settings = UnmixSettings(cfg), on_summary = (k, r) -> (summaries[k] = r))
    @test length(summaries) == n
    @test sum(r["points"] for r in values(summaries)) == n^3
    @test inv ≈ generate_inverse(model, n; threads = 1, settings = UnmixSettings(cfg))
    ft = FloatTables(n, fwd, inv)
    payload = build_model(ft, cfg, Dict("inputs" => Dict{String, String}(), "config_hash" => "x"))
    changed = FloatTables(n, copy(fwd), copy(inv))
    changed.forward[1] = 1 - changed.forward[1]
    @test model_id(build_model(changed, cfg, Dict())) != model_id(payload)
    report = error_decomposition(payload, model, ft, cfg; measured = model, samples)
    @test report["surrogate_bias_oklab"]["max"] == 0
    @test haskey(report, "concentration_quantization_channel")
    mktempdir() do dir
        cp = CheckpointStore(dir, "job", n, "inverse")
        values = inv[1:(3n^2)]
        write_slab(cp, 0, values)
        @test slab_complete(cp, 0)
        path = joinpath(dir, "inverse_k00000.bin")
        open(path, "r+") do io
            write(io, UInt8(0xff))
        end
        @test !slab_complete(cp, 0)
        @test_throws InputError CheckpointStore(dir, "job", n + 1, "inverse")
    end
end

@testset "independent inverse solver adapter" begin
    cfg, db, spectra, quad, model = synthetic_model()
    solver = OptimizationInverse(PaintMixPrecompute.Optim.BFGS())
    scratch = solver_scratch(solver)
    c = SVector(0.2, 0.3, 0.1, 0.4)
    rgb = mix_rgb(model, c)
    result = unmix_reference(model, rgb; scratch, settings = UnmixSettings{Float64}(100, 1.0e-10, 1.0e-6, 2))
    @test result.sse < 1.0e-8
    @test sum(result.c) ≈ 1
    @test result.function_evaluations > 0
end

@testset "reusable stages verify cached identities" begin
    cfg, db, spectra, quad, model = synthetic_model()
    cfg["grid"]["dev_n"] = 4
    cfg["surrogate"]["surface_divisions"] = 2
    cfg["surrogate"]["max_iterations"] = 2
    cfg["surrogate"]["alpha_halvings"] = 0
    cfg["validation"]["colors"] = 12
    cfg["validation"]["pairs"] = 32
    inputs = (db, spectra, quad)
    mktempdir() do dir
        result = run_pipeline(cfg; inputs, out_dir = dir, promote = false, log = nothing)
        @test isfile(result["payload_path"])
        @test haskey(result["validation"], "decomposition")
        validated = run_pipeline(cfg; inputs, out_dir = dir, stage = :validate, log = nothing)
        @test validated["model_id"] == result["model_id"]
        changed = deepcopy(cfg)
        changed["unmix"]["max_iterations"] += 1
        @test_throws InputError run_pipeline(changed; inputs, out_dir = dir, stage = :validate, log = nothing)
        open(result["payload_path"], "a") do io
            write(io, UInt8(0))
        end
        @test_throws InputError run_pipeline(cfg; inputs, out_dir = dir, stage = :validate, log = nothing)
    end
end

@testset "Optimization inverse with a table evaluator" begin
    model, _ = affine_lut()
    c = SVector(0.2, 0.3, 0.1, 0.4)
    rgb = eval_mix(model, c)
    scratch = solver_scratch(OptimizationInverse(PaintMixPrecompute.Optim.BFGS()))
    r = unmix_reference(model, rgb; scratch)
    @test r.sse < 1.0e-8
    @test sum(r.c) ≈ 1
end

function slab_allocations(model)
    st = TableSettings(8, :coarse, :spectral, 4, 15, 1)
    coarse = generate_inverse(model, 4; threads = 1)
    dest = Vector{Float64}(undef, 3 * st.n^3)
    ctx = (SolverScratch(), Vector{SVector{4, Float64}}(undef, st.n), Vector{SVector{4, Float64}}(undef, st.n))
    settings = UnmixSettings{Float64}(15, 1.0e-10, 1.0e-6, 1)
    report = SlabSummary()
    run! = () -> PaintMixPrecompute._inverse_strategy_slab!(Val(:coarse), dest, model, st, 0, ctx, settings, coarse, report)
    run!()
    return @allocated run!()
end

@testset "default fine slab is allocation-free" begin
    _, _, _, _, model = synthetic_model()
    @test slab_allocations(model) == 0
end

@testset "workers own scratch and preserve numerical results" begin
    if Threads.nthreads() >= 2
        owned = fill(false, 32)
        PaintMixPrecompute._run_slabs!(
            () -> objectid(current_task()),
            (owner, k) -> begin
                owned[k + 1] = owner == objectid(current_task())
                yield()
            end, 32, 2
        )
        @test all(owned)
        cfg, _, _, _, model = synthetic_model()
        for strategy in (:reference, :bulk, :coarse)
            serial = generate_tables_inverse(model, TableSettings(8, strategy, :spectral, 4, 15, 1))
            parallel = generate_tables_inverse(model, TableSettings(8, strategy, :spectral, 4, 15, 2))
            @test serial == parallel
        end
        @test generate_forward(model, 8; threads = 1) == generate_forward(model, 8; threads = 2)
    else
        @test_skip "run with at least two Julia threads for worker ownership checks"
    end
end

mutable struct StatusScratch
    function_evaluations::Int
    gradient_evaluations::Int
end
function PaintMixPrecompute._polish!(
        sc::StatusScratch, model, rgb,
        active::NTuple{K, Int}, c, settings
    ) where {K}
    if K == 4
        return SVector(0.5, 0.25, 0.25 - 1.0e-9, 1.0e-9), 1.0, 1, true
    end
    return SVector(0.5, 0.25, 0.25, 0.0), 0.9, 1, false
end
@testset "convergence belongs to the final selected face" begin
    scratch = StatusScratch(0, 0)
    seed = SVector(0.25, 0.25, 0.25, 0.25)
    result = unmix_bulk!(scratch, nothing, SVector(0.0, 0.0, 0.0), (seed,), UnmixSettings(load_config()))
    @test result.nactive == 3
    @test !result.converged
end
