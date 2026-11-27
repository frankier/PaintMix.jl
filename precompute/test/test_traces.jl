const TracePC = PaintMixPrecompute

@testset "fit trace cache and resumed export" begin
    cfg, db, spectra, quad, _ = synthetic_model()
    cfg["grid"]["dev_n"] = 3
    cfg["surrogate"]["surface_divisions"] = 2
    cfg["surrogate"]["max_iterations"] = 2
    cfg["surrogate"]["alpha_halvings"] = 2
    cfg["surrogate"]["push_tolerance"] = 0.0
    cfg["unmix"]["max_iterations"] = 5
    inputs = (db, spectra, quad)
    st = TableSettings(cfg; strategy = :bulk, evaluator = :float_lut, threads = 1)

    mktempdir() do dir
        off = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :fit, log = nothing)
        @test off["fit"].trace === nothing
        on = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :fit, log = nothing, trace_fit = true)
        fit = on["fit"]
        trace = fit.trace
        @test length(trace["snapshots"]) == 4
        @test first(trace["snapshots"])["step"] == 0
        @test first(trace["snapshots"])["Epull"] ≈ 0 atol = 1.0e-20
        @test sum(trace["boundary"]["weights"]) ≈ 1
        @test intersect(trace["boundary"]["concentrations"], trace["validation"]["concentrations"]) |> isempty
        @test all(row -> isfinite(row["sampled_max_cube_violation"]), trace["snapshots"])
        for (i, snapshot) in enumerate(trace["snapshots"])
            @test snapshot["step"] == i - 1
            replay = TracePC.replay_fit_snapshot(trace, i)
            @test collect(eachcol(replay.K)) == snapshot["absorption"]
            @test collect(eachcol(replay.S)) == snapshot["scattering"]
        end
        @test TracePC.replay_fit_snapshot(trace, 4).K == permutedims(fit.K)
        resolved = profile_config(cfg, "dev")
        repeated = fit_surrogates(resolved, spectra, quad; trace_fit = true)
        without_time(row) = Dict(k => v for (k, v) in row if k != "seconds")
        @test without_time.(repeated.trace["snapshots"]) == without_time.(trace["snapshots"])

        baseline = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :tables, log = nothing)
        result = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :tables, log = nothing, trace_inverse = true)
        @test baseline["model_id"] == result["model_id"]
        @test keys(result["traces"]) == Set(["fit", "inverse"])
        sidecar = TracePC.TOML.parsefile(result["sidecar_path"])
        @test all(r -> r["resumed"], values(sidecar["slab_summaries"]))
        exported = TracePC.read_trace(dir, sidecar["traces"]["fit"], "fit")
        @test last(exported["snapshots"])["absorption"] == sidecar["parameters"]["absorption"]
        @test last(exported["snapshots"])["scattering"] == sidecar["parameters"]["scattering"]
        inverse = TracePC.read_trace(dir, sidecar["traces"]["inverse"], "inverse")
        @test inverse["payload_sha256"] == exported["payload_sha256"] == sidecar["payload_sha256"]
        @test length(inverse["probes"]) <= inverse["probe_limit"]
        @test inverse["evaluator"] == "float_lut"
        for probe in inverse["probes"]
            runtime = probe["runtime"]
            @test sum(runtime["rounded_bytes_of_diagnostic_winner"]) == 255
            @test sum(vertex["weight"] for vertex in runtime["inverse_cell_vertices"]) ≈ 1
            @test all(v -> sum(v["stored_bytes"]) == 255, runtime["inverse_cell_vertices"])
            @test runtime["runtime_rgb"] + runtime["runtime_residual"] ≈ probe["target"]
            @test runtime["forward_interpolation_max_channel_delta"] >= 0
        end
        again = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :tables, log = nothing, trace_inverse = true)
        @test again["traces"] == result["traces"]
        ordinary = run_pipeline(cfg; inputs, settings = st, out_dir = dir, stage = :tables, log = nothing)
        @test ordinary["traces"] == result["traces"]
        path = joinpath(dir, result["traces"]["inverse"]["file"])
        open(path, "a") do io
            print(io, "\n")
        end
        @test_throws InputError TracePC.read_trace(dir, result["traces"]["inverse"], "inverse")
    end

    cfg["surrogate"]["alpha_halvings"] = 64
    @test_throws ConfigError fit_surrogates(cfg, spectra, quad; trace_fit = true)
end

@testset "traced inverse uses the normal solve" begin
    model, _ = affine_lut()
    settings = UnmixSettings{Float64}(50, 1.0e-10, 1.0e-6, 4)
    for rgb in (
            eval_mix(model, SVector(0.25, 0.25, 0.25, 0.25)),
            eval_mix(model, SVector(0.0, 1.0, 0.0, 0.0)),
            SVector(0.0, 0.0, 0.99),
        )
        seeds = TracePC.diagnostic_seeds(model, rgb)
        probe = TracePC.diagnose_inverse(model, rgb, settings; seeds)
        bulk = unmix_bulk!(SolverScratch(), model, rgb, seeds, settings)
        @test probe["concentrations"] == collect(bulk.c)
        @test probe["sse"] == bulk.sse
        @test probe["evaluations"] == bulk.function_evaluations
        @test last(probe["seeds"][probe["winning_seed"]]["faces"])["concentrations"] == probe["concentrations"]
        count = 0
        for seed in probe["seeds"], face in seed["faces"]
            last_sse = face["initial_sse"]
            for row in face["accepted"]
                c = row["concentrations"]
                @test all(>=(0), c) && sum(c) ≈ 1
                @test row["sse"] < last_sse
                last_sse = row["sse"]
                count += 1
            end
            @test isfinite(face["simplex_projected_gradient_inf"])
        end
        @test count <= probe["accepted_iterate_limit"]
        @test isfinite(probe["reference"]["sse_gap"])
    end
end
