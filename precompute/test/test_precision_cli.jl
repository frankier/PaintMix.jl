module PrecisionCLITests
    include(joinpath(@__DIR__, "..", "scripts", "precision.jl"))
end

@testset "precision CLI uses authoritative coefficients and verified slabs" begin
    cfg, db, spectra, quad, _ = synthetic_model()
    cfg["grid"]["dev_n"] = 4
    cfg["surrogate"]["surface_divisions"] = 2
    cfg["surrogate"]["max_iterations"] = 2
    cfg["surrogate"]["alpha_halvings"] = 0
    cfg["validation"]["colors"] = 12
    cfg["validation"]["pairs"] = 32
    mktempdir() do dir
        candidate, output = joinpath(dir, "candidate"), joinpath(dir, "precision")
        source = run_pipeline(cfg; inputs = (db, spectra, quad), out_dir = candidate, promote = false, log = nothing)
        report = PrecisionCLITests.main([candidate, output]; inputs = (db, spectra, quad))
        @test report["source_model_id"] == source["model_id"]
        @test length(report["runs"]) == 4
        for row in report["runs"]
            path = joinpath(output, row["variant"] * ".pmxe")
            @test model_id(read_model(path)) == row["model_id"]
            @test row["report"]["quality"]["corpus_sha256"] == report["corpus_sha256"]
            @test row["report"]["runtime"]["allocated_bytes"] == 0
        end
    end
end
