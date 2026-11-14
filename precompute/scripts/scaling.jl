#!/usr/bin/env julia
# Measure useful table parallelism before selecting SLURM resources.
using PaintMixPrecompute
using TOML

function main(args = ARGS)
    out = args[1]
    cfg = load_config()
    db, spectra, quad = pipeline_inputs(cfg)
    model = SpectralModel(spectra, quad)
    # Compilation precedes timing; multiple repetitions reduce scheduling noise.
    generate_forward(model, 8; threads = 1)
    rows = Dict{String, Any}[]
    for workers in filter(<=(Threads.nthreads()), [1, 2, 4, 8])
        times = [@elapsed generate_forward(model, 96; threads = workers) for _ in 1:3]
        st = TableSettings(16, :reference, :spectral, 2, 15, workers)
        generate_tables_inverse(model, TableSettings(2, :reference, :spectral, 2, 15, 1); settings = UnmixSettings(cfg))
        inverse_seconds = @elapsed generate_tables_inverse(model, st; settings = UnmixSettings(cfg))
        fine_settings = TableSettings(64, :coarse, :spectral, 8, 15, workers)
        generate_tables_inverse(model, TableSettings(8, :coarse, :spectral, 4, 15, workers); settings = UnmixSettings(cfg))
        fine_seconds = minimum([@elapsed generate_tables_inverse(model, fine_settings; settings = UnmixSettings(cfg)) for _ in 1:2])
        push!(rows, Dict("workers" => workers, "forward_seconds" => minimum(times), "inverse_seconds" => inverse_seconds, "fine_seconds" => fine_seconds))
    end
    mkpath(dirname(out))
    open(out, "w") do io
        TOML.print(
            io, Dict(
                "schema_version" => 1, "measurements" => rows, "cpu_name" => Sys.CPU_NAME,
                "slurm_job_id" => get(ENV, "SLURM_JOB_ID", ""), "julia_version" => string(VERSION)
            ); sorted = true
        )
    end
    return rows
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
