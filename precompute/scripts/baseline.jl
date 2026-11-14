#!/usr/bin/env julia
# Run this unchanged against the original and cleaned source trees.
using PaintMixPrecompute
using TOML
using SHA
include(joinpath(@__DIR__, "..", "test", "fixtures.jl"))

function main(args = ARGS)
    cfg, db, spectra, quad, model = synthetic_model()
    tables = Dict{String, Any}()
    for solver in (:reference, :bulk)
        values = generate_inverse(model, 8; threads = 1, solver)
        tables[string(solver)] = bytes2hex(sha256(reinterpret(UInt8, values)))
    end
    forward = generate_forward(model, 8; threads = 1)
    tables["forward"] = bytes2hex(sha256(reinterpret(UInt8, forward)))
    report = Dict(
        "schema_version" => 1, "julia_version" => string(VERSION),
        "synthetic_samples" => 6, "grid_n" => 8, "table_sha256" => tables
    )
    mkpath(dirname(args[1]))
    return open(args[1], "w") do io
        TOML.print(io, report; sorted = true)
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
