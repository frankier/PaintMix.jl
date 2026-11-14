#!/usr/bin/env julia
# CLI adapter. Numerical stages and artifact verification live in the package.
using ArgParse
using PaintMixPrecompute

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const STAGES = ("inputs", "fit", "tables", "validate", "promote", "all")

function parse_args(args)
    s = ArgParseSettings(prog = "generate.jl", description = "Generate a PaintMix payload.")
    @add_arg_table! s begin
        "--config"
        help = "configuration file"
        arg_type = String
        default = PaintMixPrecompute.DEFAULT_CONFIG_PATH
        "--profile"
        help = "dev or release"
        arg_type = String
        default = "dev"
        "--out"
        help = "artifact directory (default: precompute/output/<profile>)"
        arg_type = String
        "--stage"
        help = "inputs|fit|tables|validate|promote|all"
        arg_type = String
        default = "all"
        "--solver"
        help = "coarse|reference|bulk inverse solver"
        arg_type = String
        default = "coarse"
        "--evaluator"
        help = "spectral|float_lut inverse evaluator"
        arg_type = String
        default = "spectral"
        "--coarse"
        help = "coarse grid edge for the coarse solver (default: config)"
        arg_type = Int
        "--threads"
        help = "worker threads for table generation"
        arg_type = Int
        default = Threads.nthreads()
        "--no-promote"
        help = "never copy into data/default/"
        action = :store_false
        dest_name = "promote"
        "--force"
        help = "promote even when the provisional quality targets fail"
        action = :store_true
        "--check-inputs"
        help = "verify input checksums and exit"
        action = :store_true
        "--force-import"
        help = "re-import the DuckDB database"
        action = :store_true
    end
    opts = ArgParse.parse_args(args, s)
    opts["profile"] in ("dev", "release") ||
        error("--profile must be dev or release, got $(opts["profile"])")
    opts["stage"] in STAGES || error("--stage must be one of $(join(STAGES, ", "))")
    opts["solver"] in ("reference", "bulk", "coarse") ||
        error("--solver must be reference, bulk, or coarse")
    return opts
end

function main(args = ARGS)
    opts = parse_args(args)
    cfg = load_config(opts["config"])
    if opts["check-inputs"]
        checked, failures = check_checksums(ROOT, joinpath(ROOT, cfg["inputs"]["checksums"]))
        isempty(failures) || error(join(failures, "\n"))
        println("checked ", checked, " input files")
        return nothing
    end
    st = TableSettings(
        cfg; profile = opts["profile"], strategy = Symbol(opts["solver"]),
        evaluator = Symbol(opts["evaluator"]), coarse_n = opts["coarse"], threads = opts["threads"]
    )
    out = something(opts["out"], joinpath(ROOT, "precompute", "output", opts["profile"]))
    result = run_pipeline(
        cfg; out_dir = out, profile = opts["profile"], stage = Symbol(opts["stage"]),
        settings = st, promote = opts["promote"], force = opts["force"],
        force_import = opts["force-import"], cfg_path = opts["config"]
    )
    for (name, ok) in sort!(collect(get(result, "acceptance_gates", Dict())); by = first)
        println(name, ": ", ok ? "pass" : "FAIL")
    end
    haskey(result, "payload_path") && println("candidate: ", result["payload_path"])
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
