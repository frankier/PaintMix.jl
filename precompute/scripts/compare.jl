#!/usr/bin/env julia
# Independent candidates; never promotes or rewrites the release payload.
using PaintMixPrecompute
using TOML
using PaintMix

const VARIANTS = Dict(
    "softplus-lbfgs" => ("softplus", "lbfgs"),
    "softplus-lbfgsb" => ("softplus", "lbfgsb"),
    "direct-lbfgsb" => ("direct", "lbfgsb"),
)

const MATCHED_VARIANTS = merge(VARIANTS, Dict("gauge-lbfgsb" => ("gauge", "lbfgsb")))

function comparison_config(name; matched = false)
    variants = matched ? MATCHED_VARIANTS : VARIANTS
    haskey(variants, name) || error("unknown variant $name")
    cfg = load_config()
    cfg["surrogate"]["parameterization"], cfg["surrogate"]["optimizer"] = variants[name]
    if matched
        cfg["surrogate"]["stop_rule"] = "heldout_gamut"
        cfg["surrogate"]["gamut_tolerance"] = 1.0e-3
        cfg["surrogate"]["alpha_final"] = 1.0e-6
    end
    return cfg
end

function comparison_summary(out; matched = false)
    runs = Dict{String, Any}[]
    for name in sort!(collect(keys(matched ? MATCHED_VARIANTS : VARIANTS)))
        dir = joinpath(out, name)
        artifact = TOML.parsefile(joinpath(dir, "artifact.toml"))
        d = TOML.parsefile(joinpath(dir, "$(artifact["model_id"]).toml"))
        history = d["surrogate"]["history"]
        physical = findfirst(h -> get(h, "heldout_target_reached", false), history)
        reached = findfirst(h -> h["Epush"] <= d["surrogate"]["settings"]["push_tolerance"], history)
        push!(
            runs, Dict{String, Any}(
                "variant" => name, "model_id" => artifact["model_id"],
                "payload_sha256" => artifact["payload_sha256"], "fit_job_id" => d["fit_job_id"],
                "table_job_id" => d["table_job_id"], "table_settings" => d["table_settings"],
                "fit_seconds" => isempty(history) ? 0.0 : last(history)["seconds"],
                "generation_seconds" => d["generation_seconds"],
                "function_evaluations" => sum(get(h, "function_evaluations", 0) for h in history),
                "gradient_evaluations" => sum(get(h, "gradient_evaluations", 0) for h in history),
                "termination" => d["surrogate"]["diagnostics"]["termination"],
                "physical_target" => get(d["surrogate"]["diagnostics"], "physical_target", Dict()),
                "final_heldout_gamut_violation" => isempty(history) ? -1.0 : get(last(history), "heldout_gamut_violation", -1.0),
                "physical_target_reached" => get(d["surrogate"]["diagnostics"], "physical_target_reached", false),
                "first_physical_target_seconds" => physical === nothing ? -1.0 : history[physical]["seconds"],
                "push_target_reached" => reached !== nothing,
                "first_push_target_seconds" => reached === nothing ? -1.0 : history[reached]["seconds"],
                "diagnostics" => d["surrogate"]["diagnostics"],
                "validation" => d["validation"], "acceptance_gates" => d["acceptance_gates"],
                "environment" => d["environment"]
            )
        )
    end
    digests = [r["validation"]["quality"]["corpus_sha256"] for r in runs]
    all(==(first(digests)), digests) || error("comparison corpus mismatch")
    report = Dict("schema_version" => 2, "matched_accuracy" => matched, "corpus_sha256" => first(digests), "runs" => runs)
    open(joinpath(out, "comparison.toml"), "w") do io
        TOML.print(io, report; sorted = true)
    end
    open(joinpath(out, "comparison.md"), "w") do io
        println(io, "# Surrogate and release table comparison\n")
        println(io, "All candidates use the same measured inputs and validation corpus. No candidate is promoted.\n")
        tail = matched ? " Held-out max violation | Physical target reached |" : ""
        println(io, "| Variant | Fit seconds | Table seconds | Epush dense | Dense max violation | Mixing Oklab p99 | Encoded channel p99 |", tail)
        println(io, "| --- | ---: | ---: | ---: | ---: | ---: | ---: |", matched ? " ---: | --- |" : "")
        for r in runs
            q = r["validation"]["quality"]
            d = r["diagnostics"]
            row_values = (
                r["variant"], r["fit_seconds"], r["generation_seconds"], d["Epush_dense"],
                d["max_cube_violation"], q["oklab_error"]["p99"], q["encoded_channel_error"]["p99"],
            )
            matched && (row_values = (row_values..., r["final_heldout_gamut_violation"], r["physical_target_reached"]))
            println(io, "| ", join(row_values, " | "), " |")
        end
        matched && println(io, "\nMatched runs stop at held-out gamut violation <= 1e-3, or at the common schedule/time budget. Gauge coordinates remove one common K/S scale per wavelength.\n")
        println(io, "\nThe softplus/L-BFGS-B row isolates the solver change; the direct/L-BFGS-B row then isolates coefficient coordinates. These are local, budget-limited fits. Layer errors are reported separately and are not additive.")
    end
    return report
end

function main(args = ARGS)
    matched = "--matched" in args
    args = filter(!=("--matched"), args)
    length(args) >= 2 || error("usage: compare.jl OUT fit|tables|all|summary [VARIANT] [dev|release]")
    out, stage = abspath(args[1]), args[2]
    stage == "summary" && return comparison_summary(out; matched)
    name = args[3]
    profile = length(args) >= 4 ? args[4] : "release"
    cfg = comparison_config(name; matched)
    st = TableSettings(cfg; profile, threads = Threads.nthreads())
    return run_pipeline(
        cfg; out_dir = joinpath(out, name), profile, stage = Symbol(stage),
        settings = st, promote = false
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
