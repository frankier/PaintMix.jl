#!/usr/bin/env julia
# Generate a PaintMix payload.
#
#   julia -t auto --project=precompute scripts/generate.jl --profile dev --out /tmp/pmdev
#   julia -t auto --project=precompute scripts/generate.jl --profile release
#
# Plan steps 3-6: import the selected spectra, build the spectral reference,
# fit the surrogate pigments, solve the inverse, quantize both tables, export
# the `.pmx` payload and its provenance sidecar, validate against the
# acceptance gates, and promote a validated payload into `data/default/`.
#
# Stages:
#
#   inputs    import the spreadsheets into DuckDB and stop
#   fit       fit (or reuse) the surrogate pigments and stop
#   tables    generate both float tables and the payload
#   validate  run the quality, padding, round-trip, and quantization reports
#   promote   copy into data/default/ when every gate passes
#   all       everything above
#
# Expensive artifacts are cached under `--out` and reused only when the
# configuration, input hashes, grid, and surrogate settings match.

import PaintMix
using PaintMixPrecompute
using Printf: @printf
using SHA: sha256
using TOML: TOML

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const STAGES = ("inputs", "fit", "tables", "validate", "promote", "all")

function usage()
    return """
    usage: generate.jl [options]

      --config PATH     configuration file (default: config/default.toml)
      --profile NAME    dev or release (default: dev)
      --out DIR         artifact directory (default: precompute/output/<profile>)
      --stage NAME      inputs|fit|tables|validate|promote|all (default: all)
      --solver NAME     coarse|reference|bulk inverse solver (default: coarse)
      --coarse N        coarse grid edge for the coarse solver (default: config)
      --threads N       worker threads for table generation (default: Julia's)
      --no-promote      never copy into data/default/
      --force           promote even when the provisional quality targets fail,
                        recording which gates were overridden
      --check-inputs    verify input checksums and exit
      --force-import    re-import the DuckDB database
      --help
    """
end

function parse_args(args)
    opts = Dict{String, Any}(
        "config" => PaintMixPrecompute.DEFAULT_CONFIG_PATH,
        "profile" => "dev",
        "out" => nothing,
        "stage" => "all",
        "solver" => "coarse",
        "threads" => Threads.nthreads(),
        "coarse" => nothing,
        "promote" => true,
        "force" => false,
        "check-inputs" => false,
        "force-import" => false,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--help" || a == "-h"
            print(usage())
            exit(0)
        elseif a == "--no-promote"
            opts["promote"] = false
        elseif a == "--force"
            opts["force"] = true
        elseif a == "--check-inputs"
            opts["check-inputs"] = true
        elseif a == "--force-import"
            opts["force-import"] = true
        elseif a in ("--config", "--profile", "--out", "--stage", "--solver", "--threads", "--coarse")
            i < length(args) || error("$a needs a value")
            opts[a[3:end]] = args[i + 1]
            i += 1
        else
            error("unknown argument $a\n" * usage())
        end
        i += 1
    end
    opts["profile"] in ("dev", "release") ||
        error("--profile must be dev or release, got $(opts["profile"])")
    opts["stage"] in STAGES || error("--stage must be one of $(join(STAGES, ", "))")
    opts["solver"] in ("reference", "bulk", "coarse") ||
        error("--solver must be reference, bulk, or coarse")
    opts["threads"] = parse(Int, string(opts["threads"]))
    return opts
end

check_inputs(checksum_file) = check_checksums(ROOT, checksum_file)

"""
    profile_config(cfg, profile) -> (cfg, n)

The dev profile caps the surrogate effort so a smoke run stays minutes. Both
profiles share conventions, inputs, and the unmix settings.
"""
function profile_config(cfg::AbstractDict, profile::AbstractString)
    out = deepcopy(cfg)
    g = get(cfg, "generation", Dict{String, Any}())
    if profile == "dev"
        s = out["surrogate"]
        s["surface_divisions"] =
            min(Int(s["surface_divisions"]), Int(get(g, "dev_surface_divisions", 8)))
        s["max_iterations"] = min(Int(s["max_iterations"]), Int(get(g, "dev_max_iterations", 40)))
        s["alpha_halvings"] = min(Int(s["alpha_halvings"]), Int(get(g, "dev_alpha_halvings", 12)))
    end
    return out
end

"""
    fit_hash(cfg, db, profile) -> String

Identifies a surrogate-fitting job. Only the settings the fit actually reads
enter the hash, so changing, say, the inverse solver does not throw away a
converged fit.
"""
function fit_hash(cfg, db, profile)
    io = IOBuffer()
    print(io, "profile=", profile, "\n")
    print(io, "spectra=", repr(cfg["spectra"]), "\n")
    print(io, "color=", repr(cfg["color"]), "\n")
    print(io, "saunderson=", repr(cfg["saunderson"]), "\n")
    print(io, "pigments=", repr(cfg["pigments"]), "\n")
    print(io, "surrogate=", repr(cfg["surrogate"]), "\n")
    for f in sort(db.source_files; by = f -> f.role)
        print(io, f.role, "=", f.sha256, "\n")
    end
    print(io, "observer=", db.observer_source.sha256, "\n")
    return bytes2hex(sha256(take!(io)))
end

"""
    table_hash(fit_hash_value, cfg, profile, solver) -> String

Identifies a table-generation job: the fit plus everything that changes the
generated bytes. Checkpoints are reused only when this matches.
"""
function table_hash(fit_hash_value, cfg, profile, solver)
    io = IOBuffer()
    print(io, "fit=", fit_hash_value, "\n")
    print(io, "grid_n=", cfg["grid"]["release_n"], "\n")
    print(io, "quantization=", repr(cfg["quantization"]), "\n")
    print(io, "unmix=", repr(cfg["unmix"]), "\n")
    print(io, "generation=", repr(get(cfg, "generation", Dict{String, Any}())), "\n")
    print(io, "solver=", solver, "\n")
    return bytes2hex(sha256(take!(io)))
end

function save_surrogate(path, fit, cfg, hash, W)
    doc = Dict{String, Any}(
        "schema_version" => 1,
        "job_hash" => hash,
        "codes" => collect(pigment_codes(cfg)),
        "samples" => W,
        "theta" => collect(fit.theta),
        "absorption" => vec(fit.K),
        "scattering" => vec(fit.S),
        "history" => fit.history,
        "diagnostics" => fit.diagnostics,
    )
    open(path, "w") do io
        TOML.print(io, doc)
    end
    return path
end

# TOML returns `Vector{Any}` of `Dict{String,Any}` for an array of tables.
_as_dicts(v) = Dict{String, Any}[
    Dict{String, Any}(String(k) => x for (k, x) in pairs(d)) for d in v
]

function load_surrogate(path, cfg, hash, W)
    isfile(path) || return nothing
    doc = try
        TOML.parsefile(path)
    catch
        return nothing
    end
    get(doc, "job_hash", "") == hash || return nothing
    Int(get(doc, "samples", 0)) == W || return nothing
    theta = Float64[Float64(x) for x in doc["theta"]]
    K = reshape(Float64[Float64(x) for x in doc["absorption"]], 4, W)
    S = reshape(Float64[Float64(x) for x in doc["scattering"]], 4, W)
    return SurrogateFit(
        K, S, theta, _as_dicts(get(doc, "history", Any[])),
        Dict{String, Any}(String(k) => v for (k, v) in pairs(get(doc, "diagnostics", Dict{String, Any}()))),
        surface_quadrature(cfg),
    )
end

function main(args)
    opts = parse_args(args)
    cfg_path = opts["config"]
    cfg = load_config(cfg_path)
    profile = opts["profile"]
    pcfg = profile_config(cfg, profile)
    n = profile == "dev" ? Int(pcfg["grid"]["dev_n"]) : Int(pcfg["grid"]["release_n"])
    out = something(opts["out"], joinpath(ROOT, "precompute", "output", profile))
    mkpath(out)
    @printf("PaintMix generation\n")
    @printf("  config:   %s\n", cfg_path)
    @printf("  profile:  %s (n = %d, %d threads)\n", profile, n, opts["threads"])
    @printf("  out:      %s\n", out)

    checksum_file = joinpath(ROOT, cfg["inputs"]["checksums"])
    checked, failures = check_inputs(checksum_file)
    @printf("  inputs:   %d files checked\n", checked)
    if !isempty(failures)
        for f in failures
            @printf("    %s\n", f)
        end
        error("input checksum check failed")
    end
    opts["check-inputs"] && return nothing

    # --- step 3: inputs and the spectral reference -----------------------
    db_path = joinpath(ROOT, cfg["inputs"]["database"])
    if opts["force-import"] || !isfile(db_path)
        @printf("importing inputs into %s\n", db_path)
        db = load_spreadsheet_inputs(cfg, cfg_path; root = ROOT)
        save_database(db, db_path)
    else
        db = open_database(db_path)
    end
    spectra = load_spectra(cfg, db)
    quad = build_quadrature(cfg, db)
    base = spectral_model(spectra, quad)
    W = length(spectra.wavelength)
    @printf("  spectra:  %d wavelengths, k1 = %g, k2 = %g\n", W, spectra.k1, spectra.k2)
    opts["stage"] == "inputs" && return nothing

    hash = fit_hash(pcfg, db, profile)
    thash = table_hash(hash, pcfg, profile, opts["solver"])

    # --- step 4: surrogate fit -------------------------------------------
    surrogate_path = joinpath(out, "surrogate.toml")
    fit = nothing
    if opts["stage"] in ("all", "fit")
        fit = load_surrogate(surrogate_path, pcfg, hash, W)
        if fit === nothing
            @printf(
                "fitting surrogates (surface divisions = %d)\n",
                pcfg["surrogate"]["surface_divisions"]
            )
            fit = fit_surrogates(pcfg, spectra, quad; log = stdout)
            save_surrogate(surrogate_path, fit, pcfg, hash, W)
            @printf("  cached fit at %s\n", surrogate_path)
        else
            @printf("reusing cached surrogate fit at %s\n", surrogate_path)
        end
        @printf(
            "  diagnostics: max cube violation %.3g, max Oklab deviation %.3g\n",
            fit.diagnostics["max_cube_violation"], fit.diagnostics["max_oklab_deviation"]
        )
    end
    opts["stage"] == "fit" && return nothing
    fit === nothing && (fit = load_surrogate(surrogate_path, pcfg, hash, W))
    fit === nothing &&
        error("no cached surrogate fit at $surrogate_path; run --stage fit first")
    surrogate = with_parameters(base, fit.K, fit.S)

    result = nothing
    ft = nothing
    provenance = nothing
    payload_path = joinpath(out, "default.pmx")
    need_generate = opts["stage"] in ("all", "tables") || !isfile(payload_path)
    if need_generate
        @printf("generating forward table (n = %d)\n", n)
        t0 = time()
        forward = generate_forward(surrogate, n; threads = opts["threads"])
        @printf("  forward:  %.1f s\n", time() - t0)

        settings = unmix_settings(pcfg)
        if opts["solver"] == "coarse"
            g = get(pcfg, "generation", Dict{String, Any}())
            coarse_n = min(
                opts["coarse"] === nothing ? Int(get(g, "coarse_n", 64)) : parse(Int, string(opts["coarse"])),
                max(2, n ÷ 2),
            )
            fine_iter = Int(get(g, "fine_max_iterations", 15))
            @printf("generating inverse table coarse-to-fine (coarse n = %d)\n", coarse_n)
            thash = table_hash(hash, pcfg, profile, "coarse$(coarse_n)i$(fine_iter)")
            cp = CheckpointStore(joinpath(out, "checkpoints", thash[1:16]), thash, n, "inverse")
            done = completed_slabs(cp)
            isempty(done) || @printf("  resuming: %d/%d slabs already complete\n", length(done), n)
            t0 = time()
            inverse = generate_inverse_coarse_to_fine(
                surrogate, n; coarse_n = coarse_n, threads = opts["threads"],
                settings = UnmixSettings{Float64}(fine_iter, 1.0e-10, 1.0e-6, 1),
                coarse_settings = settings,
                resume = slab_resume_function(cp), on_slab = slab_callback_function(cp),
            )
        else
            @printf("generating inverse table with the %s solver\n", opts["solver"])
            cp = CheckpointStore(joinpath(out, "checkpoints", thash[1:16]), thash, n, "inverse")
            done = completed_slabs(cp)
            isempty(done) || @printf("  resuming: %d/%d slabs already complete\n", length(done), n)
            t0 = time()
            inverse = generate_inverse(
                surrogate, n;
                threads = opts["threads"], settings = settings, solver = Symbol(opts["solver"]),
                resume = slab_resume_function(cp), on_slab = slab_callback_function(cp),
            )
        end
        @printf("  inverse:  %.1f s\n", time() - t0)

        ft = FloatTables(n, forward, inverse)
        provenance = build_provenance(
            pcfg, cfg_path, db, fit, Dict{String, Any}(); grid_n = n
        )
        provenance["profile"] = profile
        result = export_model(ft, pcfg, provenance; out_dir = out)
        @printf("  model id: %s\n", result["model_id"])
        @printf("  payload:  %s\n", result["payload_path"])
    else
        @printf("reusing payload %s\n", payload_path)
        model = PaintMix.read_model(payload_path)
        provenance = build_provenance(
            pcfg, cfg_path, db, fit, Dict{String, Any}(); grid_n = n
        )
        provenance["profile"] = profile
        result = Dict{String, Any}(
            "model" => model,
            "model_id" => PaintMix.model_id(model),
            "payload_path" => payload_path,
            "sidecar_path" => joinpath(out, "$(PaintMix.model_id(model)).toml"),
        )
    end
    opts["stage"] == "tables" && return nothing
    result === nothing && error("run --stage tables first")

    # --- validation ------------------------------------------------------
    @printf("validating\n")
    model = result["model"]
    t0 = time()
    reports = Dict{String, Any}(
        "quality" => quality_report(model, surrogate, cfg),
        "padding" => padding_report(model, surrogate, cfg),
        "roundtrip" => roundtrip_report(model, cfg),
        "quantization" => quantization_report(model, cfg),
        "continuity" => continuity_report(model, cfg),
        "behavior" => behavior_report(model, surrogate, cfg),
    )
    gates = acceptance_gates(cfg, reports)
    @printf("  validation took %.1f s\n", time() - t0)
    for (name, ok) in sort(collect(gates))
        @printf("    %-34s %s\n", name, ok ? "pass" : "FAIL")
    end
    for (name, ok) in (
            ("blue_yellow_is_green", reports["behavior"]["blue_yellow_is_green"]),
            ("magenta_yellow_is_orange", reports["behavior"]["magenta_yellow_is_orange"]),
            ("white_tint_monotone", reports["behavior"]["white_tint_monotone"]),
        )
        @printf("    %-34s %s\n", name, ok ? "pass" : "FAIL")
        gates[name] = ok
    end
    provenance["validation"] = reports
    provenance["acceptance_gates"] = Dict{String, Any}(k => v for (k, v) in gates)
    write_sidecar(result["sidecar_path"], provenance, model, ft)
    @printf("  sidecar updated with validation results\n")
    opts["stage"] == "validate" && return nothing

    # --- promotion --------------------------------------------------------
    if !opts["promote"]
        @printf("promotion disabled (--no-promote)\n")
        return nothing
    end
    target = joinpath(ROOT, "data", "default", "default.pmx")
    failed = sort([String(k) for (k, v) in gates if !v])
    if isempty(failed)
        promote_validated(result, pcfg, reports; target = target, gates = gates)
        @printf("promoted to %s\n", target)
    elseif opts["force"]
        @printf("overriding failed gates: %s\n", join(failed, ", "))
        provenance["promoted_with_failed_gates"] = failed
        write_sidecar(result["sidecar_path"], provenance, model, ft)
        PaintMixPrecompute.promote(result["payload_path"], target; expected_id = result["model_id"])
        @printf("promoted to %s with recorded override\n", target)
    else
        @printf(
            "not promoted: failed gates %s; rerun with --force to promote and record the override\n",
            join(failed, ", ")
        )
    end
    return nothing
end

main(ARGS)
