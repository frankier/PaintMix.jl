# Stage orchestration and verified artifact reuse. CLI parsing stays in scripts/.

"""Write a versioned physical fit and a digest, atomically."""
function save_surrogate(path, fit, cfg, job_hash, W)
    doc = Dict{String, Any}(
        "schema_version" => 2, "job_hash" => job_hash,
        "codes" => collect(pigment_codes(cfg)), "samples" => W,
        "theta" => fit.theta, "absorption" => vec(fit.K), "scattering" => vec(fit.S),
        "history" => fit.history, "diagnostics" => fit.diagnostics,
        "settings" => settings_dict(FitSettings(cfg))
    )
    doc["coefficients_sha256"] = coefficient_digest(fit.K, fit.S)
    if fit.trace !== nothing
        doc["trace"] = write_trace(dirname(path), merge(fit.trace, Dict("fit_job_id" => job_hash)))
    end
    mkpath(dirname(path))
    open(path * ".tmp", "w") do io
        TOML.print(io, doc; sorted = true)
    end
    mv(path * ".tmp", path; force = true)
    return path
end

_as_dicts(v) = [Dict{String, Any}(String(k) => x for (k, x) in pairs(d)) for d in v]

"""Read only a fit with matching settings, job identity, shape, and physical digest."""
function load_surrogate(path, cfg, job_hash, W)
    isfile(path) || return nothing
    doc = TOML.parsefile(path)
    get(doc, "schema_version", 0) == 2 || return nothing
    get(doc, "job_hash", "") == job_hash || return nothing
    get(doc, "settings", nothing) == settings_dict(FitSettings(cfg)) || return nothing
    Int(get(doc, "samples", 0)) == W || return nothing
    get(doc, "codes", []) == collect(pigment_codes(cfg)) || throw(InputError("fit pigment order mismatch"))
    K = reshape(Float64.(doc["absorption"]), 4, W)
    S = reshape(Float64.(doc["scattering"]), 4, W)
    all(x -> isfinite(x) && x > 0, K) && all(x -> isfinite(x) && x > 0, S) || throw(InputError("invalid fitted coefficients"))
    coefficient_digest(K, S) == doc["coefficients_sha256"] || throw(InputError("fit coefficient digest mismatch"))
    trace = if haskey(doc, "trace")
        candidate = read_trace(dirname(path), doc["trace"], "fit")
        candidate["fit_job_id"] == job_hash && candidate["surrogate_sha256"] == doc["coefficients_sha256"] ||
            throw(InputError("fit trace identity mismatch"))
        candidate
    end
    return SurrogateFit(
        K, S, Float64.(doc["theta"]), _as_dicts(doc["history"]),
        Dict{String, Any}(doc["diagnostics"]), SurfaceQuadrature(cfg), trace
    )
end

"""Load measured inputs, importing only when requested or the database is absent."""
function pipeline_inputs(cfg; root = normpath(joinpath(@__DIR__, "..", "..")), force_import = false)
    checked, failures = check_checksums(root, joinpath(root, cfg["inputs"]["checksums"]))
    mismatches = filter(f -> startswith(f, "checksum mismatch"), failures)
    isempty(mismatches) || throw(InputError(join(mismatches, "; ")))
    path = joinpath(root, cfg["inputs"]["database"])
    if force_import || !isfile(path)
        db = load_spreadsheet_inputs(cfg, DEFAULT_CONFIG_PATH; root)
        save_database(db, path)
    else
        db = open_database(path)
    end
    return db, PigmentSpectra(cfg, db), Quadrature(cfg, db)
end

function _validation_reports(model, spectral, cfg)
    reports = Dict{String, Any}(
        "quality" => quality_report(model, spectral, cfg),
        "padding" => padding_report(model, spectral, cfg),
        "roundtrip" => roundtrip_report(model, cfg),
        "quantization" => quantization_report(model, cfg),
        "continuity" => continuity_report(model, cfg),
        "behavior" => behavior_report(model, spectral, cfg),
    )
    gates = acceptance_gates(cfg, reports)
    for key in ("blue_yellow_is_green", "magenta_yellow_is_orange", "white_tint_monotone")
        gates[key] = reports["behavior"][key]
    end
    return reports, gates
end

function _reuse_payload(out, job_hash)
    manifest = joinpath(out, "artifact.toml")
    isfile(manifest) || throw(InputError("no verified table artifact; run --stage tables"))
    doc = TOML.parsefile(manifest)
    doc["job_hash"] == job_hash || throw(InputError("table artifact belongs to different resolved settings"))
    payload = joinpath(out, "default.pmx")
    payload_digest(payload) == doc["payload_sha256"] || throw(InputError("payload digest mismatch"))
    model = PaintMix.read_model(payload)
    PaintMix.model_id(model) == doc["model_id"] || throw(InputError("payload model id mismatch"))
    sidecar = joinpath(out, "$(doc["model_id"]).toml")
    provenance = TOML.parsefile(sidecar)
    provenance["table_job_id"] == job_hash || throw(InputError("sidecar job identity mismatch"))
    provenance["payload_sha256"] == doc["payload_sha256"] || throw(InputError("sidecar payload digest mismatch"))
    return Dict{String, Any}(
            "model" => model, "model_id" => PaintMix.model_id(model),
            "payload_path" => payload, "sidecar_path" => sidecar,
            "traces" => get(provenance, "traces", Dict{String, Any}())
        ), provenance
end

"""
    run_pipeline(cfg; out_dir, profile, stage, settings, promote, force, root, log)

Execute reusable offline stages. Existing fits and tables are reused only after
identity and integrity checks. Candidates are promoted only by explicit request.
"""
function run_pipeline(
        cfg::AbstractDict; out_dir, profile = "dev", stage = :all,
        settings = nothing, inputs = nothing, fit_solver = nothing, point_solver = nothing,
        solver_source = "", promote = false, force = false, force_import = false,
        trace_fit = false, trace_inverse = false,
        root = normpath(joinpath(@__DIR__, "..", "..")), cfg_path = DEFAULT_CONFIG_PATH, log = stdout
    )
    stage in (:inputs, :fit, :tables, :validate, :promote, :all) || throw(ConfigError("unknown stage"))
    trace_inverse && !(stage in (:all, :tables)) && throw(ConfigError("inverse trace requires the tables or all stage"))
    cfg = profile_config(cfg, profile)
    st = settings === nothing ? TableSettings(cfg; profile) : settings
    mkpath(out_dir)
    db, spectra, quad = inputs === nothing ? pipeline_inputs(cfg; root, force_import) : inputs
    stage == :inputs && return Dict("database" => db)
    fit_solver = fit_solver === nothing ? fit_solver_spec(cfg) : fit_solver
    point_solver = point_solver === nothing ? inverse_solver_spec(cfg) : point_solver
    fit_id = fit_job_id(cfg, db; solver = fit_solver, source_id = solver_source)
    fit_path = joinpath(out_dir, "fits", "$fit_id.toml")
    W = length(spectra.wavelength)
    fit = load_surrogate(fit_path, cfg, fit_id, W)
    if fit === nothing || (trace_fit && fit.trace === nothing)
        stage in (:all, :fit) || throw(InputError("no matching fit with requested capture; run --stage fit"))
        fit = fit_surrogates(
            cfg, spectra, quad; log, algorithm = fit_solver.algorithm,
            adtype = fit_solver.adtype, solve_kwargs = fit_solver.solve_kwargs, trace_fit
        )
        save_surrogate(fit_path, fit, cfg, fit_id, W)
    end
    if stage == :fit
        metadata = TOML.parsefile(fit_path)
        return Dict(
            "fit" => fit, "fit_job_id" => fit_id, "fit_path" => fit_path,
            "fit_trace" => get(metadata, "trace", nothing)
        )
    end
    spectral = SpectralModel(SpectralModel(spectra, quad); K = fit.K, S = fit.S)
    job_hash = experiment_id(Dict("table" => table_job_id(cfg, fit, st; point_solver, source_id = solver_source), "fit" => fit_id))
    ft = nothing
    if stage in (:all, :tables)
        started = time()
        forward = generate_forward(spectral, st.n; threads = st.threads)
        evaluator = st.evaluator == :spectral ? spectral : forward_float_lut(st.n, forward)
        cp = CheckpointStore(joinpath(out_dir, "checkpoints", job_hash), job_hash, st.n, "inverse")
        summary = Dict{Int, Any}()
        summary_lock = ReentrantLock()
        on_summary = (k, report) -> lock(summary_lock) do
            summary[k] = report
        end
        inverse = generate_tables_inverse(
            evaluator, st; settings = UnmixSettings(cfg), point_solver,
            resume = slab_resume_function(cp), on_slab = slab_callback_function(cp), on_summary
        )
        ft = FloatTables(st.n, forward, inverse)
        provenance = build_provenance(cfg, cfg_path, db, fit, Dict{String, Any}(); grid_n = st.n)
        merge!(
            provenance, Dict(
                "fit_job_id" => fit_id, "table_job_id" => job_hash,
                "fit_solver" => solver_description(fit_solver), "point_solver" => solver_description(point_solver),
                "numerical_environment" => numerical_environment(fit_solver, point_solver), "solver_source" => solver_source,
                "resolved_config" => cfg, "resolved_config_sha256" => experiment_id(cfg), "table_settings" => settings_dict(st), "profile" => profile,
                "generation_seconds" => time() - started,
                "slab_summaries" => Dict(string(k) => v for (k, v) in summary)
            )
        )
        capture_inverse = trace_inverse ?
            ((model, doc) -> build_inverse_trace(model, spectral, evaluator, ft, doc, UnmixSettings(cfg))) : nothing
        result = export_model(ft, cfg, provenance; out_dir = out_dir, fit_trace = fit.trace, capture_inverse)
        provenance = TOML.parsefile(result["sidecar_path"])
        open(joinpath(out_dir, "artifact.toml"), "w") do io
            TOML.print(
                io, Dict(
                    "schema_version" => 2, "job_hash" => job_hash,
                    "payload_sha256" => provenance["payload_sha256"], "model_id" => result["model_id"]
                ); sorted = true
            )
        end
    else
        result, provenance = _reuse_payload(out_dir, job_hash)
    end
    stage == :tables && return result
    reports, gates = _validation_reports(result["model"], spectral, cfg)
    if ft !== nothing
        reports["decomposition"] = error_decomposition(result["model"], spectral, ft, cfg; measured = SpectralModel(spectra, quad))
    elseif haskey(get(provenance, "validation", Dict()), "decomposition")
        reports["decomposition"] = provenance["validation"]["decomposition"]
    end
    provenance["validation"] = reports
    provenance["acceptance_gates"] = gates
    write_sidecar(result["sidecar_path"], provenance, result["model"], ft)
    result["validation"] = reports
    result["acceptance_gates"] = gates
    if promote && stage != :validate
        failed = sort!([k for (k, ok) in gates if !ok])
        target = joinpath(root, "data", "default", "default.pmx")
        if isempty(failed)
            promote_validated(result, cfg, reports; target, gates)
        elseif force
            provenance["promoted_with_failed_gates"] = failed
            write_sidecar(result["sidecar_path"], provenance, result["model"], ft)
            PaintMixPrecompute.promote(result["payload_path"], target; expected_id = result["model_id"])
        else
            log === nothing || println(log, "not promoted: failed gates ", join(failed, ", "))
        end
    end
    return result
end
