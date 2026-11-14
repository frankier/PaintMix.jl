#!/usr/bin/env julia
# Reuse verified floating inverse slabs; never refit or promote a model.
using PaintMixPrecompute
using PaintMix
using SHA
using TOML

function main(args = ARGS; inputs = nothing)
    length(args) == 2 || error("usage: precision.jl CANDIDATE_DIR OUT")
    candidate, out = args
    artifact = TOML.parsefile(joinpath(candidate, "artifact.toml"))
    sidecar = TOML.parsefile(joinpath(candidate, "$(artifact["model_id"]).toml"))
    bytes2hex(sha256(read(joinpath(candidate, "default.pmx")))) == artifact["payload_sha256"] || error("candidate digest mismatch")
    sidecar["table_job_id"] == artifact["job_hash"] || error("sidecar job mismatch")
    sidecar["payload_sha256"] == artifact["payload_sha256"] || error("sidecar digest mismatch")
    cfg = sidecar["resolved_config"]
    _, spectra, quad = inputs === nothing ? pipeline_inputs(cfg) : inputs
    W = length(spectra.wavelength)
    physical = sidecar["parameters"]
    K = permutedims(hcat(physical["absorption"]...))
    S = permutedims(hcat(physical["scattering"]...))
    size(K) == size(S) == (4, W) || error("physical coefficient shape mismatch")
    spectral = SpectralModel(SpectralModel(spectra, quad); K, S)
    n = Int(sidecar["table_settings"]["n"])
    cp = CheckpointStore(joinpath(candidate, "checkpoints", artifact["job_hash"]), artifact["job_hash"], n, "inverse")
    inverse = Vector{Float64}(undef, 3n^3)
    for k in 0:(n - 1)
        slab = read_slab(cp, k)
        slab === nothing && error("missing verified floating inverse slab $k")
        copyto!(inverse, 3k * n^2 + 1, slab, 1, length(slab))
    end
    forward = generate_forward(spectral, n)
    cube = ByteLUT(n, round.(UInt8, 255 .* clamp.(forward, 0.0, 1.0)))
    simplex = generate_simplex_forward(spectral, n)
    samples = validation_corpus()
    runs = Dict{String, Any}[]
    mkpath(out)
    for U in (UInt8, UInt16)
        q = quantize_inverse(inverse, n, U)
        inv = U == UInt8 ? ByteLUT(n, q) : WordLUT(n, q)
        for (mode, fwd) in (("cube", cube), ("simplex", simplex))
            name = "$(8sizeof(U))-$mode"
            model = experimental_model(inv, fwd)
            payload = joinpath(out, "$name.pmxe")
            write_model(payload, model)
            loaded = read_model(payload)
            @assert model_id(loaded) == model_id(model)
            report = precision_report(loaded, spectral, cfg; samples)
            push!(runs, Dict("variant" => name, "model_id" => model_id(model), "payload_bytes" => filesize(payload), "payload_sha256" => bytes2hex(sha256(read(payload))), "report" => report))
            println(name, " Oklab p99=", report["quality"]["oklab_error"]["p99"])
            flush(stdout)
        end
    end
    doc = Dict("schema_version" => 1, "grid_n" => n, "source_model_id" => artifact["model_id"], "source_payload_sha256" => artifact["payload_sha256"], "source_fit_job_id" => sidecar["fit_job_id"], "corpus_sha256" => samples.digest, "runs" => runs, "implementation_sha256" => PaintMixPrecompute.implementation_hash(["precision.jl", "tables.jl", "export.jl"]), "script_sha256" => bytes2hex(sha256(read(@__FILE__))), "environment" => numerical_environment())
    open(joinpath(out, "precision.toml"), "w") do io
        TOML.print(io, doc; sorted = true)
    end
    return doc
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
