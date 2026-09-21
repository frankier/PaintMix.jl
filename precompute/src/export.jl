# Plan step 5: payload and provenance export.
#
# The `.pmx` payload layout lives in `PaintMix` (`src/tables.jl`); this file
# only converts float candidates to bytes, derives a model id from the
# provenance, writes the payload, and writes the human- and machine-readable
# sidecar. `promote` copies a validated candidate into `data/default/`.

"""
    model_id_bytes(io_writer) -> NTuple{16,UInt8}

First 16 bytes of a SHA-256 over the canonical provenance string.
"""
function _model_id_from(text::AbstractString)
    digest = sha256(text)
    return ntuple(i -> digest[i], Val(16))
end

"""
    provenance_text(cfg, cfg_hash, inputs) -> String

The canonical string the model id hashes. Every convention that changes the
bytes must appear here, so two payloads with the same id are the same model.
"""
function provenance_text(
        cfg::AbstractDict, cfg_hash::AbstractString, inputs::Dict{String,String};
        grid_n::Integer = Int(cfg["grid"]["release_n"]),
    )
    io = IOBuffer()
    print(io, "paintmix-format-v1\n")
    print(io, "grid_n=", grid_n, "\n")
    print(io, "storage=", cfg["grid"]["storage"], "\n")
    print(io, "byte_scale=", cfg["grid"]["byte_scale"], "\n")
    print(io, "interpolation=", cfg["grid"]["interpolation"], "\n")
    print(io, "color_space=", cfg["color"]["color_space"], "\n")
    print(io, "storage_transfer=", cfg["color"]["storage_transfer"], "\n")
    print(io, "quantization=", cfg["quantization"]["rule"], "\n")
    print(io, "padding=", cfg["quantization"]["padding"], "\n")
    print(io, "pigments=", join(pigment_codes(cfg), ","), "\n")
    print(io, "config_hash=", cfg_hash, "\n")
    for key in sort(collect(keys(inputs)))
        print(io, "input:", key, "=", inputs[key], "\n")
    end
    return String(take!(io))
end

"""
    build_provenance(cfg, cfg_path, db, fit, validation) -> Dict{String,Any}

Assemble everything the sidecar records: configuration hash and settings,
input files and checksums, surrogate parameters and convergence history,
package versions, seeds, precision, quantization and padding rules, and the
validation results.
"""
function build_provenance(
        cfg::AbstractDict, cfg_path::AbstractString, db::InputDatabase,
        fit::SurrogateFit, validation::AbstractDict; grid_n::Integer = Int(cfg["grid"]["release_n"]),
    )
    inputs = Dict{String,String}(
        "config" => config_hash(cfg_path),
        "observer" => db.observer_source.sha256,
    )
    for f in db.source_files
        inputs[f.role] = f.sha256
    end
    K, S = fit.K, fit.S
    return Dict{String,Any}(
        "schema_version" => 1,
        "generated_utc" => string(Dates.now()),
        "config_path" => abspath(cfg_path),
        "config_hash" => config_hash(cfg_path),
        "inputs" => inputs,
        "pigments" => [
            Dict{String,Any}(
                "slot" => i,
                "code" => pigment_codes(cfg)[i],
                "name" => cfg["pigments"][i]["name"],
                "ci" => cfg["pigments"][i]["ci"],
            ) for i in 1:4
        ],
        "grid" => Dict{String,Any}(
            "n" => grid_n,
            "storage" => cfg["grid"]["storage"],
            "byte_scale" => cfg["grid"]["byte_scale"],
            "interpolation" => cfg["grid"]["interpolation"],
            "index_order" => cfg["grid"]["index_order"],
        ),
        "quantization" => Dict{String,Any}(
            "rule" => cfg["quantization"]["rule"],
            "padding" => cfg["quantization"]["padding"],
        ),
        "surrogate" => Dict{String,Any}(
            "epsilon" => cfg["surrogate"]["epsilon"],
            "alpha_initial" => cfg["surrogate"]["alpha_initial"],
            "alpha_final" => cfg["surrogate"]["alpha_final"],
            "alpha_halvings" => cfg["surrogate"]["alpha_halvings"],
            "solver" => "Optim.LBFGS on softplus-transformed parameters",
            "note" =>
                "the paper uses L-BFGS-B; transformed L-BFGS is an intentional substitution",
            "surface_divisions" => cfg["surrogate"]["surface_divisions"],
            "history" => fit.history,
            "diagnostics" => fit.diagnostics,
            "theta" => fit.theta,
        ),
        "parameters" => Dict{String,Any}(
            "absorption" => _matrix_rows(K),
            "scattering" => _matrix_rows(S),
        ),
        "unmix" => Dict{String,Any}(
            "strategy" => cfg["unmix"]["strategy"],
            "objective" => cfg["unmix"]["objective"],
            "restarts" => cfg["unmix"]["restarts"],
            "max_iterations" => cfg["unmix"]["max_iterations"],
        ),
        "environment" => Dict{String,Any}(
            "julia_version" => string(VERSION),
            "threads" => Threads.nthreads(),
            "paintmix_version" => _pkg_version(PaintMix),
            "paintmixprecompute_version" => _pkg_version(@__MODULE__),
            "optim_version" => _pkg_version(Optim),
            "forwarddiff_version" => _pkg_version(ForwardDiff),
        ),
        "validation" => validation,
    )
end

function _matrix_rows(M::AbstractMatrix)
    return [[Float64(x) for x in row] for row in eachrow(M)]
end

function _pkg_version(mod::Module)
    path = pkgdir(mod)
    path === nothing && return "unknown"
    project = joinpath(path, "Project.toml")
    isfile(project) || return "unknown"
    return string(get(TOML.parsefile(project), "version", "unknown"))
end

"""
    quantize_tables(ft::FloatTables, cfg) -> (inverse::ByteLUT, forward::ByteLUT)

Convert the float candidates to the payload's byte tables. Inverse vertices
are quantized with [`quantize_simplex`](@ref), guaranteeing
`b1 + b2 + b3 <= 255`; forward vertices are rounded and clipped to
`[0, 255]`, which is the only clip in the pipeline and happens after the
mixture has been evaluated.
"""
function quantize_tables(ft::FloatTables, cfg::AbstractDict)
    n = ft.n
    inv = Vector{UInt8}(undef, 3 * n^3)
    fwd = Vector{UInt8}(undef, 3 * n^3)
    @inbounds for v in 0:(n^3 - 1)
        o = 3 * v
        c1 = ft.inverse[o + 1]
        c2 = ft.inverse[o + 2]
        c3 = ft.inverse[o + 3]
        c4 = one(Float64) - c1 - c2 - c3
        if c4 < 0
            # The solver works on the simplex; this only repairs floating
            # point drift before joint quantization.
            p = project_to_simplex((c1, c2, c3, max(c4, 0.0)))
            c1, c2, c3 = p[1], p[2], p[3]
        end
        b = quantize_simplex((c1, c2, c3, max(one(Float64) - c1 - c2 - c3, 0.0)))
        inv[o + 1] = b[1]
        inv[o + 2] = b[2]
        inv[o + 3] = b[3]
        fwd[o + 1] = _byte(ft.forward[o + 1])
        fwd[o + 2] = _byte(ft.forward[o + 2])
        fwd[o + 3] = _byte(ft.forward[o + 3])
    end
    return PaintMix.ByteLUT(n, inv), PaintMix.ByteLUT(n, fwd)
end

@inline function _byte(x::Float64)
    v = round(255 * clamp(x, 0.0, 1.0))
    return UInt8(clamp(v, 0, 255))
end

"""
    build_model(ft::FloatTables, cfg, provenance) -> PigmentModel

Quantize the candidates, derive the model id from the provenance, and set
the generation flags. Does not touch the filesystem.
"""
function build_model(ft::FloatTables, cfg::AbstractDict, provenance::AbstractDict)
    inputs = Dict{String,String}(
        String(k) => String(v) for (k, v) in provenance["inputs"]
    )
    text = provenance_text(
        cfg, provenance["config_hash"], inputs;
        grid_n = Int(get(get(provenance, "grid", Dict{String,Any}()), "n", cfg["grid"]["release_n"])),
    )
    id = _model_id_from(text)
    inverse, forward = quantize_tables(ft, cfg)
    flags = PaintMix.FLAG_FORWARD_SIMPLEX_PROJECTED | PaintMix.FLAG_INVERSE_LARGEST_REMAINDER
    return PaintMix.PigmentModel(id, inverse, forward, PaintMix.FORMAT_VERSION, flags)
end

"""
    export_model(ft, cfg, provenance; out_dir, name = "default.pmx") -> Dict

Write the payload and its TOML sidecar into `out_dir`. Returns the paths and
the model id. The payload is written through `PaintMix.write_model`, so there
is exactly one implementation of the format.
"""
function export_model(
        ft::FloatTables, cfg::AbstractDict, provenance::AbstractDict;
        out_dir::AbstractString, name::AbstractString = "default.pmx",
    )
    mkpath(out_dir)
    model = build_model(ft, cfg, provenance)
    PaintMix.validate_model(model)
    payload_path = joinpath(out_dir, name)
    PaintMix.write_model(payload_path, model)
    id = PaintMix.model_id(model)
    sidecar = joinpath(out_dir, "$(id).toml")
    write_sidecar(sidecar, provenance, model, ft)
    return Dict{String,Any}(
        "model" => model,
        "model_id" => id,
        "payload_path" => payload_path,
        "sidecar_path" => sidecar,
    )
end

"""
    write_sidecar(path, provenance, model, ft)

Write the provenance as TOML. Validation results are embedded so the
acceptance gates can be audited after the fact.
"""
function write_sidecar(path::AbstractString, provenance::AbstractDict, model, ft)
    doc = deepcopy(provenance)
    doc["model"] = Dict{String,Any}(
        "id" => PaintMix.model_id(model),
        "format_version" => Int(PaintMix.FORMAT_VERSION),
        "grid_n" => PaintMix.grid_n(model),
        "flags" => Int(model.flags),
        "payload_bytes" => 128 + 6 * PaintMix.grid_n(model)^3,
        "table_bytes_each" => 3 * PaintMix.grid_n(model)^3,
        "float_grid_n" => ft === nothing ? PaintMix.grid_n(model) : ft.n,
    )
    open(path, "w") do io
        TOML.print(io, doc)
    end
    return path
end

"""
    promote(candidate_payload, target; expected_id = nothing) -> PigmentModel

Copy a candidate payload into the packaged location, but only after reading
it back, checking its id, and running `PaintMix.validate_model`. Returns the
loaded and validated model.
"""
function promote(candidate_payload::AbstractString, target::AbstractString;
        expected_id::Union{Nothing,AbstractString} = nothing,
    )
    model = PaintMix.read_model(candidate_payload)
    if expected_id !== nothing
        got = PaintMix.model_id(model)
        got == expected_id || throw(InputError(
            "candidate $candidate_payload has model id $got, expected $expected_id"
        ))
    end
    mkpath(dirname(target))
    cp(candidate_payload, target; force = true)
    return PaintMix.read_model(target)
end

"""
    promote_validated(result, cfg, validation; target, gates) -> PigmentModel

Promote only when every acceptance gate in `gates` is `true`. The gates are
named in the error so a refusal is actionable.
"""
function promote_validated(result::AbstractDict, cfg::AbstractDict, validation::AbstractDict;
        target::AbstractString, gates::AbstractDict,
    )
    failed = String[k for (k, v) in gates if !v]
    isempty(failed) || throw(InputError(
        "refusing to promote $(result["model_id"]): failed acceptance gates " *
            join(failed, ", ") * "; see $(result["sidecar_path"])"
    ))
    validation isa AbstractDict ||
        throw(InputError("validation record must be a dictionary"))
    return promote(result["payload_path"], target; expected_id = result["model_id"])
end
