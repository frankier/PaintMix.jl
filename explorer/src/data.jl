# ExplorerData assembly.
#
# One value is assembled at startup and cached for the process. Every source
# degrades rather than throwing: a missing database, payload, or sidecar adds a
# note and the pages show a "not available" state. `ExplorerData` also has an
# injectable keyword constructor, so tests can swap in a synthetic model
# without touching disk.

# The repository root, resolved from this file rather than the working
# directory, exactly as `precompute/scripts/generate.jl` does.
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

const TEMPLATES = normpath(joinpath(@__DIR__, "..", "templates"))
const STATIC = normpath(joinpath(@__DIR__, "..", "static"))
const OUTPUT = joinpath(ROOT, "precompute", "output")

"""
    ExplorerData

Everything the viewer reads, loaded once.

`cfg`, `db`, `spectra`, `quad`, `spectral`, `fitted`, and `sidecar` are
optional so the viewer can run against synthetic tables in tests and against a
partial checkout. `notes` records every degradation for the provenance banner.

Construct directly for tests:

    ExplorerData(; model = synthetic_model(), source = "synthetic")
"""
struct ExplorerData
    cfg::Dict{String, Any}
    db::Union{Nothing, InputDatabase}
    spectra::Union{Nothing, PigmentSpectra{Float64}}
    quad::Union{Nothing, Quadrature{Float64}}
    spectral::Union{Nothing, SpectralModel{Float64}}
    fitted::Union{Nothing, SpectralModel{Float64}}
    model::PigmentModel
    source::String
    sidecar::Union{Nothing, Dict{String, Any}}
    derived::DerivedCurves
    notes::Vector{String}
end

function ExplorerData(;
        cfg::Dict{String, Any} = Dict{String, Any}(),
        db::Union{Nothing, InputDatabase} = nothing,
        spectra::Union{Nothing, PigmentSpectra{Float64}} = nothing,
        quad::Union{Nothing, Quadrature{Float64}} = nothing,
        spectral::Union{Nothing, SpectralModel{Float64}} = nothing,
        fitted::Union{Nothing, SpectralModel{Float64}} = nothing,
        model::PigmentModel,
        source::AbstractString = "synthetic",
        sidecar::Union{Nothing, Dict{String, Any}} = nothing,
        derived::DerivedCurves = DerivedCurves(),
        notes::Vector{String} = String[],
    )
    return ExplorerData(
        cfg, db, spectra, quad, spectral, fitted, model, String(source),
        sidecar, derived, notes
    )
end

"""
    synthetic_model(; n = 2) -> PigmentModel

A tiny valid payload used only when no `.pmx` exists. The forward table samples
the identity on the cube; the inverse table projects to the simplex and
quantizes with largest-remainder rounding. It carries the all-zero id, so it
cannot be mistaken for a generated model.
"""
function synthetic_model(; n::Integer = 2)
    m = Int(n)
    payload = 3 * m^3
    fwd = Vector{UInt8}(undef, payload)
    inv = Vector{UInt8}(undef, payload)
    for k in 0:(m - 1), j in 0:(m - 1), i in 0:(m - 1)
        base = 3 * (i + m * (j + m * k))
        rgb = m == 1 ? (0.0, 0.0, 0.0) : (i / (m - 1), j / (m - 1), k / (m - 1))
        for ch in 1:3
            fwd[base + ch] = UInt8(clamp(round(Int, 255 * rgb[ch]), 0, 255))
        end
        q = _simplex_bytes(rgb)
        inv[base + 1] = UInt8(q[1])
        inv[base + 2] = UInt8(q[2])
        inv[base + 3] = UInt8(q[3])
    end
    model = PigmentModel(
        ntuple(_ -> 0x00, Val(16)), ByteLUT(m, inv), ByteLUT(m, fwd), FORMAT_VERSION,
        FLAG_FORWARD_SIMPLEX_PROJECTED | FLAG_INVERSE_LARGEST_REMAINDER,
    )
    return validate_model(model)
end

# Project to the concentration simplex and quantize the first three slots so
# that `b1 + b2 + b3 <= 255` holds.
function _simplex_bytes(rgb::NTuple{3, Float64})
    s = rgb[1] + rgb[2] + rgb[3]
    q = s > 1 ? rgb ./ s : rgb
    b = [round(Int, 255 * q[i]) for i in 1:3]
    while sum(b) > 255
        b[argmax(b)] -= 1
    end
    return (b[1], b[2], b[3])
end

# Payload candidates, first found wins. The release payload is a 96 MiB
# untracked artifact, so its absence is normal.
function payload_candidates()
    return (
        (joinpath(ROOT, "data", "default", "default.pmx"), "release (data/default)"),
        (joinpath(OUTPUT, "release", "default.pmx"), "release (output/release)"),
        (joinpath(OUTPUT, "dev", "default.pmx"), "dev (output/dev)"),
        (joinpath(ROOT, "build", "data", "payload.pmx"), "build dummy payload"),
    )
end

function load_payload(override::Union{Nothing, AbstractString}, notes::Vector{String})
    candidates = override === nothing ? payload_candidates() :
        ((abspath(override), "explicit ($override)"),)
    for (path, label) in candidates
        isfile(path) || continue
        try
            return read_model(path), label, path
        catch err
            push!(notes, "payload at $path unreadable: $(sprint(showerror, err))")
        end
    end
    push!(notes, "no payload found; using synthetic in-memory tables")
    return synthetic_model(), "synthetic", nothing
end

# The sidecar sits next to the payload or in the release artifact directory.
function find_sidecar(model::PigmentModel, payload_path::Union{Nothing, AbstractString})
    id = PaintMix.model_id(model)
    if payload_path !== nothing
        path = joinpath(dirname(payload_path), "$id.toml")
        isfile(path) && return TOML.parsefile(path)
    end
    path = joinpath(OUTPUT, "release", "$id.toml")
    isfile(path) && return TOML.parsefile(path)
    return nothing
end

function build_fitted(sidecar::Dict{String, Any}, spectral::SpectralModel)
    surrogate = Base.get(sidecar, "surrogate", nothing)
    theta = surrogate isa AbstractDict ? Base.get(surrogate, "theta", nothing) : nothing
    theta === nothing && return nothing
    epsilon = Base.get(surrogate, "epsilon", 1.0e-6)
    W = length(spectral.quad.wavelength)
    K, S = theta_parameters(Float64.(theta), W, epsilon)
    return SpectralModel(spectral; K = K, S = S)
end

"""
    load_explorer_data(; with_sidecar = true, payload_path = nothing) -> ExplorerData

Assemble the viewer's data from the repository, degrading at every step. Pass
`with_sidecar = false` to exercise the missing-sidecar path; pass
`payload_path` to pin a specific `.pmx`.
"""
function load_explorer_data(;
        with_sidecar::Bool = true, payload_path::Union{Nothing, AbstractString} = nothing
    )
    notes = String[]

    cfg = try
        load_config()
    catch err
        push!(notes, "configuration unavailable: $(sprint(showerror, err))")
        Dict{String, Any}()
    end

    db = nothing
    if !isempty(cfg)
        db_path = joinpath(ROOT, cfg["inputs"]["database"])
        if isfile(db_path)
            db = try
                open_database(db_path)
            catch err
                push!(notes, "input database unreadable: $(sprint(showerror, err))")
                nothing
            end
        else
            push!(notes, "input database absent at $db_path")
        end
    end

    spectra = quad = spectral = nothing
    if db !== nothing
        try
            spectra = PigmentSpectra(cfg, db)
            quad = Quadrature(cfg, db)
            spectral = SpectralModel(spectra, quad)
        catch err
            push!(notes, "spectral reference unavailable: $(sprint(showerror, err))")
            spectra = quad = spectral = nothing
        end
    end

    model, source, used_path = load_payload(payload_path, notes)

    sidecar = nothing
    if with_sidecar
        sidecar = find_sidecar(model, used_path)
        if sidecar === nothing
            push!(notes, "provenance sidecar absent; fitted curves unavailable")
        end
    else
        push!(notes, "sidecar disabled (--no-sidecar)")
    end

    fitted = nothing
    if sidecar !== nothing && spectral !== nothing
        try
            fitted = build_fitted(sidecar, spectral)
        catch err
            push!(notes, "surrogate fit unavailable: $(sprint(showerror, err))")
            fitted = nothing
        end
    end

    derived = if spectra !== nothing && quad !== nothing && spectral !== nothing
        derive_curves(spectra, quad, spectral, fitted, model, db)
    else
        DerivedCurves()
    end

    return ExplorerData(;
        cfg = cfg, db = db, spectra = spectra, quad = quad, spectral = spectral,
        fitted = fitted, model = model, source = source, sidecar = sidecar,
        derived = derived, notes = notes,
    )
end
