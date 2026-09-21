"""
    PaintMixPrecompute

Offline generation of the lookup tables that `PaintMix` evaluates.

The pipeline is, in order:

  1. import the selected spreadsheet ranges into DuckDB and load validated
     absorption and scattering spectra plus the Saunderson constants
     (`inputdb.jl`, `spectra.jl`, plan step 3),
  2. evaluate equations (1)-(7) in generic floating-point arithmetic
     (`kubelka_munk.jl`),
  3. fit surrogate pigments in gamut (`surrogates.jl`, equations (15)-(17)),
  4. solve the constrained RGB-to-concentration inverse (`unmix.jl`,
     equation (9)),
  5. generate the forward and inverse tables, quantize them jointly, pad the
     off-simplex cube vertices, and write the `.pmx` payload together with a
     provenance sidecar (`tables.jl`, `export.jl`), and
  6. validate the result against the acceptance gates and promote it
     (`validate.jl`).

Only the promoted payload is used at run time; nothing in this package is
reachable from the `PaintMix` mixing path.
"""
module PaintMixPrecompute

import PaintMix
using DataFrames: DataFrame, eachrow, nrow, rename!, select, stack
using Dates: Dates
using DuckDB
using ForwardDiff: ForwardDiff
using Optim: Optim
using Printf: @sprintf
using Random: AbstractRNG, Xoshiro
using SHA: sha256
using StaticArrays: SMatrix, SVector
using TOML: TOML

include("inputdb.jl")
include("spectra.jl")
include("kubelka_munk.jl")
include("surrogates.jl")
include("unmix.jl")
include("tables.jl")
include("validate.jl")
include("export.jl")

export ConfigError,
    InputError,
    InputDatabase,
    SourceFile,
    PigmentRecord,
    SpectrumRecord,
    SaundersonRecord,
    ObserverRecord,
    open_database,
    save_database,
    load_spreadsheet_inputs,
    check_checksums,
    validate_database,
    load_config,
    validate_config,
    config_hash,
    pigment_codes,
    wavelength_grid,
    DEFAULT_CONFIG_PATH,
    PigmentSpectra,
    Quadrature,
    SpectralModel,
    load_spectra,
    build_quadrature,
    spectral_model,
    with_parameters,
    pigment_parameters,
    mix_rgb,
    mix_rgb_params,
    mix_rgb_simplex,
    mix_rgb_jacobian!,
    km_reflectance,
    saunderson,
    linear_srgb_to_oklab,
    oklab_distance_squared,
    cube_signed_distance,
    cube_outside_penalty,
    softplus,
    inv_softplus,
    SurfaceQuadrature,
    surface_quadrature,
    surface_targets,
    Epush,
    Epull,
    rgb_surface_weights,
    quadrature_report,
    SurrogateFit,
    initial_theta,
    theta_parameters,
    theta_parameters_wl,
    alpha_schedule,
    fit_surrogates,
    surrogate_diagnostics,
    UnmixSettings,
    unmix_settings,
    UnmixResult,
    SolverScratch,
    unmix_reference,
    unmix_bulk!,
    lm_active!,
    project_to_simplex,
    quantize_simplex,
    FloatTables,
    ForwardFloatLUT,
    _scalar_type,
    forward_float_lut,
    forward_vertex!,
    generate_forward,
    generate_inverse,
    inverse_slab!,
    reference_slab!,
    coarse_seed,
    coarse_to_fine_slab!,
    generate_inverse_coarse_to_fine,
    CheckpointStore,
    slab_complete,
    write_slab,
    read_slab,
    slab_resume_function,
    slab_callback_function,
    completed_slabs,
    provenance_text,
    build_provenance,
    quantize_tables,
    build_model,
    export_model,
    write_sidecar,
    promote,
    promote_validated,
    spectral_encode,
    spectral_mix,
    quality_report,
    padding_report,
    roundtrip_report,
    quantization_report,
    continuity_report,
    behavior_report,
    acceptance_gates

"""
    DEFAULT_CONFIG_PATH

Absolute path of `precompute/config/default.toml`.
"""
const DEFAULT_CONFIG_PATH = normpath(joinpath(@__DIR__, "..", "config", "default.toml"))

struct ConfigError <: Exception
    msg::String
end

Base.showerror(io::IO, e::ConfigError) = print(io, "ConfigError: ", e.msg)

"""
    load_config(path = DEFAULT_CONFIG_PATH) -> Dict{String,Any}

Parse a precompute configuration and check it with [`validate_config`](@ref).
"""
function load_config(path::AbstractString = DEFAULT_CONFIG_PATH)
    isfile(path) || throw(ConfigError("no configuration at $path"))
    cfg = TOML.parsefile(path)
    validate_config(cfg)
    return cfg
end

"""
    validate_config(cfg) -> cfg

Check the invariants the pipeline depends on: a consistent wavelength grid, a
fixed four-pigment order with unique codes, a 3x3 color matrix, a storage
precision and grid size the runtime understands, and a positive surrogate
floor.

Throws [`ConfigError`](@ref) naming the first problem. The runtime's own
payload validation is deliberately separate: this checks a *plan*, while
`PaintMix.model_from_bytes` checks bytes.
"""
function validate_config(cfg::AbstractDict)
    get(cfg, "schema_version", nothing) == 1 ||
        throw(ConfigError("schema_version must be 1, got $(repr(get(cfg, "schema_version", nothing)))"))

    color = _section(cfg, "color")
    matrix = get(color, "xyz_to_rgb", nothing)
    matrix isa AbstractVector && length(matrix) == 3 ||
        throw(ConfigError("color.xyz_to_rgb must have 3 rows"))
    for row in matrix
        row isa AbstractVector && length(row) == 3 ||
            throw(ConfigError("every color.xyz_to_rgb row needs 3 entries"))
        all(v -> v isa Real && isfinite(v), row) ||
            throw(ConfigError("color.xyz_to_rgb entries must be finite numbers"))
    end
    color["storage_transfer"] == "linear" || throw(
        ConfigError(
            "color.storage_transfer must be \"linear\": encoded sRGB must never be " *
                "combined with linear-light residuals"
        )
    )

    spectra = _section(cfg, "spectra")
    lo = _positive(spectra, "wavelength_min_nm")
    hi = _positive(spectra, "wavelength_max_nm")
    step = _positive(spectra, "wavelength_step_nm")
    hi > lo || throw(ConfigError("spectra.wavelength_max_nm must exceed the minimum"))
    count = (hi - lo) / step + 1
    isinteger(count) || throw(
        ConfigError(
            "spectra grid ($lo:$step:$hi) does not contain an integer number of samples"
        )
    )
    count >= 4 || throw(ConfigError("spectra grid needs at least 4 samples"))

    saunderson = _section(cfg, "saunderson")
    for key in ("k1", "k2")
        v = get(saunderson, key, nothing)
        v isa Real && isfinite(v) || throw(ConfigError("saunderson.$key must be a finite number"))
    end
    0 <= saunderson["k1"] < 1 || throw(ConfigError("saunderson.k1 must be in [0, 1)"))
    0 <= saunderson["k2"] < 1 || throw(ConfigError("saunderson.k2 must be in [0, 1)"))

    pigments = get(cfg, "pigments", nothing)
    pigments isa AbstractVector && length(pigments) == 4 || throw(
        ConfigError(
            "exactly four pigments are required: the latent vector has four concentrations"
        )
    )
    codes = String[]
    for p in pigments
        p isa AbstractDict || throw(ConfigError("each [[pigments]] entry must be a table"))
        for key in ("code", "name", "ci", "column")
            haskey(p, key) && p[key] isa AbstractString ||
                throw(ConfigError("pigment entry is missing the string field \"$key\""))
        end
        push!(codes, p["code"])
    end
    length(unique(codes)) == 4 ||
        throw(ConfigError("pigment codes must be unique, got $(join(codes, ", "))"))

    grid = _section(cfg, "grid")
    for key in ("release_n", "dev_n")
        n = get(grid, key, nothing)
        n isa Integer && n >= 2 || throw(ConfigError("grid.$key must be an integer >= 2"))
    end
    grid["storage"] == "u8" || throw(ConfigError("grid.storage must be \"u8\" for the v1 format"))
    grid["byte_scale"] == "b/255" ||
        throw(ConfigError("grid.byte_scale must be \"b/255\" for the v1 format"))

    surrogate = _section(cfg, "surrogate")
    epsilon = get(surrogate, "epsilon", nothing)
    epsilon isa Real && epsilon > 0 ||
        throw(ConfigError("surrogate.epsilon must be a positive floor"))
    alpha = get(surrogate, "alpha_initial", nothing)
    alpha isa Real && alpha > 0 || throw(ConfigError("surrogate.alpha_initial must be positive"))
    divisions = get(surrogate, "surface_divisions", nothing)
    divisions isa Integer && divisions >= 1 || throw(
        ConfigError(
            "surrogate.surface_divisions must be an integer >= 1"
        )
    )

    inputs = _section(cfg, "inputs")
    for key in ("primary", "observer_file", "database", "checksums")
        get(inputs, key, nothing) isa AbstractString ||
            throw(ConfigError("inputs.$key must be a path string"))
    end

    unmix = _section(cfg, "unmix")
    unmix["objective"] == "rgb-least-squares" || throw(
        ConfigError(
            "unmix.objective must remain \"rgb-least-squares\"; changing it is a model variant"
        )
    )
    tol = get(unmix, "tolerance", nothing)
    tol isa Real && tol > 0 || throw(ConfigError("unmix.tolerance must be positive"))
    face = get(unmix, "face_threshold", nothing)
    face isa Real && face >= 0 || throw(ConfigError("unmix.face_threshold must be non-negative"))

    return cfg
end

function _section(cfg::AbstractDict, name::AbstractString)
    s = get(cfg, name, nothing)
    s isa AbstractDict || throw(ConfigError("missing [$(name)] section"))
    return s
end

function _positive(section::AbstractDict, key::AbstractString)
    v = get(section, key, nothing)
    v isa Real && isfinite(v) && v > 0 ||
        throw(ConfigError("$(key) must be a positive finite number"))
    return v
end

"""
    pigment_codes(cfg) -> NTuple{4,String}

The fixed pigment order recorded in the configuration.
"""
pigment_codes(cfg::AbstractDict) = ntuple(i -> String(cfg["pigments"][i]["code"]), Val(4))

"""
    wavelength_grid(cfg) -> StepRangeLen

The validated wavelength grid in nanometres, inclusive of both ends.
"""
function wavelength_grid(cfg::AbstractDict)
    validate_config(cfg)
    s = cfg["spectra"]
    return s["wavelength_min_nm"]:s["wavelength_step_nm"]:s["wavelength_max_nm"]
end

"""
    config_hash(path = DEFAULT_CONFIG_PATH) -> String

SHA-256 of the configuration file, recorded in the provenance sidecar so a
payload can be traced back to the exact settings that produced it.
"""
function config_hash(path::AbstractString = DEFAULT_CONFIG_PATH)
    return bytes2hex(sha256(read(path)))
end

end # module PaintMixPrecompute
