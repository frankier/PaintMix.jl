"""
    PaintMix

Practical pigment mixing in RGB, following Sochorová and Jamriška,
*Practical Pigment Mixing for Digital Painting* (2021).

The package evaluates a pair of precomputed lookup tables:

  * an inverse table `U`, which maps a linear-light sRGB color to four pigment
    concentrations, and
  * a forward table `M`, which maps those concentrations back to linear-light
    sRGB.

Encoding stores four concentrations plus a signed RGB residual, so that
linear operations on encoded colors behave like paint mixtures while still
reproducing colors that no mixture of the four pigments can produce:

    encode(x)       = Latent(U(x), x - M(U(x)))
    decode(z)       = M(z.c) + z.r
    mix(a, b, t)    = decode((1 - t) * encode(a) + t * encode(b))

Colors are linear-light sRGB with sRGB primaries and the D65 white point,
matching equation (7) of the paper. Use [`linear_from_srgb8`](@ref) and
[`srgb8_from_linear`](@ref) for byte-oriented image data.

The module is a pure evaluator. It owns no optimizer, no spectral model, and
no file-format dependencies beyond the packaged `.pmx` payload. Table
generation lives in the separate `PaintMixPrecompute` package.
"""
module PaintMix

export ByteLUT,
    PigmentModel,
    Latent,
    ModelHeader,
    InvalidPayload,
    default_model,
    default_payload_path,
    model_id,
    grid_n,
    encode,
    decode,
    mix,
    mix!,
    bulk_mix!,
    weighted_mix,
    weighted_mix!,
    linear_from_srgb,
    srgb_from_linear,
    linear_from_srgb8,
    srgb8_from_linear,
    read_model,
    write_model,
    model_from_bytes,
    model_to_bytes,
    validate_model,
    concentrations,
    residual,
    error_message,
    FORMAT_VERSION,
    HEADER_BYTES,
    FLAG_FORWARD_SIMPLEX_PROJECTED,
    FLAG_INVERSE_LARGEST_REMAINDER,
    PM_OK,
    PM_ERR_NULL,
    PM_ERR_LENGTH,
    PM_ERR_NONFINITE,
    PM_ERR_WEIGHT,
    PM_ERR_TOTAL

include("types.jl")
include("srgb.jl")
include("lookup.jl")
include("mixing.jl")
include("tables.jl")

end # module PaintMix
