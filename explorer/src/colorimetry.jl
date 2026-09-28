# Derived per-pigment curves and the linear/encoded-sRGB adapters every page
# shares.
#
# All colorimetry reuses `PaintMixPrecompute` (`mix_rgb`, `km_reflectance`,
# `saunderson`, `pigment_parameters`); nothing here re-derives the observer or
# the color matrix. The two conversion helpers exist so linear and encoded
# values are never mixed up: `encode` and `mix` take linear light, byte colors
# are encoded sRGB.

"""
    encoded(c) -> RGB8

Encode a linear-light sRGB triple as three sRGB bytes.
"""
encoded(c) = PaintMix.srgb8_from_linear(SVector{3, Float32}(c...))

"""
    linear(b) -> RGB{Float32}

Decode an encoded-sRGB byte triple to linear light.
"""
linear(b) = PaintMix.linear_from_srgb8(b)

"""
    xy_of_linear(quad, c) -> (x, y)

Chromaticity of a linear-light sRGB triple. Maps back to XYZ with the pinned
matrix and sum-normalizes, so it is the exact inverse of the transform
`mix_rgb` applies.
"""
function xy_of_linear(quad::Quadrature, c)
    xyz = inv(quad.xyz_to_rgb) * collect(Float64, c)
    s = sum(xyz)
    return (xyz[1] / s, xyz[2] / s)
end

# The four pure-pigment concentrations, in the model's fixed slot order.
const _PIGMENT_UNITS = ntuple(
    i -> SVector{4, Float64}(ntuple(j -> i == j ? 1.0 : 0.0, Val(4))), Val(4)
)

const _ZERO_RGB3 = SVector(0.0, 0.0, 0.0)
const _ZERO_XY = (0.0, 0.0)

# The spectral model evaluated at each pure pigment.
pure_rgb(m::SpectralModel) =
    ntuple(i -> SVector{3, Float64}(mix_rgb(m, _PIGMENT_UNITS[i])), Val(4))

"""
    DerivedCurves

Every curve and chromaticity the viewer shows, computed once per process from
the measured spectra, the fitted surrogate, and the runtime payload.

Fields are `4 x W` pigment-major for the measured and fitted spectra (the
layout of [`PigmentSpectra`](@ref PaintMixPrecompute.PigmentSpectra)), and
length-four tuples for the per-pigment color points. `K_cells` and `S_cells`
hold the source spreadsheet cell of each sample for the readout table.
"""
Base.@kwdef struct DerivedCurves
    wavelength::Vector{Float64} = Float64[]
    codes::NTuple{4, String} = ntuple(_ -> "", Val(4))
    K::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    S::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    KS::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    Rinf::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    Rprime::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    K_cells::Matrix{String} = Matrix{String}(undef, 4, 0)
    S_cells::Matrix{String} = Matrix{String}(undef, 4, 0)
    fitted_K::Union{Nothing, Matrix{Float64}} = nothing
    fitted_S::Union{Nothing, Matrix{Float64}} = nothing
    fitted_Rprime::Union{Nothing, Matrix{Float64}} = nothing
    measured_rgb::NTuple{4, SVector{3, Float64}} = ntuple(_ -> _ZERO_RGB3, Val(4))
    measured_xy::NTuple{4, Tuple{Float64, Float64}} = ntuple(_ -> _ZERO_XY, Val(4))
    fitted_rgb::Union{Nothing, NTuple{4, SVector{3, Float64}}} = nothing
    fitted_xy::Union{Nothing, NTuple{4, Tuple{Float64, Float64}}} = nothing
    runtime_rgb::NTuple{4, SVector{3, Float64}} = ntuple(_ -> _ZERO_RGB3, Val(4))
    runtime_xy::NTuple{4, Tuple{Float64, Float64}} = ntuple(_ -> _ZERO_XY, Val(4))
end

"""
    derive_curves(spectra, quad, spectral, fitted, model, db) -> DerivedCurves

Compute the measured `K`, `S`, `K/S`, `R∞`, `R′`, the source cell map, and the
three chromaticity layers for every pigment:

  * `measured_*` — the raw K/S spectral model at `c = eᵢ`,
  * `fitted_*` — the surrogate pigments, when a sidecar supplied `theta`,
  * `runtime_*` — the 8-bit forward table at `c = eᵢ`.
"""
function derive_curves(
        spectra::PigmentSpectra, quad::Quadrature, spectral::SpectralModel,
        fitted::Union{Nothing, SpectralModel}, model::PigmentModel,
        db::Union{Nothing, InputDatabase},
    )
    codes = spectra.codes
    W = length(spectra.wavelength)
    K = copy(spectra.K)
    S = copy(spectra.S)
    KS = K ./ S
    Rinf = km_reflectance.(KS)
    Rprime = saunderson.(Rinf, spectra.k1, spectra.k2)

    K_cells = fill("", 4, W)
    S_cells = fill("", 4, W)
    if db !== nothing
        slots = Dict(code => i for (i, code) in enumerate(codes))
        for r in db.spectra
            i = Base.get(slots, r.code, 0)
            i == 0 && continue
            j = findfirst(==(r.wavelength_nm), spectra.wavelength)
            j === nothing && continue
            (r.quantity == "K" ? K_cells : S_cells)[i, j] = r.cell_range
        end
    end

    measured_rgb = pure_rgb(spectral)
    measured_xy = ntuple(i -> xy_of_linear(quad, measured_rgb[i]), Val(4))
    runtime_rgb = ntuple(
        i -> SVector{3, Float64}(PaintMix.forward_rgb(model, _PIGMENT_UNITS[i])), Val(4)
    )
    runtime_xy = ntuple(i -> xy_of_linear(quad, runtime_rgb[i]), Val(4))

    fitted_K = fitted_S = fitted_Rprime = nothing
    fitted_rgb = fitted_xy = nothing
    if fitted !== nothing
        fitted_K, fitted_S = pigment_parameters(fitted)
        fitted_Rprime = saunderson.(
            km_reflectance.(fitted_K ./ fitted_S), fitted.k1, fitted.k2
        )
        fitted_rgb = pure_rgb(fitted)
        fitted_xy = ntuple(i -> xy_of_linear(quad, fitted_rgb[i]), Val(4))
    end

    return DerivedCurves(;
        wavelength = spectra.wavelength, codes = codes, K = K, S = S, KS = KS,
        Rinf = Rinf, Rprime = Rprime, K_cells = K_cells, S_cells = S_cells,
        fitted_K = fitted_K, fitted_S = fitted_S, fitted_Rprime = fitted_Rprime,
        measured_rgb = measured_rgb, measured_xy = measured_xy,
        fitted_rgb = fitted_rgb, fitted_xy = fitted_xy,
        runtime_rgb = runtime_rgb, runtime_xy = runtime_xy,
    )
end
