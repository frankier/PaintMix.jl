# Plan step 3: the validated spectral reference.
#
#   * `PigmentSpectra` — the four absorption and scattering curves on the
#     configured grid, plus the Saunderson constants.
#   * `load_spectra` — turn an [`InputDatabase`](@ref) into checked matrices.
#   * `Quadrature` — trapezoidal weights and the D65/CIE 1931 tables, plus
#     the pinned XYZ-to-RGB matrix and the `1 / Y_D65` factor.
#   * `SpectralModel` — spectra plus quadrature, the object equations (1)-(7)
#     are evaluated through.
#
# The source workbook's `kins` field is contradictory (1.0 in `k and s data`,
# 0.0 in `Details`) and is not used: the paper's convention has no added
# specular term. See `inputs/README.md`.

"""
    PigmentSpectra{T}

The four pigments' absorption and scattering coefficients on the configured
wavelength grid.

  * `codes`, `names`: the fixed pigment order from the configuration.
  * `K`, `S`: `4 x W` matrices, pigment-major, in the units the source
    measured (the K-M equations are invariant to a common scale).
  * `wavelength`: the `W` nanometre samples shared by every curve.
  * `k1`, `k2`: the Saunderson constants of equation (6).
  * `provenance`: free-form record of the source file, sheet, and hashes.
"""
struct PigmentSpectra{T <: AbstractFloat}
    codes::NTuple{4, String}
    names::NTuple{4, String}
    wavelength::Vector{T}
    K::Matrix{T}
    S::Matrix{T}
    k1::T
    k2::T
    provenance::Dict{String, Any}
end

"""
    Quadrature{T}

Everything equation (3)-(7) integration needs, cached once per grid:
trapezoidal weights, the CIE 1931 2-degree observer functions, D65, the
pinned XYZ-to-linear-sRGB matrix, and `1 / Y_D65`.
"""
struct Quadrature{T <: AbstractFloat}
    wavelength::Vector{T}
    weights::Vector{T}
    x_bar::Vector{T}
    y_bar::Vector{T}
    z_bar::Vector{T}
    d65::Vector{T}
    xyz_to_rgb::NTuple{9, T}
    norm::T
end

"""
    SpectralModel{T}

`PigmentSpectra` and `Quadrature` fused into the layout the hot loops want:
`K` and `S` are `W x 4` (wavelength-major) so a wavelength's four pigment
coefficients are adjacent.

A surrogate fit produces one of these from fitted `K`/`S` while reusing the
same quadrature, which is why the two are separate fields rather than one
merged struct.
"""
struct SpectralModel{T <: AbstractFloat}
    K::Matrix{T}
    S::Matrix{T}
    quad::Quadrature{T}
    k1::T
    k2::T
end

"""
    load_spectra(cfg, db) -> PigmentSpectra

Build the four spectral curves from the validated input database.

Checks that every pigment/quantity pair has exactly the configuration's
wavelength grid, in order, and finite positive values, and that the source
Saunderson constants agree with the configuration. Throws
[`InputError`](@ref) on the first problem.
"""
function load_spectra(cfg::AbstractDict, db::InputDatabase)
    validate_database(db)
    validate_config(cfg)
    grid = collect(Float64, wavelength_grid(cfg))
    expected = length(grid)
    codes = pigment_codes(cfg)
    names = ntuple(i -> String(cfg["pigments"][i]["name"]), Val(4))

    K = Matrix{Float64}(undef, 4, expected)
    S = Matrix{Float64}(undef, 4, expected)
    for (i, code) in enumerate(codes)
        for (quantity, out) in (("K", K), ("S", S))
            rows = sort(
                [r for r in db.spectra if r.code == code && r.quantity == quantity];
                by = r -> r.wavelength_nm,
            )
            isempty(rows) && throw(InputError("no $quantity spectrum for pigment $code"))
            length(rows) == expected || throw(
                InputError(
                    "$code $quantity has $(length(rows)) samples, the configuration grid has $expected"
                )
            )
            for (j, r) in enumerate(rows)
                r.wavelength_nm == grid[j] || throw(
                    InputError(
                        "$code $quantity sample $j is at $(r.wavelength_nm) nm, expected $(grid[j])"
                    )
                )
                out[i, j] = r.value
            end
        end
    end

    # The source states k1/k2 once for all pigments; the configuration pins
    # the same values. A mismatch is a different model, so it is an error.
    k1 = Float64(cfg["saunderson"]["k1"])
    k2 = Float64(cfg["saunderson"]["k2"])
    for (name, want, got) in (("k1", k1, db.saunderson.k1), ("k2", k2, db.saunderson.k2))
        isapprox(want, got; atol = 1.0e-9) || throw(
            InputError(
                "configuration Saunderson $name = $want, source has $got"
            )
        )
    end

    provenance = Dict{String, Any}(
        "source_file" => db.saunderson.source_file,
        "sheet" => db.saunderson.sheet,
        "wavelength_min_nm" => grid[1],
        "wavelength_max_nm" => grid[end],
        "wavelength_step_nm" => grid[2] - grid[1],
        "samples" => expected,
        "codes" => collect(codes),
        "saunderson_k1" => k1,
        "saunderson_k2" => k2,
        "kins_note" =>
            "kins is ignored: the source lists 1.0 and 0.0, and the paper's " *
            "equation (6) has no added specular term",
    )
    return PigmentSpectra(codes, names, grid, K, S, k1, k2, provenance)
end

"""
    build_quadrature(cfg, db) -> Quadrature

Validate the observer table against the configuration grid and precompute
trapezoidal weights. `norm = 1 / sum(w * y_bar * D65)` is the `Y_D65`
normalizer of equation (7).

The observer grid must be exactly the configuration's grid: resampling would
be a different model, not a runtime option.
"""
function build_quadrature(cfg::AbstractDict, db::InputDatabase)
    grid = collect(Float64, wavelength_grid(cfg))
    obs = sort(db.observer; by = r -> r.wavelength_nm)
    length(obs) == length(grid) || throw(
        InputError(
            "observer table has $(length(obs)) samples, the configuration grid has $(length(grid))"
        )
    )
    W = length(grid)
    xbar = Vector{Float64}(undef, W)
    ybar = Vector{Float64}(undef, W)
    zbar = Vector{Float64}(undef, W)
    d65 = Vector{Float64}(undef, W)
    for (j, r) in enumerate(obs)
        r.wavelength_nm == grid[j] || throw(
            InputError(
                "observer sample $j is at $(r.wavelength_nm) nm, expected $(grid[j])"
            )
        )
        xbar[j] = r.x_bar
        ybar[j] = r.y_bar
        zbar[j] = r.z_bar
        d65[j] = r.d65
    end
    integ = String(get(cfg["spectra"], "integration", "trapezoidal"))
    integ == "trapezoidal" || throw(
        ConfigError(
            "spectra.integration must remain \"trapezoidal\", got $(repr(integ))"
        )
    )
    weights = _trapezoid_weights(grid)
    ynorm = sum(weights[j] * ybar[j] * d65[j] for j in 1:W)
    ynorm > 0 || throw(InputError("integral of y_bar * D65 is not positive"))
    matrix = _color_matrix(cfg)
    return Quadrature(
        grid, weights, xbar, ybar, zbar, d65, matrix, 1 / ynorm,
    )
end

function _trapezoid_weights(grid::AbstractVector{T}) where {T}
    W = length(grid)
    w = Vector{T}(undef, W)
    W == 1 && return fill(one(T), 1)
    @inbounds for j in 1:W
        if j == 1
            w[j] = (grid[2] - grid[1]) / 2
        elseif j == W
            w[j] = (grid[W] - grid[W - 1]) / 2
        else
            w[j] = (grid[j + 1] - grid[j - 1]) / 2
        end
    end
    return w
end

function _color_matrix(cfg::AbstractDict)
    rows = cfg["color"]["xyz_to_rgb"]
    return ntuple(9) do i
        row = rows[(i - 1) ÷ 3 + 1]
        Float64(row[(i - 1) % 3 + 1])
    end
end

"""
    spectral_model(spectra::PigmentSpectra, quad::Quadrature) -> SpectralModel

Transpose the pigment-major spectra into the wavelength-major layout the
kernels read.
"""
function spectral_model(spectra::PigmentSpectra{T}, quad::Quadrature{T}) where {T}
    W = length(quad.wavelength)
    length(spectra.wavelength) == W || throw(
        InputError(
            "spectra have $(length(spectra.wavelength)) samples, quadrature has $W"
        )
    )
    K = Matrix{T}(undef, W, 4)
    S = Matrix{T}(undef, W, 4)
    @inbounds for i in 1:4, j in 1:W
        K[j, i] = spectra.K[i, j]
        S[j, i] = spectra.S[i, j]
    end
    return SpectralModel(K, S, quad, spectra.k1, spectra.k2)
end

"""
    with_parameters(model::SpectralModel, K, S) -> SpectralModel

A copy of `model` with new `4 x W` pigment-major absorption and scattering
matrices. Used by the surrogate fit to evaluate candidate pigments against
the fixed quadrature.
"""
function with_parameters(model::SpectralModel{T}, K::AbstractMatrix, S::AbstractMatrix) where {T}
    W, _ = size(model.K)
    size(K) == (4, W) || throw(InputError("K must be 4 x $W, got $(size(K))"))
    size(S) == (4, W) || throw(InputError("S must be 4 x $W, got $(size(S))"))
    Kt = Matrix{T}(undef, W, 4)
    St = Matrix{T}(undef, W, 4)
    @inbounds for i in 1:4, j in 1:W
        Kt[j, i] = T(K[i, j])
        St[j, i] = T(S[i, j])
    end
    return SpectralModel(Kt, St, model.quad, model.k1, model.k2)
end

"""
    pigment_parameters(model) -> (K, S)

The surrogate pigment parameters in the `4 x W` pigment-major layout, the
form the configuration and the provenance sidecar record.
"""
function pigment_parameters(model::SpectralModel{T}) where {T}
    W = size(model.K, 1)
    K = Matrix{T}(undef, 4, W)
    S = Matrix{T}(undef, 4, W)
    @inbounds for i in 1:4, j in 1:W
        K[i, j] = model.K[j, i]
        S[i, j] = model.S[j, i]
    end
    return K, S
end
