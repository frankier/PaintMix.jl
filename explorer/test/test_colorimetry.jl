# Colorimetry unit tests. The sRGB primaries and D65 come from the pinned
# matrix alone. The pigment chromaticities need the input database; when it is
# absent the testset skips, as the measurement workbooks are untracked.

using Test
using StaticArrays: SMatrix, SVector
using PaintMix
using PaintMixPrecompute
using Explorer

function config_matrix(cfg)
    r = cfg["color"]["xyz_to_rgb"]
    # `SMatrix` fills column-major; the configuration stores rows.
    return SMatrix{3, 3}(
        Float64(r[1][1]), Float64(r[2][1]), Float64(r[3][1]),
        Float64(r[1][2]), Float64(r[2][2]), Float64(r[3][2]),
        Float64(r[1][3]), Float64(r[2][3]), Float64(r[3][3]),
    )
end

function spectral_reference()
    cfg = load_config()
    db_path = joinpath(Explorer.ROOT, cfg["inputs"]["database"])
    isfile(db_path) || return nothing
    db = open_database(db_path)
    spectra = PigmentSpectra(cfg, db)
    quad = Quadrature(cfg, db)
    return (; cfg, db, spectra, quad, spectral = SpectralModel(spectra, quad))
end

# The independent route to the same chromaticity: integrate R' * D65 against
# the observer directly, without going through `mix_rgb` or the color matrix.
function xy_from_reflectance(spectra, quad, i::Integer)
    X = Y = Z = 0.0
    for l in eachindex(quad.wavelength)
        R = km_reflectance(spectra.K[i, l] / spectra.S[i, l])
        Rc = saunderson(R, spectra.k1, spectra.k2)
        w = quad.weights[l] * quad.d65[l] * Rc
        X += w * quad.x_bar[l]
        Y += w * quad.y_bar[l]
        Z += w * quad.z_bar[l]
    end
    s = X + Y + Z
    return (X / s, Y / s)
end

@testset "colorimetry" begin
    @testset "sRGB primaries and D65 from the pinned matrix" begin
        cfg = load_config()
        m_inv = inv(config_matrix(cfg))
        primaries = ((0.6401, 0.33), (0.3, 0.6), (0.15, 0.06))
        for (i, want) in enumerate(primaries)
            col = m_inv[:, i]
            s = sum(col)
            @test isapprox(col[1] / s, want[1]; atol = 1.0e-3)
            @test isapprox(col[2] / s, want[2]; atol = 1.0e-3)
        end
        w = m_inv * SVector(1.0, 1.0, 1.0)
        s = sum(w)
        @test isapprox(w[1] / s, 0.3127; atol = 1.0e-3)
        @test isapprox(w[2] / s, 0.329; atol = 1.0e-3)
    end

    ref = spectral_reference()
    if ref === nothing
        @test_skip "input database absent; spectral chromaticity not checked"
    else
        @testset "raw K/S chromaticity matches the plan" begin
            # The plan's verified table: raw K/S chromaticity (xy).
            want = (
                (0.1972, 0.1186), (0.5292, 0.2731), (0.4747, 0.4929), (0.3124, 0.3301),
            )
            for i in 1:4
                rgb = mix_rgb(ref.spectral, Explorer._PIGMENT_UNITS[i])
                xy = xy_of_linear(ref.quad, rgb)
                @test isapprox(xy[1], want[i][1]; atol = 5.0e-4)
                @test isapprox(xy[2], want[i][2]; atol = 5.0e-4)
            end
        end

        @testset "R'·D65 integral equals mix_rgb chromaticity" begin
            # Two independent routes to the same point: the library's
            # `mix_rgb` and a direct integration here.
            for i in 1:4
                direct = xy_from_reflectance(ref.spectra, ref.quad, i)
                rgb = mix_rgb(ref.spectral, Explorer._PIGMENT_UNITS[i])
                library = xy_of_linear(ref.quad, rgb)
                @test isapprox(direct[1], library[1]; atol = 1.0e-10)
                @test isapprox(direct[2], library[2]; atol = 1.0e-10)
            end
        end
    end
end
