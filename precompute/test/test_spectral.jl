@testset "Kubelka-Munk limits" begin
    @test km_reflectance(0.0) == 1.0
    @test km_reflectance(1.0e12) < 1.0e-6
    @test km_reflectance(1.0) ≈ 1 / (2 + sqrt(3))
    # Equation (2) and its cancellation-free form agree.
    q = 0.3
    @test km_reflectance(q) ≈ 1 + q - sqrt(q * q + 2q)
    @test saunderson(1.0, 0.03, 0.65) ≈ (0.97 * 0.35) / 0.35
    @test saunderson(0.0, 0.03, 0.65) == 0.0
end


@testset "Kubelka-Munk is invariant to a common K and S scale" begin
    _, _, spectra, quad, model = synthetic_model()
    scaled = deepcopy(spectra)
    scaled.K .= 37.0 .* spectra.K
    scaled.S .= 37.0 .* spectra.S
    other = SpectralModel(scaled, quad)
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(0.7, 0.1, 0.1, 0.1),
            SVector(0.0, 0.5, 0.5, 0.0),
        )
        @test maximum(
            abs.(collect(mix_rgb(model, c)) .- collect(mix_rgb(other, c)))
        ) < 1.0e-12
    end
end


@testset "spectral Jacobian matches finite differences" begin
    _, _, _, _, model = synthetic_model()
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(0.7, 0.1, 0.1, 0.1),
            SVector(0.0, 0.4, 0.3, 0.3),
        )
        J = zeros(3, 4)
        mix_rgb_jacobian!(J, model, c)
        h = 1.0e-6
        for i in 1:4
            cp = SVector{4}(ntuple(j -> j == i ? c[j] + h : c[j], 4))
            cm = SVector{4}(ntuple(j -> j == i ? c[j] - h : c[j], 4))
            rp = mix_rgb(model, cp)
            rm = mix_rgb(model, cm)
            for r in 1:3
                @test J[r, i] ≈ (rp[r] - rm[r]) / (2h) rtol = 1.0e-5 atol = 1.0e-7
            end
        end
    end
end


@testset "Oklab conventions" begin
    black = linear_srgb_to_oklab(SVector(0.0, 0.0, 0.0))
    white = linear_srgb_to_oklab(SVector(1.0, 1.0, 1.0))
    @test maximum(abs, black) < 1.0e-12
    @test white[1] ≈ 1.0 atol = 1.0e-6
    @test abs(white[2]) < 1.0e-6
    @test abs(white[3]) < 1.0e-6
    # Signed cube root keeps out-of-gamut mixtures finite.
    o = linear_srgb_to_oklab(SVector(1.2, -0.1, 0.5))
    @test all(isfinite, o)
    @test oklab_distance_squared(white, white) == 0.0
end


@testset "cube penalty equals the signed distance squared" begin
    for p in (
            SVector(0.5, 0.5, 0.5), SVector(1.2, 0.5, 0.5),
            SVector(-0.1, -0.2, 0.5), SVector(1.1, 1.2, 1.3),
        )
        d = cube_signed_distance(p)
        want = d > 0 ? d^2 : 0.0
        @test cube_outside_penalty(p) ≈ want atol = 1.0e-12
    end
end
