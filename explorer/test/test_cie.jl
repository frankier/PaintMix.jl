# CIE-page tests. The geometry helpers (`convex_hull`, `clip_convex`) and the
# encoding helper run on hand-built inputs. The locus, primaries, pigment
# layers, mixture gamut, and probe are checked against the input database and
# skip when it is absent.

using Test
using StaticArrays: SVector
using PaintMix
using PaintMixPrecompute
using Explorer

@testset "cie" begin
    @testset "hex_from_linear clips and encodes" begin
        @test Explorer.hex_from_linear((1.0, 1.0, 1.0)) == "#ffffff"
        @test Explorer.hex_from_linear((0.0, 0.0, 0.0)) == "#000000"
        @test Explorer.hex_from_linear((-1.0, 2.0, 0.0)) == "#00ff00"
    end

    @testset "convex_hull drops interior and collinear points" begin
        pts = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.5, 0.5), (0.5, 0.0)]
        hull = Explorer.convex_hull(pts)
        @test length(hull) == 4
        @test all(p -> p in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)), hull)
        @test Explorer.convex_hull([(0.0, 0.0)]) == [(0.0, 0.0)]
    end

    @testset "clip_convex keeps a polygon inside the triangle" begin
        tri = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        square = [(0.25, 0.25), (2.0, 0.25), (2.0, 2.0), (0.25, 2.0)]
        clipped = Explorer.clip_convex(square, tri)
        @test !isempty(clipped)
        @test all(
            p -> p[1] >= -1.0e-12 && p[2] >= -1.0e-12 && p[1] + p[2] <= 1.0 + 1.0e-12,
            clipped,
        )
        # A clockwise clip polygon is normalized to counterclockwise.
        cw = [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)]
        @test !isempty(Explorer.clip_convex(square, cw))
    end

    ref = spectral_reference()
    if ref === nothing
        @test_skip "input database absent; CIE colorimetry not checked"
    else
        @testset "locus, primaries, and D65" begin
            locus = Explorer.spectral_locus(ref.quad)
            @test length(locus) == length(ref.quad.wavelength)
            @test all(
                p ->
                p[1] >= -1.0e-12 && p[2] >= -1.0e-12 && p[1] + p[2] <= 1.0 + 1.0e-12,
                locus,
            )
            primaries = Explorer.srgb_primaries(ref.quad)
            want = ((0.6401, 0.33), (0.3, 0.6), (0.15, 0.06))
            for i in 1:3
                @test isapprox(primaries[i][1], want[i][1]; atol = 1.0e-3)
                @test isapprox(primaries[i][2], want[i][2]; atol = 1.0e-3)
            end
            d65 = Explorer.d65_xy(ref.quad)
            @test isapprox(d65[1], 0.3127; atol = 1.0e-3)
            @test isapprox(d65[2], 0.329; atol = 1.0e-3)
        end

        @testset "wavelength fill is displayable" begin
            for j in eachindex(ref.quad.wavelength)
                c = Explorer.wavelength_linear(ref.quad, j)
                @test all(isfinite, c)
                @test all(x -> 0.0 <= x <= 1.0, c)
            end
        end

        @testset "pigment layers match the plan" begin
            curves = real_curves()
            @test curves !== nothing
            layers = Explorer.pigment_chromaticities(curves)
            # The plan's verified table: raw K/S chromaticity.
            measured = (
                (0.1972, 0.1186), (0.5292, 0.2731), (0.4747, 0.4929), (0.3124, 0.3301),
            )
            for i in 1:4
                @test isapprox(layers[i].measured[1], measured[i][1]; atol = 5.0e-4)
                @test isapprox(layers[i].measured[2], measured[i][2]; atol = 5.0e-4)
            end
            # The runtime forward table needs the untracked release payload.
            payload = joinpath(Explorer.ROOT, "data", "default", "default.pmx")
            if !isfile(payload)
                @test_skip "release payload absent; runtime chromaticity not checked"
            else
                model = PaintMix.read_model(payload)
                full = Explorer.derive_curves(
                    ref.spectra, ref.quad, ref.spectral, nothing, model, ref.db
                )
                runtime = (
                    (0.1954, 0.107), (0.5416, 0.2757), (0.4533, 0.4761), (0.316, 0.3344),
                )
                for i in 1:4
                    @test isapprox(full.runtime_xy[i][1], runtime[i][1]; atol = 5.0e-4)
                    @test isapprox(full.runtime_xy[i][2], runtime[i][2]; atol = 5.0e-4)
                end
            end
        end

        @testset "mixture gamut sits inside the sRGB triangle" begin
            tri = collect(Explorer.srgb_primaries(ref.quad))
            spectral_hull = Explorer.clip_convex(
                Explorer.convex_hull(Explorer.mixture_gamut(ref.spectral, ref.quad; d = 6)),
                tri,
            )
            @test length(spectral_hull) >= 3
            @test all(
                p ->
                p[1] >= -1.0e-9 && p[2] >= -1.0e-9 && p[1] + p[2] <= 1.0 + 1.0e-9,
                spectral_hull,
            )
        end

        @testset "probe runs the inverse lookup" begin
            curves = real_curves()
            model = Explorer.synthetic_model()
            p = Explorer.probe_xy(model, ref.quad, 0.3127, 0.329; derived = curves)
            @test p isa Explorer.CieProbe
            @test isapprox(sum(p.concentrations), 1.0; atol = 1.0e-4)
            @test p.nearest_pigment in curves.codes
            # Encoding then decoding cancels the residual exactly.
            @test all(
                isapprox(p.reconstructed_rgb[k], p.target_rgb[k]; atol = 1.0e-6) for k in 1:3
            )
            @test Explorer.probe_xy(model, ref.quad, 0.3, 0.0) === nothing
        end

        @testset "caption lists the reference points and nominal colors" begin
            ref_data = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
                derived = real_curves(),
            )
            caption = Explorer.cie_caption(ref_data)
            @test occursin("R primary", caption)
            @test occursin("D65 white", caption)
            @test occursin("Pigment chromaticity layers", caption)
            @test occursin("Nominal paper colors", caption)
            @test occursin("#2755ad", caption)
        end
    end
end
