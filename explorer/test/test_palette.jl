# Palette-page tests. The ramp, mixer, and nominal-table helpers run on the
# synthetic model, which needs no disk. The pure-pigment spectral RGB values
# are checked against the plan's verified table and skip when the input
# database is absent.

using Test
using StaticArrays: SVector
using PaintMix
using PaintMixPrecompute
using Explorer

@testset "palette" begin
    model = Explorer.synthetic_model()
    data = Explorer.ExplorerData(; model = model)

    @testset "pure-pigment vertices" begin
        for i in 1:4
            @test Explorer.pigment_rgb(model, i) == SVector{3, Float64}(
                ntuple(j -> i == j ? 1.0 : 0.0, Val(3))
            )
        end
        @test length(Explorer.pigment_rgbs(model)) == 4
        # Without measured spectra the labels fall back to the nominal names.
        @test Explorer.pigment_label(data, 1) == "blue"
        @test Explorer.pigment_label(data, 4) == "white"
        @test Explorer.nominal_rgb(1) == SVector(0.02, 0.09, 0.42)
        @test Explorer.nominal_rgb(4) == SVector(1.0, 1.0, 1.0)
    end

    @testset "mix_curve hits the endpoints and the naive midpoint" begin
        a = (0.1, 0.2, 0.3)
        b = (0.9, 0.4, 0.1)
        curve = Explorer.mix_curve(model, a, b; n = 11)
        @test length(curve.t) == 11
        @test size(curve.paint) == (3, 11)
        @test size(curve.naive) == (3, 11)
        @test curve.t[1] == 0.0
        @test curve.t[end] == 1.0
        # Mixing runs in the runtime scalar type, Float32, so the endpoints
        # and the naive line agree to Float32 rounding, not exactly.
        @test all(isapprox.(curve.paint[:, 1], collect(a); atol = 1.0e-6))
        @test all(isapprox.(curve.paint[:, end], collect(b); atol = 1.0e-6))
        @test all(isapprox.(curve.naive[:, 1], collect(a); atol = 1.0e-6))
        @test all(isapprox.(curve.naive[:, end], collect(b); atol = 1.0e-6))
        # The naive line is linear interpolation.
        mid = 6
        for k in 1:3
            @test isapprox(curve.naive[k, mid], (a[k] + b[k]) / 2; atol = 1.0e-6)
        end
        @test_throws ArgumentError Explorer.mix_curve(model, a, b; n = 1)
    end

    @testset "ramp hex and gradient" begin
        a = (0.0, 0.0, 0.0)
        b = (1.0, 1.0, 1.0)
        hexes = Explorer.ramp_hex(model, a, b; n = 8)
        @test length(hexes) == 8
        @test first(hexes) == "#000000"
        @test last(hexes) == "#ffffff"
        gradient = Explorer.ramp_gradient(hexes)
        @test startswith(gradient, "linear-gradient(to right, ")
        @test occursin("#000000", gradient)
        @test occursin("#ffffff", gradient)
    end

    @testset "pairwise ramps are a 4x4 matrix" begin
        ramps = Explorer.pairwise_ramps(model; n = 6)
        @test size(ramps) == (4, 4)
        for i in 1:4
            # A pigment mixed with itself stays that pigment.
            @test all(==(ramps[i, i][1]), ramps[i, i])
        end
        @test ramps[1, 1][1] == Explorer.hex_from_linear((1.0, 0.0, 0.0))
        @test ramps[1, 2][end] == Explorer.hex_from_linear((0.0, 1.0, 0.0))
    end

    @testset "mixer result and ramp" begin
        r = Explorer.mixer_result(model, (1.0, 0.0, 0.0, 0.0); t = 1.0)
        @test r isa Explorer.MixerResult
        @test r.hex == Explorer.hex_from_linear((1.0, 0.0, 0.0))
        @test isapprox(sum(r.concentrations), 1.0; atol = 1.0e-4)

        blank = Explorer.mixer_result(model, (0.0, 0.0, 0.0, 0.0))
        @test blank.hex == "#ffffff"

        ramp = Explorer.mixer_ramp(model, (1.0, 0.0, 0.0, 0.0); n = 5)
        @test length(ramp) == 5
        @test first(ramp) == "#ffffff"
        @test last(ramp) == r.hex
    end

    @testset "nominal rows" begin
        rows = Explorer.nominal_rows(data)
        @test length(rows) == 4
        @test rows[1].name == "blue"
        @test rows[1].nominal_hex == Explorer.hex_from_linear(Explorer.nominal_rgb(1))
        @test rows[1].runtime_hex == Explorer.hex_from_linear(
            Explorer.pigment_rgb(model, 1)
        )
        @test rows[1].residual == Explorer.pigment_rgb(model, 1) - Explorer.nominal_rgb(1)
        # No spectral reference, so no chromaticity column.
        @test all(row -> row.xy === nothing, rows)
    end

    @testset "server-rendered palette fragments" begin
        caption = Explorer.palette_caption(data)
        @test occursin("Nominal paper colors", caption)
        @test occursin("#2755ad", caption)
        @test occursin("other twenty pigments", caption)

        matrix = Explorer.palette_ramp_matrix(data)
        @test occursin("Pairwise mixture ramps", matrix)
        @test occursin("linear-gradient", matrix)
        @test length(findall("ramp", matrix)) >= 16
    end

    ref = spectral_reference()
    if ref === nothing
        @test_skip "input database absent; spectral pure RGB not checked"
    else
        @testset "spectral pure RGB matches the plan" begin
            # The plan's verified table: spectral pure c = eᵢ in linear sRGB.
            want = (
                (0.0115, 0.0059, 0.0706),
                (0.2048, 0.0013, 0.0313),
                (1.0389, 0.6332, -0.0542),
                (0.945, 0.9578, 0.9483),
            )
            for i in 1:4
                rgb = mix_rgb(ref.spectral, Explorer._PIGMENT_UNITS[i])
                for k in 1:3
                    @test isapprox(rgb[k], want[i][k]; atol = 5.0e-4)
                end
            end
        end

        @testset "real palette caption gets the chromaticity column" begin
            curves = real_curves()
            real = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
                derived = curves,
            )
            caption = Explorer.palette_caption(real)
            @test occursin("runtime xy", caption)
            @test occursin("sRGB reference", caption)
            rows = Explorer.nominal_rows(real)
            @test all(row -> row.xy !== nothing, rows)
        end
    end
end
