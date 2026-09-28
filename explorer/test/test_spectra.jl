# Spectra-page tests. The pure readout and quantity-selection functions run on a
# hand-built `DerivedCurves`, so they never touch disk. The caption and the real
# figure are checked against the input database and skip when it is absent.

using Test
using PaintMix
using PaintMixPrecompute
using Explorer

function synthetic_curves()
    w = collect(380.0:10.0:750.0)
    W = length(w)
    K = [Float64(10i + j) for i in 1:4, j in 1:W]
    S = fill(2.0, 4, W)
    KS = K ./ S
    Rinf = km_reflectance.(KS)
    Rprime = saunderson.(Rinf, 0.03, 0.65)
    cells = [string("K", j) for _ in 1:4, j in 1:W]
    scells = [string("S", j) for _ in 1:4, j in 1:W]
    return Explorer.DerivedCurves(;
        wavelength = w, codes = ("A", "B", "C", "D"),
        K = K, S = S, KS = KS, Rinf = Rinf, Rprime = Rprime,
        K_cells = cells, S_cells = scells,
        fitted_K = K ./ 2, fitted_S = S, fitted_KS = (K ./ 2) ./ S,
        fitted_Rinf = km_reflectance.((K ./ 2) ./ S),
        fitted_Rprime = saunderson.(km_reflectance.((K ./ 2) ./ S), 0.03, 0.65),
    )
end

function real_curves()
    ref = spectral_reference()
    ref === nothing && return nothing
    curves = Explorer.derive_curves(
        ref.spectra, ref.quad, ref.spectral, nothing, Explorer.synthetic_model(), ref.db
    )
    return curves
end

@testset "spectra" begin
    c = synthetic_curves()

    @testset "quantity metadata" begin
        @test Explorer.SPECTRA_QUANTITY_KEYS == (:K, :S, :KS, :Rinf, :Rprime)
        @test Explorer.quantity_label(:K) == "K"
        @test Explorer.quantity_label(:KS) == "K/S"
        @test Explorer.quantity_label(:Rprime) == "R′"
        @test Explorer.quantity_log(:K)
        @test Explorer.quantity_log(:S)
        @test Explorer.quantity_log(:KS)
        @test !Explorer.quantity_log(:Rinf)
        @test !Explorer.quantity_log(:Rprime)
    end

    @testset "quantity matrices dispatch on the quantity" begin
        @test Explorer.quantity_matrix(c, :K) === c.K
        @test Explorer.quantity_matrix(c, :Rprime) === c.Rprime
        @test Explorer.fitted_matrix(c, :KS) === c.fitted_KS
        @test Explorer.quantity_cell(c, 1, 1, :K) == "K1"
        @test Explorer.quantity_cell(c, 1, 1, :KS) == "K1 / S1"
    end

    @testset "wavelength_index picks the nearest sample" begin
        j550 = findfirst(==(550.0), c.wavelength)
        @test Explorer.wavelength_index(c, 550.0) == j550
        @test Explorer.wavelength_index(c, 553.0) == j550
        @test Explorer.wavelength_index(c, 557.0) == j550 + 1
        @test Explorer.wavelength_index(c, 1.0) == 1
        @test Explorer.wavelength_index(c, 10_000.0) == length(c.wavelength)
        @test Explorer.wavelength_index(Explorer.DerivedCurves(), 550.0) == 0
    end

    @testset "readout samples every enabled pigment and quantity" begin
        # The fitted matrices exist, so fitted rows are on by default.
        r = Explorer.spectra_readout(c, 550.0)
        @test r.wavelength == 550.0
        @test r.quantities == collect(Explorer.SPECTRA_QUANTITY_KEYS)
        @test length(r.rows) == 8
        @test count(row -> row.layer == :fitted, r.rows) == 4

        measured = Explorer.spectra_readout(c, 550.0; fitted = false)
        @test length(measured.rows) == 4
        @test all(row -> row.layer == :measured, measured.rows)
        @test measured.rows[1].pigment == "A"
        j = Explorer.wavelength_index(c, 550.0)
        @test measured.rows[1].values[1] == c.K[1, j]
        @test measured.rows[1].values[3] == c.KS[1, j]
        @test measured.rows[1].cells == "K$j / S$j"

        subset = Explorer.spectra_readout(
            c, 550.0;
            pigments = (true, false, false, false),
            quantities = (true, false, false, false, false), fitted = false,
        )
        @test subset.quantities == [:K]
        @test length(subset.rows) == 1
        @test subset.rows[1].values == [c.K[1, j]]

        @test r.rows[2].values[1] == c.fitted_K[1, j]
    end

    @testset "real data: caption and readout" begin
        curves = real_curves()
        if curves === nothing
            @test_skip "input database absent; spectra caption not checked"
        else
            ref = spectral_reference()
            d = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = Explorer.synthetic_model(),
                derived = curves,
            )
            caption = Explorer.spectra_caption(d)
            @test occursin("380–750 nm", caption)
            @test occursin("Saunderson k1", caption)
            @test occursin("kins", caption)

            r = Explorer.spectra_readout(curves, 550.0)
            @test all(row -> !isempty(row.cells), r.rows)
            @test all(row -> all(!ismissing, row.values), r.rows)
        end
    end
end
