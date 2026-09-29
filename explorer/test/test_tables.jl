# Table-slice tests. The slice maps run on the synthetic model, which needs no
# disk. The spectral-error slice needs the input database and skips without it.

using Test
using StaticArrays: SVector
using PaintMix
using PaintMixPrecompute
using Explorer

@testset "tables" begin
    model = Explorer.synthetic_model()
    data = Explorer.ExplorerData(; model = model)

    @testset "slice index, level, and labels" begin
        @test Explorer.TABLE_AXES == (:r, :g, :b)
        @test Explorer.slice_index(256, 0.0) == 0
        @test Explorer.slice_index(256, 1.0) == 255
        @test Explorer.slice_index(256, 0.5) == 128
        @test Explorer.slice_index(2, 0.5) == 0
        @test Explorer.slice_index(1, 0.5) == 0
        @test Explorer.slice_level_value(256, 0.5) == 128 / 255
        @test Explorer.slice_axis_labels(:b) == ("R", "G")
        @test Explorer.slice_axis_labels(:g) == ("R", "B")
        @test Explorer.slice_axis_labels(:r) == ("G", "B")
    end

    @testset "concentration fields form a simplex" begin
        fields = Explorer.concentration_fields(model, :b, 0.5)
        @test length(fields) == 4
        @test all(f -> size(f) == (2, 2), fields)
        for iy in 1:2, ix in 1:2
            @test isapprox(sum(fields[slot][iy, ix] for slot in 1:4), 1.0; atol = 1.0e-12)
            @test all(f -> 0 <= f[iy, ix] <= 1, fields)
        end
        # The (0, 0, 0) vertex stores no pigment, so the white slot holds all.
        @test fields[1][1, 1] == 0.0
        @test fields[4][1, 1] == 1.0
    end

    @testset "basin map" begin
        fields = Explorer.concentration_fields(model, :b, 0.5)
        basin = Explorer.basin_map(fields)
        @test size(basin) == (2, 2)
        @test all(x -> 1 <= x <= 4, basin)
    end

    @testset "residual field" begin
        residual = Explorer.residual_field(model, :b, 0.5)
        @test size(residual) == (2, 2)
        @test all(>=(0), residual)
        # The (0, 0, 0) vertex is its own concentrations, so the residual is zero.
        @test residual[1, 1] == 0.0
        # The (1, 1, 0) vertex projects onto the simplex, so it is not.
        @test residual[2, 2] > 0
    end

    @testset "forward plane" begin
        plane = Explorer.forward_plane(model, 0.0; n = 3)
        @test plane isa Explorer.ForwardPlane
        @test plane.c4 == 0.0
        @test size(plane.rgb) == (3, 3)
        @test plane.valid[1, 1]
        @test !plane.valid[3, 3]
        @test plane.rgb[1, 1] == SVector(0.0, 0.0, 1.0)
        @test all(isnan, plane.rgb[3, 3])
        img = Explorer.forward_image(plane)
        @test size(img) == (3, 3)
        @test img[1, 1].alpha == 1.0
        @test img[3, 3].alpha == 0.0
    end

    @testset "slice data and caption without a sidecar" begin
        slices = Explorer.tables_slice_data(data, :b, 0.5)
        @test slices.spectral === nothing
        @test size(slices.basin) == (2, 2)
        caption = Explorer.tables_caption(data)
        @test occursin("Acceptance gates", caption)
        @test occursin("Not available: no sidecar.", caption)
        @test isempty(Explorer.acceptance_gate_rows(data))
        @test isempty(Explorer.promoted_gate_names(data))
    end

    @testset "sidecar validation accessors" begin
        sidecar = Dict{String, Any}(
            "acceptance_gates" => Dict{String, Any}("a" => true, "b" => false),
            "promoted_with_failed_gates" => ["b"],
            "validation" => Dict{String, Any}(
                "quality" => Dict{String, Any}(
                    "pairs" => 10,
                    "encoded_channel_error" => Dict{String, Any}("max" => 0.5),
                ),
            ),
        )
        d = Explorer.ExplorerData(; model = model, sidecar = sidecar)
        @test Explorer.acceptance_gate_rows(d) == [("a", true), ("b", false)]
        @test Explorer.promoted_gate_names(d) == ["b"]
        @test Explorer.validation_report(d, "quality")["pairs"] == 10
        @test Explorer.validation_report(d, "nope") === nothing
        caption = Explorer.tables_caption(d)
        @test occursin("overridden on promotion", caption)
        @test occursin("FAIL", caption)
    end

    @testset "number formatting" begin
        @test Explorer.number(1.23456) == "1.235"
        @test Explorer.number(nothing) == "—"
    end

    ref = spectral_reference()
    if ref === nothing
        @test_skip "input database absent; spectral-error slice not checked"
    else
        @testset "spectral-error slice" begin
            d = Explorer.ExplorerData(;
                cfg = ref.cfg, db = ref.db, spectra = ref.spectra, quad = ref.quad,
                spectral = ref.spectral, model = model,
            )
            field = Explorer.spectral_error_field(model, ref.spectral, :b, 0.5)
            @test size(field) == (Explorer.TABLE_SPECTRAL_N, Explorer.TABLE_SPECTRAL_N)
            @test all(isfinite, field)
            @test all(>=(0), field)
        end
    end
end
