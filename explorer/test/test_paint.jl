# Paint-page tests. The canvas, brush, and dab run on the synthetic model,
# which needs no disk and no browser.

using Test
using StaticArrays: SVector
using PaintMix
using Explorer

@testset "paint" begin
    model = Explorer.synthetic_model()
    data = Explorer.ExplorerData(; model = model)

    @testset "canvas starts as paper white" begin
        canvas = Explorer.paint_canvas(model; width = 4, height = 3)
        @test size(canvas.latent) == (3, 4)
        white = PaintMix.encode(model, SVector{3, Float32}(1, 1, 1))
        @test all(==(white), canvas.latent)
        @test all(c -> c == SVector{3, Float32}(1, 1, 1), canvas.linear)
        @test all(c -> c == Explorer.linear_color(SVector{3, Float32}(1, 1, 1)), canvas.display)
        @test_throws ArgumentError Explorer.paint_canvas(model; width = 0, height = 4)
    end

    @testset "brush from weights" begin
        for i in 1:4
            weights = ntuple(j -> j == i ? 1.0 : 0.0, Val(4))
            brush = Explorer.brush_from_weights(model, weights)
            @test brush.linear == SVector{3, Float32}(Explorer.pigment_rgb(model, i))
            @test brush.latent == PaintMix.encode(model, brush.linear)
        end
        blank = Explorer.brush_from_weights(model, (0.0, 0.0, 0.0, 0.0))
        @test blank.linear == SVector{3, Float32}(1, 1, 1)
    end

    @testset "brush from hex and picked result" begin
        red = Explorer.brush_from_hex(model, "#ff0000")
        @test red.linear == PaintMix.linear_from_srgb8(SVector{3, UInt8}(255, 0, 0))
        # An unparsable color degrades to paper white, not an error.
        @test Explorer.brush_from_hex(model, "nope").linear ==
            SVector{3, Float32}(1, 1, 1)
        @test Explorer.brush_from_hex(model, "#12345").linear ==
            SVector{3, Float32}(1, 1, 1)

        r = Explorer.picked_result(model, "#c42761")
        @test r isa Explorer.PickResult
        @test r.hex == "#c42761"
        @test isapprox(sum(r.concentrations), 1.0; atol = 1.0e-4)
        z = PaintMix.encode(model, r.linear)
        @test r.concentrations == PaintMix.concentrations(z)
        @test r.residual == PaintMix.residual(z)
        @test Explorer.picked_result(model, "bogus") === nothing
    end

    @testset "dab hits the disc and only the disc" begin
        canvas = Explorer.paint_canvas(model; width = 8, height = 8)
        brush = Explorer.brush_from_weights(model, (1.0, 0.0, 0.0, 0.0))
        rect = Explorer.dab!(canvas, model, brush, 4, 4, 1.5, 1.0, Val(:paint))
        @test rect == (2, 6, 2, 6)
        touched = 0
        for row in 1:8, col in 1:8
            inside = (row - 4)^2 + (col - 4)^2 <= 1.5^2
            if inside
                touched += 1
                @test canvas.latent[row, col] == brush.latent
                @test canvas.linear[row, col] == brush.linear
            else
                @test canvas.latent[row, col] ==
                    PaintMix.encode(model, SVector{3, Float32}(1, 1, 1))
            end
        end
        @test touched == 9
    end

    @testset "dab clamps to the canvas and misses cleanly" begin
        canvas = Explorer.paint_canvas(model; width = 6, height = 6)
        brush = Explorer.brush_from_weights(model, (0.0, 1.0, 0.0, 0.0))
        rect = Explorer.dab!(canvas, model, brush, 1, 1, 2.0, 1.0, Val(:paint))
        @test rect == (1, 3, 1, 3)
        @test Explorer.dab!(canvas, model, brush, 50, 50, 2.0, 1.0, Val(:paint)) ===
            nothing
    end

    @testset "alpha zero is a no-op" begin
        canvas = Explorer.paint_canvas(model; width = 4, height = 4)
        before = copy(canvas.latent)
        brush = Explorer.brush_from_weights(model, (1.0, 0.0, 0.0, 0.0))
        Explorer.dab!(canvas, model, brush, 2, 2, 2.0, 0.0, Val(:paint))
        @test canvas.latent == before
    end

    @testset "RGB mode blends the linear display" begin
        canvas = Explorer.paint_canvas(model; width = 4, height = 4)
        brush = Explorer.brush_from_weights(model, (0.0, 0.0, 1.0, 0.0))
        α = 0.25f0
        Explorer.dab!(canvas, model, brush, 2, 2, 0.5, α, Val(:rgb))
        want = (1 - α) * SVector{3, Float32}(1, 1, 1) + α * brush.linear
        @test canvas.linear[2, 2] == want
        @test canvas.latent[2, 2] == PaintMix.encode(model, want)
        @test canvas.display[2, 2] == Explorer.linear_color(want)
    end

    @testset "clear resets the canvas" begin
        canvas = Explorer.paint_canvas(model; width = 5, height = 5)
        brush = Explorer.brush_from_weights(model, (1.0, 0.0, 0.0, 0.0))
        Explorer.dab!(canvas, model, brush, 3, 3, 2.0, 1.0, Val(:paint))
        rect = Explorer.clear!(canvas, model)
        @test rect == (1, 5, 1, 5)
        @test all(c -> c == SVector{3, Float32}(1, 1, 1), canvas.linear)
    end

    @testset "paint_pixel maps the image axis" begin
        # Row 1 is at the top, so the highest data y maps to row 1.
        @test Explorer.paint_pixel((3.0, 8.0), 8) == (3, 1)
        @test Explorer.paint_pixel((3.0, 1.0), 8) == (3, 8)
        @test Explorer.paint_pixel((2.6, 4.4), 8) == (3, 5)
    end

    @testset "canvas_image is a copy" begin
        canvas = Explorer.paint_canvas(model; width = 3, height = 3)
        image = Explorer.canvas_image(canvas)
        brush = Explorer.brush_from_weights(model, (1.0, 0.0, 0.0, 0.0))
        Explorer.dab!(canvas, model, brush, 2, 2, 1.0, 1.0, Val(:paint))
        @test image[2, 2] != canvas.display[2, 2]
    end

    @testset "caption and default size" begin
        caption = Explorer.paint_caption(data)
        @test occursin("How the dab works", caption)
        @test occursin("decode", caption)
        @test occursin("$(Explorer.PAINT_DEFAULT_SIZE)", caption)
    end

    @testset "a default frame is 192 KiB" begin
        canvas = Explorer.paint_canvas(
            model; width = Explorer.PAINT_DEFAULT_SIZE, height = Explorer.PAINT_DEFAULT_SIZE
        )
        image = Explorer.canvas_image(canvas)
        # Three Float32 channels per pixel; this is the buffer that crosses the
        # websocket on every throttled frame.
        @test sizeof(image) == 3 * 4 * Explorer.PAINT_DEFAULT_SIZE^2
    end
end
