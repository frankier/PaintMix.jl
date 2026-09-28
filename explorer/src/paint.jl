# The `/paint` canvas: a coarse bitmap of `Latent{Float32}` pixels plus the
# linear and display caches the figure streams.
#
# The dab is exactly the runtime blend: lerp the latents, then `decode`. The
# residual rides along, so overlapping strokes compose the way `mix` does and
# no stroke re-encodes. The naive-RGB mode lerps the display instead and
# re-encodes, which is the comparison the page is for. Everything here needs
# only the runtime payload and is free of Bonito, so it is testable on its own.

"""
    Brush

A dab color in both representations: the `Latent` the paint mode lerps and
the linear-light `RGB` the naive mode lerps.
"""
struct Brush
    latent::Latent{Float32}
    linear::SVector{3, Float32}
end

"""
    PaintCanvas

The paint state. `latent` is authoritative; `linear` and `display` are caches
kept in step over the dirty rectangle so a frame costs a copy, not a decode.

`display` holds sRGB-encoded colors, the form WGLMakie's `image!` expects.
`linear` holds linear light, the form the model and the naive blend use.
"""
struct PaintCanvas
    latent::Matrix{Latent{Float32}}
    linear::Matrix{SVector{3, Float32}}
    display::Matrix{Makie.RGBf}
end

const _PAPER_WHITE = SVector{3, Float32}(1, 1, 1)

"""
    paint_canvas(model; width = 128, height = 128) -> PaintCanvas

A canvas of paper white, `width` columns by `height` rows. The default edge is
128 because the whole display matrix crosses the websocket on every frame.
"""
function paint_canvas(
        model::PigmentModel; width::Integer = 128, height::Integer = 128
    )
    w = Int(width)
    h = Int(height)
    (w > 0 && h > 0) || throw(
        ArgumentError("canvas size must be positive, got $(w) x $(h)")
    )
    white = PaintMix.encode(model, _PAPER_WHITE)
    return PaintCanvas(
        fill(white, h, w), fill(_PAPER_WHITE, h, w), fill(linear_color(_PAPER_WHITE), h, w)
    )
end

"""
    brush_from_weights(model, weights) -> Brush

A brush from the four pigment weights: a `weighted_mix` of the pure pigments,
or paper white when every weight is zero.
"""
function brush_from_weights(model::PigmentModel, weights)
    c = _mixer_base(model, weights)
    return Brush(PaintMix.encode(model, c), c)
end

"""
    brush_from_linear(model, c) -> Brush

A brush from a linear-light color: its `encode` latent and the color itself.
"""
function brush_from_linear(model::PigmentModel, c)
    f = SVector{3, Float32}(Float32(c[1]), Float32(c[2]), Float32(c[3]))
    return Brush(PaintMix.encode(model, f), f)
end

# Parse `#rrggbb` or `rrggbb` into three bytes, or `nothing`.
function _parse_hex(hex::AbstractString)
    s = lstrip(hex, '#')
    length(s) == 6 || return nothing
    try
        return SVector{3, UInt8}(
            parse(UInt8, s[1:2]; base = 16), parse(UInt8, s[3:4]; base = 16),
            parse(UInt8, s[5:6]; base = 16),
        )
    catch
        return nothing
    end
end

"""
    brush_from_hex(model, hex) -> Brush

A brush from an encoded-sRGB `#rrggbb` string. An unparsable string gives
paper white, so the color picker never errors the session.
"""
function brush_from_hex(model::PigmentModel, hex::AbstractString)
    bytes = _parse_hex(hex)
    bytes === nothing && return brush_from_weights(model, (0.0, 0.0, 0.0, 0.0))
    return brush_from_linear(model, linear(bytes))
end

"""
    PickResult

The inverse lookup of a picked color: the encoded input, its linear light,
and the concentrations and residual `encode` recovers.
"""
struct PickResult
    hex::String
    linear::SVector{3, Float32}
    concentrations::SVector{4, Float32}
    residual::SVector{3, Float32}
end

"""
    picked_result(model, hex) -> Union{Nothing, PickResult}

The `encode` of a picked `#rrggbb` color, or `nothing` when the string is not
a color.
"""
function picked_result(model::PigmentModel, hex::AbstractString)
    bytes = _parse_hex(hex)
    bytes === nothing && return nothing
    c = SVector{3, Float32}(linear(bytes))
    z = PaintMix.encode(model, c)
    return PickResult(String(hex), c, PaintMix.concentrations(z), PaintMix.residual(z))
end

"""
    paint_pixel(pos, height) -> (col, row)

Map a mouse data position on an `image!` axis to canvas indices. `image!`
puts row 1 at the top, so the data y axis runs upward.
"""
paint_pixel(pos, height::Integer) = (round(Int, pos[1]), round(Int, Int(height) - pos[2] + 1))

# The bounding box of the dab disc, clamped to the canvas, or `nothing` when
# the disc lies entirely outside.
function _disc_rect(canvas::PaintCanvas, cx::Real, cy::Real, radius::Real)
    h, w = size(canvas.latent)
    r = Float64(radius)
    i0 = max(1, floor(Int, cy - r))
    i1 = min(h, ceil(Int, cy + r))
    j0 = max(1, floor(Int, cx - r))
    j1 = min(w, ceil(Int, cx + r))
    (i0 > i1 || j0 > j1) && return nothing
    return (i0, i1, j0, j1)
end

"""
    dab!(canvas, model, brush, cx, cy, radius, alpha, mode) -> rect

Stamp a hard-edged disc of radius `radius` centered at column `cx`, row `cy`,
blending each covered pixel toward `brush` by `alpha`. Returns the dirty
rectangle `(i0, i1, j0, j1)` in row-major bounds, or `nothing` when the dab
misses the canvas.

`mode` is `Val(:paint)` to lerp the latents, or `Val(:rgb)` to lerp the
linear display and re-encode.
"""
function dab!(
        canvas::PaintCanvas, model::PigmentModel, brush::Brush, cx::Real, cy::Real,
        radius::Real, alpha::Real, mode::Val
    )
    rect = _disc_rect(canvas, cx, cy, radius)
    rect === nothing && return nothing
    i0, i1, j0, j1 = rect
    α = Float32(clamp(Float64(alpha), 0.0, 1.0))
    r2 = Float64(radius)^2
    @inbounds for row in i0:i1, col in j0:j1
        (row - cy)^2 + (col - cx)^2 <= r2 || continue
        _dab_pixel!(canvas, model, brush, row, col, α, mode)
    end
    return rect
end

function _dab_pixel!(
        canvas::PaintCanvas, model::PigmentModel, brush::Brush, row::Int, col::Int,
        α::Float32, ::Val{:paint}
    )
    z = canvas.latent[row, col]
    znew = Latent(
        (1 - α) * z.c + α * brush.latent.c, (1 - α) * z.r + α * brush.latent.r
    )
    lin = PaintMix.decode(model, znew)
    canvas.latent[row, col] = znew
    canvas.linear[row, col] = lin
    canvas.display[row, col] = linear_color(lin)
    return
end

function _dab_pixel!(
        canvas::PaintCanvas, model::PigmentModel, brush::Brush, row::Int, col::Int,
        α::Float32, ::Val{:rgb}
    )
    lin = (1 - α) * canvas.linear[row, col] + α * brush.linear
    canvas.linear[row, col] = lin
    canvas.latent[row, col] = PaintMix.encode(model, lin)
    canvas.display[row, col] = linear_color(lin)
    return
end

"""
    clear!(canvas, model) -> rect

Reset every pixel to paper white and return the full-canvas dirty rectangle.
"""
function clear!(canvas::PaintCanvas, model::PigmentModel)
    white = PaintMix.encode(model, _PAPER_WHITE)
    disp = linear_color(_PAPER_WHITE)
    fill!(canvas.latent, white)
    fill!(canvas.linear, _PAPER_WHITE)
    fill!(canvas.display, disp)
    h, w = size(canvas.latent)
    return (1, h, 1, w)
end

"""
    canvas_image(canvas) -> Matrix{Makie.RGBf}

A copy of the display matrix, ready to hand to an `Observable` and `image!`.
"""
canvas_image(canvas::PaintCanvas) = copy(canvas.display)
