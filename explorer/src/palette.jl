# Palette data: pure-pigment vertices, nominal colors, mixing ramps, and the
# mixer that `/palette` and `/paint` share.
#
# Everything here needs only the runtime payload, so the page still works when
# the input database (and therefore the spectral reference) is absent. Mixing
# goes through `PaintMix.mix` and `PaintMix.weighted_mix`, never through a
# hand-rolled latent lerp, so the viewer shows exactly what the runtime does.

"""
    pigment_rgb(model, i) -> SVector{3, Float64}

The runtime forward-table color of pure pigment `i`, evaluated at `c = eᵢ`.
"""
pigment_rgb(model::PigmentModel, i::Integer) =
    SVector{3, Float64}(PaintMix.forward_rgb(model, _PIGMENT_UNITS[i]))

"""
    pigment_rgbs(model) -> NTuple{4, SVector{3, Float64}}

The four pure-pigment runtime colors, in slot order.
"""
pigment_rgbs(model::PigmentModel) = ntuple(i -> pigment_rgb(model, i), Val(4))

"""
    pigment_label(d, i) -> String

The measured pigment code, or the nominal paper name when the input database
is absent and no codes were derived.
"""
pigment_label(d::ExplorerData, i::Integer) =
    isempty(d.derived.codes[i]) ? NOMINAL_COLORS[i][1] : d.derived.codes[i]

"""
    nominal_rgb(i) -> SVector{3, Float64}

The nominal paper color of slot `i` in linear light.
"""
nominal_rgb(i::Integer) = SVector{3, Float64}(NOMINAL_COLORS[i][2]...)

"""
    MixCurve

A mixing path between two linear-light colors, sampled at `n` fractions: the
paint path through the model (`paint`) and the naive linear-RGB line
(`naive`), each as a `3 x n` matrix. `t` holds the fractions.
"""
struct MixCurve
    t::Vector{Float64}
    paint::Matrix{Float64}
    naive::Matrix{Float64}
end

function _rgb32(c)
    return SVector{3, Float32}(Float32(c[1]), Float32(c[2]), Float32(c[3]))
end

"""
    mix_curve(model, a, b; n = 64) -> MixCurve

Sample the paint mixture and the naive linear-RGB line between two
linear-light colors. `t = 0` returns `a` and `t = 1` returns `b` exactly.
"""
function mix_curve(model::PigmentModel, a, b; n::Integer = 64)
    n >= 2 || throw(ArgumentError("mix_curve needs at least two samples, got $n"))
    af = _rgb32(a)
    bf = _rgb32(b)
    ts = collect(range(0.0, 1.0; length = n))
    paint = Matrix{Float64}(undef, 3, n)
    naive = Matrix{Float64}(undef, 3, n)
    for (j, t) in enumerate(ts)
        p = PaintMix.mix(model, af, bf, Float32(t))
        q = af + Float32(t) * (bf - af)
        for k in 1:3
            paint[k, j] = p[k]
            naive[k, j] = q[k]
        end
    end
    return MixCurve(ts, paint, naive)
end

"""
    ramp_hex(model, a, b; n = 12) -> Vector{String}

The encoded-sRGB hex colors along the paint mixture between `a` and `b`, in
fraction order. Used to build the CSS ramp strips.
"""
function ramp_hex(model::PigmentModel, a, b; n::Integer = 12)
    curve = mix_curve(model, a, b; n = n)
    return [hex_from_linear(@view curve.paint[:, j]) for j in 1:n]
end

"""
    ramp_gradient(hexes) -> String

A CSS `linear-gradient` with one evenly spaced stop per hex color. The stops
are the sampled mixture, not a straight interpolation, so the strip bends
with the paint model.
"""
ramp_gradient(hexes) = "linear-gradient(to right, " * join(hexes, ", ") * ")"

"""
    pairwise_ramps(model; n = 12) -> Matrix{Vector{String}}

The 4×4 matrix of pure-pigment mixture ramps: entry `(i, j)` runs from
pigment `i` to pigment `j`. The diagonal is a constant pure color.
"""
function pairwise_ramps(model::PigmentModel; n::Integer = 12)
    colors = pigment_rgbs(model)
    return [ramp_hex(model, colors[i], colors[j]; n = n) for i in 1:4, j in 1:4]
end

"""
    MixerResult

The mixer's output: the ramped color in linear light, its encoded hex, and
the concentrations the inverse lookup recovers for it.
"""
struct MixerResult
    rgb::SVector{3, Float64}
    hex::String
    concentrations::SVector{4, Float32}
end

# Weighted average of the four pure pigments, or paper white when every weight
# is zero.
function _mixer_base(model::PigmentModel, weights)
    colors = [SVector{3, Float32}(pigment_rgb(model, i)) for i in 1:4]
    w = [Float32(max(0.0, Float64(x))) for x in weights]
    sum(w) <= 0 && return SVector{3, Float32}(1, 1, 1)
    return SVector{3, Float32}(PaintMix.weighted_mix(model, colors, w))
end

"""
    mixer_result(model, weights; t = 1.0) -> MixerResult

Blend the four pure pigments with `weights`, then ramp from paper white to
that mixture by `t`. Concentrations come from `PaintMix.encode`, so the swatch
is exactly what the model reproduces at its own latent.
"""
function mixer_result(model::PigmentModel, weights; t::Real = 1.0)
    base = _mixer_base(model, weights)
    paper = SVector{3, Float32}(1, 1, 1)
    rgb = PaintMix.mix(model, paper, base, Float32(clamp(Float64(t), 0.0, 1.0)))
    z = PaintMix.encode(model, rgb)
    return MixerResult(
        SVector{3, Float64}(rgb[1], rgb[2], rgb[3]), hex_from_linear(rgb),
        PaintMix.concentrations(z),
    )
end

"""
    mixer_ramp(model, weights; n = 16) -> Vector{String}

The encoded-sRGB hex colors from paper white to the weighted pigment mixture.
"""
function mixer_ramp(model::PigmentModel, weights; n::Integer = 16)
    base = _mixer_base(model, weights)
    paper = SVector{3, Float32}(1, 1, 1)
    return [
        hex_from_linear(PaintMix.mix(model, paper, base, Float32(t))) for
            t in range(0.0, 1.0; length = n)
    ]
end

"""
    NominalRow

One slot's nominal paper color against the runtime forward-table vertex: the
linear and encoded values of each, the per-channel residual, and the runtime
chromaticity when the spectral reference is available.
"""
struct NominalRow
    slot::Int
    name::String
    code::String
    nominal::SVector{3, Float64}
    nominal_hex::String
    runtime::SVector{3, Float64}
    runtime_hex::String
    residual::SVector{3, Float64}
    xy::Union{Nothing, Tuple{Float64, Float64}}
end

"""
    nominal_rows(d) -> Vector{NominalRow}

The nominal table: every slot's design color and its runtime vertex. The
residual is the displacement the fit and quantization introduced.
"""
function nominal_rows(d::ExplorerData)
    model = d.model
    quad = d.quad
    return [
        begin
            nominal = nominal_rgb(i)
            runtime = pigment_rgb(model, i)
            NominalRow(
                i, NOMINAL_COLORS[i][1], pigment_label(d, i), nominal,
                hex_from_linear(nominal), runtime, hex_from_linear(runtime),
                runtime - nominal,
                quad === nothing ? nothing : xy_of_linear(quad, runtime),
            )
        end for i in 1:4
    ]
end
