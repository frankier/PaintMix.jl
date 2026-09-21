# Plan step 6 gate: validation and independent quality reports.
#
# The runtime tests cover algebra (exact endpoints, round trips, allocations).
# This file produces the evidence the plan's acceptance table asks for on the
# *generated* payload: mixing quality against the spectral reference, the
# padding rule's error near the simplex boundary, and quantization sanity.
#
# Reports are plain dictionaries so they serialize into the provenance
# sidecar and can be compared between runs.

# Sampling uses `Random.Xoshiro` from the standard library. The stream is not
# promised across Julia versions, but it is stable in practice, which is all
# the validation reports need.

@inline rand_rgb(rng::AbstractRNG) = (rand(rng), rand(rng), rand(rng))

# --- spectral reference implementation of encode/mix ----------------------

"""
    spectral_encode(model, rgb; settings) -> (c, r)

The paper's encoder `F` evaluated directly on the spectral model: solve
equation (9) with the reference solver, then store the residual against the
*same* spectral model. This is the reference the byte tables approximate.
"""
function spectral_encode(
        model::SpectralModel{T}, rgb::NTuple{3, T};
        settings::UnmixSettings{T} = UnmixSettings{T}(100, T(1.0e-10), T(1.0e-6), 4),
        scratch::SolverScratch = SolverScratch(),
    ) where {T}
    res = unmix_reference(model, rgb; settings = settings, scratch = scratch)
    m = mix_rgb(model, res.c)
    r = (rgb[1] - m[1], rgb[2] - m[2], rgb[3] - m[3])
    return res.c, r
end

"""
    spectral_mix(model, a, ca, ra, b, cb, rb, t) -> RGB

Reference `kmerp`: lerp the spectral latents and decode through `mix` of the
spectral model. Callers cache the encodings of `a` and `b`.
"""
@inline function spectral_mix(
        model::SpectralModel{T}, ca::NTuple{4, T}, ra::NTuple{3, T},
        cb::NTuple{4, T}, rb::NTuple{3, T}, t::T,
    ) where {T}
    s = one(T) - t
    c = ntuple(i -> s * ca[i] + t * cb[i], Val(4))
    r = ntuple(i -> s * ra[i] + t * rb[i], Val(3))
    m = mix_rgb(model, c)
    return (m[1] + r[1], m[2] + r[2], m[3] + r[3])
end

# --- reports ---------------------------------------------------------------

"""
    corpus(rng, count) -> Vector{RGB}

A deterministic validation corpus: the eight RGB cube corners, the four
pigment vertices and white, a spread of random colors, and some near-black
and near-white samples where residuals dominate.
"""
function corpus(rng::AbstractRNG, count::Integer)
    colors = NTuple{3, Float64}[]
    for r in (0.0, 1.0), g in (0.0, 1.0), b in (0.0, 1.0)
        push!(colors, (r, g, b))
    end
    push!(colors, (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (0.02, 0.02, 0.02), (0.98, 0.98, 0.98))
    while length(colors) < count
        push!(colors, rand_rgb(rng))
    end
    return colors
end

_pctl(v, p) = isempty(v) ? NaN : sort(v)[clamp(ceil(Int, p * length(v)), 1, length(v))]

"""
    quality_report(model, spectral, cfg; colors, pairs, seed) -> Dict

Compare the byte-table model against the spectral reference on a fixed
corpus: for each sampled pair and fraction, the encoded spectrum is taken
from the spectral model and compared with `PaintMix.mix` on the tables.

Reports encoded-sRGB channel errors and Oklab distances, so the two error
notions the plan separates (interpolation plus quantization, and perceptual
bias) stay separate.
"""
function quality_report(
        model::PaintMix.PigmentModel, spectral::SpectralModel{T}, cfg::AbstractDict;
        colors::Integer = 512, pairs::Integer = 20000, seed::Integer = 20210913,
    ) where {T}
    settings = unmix_settings(cfg)
    scratch = SolverScratch()
    rng = Xoshiro(seed)
    points = corpus(rng, colors)
    latents = Vector{Tuple{NTuple{4, T}, NTuple{3, T}}}(undef, length(points))
    for (i, rgb) in enumerate(points)
        latents[i] = spectral_encode(spectral, rgb; settings = settings, scratch = scratch)
    end

    ch = Float64[]
    enc = Float64[]
    ok = Float64[]
    worst = (0.0, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0), 0.0)
    @inbounds for _ in 1:pairs
        ia = 1 + floor(Int, rand(rng) * (length(points) - 1))
        ib = 1 + floor(Int, rand(rng) * (length(points) - 1))
        t = rand(rng)
        a = points[ia]
        b = points[ib]
        ca, ra = latents[ia]
        cb, rb = latents[ib]
        ref = spectral_mix(spectral, ca, ra, cb, rb, T(t))
        got = PaintMix.mix(model, a, b, Float32(t))
        e = (
            abs(Float64(got[1]) - ref[1]), abs(Float64(got[2]) - ref[2]),
            abs(Float64(got[3]) - ref[3]),
        )
        m = maximum(e)
        push!(ch, m)
        # The plan's quality target is on *encoded* sRGB, so report the
        # displayable-range error separately from the linear-light one.
        ga = PaintMix.srgb_from_linear(ntuple(i -> clamp(Float64(got[i]), 0.0, 1.0), 3))
        ra = PaintMix.srgb_from_linear(ntuple(i -> clamp(ref[i], 0.0, 1.0), 3))
        push!(enc, maximum(abs.(collect(ga) .- collect(ra))))
        d = oklab_distance_squared(
            linear_srgb_to_oklab((Float64(got[1]), Float64(got[2]), Float64(got[3]))),
            linear_srgb_to_oklab(ref)
        )
        push!(ok, sqrt(d))
        if m > worst[1]
            worst = (m, a, b, t)
        end
    end
    return Dict{String, Any}(
        "pairs" => pairs,
        "corpus_colors" => length(points),
        "channel_error" => Dict{String, Any}(
            "mean" => sum(ch) / length(ch),
            "p99" => _pctl(ch, 0.99),
            "max" => maximum(ch),
        ),
        "encoded_channel_error" => Dict{String, Any}(
            "mean" => sum(enc) / length(enc),
            "p99" => _pctl(enc, 0.99),
            "max" => maximum(enc),
        ),
        "oklab_error" => Dict{String, Any}(
            "mean" => sum(ok) / length(ok),
            "p99" => _pctl(ok, 0.99),
            "max" => maximum(ok),
        ),
        "worst_channel_sample" => Dict{String, Any}(
            "error" => worst[1],
            "a" => collect(worst[2]),
            "b" => collect(worst[3]),
            "t" => worst[4],
        ),
    )
end

"""
    padding_report(model, spectral, cfg; face_samples, inside_depth) -> Dict

Quantify the forward table's simplex padding. Samples points on each simplex
face and just inside it, and compares the byte-table forward lookup with the
spectral model evaluated at the same concentrations.

Off-simplex vertices are never a valid query, so the reported error measures
what the interpolation knots do to values *on* and just inside the domain of
interest, which is what the padding rule can affect.
"""
function padding_report(
        model::PaintMix.PigmentModel, spectral::SpectralModel{T}, cfg::AbstractDict;
        face_samples::Integer = 60, inside_depths = (0.0, 1.0e-3, 5.0e-3, 2.0e-2),
    ) where {T}
    n = PaintMix.grid_n(model)
    errs = Dict{Float64, Vector{Float64}}(d => Float64[] for d in inside_depths)
    worst = (0.0, (0.0, 0.0, 0.0, 0.0), 0.0)
    for face in 1:4
        free = [i for i in 1:4 if i != face]
        for a in 0:face_samples, b in 0:(face_samples - a)
            base = fill(0.0, 4)
            base[free[1]] = a / face_samples
            base[free[2]] = b / face_samples
            base[free[3]] = (face_samples - a - b) / face_samples
            for depth in inside_depths
                # Pull the point `depth` into the simplex along the face normal
                # so it is strictly inside when depth > 0.
                eps = depth / 4
                c = (
                    base[1] + eps, base[2] + eps, base[3] + eps, base[4] + eps,
                )
                s = c[1] + c[2] + c[3] + c[4]
                c = (c[1] / s, c[2] / s, c[3] / s, c[4] / s)
                ref = mix_rgb(spectral, c)
                got = PaintMix.forward_rgb(model, T(c[1]), T(c[2]), T(c[3]))
                e = maximum(
                    (
                        abs(Float64(got[1]) - ref[1]), abs(Float64(got[2]) - ref[2]),
                        abs(Float64(got[3]) - ref[3]),
                    )
                )
                push!(errs[depth], e)
                if e > worst[1]
                    worst = (e, c, depth)
                end
            end
        end
    end
    return Dict{String, Any}(
        "grid_n" => n,
        "padding_rule" => cfg["quantization"]["padding"],
        "samples_per_depth" => 4 * (face_samples + 1) * (face_samples + 2) ÷ 2,
        "by_depth" => Dict{String, Any}(
            string(d) => Dict{String, Any}(
                "max" => maximum(v),
                "p99" => _pctl(v, 0.99),
                "mean" => sum(v) / length(v),
            ) for (d, v) in errs
        ),
        "worst" => Dict{String, Any}(
            "error" => worst[1],
            "c" => collect(worst[2]),
            "depth" => worst[3],
        ),
    )
end

"""
    roundtrip_report(model, cfg; colors) -> Dict

Encoder/decoder round-trip error of the byte tables, at `Float32` (the
runtime scalar type) and `Float64`, over the deterministic corpus.
"""
function roundtrip_report(
        model::PaintMix.PigmentModel, cfg::AbstractDict; colors::Integer = 2000,
        seed::Integer = 7,
    )
    rng = Xoshiro(seed)
    points = corpus(rng, colors)
    err32 = Float64[]
    err64 = Float64[]
    for rgb in points
        a32 = (Float32(rgb[1]), Float32(rgb[2]), Float32(rgb[3]))
        z32 = PaintMix.encode(model, a32)
        back32 = PaintMix.decode(model, z32)
        push!(
            err32, maximum(
                abs.(
                    (
                        Float64(back32[1]) - rgb[1], Float64(back32[2]) - rgb[2],
                        Float64(back32[3]) - rgb[3],
                    )
                )
            )
        )
        z64 = PaintMix.encode(model, rgb)
        back64 = PaintMix.decode(model, z64)
        push!(
            err64, maximum(
                abs.(
                    (
                        back64[1] - rgb[1], back64[2] - rgb[2], back64[3] - rgb[3],
                    )
                )
            )
        )
    end
    return Dict{String, Any}(
        "samples" => length(points),
        "float32" => Dict{String, Any}("max" => maximum(err32), "mean" => sum(err32) / length(err32)),
        "float64" => Dict{String, Any}("max" => maximum(err64), "mean" => sum(err64) / length(err64)),
    )
end

"""
    quantization_report(model, cfg) -> Dict

Check the invariant the runtime depends on: every inverse-table vertex stores
three bytes summing to at most 255, so the implied fourth concentration is
non-negative. Also reports the distribution of the stored sums, which shows
how often the largest-remainder repair matters.
"""
function quantization_report(model::PaintMix.PigmentModel, cfg::AbstractDict)
    n = PaintMix.grid_n(model)
    d = model.inverse.data
    sums = zeros(Int, 256)
    bad = 0
    @inbounds for v in 0:(n^3 - 1)
        o = 3v
        s = Int(d[o + 1]) + Int(d[o + 2]) + Int(d[o + 3])
        s <= 255 || (bad += 1)
        sums[s + 1] += 1
    end
    n255 = count(==(255), sums)
    return Dict{String, Any}(
        "vertices" => n^3,
        "invalid_vertices" => bad,
        "sum_255_fraction" => n255 / n^3,
        "rule" => cfg["quantization"]["rule"],
        "note" =>
            "sum == 255 means the fourth concentration is zero; a sum below 255 " *
            "carries positive white pigment",
    )
end

"""
    acceptance_gates(cfg, reports) -> Dict{String,Bool}

Evaluate the provisional targets the plan sets for LUT-versus-spectral
mixing, plus the hard invariants. `reports` is the dictionary assembled by
`generate.jl`.
"""
function acceptance_gates(cfg::AbstractDict, reports::AbstractDict)
    q = reports["quality"]
    p = reports["padding"]
    rt = reports["roundtrip"]
    quant = reports["quantization"]
    targets = get(cfg, "validation", Dict{String, Any}())
    qe = get(q, "encoded_channel_error", q["channel_error"])
    p99_target = get(targets, "lut_channel_p99", 2 / 255)
    max_target = get(targets, "lut_channel_max", 8 / 255)
    oklab_target = get(targets, "lut_oklab_p99", 0.01)
    rt32 = get(targets, "roundtrip_f32", 2.0e-6)
    rt64 = get(targets, "roundtrip_f64", 1.0e-12)
    return Dict{String, Bool}(
        "lut_vs_spectral_channel_p99" => qe["p99"] <= p99_target,
        "lut_vs_spectral_channel_max" => qe["max"] <= max_target,
        "lut_vs_spectral_oklab_p99" => q["oklab_error"]["p99"] <= oklab_target,
        "padding_on_simplex" => p["by_depth"]["0.0"]["p99"] <= p99_target,
        "quantization_simplex" => quant["invalid_vertices"] == 0,
        "roundtrip_float32" => rt["float32"]["max"] <= rt32,
        "roundtrip_float64" => rt["float64"]["max"] <= rt64,
    )
end

"""
    continuity_report(model, cfg; lines, points, seed) -> Dict

Look for discontinuities between neighbouring RGB samples, which is what a
badly behaved inverse shows up as: two nearby colors whose concentrations
differ by a jump rather than a small step.

`lines` random straight segments are sampled at `points` positions each; the
reported statistic is the largest concentration jump and the largest decode
error between consecutive samples. A large jump with a small color change is
a basin change, not a modelling feature.
"""
function continuity_report(
        model::PaintMix.PigmentModel, cfg::AbstractDict;
        lines::Integer = 256, points::Integer = 65, seed::Integer = 424242,
    )
    rng = Xoshiro(seed)
    jumps = Float64[]
    decode_err = Float64[]
    worst_jump = (0.0, (0.0, 0.0, 0.0))
    for _ in 1:lines
        a = rand_rgb(rng)
        b = rand_rgb(rng)
        prev = nothing
        for i in 0:(points - 1)
            t = Float32(i / (points - 1))
            rgb = (
                Float32(a[1] + t * (b[1] - a[1])),
                Float32(a[2] + t * (b[2] - a[2])),
                Float32(a[3] + t * (b[3] - a[3])),
            )
            z = PaintMix.encode(model, rgb)
            if prev !== nothing
                c1, c2, c3, c4 = PaintMix.concentrations(prev)
                d1, d2, d3, d4 = PaintMix.concentrations(z)
                jump = max(abs(d1 - c1), abs(d2 - c2), abs(d3 - c3), abs(d4 - c4))
                push!(jumps, jump)
                if jump > worst_jump[1]
                    worst_jump = (jump, (Float64(rgb[1]), Float64(rgb[2]), Float64(rgb[3])))
                end
                back = PaintMix.decode(model, z)
                push!(
                    decode_err, maximum(
                        abs.(
                            (Float64(back[1]) - rgb[1], Float64(back[2]) - rgb[2], Float64(back[3]) - rgb[3]),
                        )
                    )
                )
            end
            prev = z
        end
    end
    return Dict{String, Any}(
        "lines" => lines,
        "points_per_line" => points,
        "concentration_jump" => Dict{String, Any}(
            "mean" => sum(jumps) / length(jumps),
            "p99" => _pctl(jumps, 0.99),
            "max" => maximum(jumps),
        ),
        "decode_error" => Dict{String, Any}(
            "p99" => _pctl(decode_err, 0.99),
            "max" => maximum(decode_err),
        ),
        "worst_jump_at" => collect(worst_jump[2]),
        "note" =>
            "runtime-table encode along random RGB segments; a concentration jump " *
            "much larger than 1/255 with a small color step is a basin change",
    )
end

"""
    behavior_report(model, spectral, cfg) -> Dict

The plan's visual-behavior row: blue and yellow must make a green, magenta
and yellow an orange, and adding white to a chromatic paint must raise
lightness. The qualitative flags are evaluated on the *spectral* reference,
because they are properties of the pigment data rather than of the
approximation; the runtime model is reported next to it, and the channel
difference between the two is recorded.

Also reports the same-color mixing residual and the reversal-symmetry
residual of the runtime model.
"""
function behavior_report(
        model::PaintMix.PigmentModel, spectral::SpectralModel{T}, cfg::AbstractDict,
    ) where {T}
    Tf = Float64
    blue = (0.02, 0.09, 0.42)
    yellow = (0.71, 0.62, 0.02)
    magenta = (0.55, 0.02, 0.12)
    white = (1.0, 1.0, 1.0)
    settings = unmix_settings(cfg)
    scratch = SolverScratch()
    encode_spec = rgb -> spectral_encode(
        spectral, (T(rgb[1]), T(rgb[2]), T(rgb[3])); settings = settings, scratch = scratch
    )
    cblue, rblue = encode_spec(blue)
    cyellow, ryellow = encode_spec(yellow)
    cmagenta, rmagenta = encode_spec(magenta)
    spec_mix = (a, b, t) -> spectral_mix(spectral, a[1], a[2], b[1], b[2], T(t))

    green_spec = spec_mix((cblue, rblue), (cyellow, ryellow), 0.5)
    orange_spec = spec_mix((cmagenta, rmagenta), (cyellow, ryellow), 0.5)
    green = PaintMix.mix(model, Tf.(blue), Tf.(yellow), 0.5)
    orange = PaintMix.mix(model, Tf.(magenta), Tf.(yellow), 0.5)

    luma(c) = 0.2126 * c[1] + 0.7152 * c[2] + 0.0722 * c[3]
    tints = [PaintMix.mix(model, Tf.(blue), Tf.(white), Float32(t)) for t in (0.0, 0.25, 0.5, 0.75, 1.0)]
    lightness = [luma(c) for c in tints]
    same = maximum(abs.(collect(PaintMix.mix(model, Tf.(blue), Tf.(blue), 0.37)) .- blue))
    reversal = maximum(
        abs.(
            collect(PaintMix.mix(model, Tf.(blue), Tf.(yellow), 0.3)) .-
                collect(PaintMix.mix(model, Tf.(yellow), Tf.(blue), 0.7)),
        )
    )
    return Dict{String, Any}(
        "spectral_blue_yellow_50_50" => collect(green_spec),
        "runtime_blue_yellow_50_50" => collect(Float64.(green)),
        "blue_yellow_channel_error" => maximum(abs.(collect(Float64.(green)) .- collect(green_spec))),
        "blue_yellow_is_green" =>
            green_spec[2] > green_spec[1] && green_spec[2] >= green_spec[3] - 0.05,
        "spectral_magenta_yellow_50_50" => collect(orange_spec),
        "runtime_magenta_yellow_50_50" => collect(Float64.(orange)),
        "magenta_yellow_channel_error" =>
            maximum(abs.(collect(Float64.(orange)) .- collect(orange_spec))),
        "magenta_yellow_is_orange" =>
            orange_spec[1] > orange_spec[2] && orange_spec[2] >= orange_spec[3] - 0.05,
        "white_tint_lightness" => collect(lightness),
        "white_tint_monotone" => all(diff(lightness) .> 0),
        "same_color_max_error" => same,
        "reversal_symmetry_error" => reversal,
        "note" =>
            "qualitative flags are evaluated on the spectral reference; the " *
            "runtime values and their channel difference show the approximation",
    )
end
