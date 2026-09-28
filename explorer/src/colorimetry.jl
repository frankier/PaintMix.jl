# Derived per-pigment curves and the linear/encoded-sRGB adapters every page
# shares.
#
# All colorimetry reuses `PaintMixPrecompute` (`mix_rgb`, `km_reflectance`,
# `saunderson`, `pigment_parameters`); nothing here re-derives the observer or
# the color matrix. The two conversion helpers exist so linear and encoded
# values are never mixed up: `encode` and `mix` take linear light, byte colors
# are encoded sRGB.

"""
    encoded(c) -> RGB8

Encode a linear-light sRGB triple as three sRGB bytes.
"""
encoded(c) = PaintMix.srgb8_from_linear(SVector{3, Float32}(c...))

"""
    linear(b) -> RGB{Float32}

Decode an encoded-sRGB byte triple to linear light.
"""
linear(b) = PaintMix.linear_from_srgb8(b)

"""
    xy_of_linear(quad, c) -> (x, y)

Chromaticity of a linear-light sRGB triple. Maps back to XYZ with the pinned
matrix and sum-normalizes, so it is the exact inverse of the transform
`mix_rgb` applies.
"""
function xy_of_linear(quad::Quadrature, c)
    xyz = inv(quad.xyz_to_rgb) * collect(Float64, c)
    s = sum(xyz)
    return (xyz[1] / s, xyz[2] / s)
end

# The four pure-pigment concentrations, in the model's fixed slot order.
const _PIGMENT_UNITS = ntuple(
    i -> SVector{4, Float64}(ntuple(j -> i == j ? 1.0 : 0.0, Val(4))), Val(4)
)

const _ZERO_RGB3 = SVector(0.0, 0.0, 0.0)
const _ZERO_XY = (0.0, 0.0)

# The spectral model evaluated at each pure pigment.
pure_rgb(m::SpectralModel) =
    ntuple(i -> SVector{3, Float64}(mix_rgb(m, _PIGMENT_UNITS[i])), Val(4))

"""
    DerivedCurves

Every curve and chromaticity the viewer shows, computed once per process from
the measured spectra, the fitted surrogate, and the runtime payload.

Fields are `4 x W` pigment-major for the measured and fitted spectra (the
layout of [`PigmentSpectra`](@ref PaintMixPrecompute.PigmentSpectra)), and
length-four tuples for the per-pigment color points. `K_cells` and `S_cells`
hold the source spreadsheet cell of each sample for the readout table.
"""
Base.@kwdef struct DerivedCurves
    wavelength::Vector{Float64} = Float64[]
    codes::NTuple{4, String} = ntuple(_ -> "", Val(4))
    K::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    S::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    KS::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    Rinf::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    Rprime::Matrix{Float64} = Matrix{Float64}(undef, 4, 0)
    K_cells::Matrix{String} = Matrix{String}(undef, 4, 0)
    S_cells::Matrix{String} = Matrix{String}(undef, 4, 0)
    fitted_K::Union{Nothing, Matrix{Float64}} = nothing
    fitted_S::Union{Nothing, Matrix{Float64}} = nothing
    fitted_KS::Union{Nothing, Matrix{Float64}} = nothing
    fitted_Rinf::Union{Nothing, Matrix{Float64}} = nothing
    fitted_Rprime::Union{Nothing, Matrix{Float64}} = nothing
    measured_rgb::NTuple{4, SVector{3, Float64}} = ntuple(_ -> _ZERO_RGB3, Val(4))
    measured_xy::NTuple{4, Tuple{Float64, Float64}} = ntuple(_ -> _ZERO_XY, Val(4))
    fitted_rgb::Union{Nothing, NTuple{4, SVector{3, Float64}}} = nothing
    fitted_xy::Union{Nothing, NTuple{4, Tuple{Float64, Float64}}} = nothing
    runtime_rgb::NTuple{4, SVector{3, Float64}} = ntuple(_ -> _ZERO_RGB3, Val(4))
    runtime_xy::NTuple{4, Tuple{Float64, Float64}} = ntuple(_ -> _ZERO_XY, Val(4))
end

"""
    derive_curves(spectra, quad, spectral, fitted, model, db) -> DerivedCurves

Compute the measured `K`, `S`, `K/S`, `R∞`, `R′`, the source cell map, and the
three chromaticity layers for every pigment:

  * `measured_*` — the raw K/S spectral model at `c = eᵢ`,
  * `fitted_*` — the surrogate pigments, when a sidecar supplied `theta`,
  * `runtime_*` — the 8-bit forward table at `c = eᵢ`.
"""
function derive_curves(
        spectra::PigmentSpectra, quad::Quadrature, spectral::SpectralModel,
        fitted::Union{Nothing, SpectralModel}, model::PigmentModel,
        db::Union{Nothing, InputDatabase},
    )
    codes = spectra.codes
    W = length(spectra.wavelength)
    K = copy(spectra.K)
    S = copy(spectra.S)
    KS = K ./ S
    Rinf = km_reflectance.(KS)
    Rprime = saunderson.(Rinf, spectra.k1, spectra.k2)

    K_cells = fill("", 4, W)
    S_cells = fill("", 4, W)
    if db !== nothing
        slots = Dict(code => i for (i, code) in enumerate(codes))
        for r in db.spectra
            i = Base.get(slots, r.code, 0)
            i == 0 && continue
            j = findfirst(==(r.wavelength_nm), spectra.wavelength)
            j === nothing && continue
            (r.quantity == "K" ? K_cells : S_cells)[i, j] = r.cell_range
        end
    end

    measured_rgb = pure_rgb(spectral)
    measured_xy = ntuple(i -> xy_of_linear(quad, measured_rgb[i]), Val(4))
    runtime_rgb = ntuple(
        i -> SVector{3, Float64}(PaintMix.forward_rgb(model, _PIGMENT_UNITS[i])), Val(4)
    )
    runtime_xy = ntuple(i -> xy_of_linear(quad, runtime_rgb[i]), Val(4))

    fitted_K = fitted_S = fitted_KS = fitted_Rinf = fitted_Rprime = nothing
    fitted_rgb = fitted_xy = nothing
    if fitted !== nothing
        fitted_K, fitted_S = pigment_parameters(fitted)
        fitted_KS = fitted_K ./ fitted_S
        fitted_Rinf = km_reflectance.(fitted_KS)
        fitted_Rprime = saunderson.(fitted_Rinf, fitted.k1, fitted.k2)
        fitted_rgb = pure_rgb(fitted)
        fitted_xy = ntuple(i -> xy_of_linear(quad, fitted_rgb[i]), Val(4))
    end

    return DerivedCurves(;
        wavelength = spectra.wavelength, codes = codes, K = K, S = S, KS = KS,
        Rinf = Rinf, Rprime = Rprime, K_cells = K_cells, S_cells = S_cells,
        fitted_K = fitted_K, fitted_S = fitted_S, fitted_KS = fitted_KS,
        fitted_Rinf = fitted_Rinf, fitted_Rprime = fitted_Rprime,
        measured_rgb = measured_rgb, measured_xy = measured_xy,
        fitted_rgb = fitted_rgb, fitted_xy = fitted_xy,
        runtime_rgb = runtime_rgb, runtime_xy = runtime_xy,
    )
end

# --- spectral quantity accessors ------------------------------------------

"""
    SpectraReadoutRow

One pigment at one layer, with `values` aligned to `SpectraReadout.quantities`
and `cells` listing the source spreadsheet cells of the underlying K and S
samples. A `missing` entry means the fitted curve is not available.
"""
struct SpectraReadoutRow
    pigment::String
    layer::Symbol
    cells::String
    values::Vector{Union{Float64, Missing}}
end

"""
    SpectraReadout

The values of a [`spectra_readout`](@ref) at one grid wavelength: the column
keys, and one row per pigment and layer (measured, and fitted when present).
"""
struct SpectraReadout
    wavelength::Float64
    quantities::Vector{Symbol}
    rows::Vector{SpectraReadoutRow}
end

"""
    SPECTRA_QUANTITY_KEYS

The five spectral quantities the viewer plots, in display order. The names are
symbols so a readout can carry them without a type parameter; each is wrapped
in `Val` at the dispatch site.
"""
const SPECTRA_QUANTITY_KEYS = (:K, :S, :KS, :Rinf, :Rprime)

quantity_label(::Val{:K}) = "K"
quantity_label(::Val{:S}) = "S"
quantity_label(::Val{:KS}) = "K/S"
quantity_label(::Val{:Rinf}) = "R∞"
quantity_label(::Val{:Rprime}) = "R′"

"""
    quantity_label(q) -> String

The short axis/table label of a spectral quantity, e.g. `"K/S"` or `"R∞"`.
"""
quantity_label(q::Symbol) = quantity_label(Val(q))

quantity_log(::Val{:K}) = true
quantity_log(::Val{:S}) = true
quantity_log(::Val{:KS}) = true
quantity_log(::Val{:Rinf}) = false
quantity_log(::Val{:Rprime}) = false

"""
    quantity_log(q) -> Bool

Whether the quantity spans enough orders of magnitude to default to a log
axis.
"""
quantity_log(q::Symbol) = quantity_log(Val(q))

quantity_matrix(d::DerivedCurves, ::Val{:K}) = d.K
quantity_matrix(d::DerivedCurves, ::Val{:S}) = d.S
quantity_matrix(d::DerivedCurves, ::Val{:KS}) = d.KS
quantity_matrix(d::DerivedCurves, ::Val{:Rinf}) = d.Rinf
quantity_matrix(d::DerivedCurves, ::Val{:Rprime}) = d.Rprime

"""
    quantity_matrix(d, q) -> Matrix{Float64}

The measured `4 x W` curve of quantity `q`, pigment-major.
"""
quantity_matrix(d::DerivedCurves, q::Symbol) = quantity_matrix(d, Val(q))

fitted_matrix(d::DerivedCurves, ::Val{:K}) = d.fitted_K
fitted_matrix(d::DerivedCurves, ::Val{:S}) = d.fitted_S
fitted_matrix(d::DerivedCurves, ::Val{:KS}) = d.fitted_KS
fitted_matrix(d::DerivedCurves, ::Val{:Rinf}) = d.fitted_Rinf
fitted_matrix(d::DerivedCurves, ::Val{:Rprime}) = d.fitted_Rprime

"""
    fitted_matrix(d, q) -> Union{Nothing, Matrix{Float64}}

The fitted `4 x W` curve of quantity `q`, or `nothing` when no sidecar supplied
surrogate parameters.
"""
fitted_matrix(d::DerivedCurves, q::Symbol) = fitted_matrix(d, Val(q))

quantity_cell(d::DerivedCurves, i::Integer, j::Integer, ::Val{:K}) = d.K_cells[i, j]
quantity_cell(d::DerivedCurves, i::Integer, j::Integer, ::Val{:S}) = d.S_cells[i, j]
function quantity_cell(d::DerivedCurves, i::Integer, j::Integer, ::Val)
    # K/S, R∞, and R′ all consume both source curves.
    return "$(d.K_cells[i, j]) / $(d.S_cells[i, j])"
end
quantity_cell(d::DerivedCurves, i::Integer, j::Integer, q::Symbol) =
    quantity_cell(d, i, j, Val(q))

"""
    wavelength_index(d, λ) -> Int

The index of the grid sample nearest `λ`, or `0` when the grid is empty.
"""
function wavelength_index(d::DerivedCurves, λ::Real)
    w = d.wavelength
    isempty(w) && return 0
    j = searchsortedfirst(w, Float64(λ))
    j = clamp(j, 1, length(w))
    if j > 1 && abs(w[j] - λ) > abs(w[j - 1] - λ)
        j -= 1
    end
    return j
end

"""
    spectra_readout(d, λ; pigments, quantities, fitted) -> SpectraReadout

Sample every enabled pigment/quantity pair at the grid wavelength nearest `λ`,
carrying the source spreadsheet cell of each value.
"""
function spectra_readout(
        d::DerivedCurves, λ::Real;
        pigments = ntuple(_ -> true, 4),
        quantities = ntuple(_ -> true, length(SPECTRA_QUANTITY_KEYS)),
        fitted::Bool = d.fitted_K !== nothing,
    )
    keys = Symbol[q for (qi, q) in enumerate(SPECTRA_QUANTITY_KEYS) if quantities[qi]]
    j = wavelength_index(d, λ)
    rows = SpectraReadoutRow[]
    if j != 0 && !isempty(keys)
        for i in 1:4
            pigments[i] || continue
            for layer in (:measured, :fitted)
                layer == :fitted && !fitted && continue
                values = Vector{Union{Float64, Missing}}(undef, length(keys))
                for (ci, q) in enumerate(keys)
                    m = layer == :fitted ? fitted_matrix(d, q) : quantity_matrix(d, q)
                    values[ci] = m === nothing ? missing : m[i, j]
                end
                push!(
                    rows, SpectraReadoutRow(
                        d.codes[i], layer,
                        join((quantity_cell(d, i, j, q) for q in (:K, :S)), " / "),
                        values,
                    )
                )
            end
        end
    end
    wavelength = j == 0 ? NaN : d.wavelength[j]
    return SpectraReadout(wavelength, keys, rows)
end

# --- CIE 1931 chromaticity diagram ----------------------------------------
#
# Everything on `/cie` is derived here: the spectral locus, the sRGB
# chromaticity triangle, the displayable color of each grid wavelength, the
# mixture gamut, and the click-to-probe inverse lookup. Chromaticities are
# plain `(x, y)` tuples; the diagram is a Makie 2D axis.

"""
    spectral_locus(quad) -> Vector{Tuple{Float64, Float64}}

The CIE 1931 locus at every grid wavelength. A monochromatic stimulus has
`XYZ ∝ (x̄, ȳ, z̄)`, so the chromaticity is the observer row sum-normalized.
The caller closes the curve with the line of purples.
"""
function spectral_locus(quad::Quadrature)
    n = length(quad.wavelength)
    points = Vector{Tuple{Float64, Float64}}(undef, n)
    @inbounds for j in 1:n
        s = quad.x_bar[j] + quad.y_bar[j] + quad.z_bar[j]
        points[j] = (Float64(quad.x_bar[j] / s), Float64(quad.y_bar[j] / s))
    end
    return points
end

"""
    srgb_primaries(quad) -> NTuple{3, Tuple{Float64, Float64}}

The chromaticity of the sRGB primaries: the columns of `inv(xyz_to_rgb)`,
sum-normalized. The `[0, 1]³` cube projects to exactly the triangle these
three points span.
"""
function srgb_primaries(quad::Quadrature)
    m = inv(quad.xyz_to_rgb)
    return ntuple(Val(3)) do i
        col = m[:, i]
        s = sum(col)
        (Float64(col[1] / s), Float64(col[2] / s))
    end
end

"""
    d65_xy(quad) -> (x, y)

The chromaticity of the D65 white point, `(1, 1, 1)` in linear light.
"""
d65_xy(quad::Quadrature) = xy_of_linear(quad, SVector(1.0, 1.0, 1.0))

"""
    wavelength_linear(quad, j) -> SVector{3, Float64}

The displayable linear-light sRGB color of grid wavelength `j`: the
monochromatic `XYZ` at unit luminance, desaturated toward white until it lies
in the sRGB cube. The locus point itself is the exact, unclamped
chromaticity; only this fill color is adjusted for display.
"""
function wavelength_linear(quad::Quadrature, j::Integer)
    y = Float64(quad.y_bar[j])
    y <= 0 && return _ZERO_RGB3
    xyz = SVector(Float64(quad.x_bar[j]) / y, 1.0, Float64(quad.z_bar[j]) / y)
    rgb = quad.xyz_to_rgb * xyz
    m = min(rgb[1], rgb[2], rgb[3])
    if m < 0
        # Add the smallest amount of white that lifts the worst channel to 0.
        t = -m / (1 - m)
        rgb = rgb .+ t .* (1 .- rgb)
    end
    return SVector{3, Float64}(
        clamp(rgb[1], 0.0, 1.0), clamp(rgb[2], 0.0, 1.0), clamp(rgb[3], 0.0, 1.0)
    )
end

"""
    hex_from_linear(c) -> String

The `#rrggbb` encoded-sRGB string of a linear-light triple, clipped to the
cube exactly as `srgb8_from_linear` clips.
"""
function hex_from_linear(c)
    b = encoded(c)
    return @sprintf("#%02x%02x%02x", b[1], b[2], b[3])
end

"""
    linear_color(c) -> Makie.RGBf

The displayable Makie color of a linear-light sRGB triple, clipped and
gamma-encoded for an sRGB canvas.
"""
function linear_color(c)
    s = PaintMix.srgb_from_linear(
        SVector{3, Float32}(Float32(c[1]), Float32(c[2]), Float32(c[3]))
    )
    return Makie.RGBf(
        clamp(s[1], 0.0f0, 1.0f0), clamp(s[2], 0.0f0, 1.0f0), clamp(s[3], 0.0f0, 1.0f0)
    )
end

"""
    PigmentChromaticity

One pigment's three chromaticity layers: the raw K/S spectral model
(`measured`), the fitted surrogate (`fitted`, or `nothing` without a
sidecar), and the 8-bit runtime forward table (`runtime`).
"""
struct PigmentChromaticity
    code::String
    measured::Tuple{Float64, Float64}
    fitted::Union{Nothing, Tuple{Float64, Float64}}
    runtime::Tuple{Float64, Float64}
end

"""
    pigment_chromaticities(d) -> Vector{PigmentChromaticity}

The three chromaticity layers of every pigment, in slot order.
"""
function pigment_chromaticities(d::DerivedCurves)
    return [
        PigmentChromaticity(
            d.codes[i], d.measured_xy[i],
            d.fitted_xy === nothing ? nothing : d.fitted_xy[i], d.runtime_xy[i],
        ) for i in 1:4
    ]
end

"""
    NOMINAL_COLORS

The nominal paper colors of the four pigments, in the same slot order as the
model, as `(name, linear RGB)` pairs. These are design targets, not
measurements.
"""
const NOMINAL_COLORS = (
    ("blue", (0.02, 0.09, 0.42)),
    ("magenta", (0.55, 0.02, 0.12)),
    ("yellow", (0.71, 0.62, 0.02)),
    ("white", (1.0, 1.0, 1.0)),
)

"""
    nominal_xy(quad) -> NTuple{4, Tuple{Float64, Float64}}

The chromaticity of the four nominal paper colors.
"""
function nominal_xy(quad::Quadrature)
    return ntuple(Val(4)) do i
        xy_of_linear(quad, SVector{3, Float64}(NOMINAL_COLORS[i][2]...))
    end
end

"""
    mixture_gamut(spectral, quad; d = 20) -> Vector{Tuple{Float64, Float64}}

The chromaticity of the spectral model at every point of the concentration
simplex surface with `d` divisions per axis. The convex hull of these points
is the reachable chromaticity region.
"""
function mixture_gamut(spectral::SpectralModel, quad::Quadrature; d::Integer = 20)
    sq = SurfaceQuadrature(d)
    return [xy_of_linear(quad, mix_rgb(spectral, c)) for c in sq.points]
end

"""
    mixture_gamut(model, quad; d = 20) -> Vector{Tuple{Float64, Float64}}

The same region for the runtime payload, evaluated on the forward table.
"""
function mixture_gamut(model::PigmentModel, quad::Quadrature; d::Integer = 20)
    sq = SurfaceQuadrature(d)
    return [xy_of_linear(quad, PaintMix.forward_rgb(model, c)) for c in sq.points]
end

"""
    convex_hull(points) -> Vector{Tuple{Float64, Float64}}

The convex hull of 2D points, counterclockwise, by Andrew's monotone chain.
Collinear points on the hull edges are dropped.
"""
function convex_hull(points)
    pts = sort!(unique!([(Float64(p[1]), Float64(p[2])) for p in points]))
    length(pts) <= 2 && return pts
    cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
    function chain(iter)
        h = Tuple{Float64, Float64}[]
        for p in iter
            while length(h) >= 2 && cross(h[end - 1], h[end], p) <= 0
                pop!(h)
            end
            push!(h, p)
        end
        return h
    end
    lower = chain(pts)
    upper = chain(Iterators.reverse(pts))
    return vcat(lower[1:(end - 1)], upper[1:(end - 1)])
end

# Signed area of a polygon; positive when counterclockwise.
function _signed_area(poly)
    a = 0.0
    n = length(poly)
    for i in 1:n
        p, q = poly[i], poly[mod1(i + 1, n)]
        a += p[1] * q[2] - q[1] * p[2]
    end
    return a / 2
end

# Is `p` on the left of the directed line `a -> b`?
_left_of(a, b, p) = (b[1] - a[1]) * (p[2] - a[2]) - (b[2] - a[2]) * (p[1] - a[1]) >= 0

# Intersection of segments `p -> q` and `a -> b`, used by the clipper.
function _intersect(p, q, a, b)
    rx, ry = q[1] - p[1], q[2] - p[2]
    sx, sy = b[1] - a[1], b[2] - a[2]
    denom = rx * sy - ry * sx
    denom == 0 && return q
    t = ((a[1] - p[1]) * sy - (a[2] - p[2]) * sx) / denom
    return (p[1] + t * rx, p[2] + t * ry)
end

"""
    clip_convex(subject, clip) -> Vector{Tuple{Float64, Float64}}

Sutherland–Hodgman clip of a polygon to a convex clip polygon. Used to keep
the mixture gamut inside the sRGB triangle for display.
"""
function clip_convex(subject, clip)
    isempty(subject) && return Tuple{Float64, Float64}[]
    poly = collect(clip)
    _signed_area(poly) < 0 && reverse!(poly)
    out = [(Float64(p[1]), Float64(p[2])) for p in subject]
    n = length(poly)
    for i in 1:n
        a, b = poly[i], poly[mod1(i + 1, n)]
        input = out
        out = Tuple{Float64, Float64}[]
        isempty(input) && break
        for j in eachindex(input)
            cur = input[j]
            prev = input[mod1(j - 1, length(input))]
            cur_in = _left_of(a, b, cur)
            prev_in = _left_of(a, b, prev)
            if cur_in
                prev_in || push!(out, _intersect(prev, cur, a, b))
                push!(out, cur)
            elseif prev_in
                push!(out, _intersect(prev, cur, a, b))
            end
        end
    end
    return out
end

"""
    CieProbe

The result of a click on the chromaticity diagram: the target point, whether
it was inside the sRGB gamut, the clamped displayable target color, the
inverse lookup's four concentrations, the reconstructed color, the residual,
and the nearest pigment.
"""
struct CieProbe
    xy::Tuple{Float64, Float64}
    in_gamut::Bool
    target_rgb::SVector{3, Float64}
    target_hex::String
    concentrations::SVector{4, Float32}
    reconstructed_rgb::SVector{3, Float64}
    reconstructed_hex::String
    residual::SVector{3, Float64}
    nearest_pigment::Union{Nothing, String}
    nearest_distance::Float64
end

"""
    probe_xy(model, quad, x, y; derived = nothing) -> Union{Nothing, CieProbe}

Run the inverse lookup at a clicked chromaticity. `y = 1` fixes luminance and
`z` follows from `x + y + z = 1`; the resulting linear-light color is clamped
to the cube for display and then encoded. Returns `nothing` for `y ≈ 0`,
where the chromaticity is not a color.
"""
function probe_xy(
        model::PigmentModel, quad::Quadrature, x::Real, y::Real;
        derived::Union{Nothing, DerivedCurves} = nothing,
    )
    xf, yf = Float64(x), Float64(y)
    yf <= 1.0e-6 && return nothing
    xyz = SVector(xf / yf, 1.0, (1.0 - xf - yf) / yf)
    raw = quad.xyz_to_rgb * xyz
    in_gamut = all(c -> -1.0e-6 <= c <= 1.0 + 1.0e-6, raw)
    lin = SVector{3, Float64}(
        clamp(raw[1], 0.0, 1.0), clamp(raw[2], 0.0, 1.0), clamp(raw[3], 0.0, 1.0)
    )
    z = PaintMix.encode(model, SVector{3, Float32}(lin[1], lin[2], lin[3]))
    rec = PaintMix.decode(model, z)
    nearest = nothing
    distance = NaN
    if derived !== nothing
        best = Inf
        for i in 1:4
            m = derived.measured_xy[i]
            dd = hypot(m[1] - xf, m[2] - yf)
            if dd < best
                best = dd
                nearest = derived.codes[i]
            end
        end
        distance = best
    end
    c = PaintMix.concentrations(z)
    r = PaintMix.residual(z)
    return CieProbe(
        (xf, yf), in_gamut, lin, hex_from_linear(lin), c,
        SVector{3, Float64}(rec[1], rec[2], rec[3]), hex_from_linear(rec),
        SVector{3, Float64}(r[1], r[2], r[3]), nearest, distance,
    )
end
