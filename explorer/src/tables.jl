# `/tables` data: 2D slices through the inverse and forward lookup tables, and
# the validation-report accessors. Pure data, so the app file stays about
# figures and the tests can check the slices without a browser.
#
# A slice fixes one sRGB channel and varies the other two. The inverse table is
# read at its stored vertices, so a slice shows exactly what the payload holds;
# the spectral-error slice is the one exception, and it samples the runtime
# `encode` path at a reduced resolution because it evaluates the spectral model
# per sample.

"""
    TABLE_AXES

The three slice axes, one per fixed sRGB channel.
"""
const TABLE_AXES = (:r, :g, :b)

"""
    TABLE_SPECTRAL_N

Edge of the reduced-resolution grid the spectral-error slice uses. Evaluating
the spectral model per sample is the cost, so this stays below the payload
grid.
"""
const TABLE_SPECTRAL_N = 96

"""
    slice_index(n, level) -> Int

The zero-based grid index nearest the fraction `level` of the edge `n`.
"""
function slice_index(n::Integer, level::Real)
    n <= 1 && return 0
    return clamp(round(Int, Float64(level) * (n - 1)), 0, n - 1)
end

# The fixed-channel value of the nearest grid plane. The vertex reads use the
# quantized index, so every other map must use the same value to stay aligned.
function slice_level_value(n::Integer, level::Real)
    n <= 1 && return 0.0
    return slice_index(n, level) / (n - 1)
end

"""
    slice_axis_labels(axis) -> (xlabel, ylabel)

The free-channel labels of a slice, in plot order.
"""
slice_axis_labels(::Val{:b}) = ("R", "G")
slice_axis_labels(::Val{:g}) = ("R", "B")
slice_axis_labels(::Val{:r}) = ("G", "B")
slice_axis_labels(axis::Symbol) = slice_axis_labels(Val(axis))

# Zero-based (i, j, k) of the vertex at free indices `ix`, `iy` on the slice
# with fixed index `k`. The free axes match `slice_axis_labels`.
_slice_vertex(::Val{:b}, ix::Integer, iy::Integer, k::Integer) = (ix, iy, k)
_slice_vertex(::Val{:g}, ix::Integer, iy::Integer, k::Integer) = (ix, k, iy)
_slice_vertex(::Val{:r}, ix::Integer, iy::Integer, k::Integer) = (k, ix, iy)
_slice_vertex(axis::Symbol, ix::Integer, iy::Integer, k::Integer) =
    _slice_vertex(Val(axis), ix, iy, k)

# Linear-light RGB of a point on the slice, `x` and `y` the free coordinates.
_slice_rgb(::Val{:b}, level::Float64, x::Float64, y::Float64) = (x, y, level)
_slice_rgb(::Val{:g}, level::Float64, x::Float64, y::Float64) = (x, level, y)
_slice_rgb(::Val{:r}, level::Float64, x::Float64, y::Float64) = (level, x, y)
_slice_rgb(axis::Symbol, level::Float64, x::Float64, y::Float64) =
    _slice_rgb(Val(axis), level, x, y)

# The four concentrations stored at an inverse-table vertex. The fourth is
# implied, exactly as `decode` reconstructs it.
function _inverse_concentrations(model::PigmentModel, i::Integer, j::Integer, k::Integer)
    b = PaintMix.vertex(model.inverse, i, j, k)
    c1 = Float64(b[1]) / 255
    c2 = Float64(b[2]) / 255
    c3 = Float64(b[3]) / 255
    return (c1, c2, c3, max(0.0, 1.0 - c1 - c2 - c3))
end

"""
    concentration_fields(model, axis, level) -> Vector{Matrix{Float64}}

One `n x n` map per pigment of the inverse table's stored concentration on the
slice with `axis` fixed at `level`. `n` is the payload grid edge; the maps are
indexed `[y, x]`, matching the heatmap row/column order.
"""
function concentration_fields(model::PigmentModel, axis::Symbol, level::Real)
    n = PaintMix.grid_n(model)
    k = slice_index(n, level)
    fields = [Matrix{Float64}(undef, n, n) for _ in 1:4]
    for iy in 1:n, ix in 1:n
        c = _inverse_concentrations(model, _slice_vertex(axis, ix - 1, iy - 1, k)...)
        for slot in 1:4
            fields[slot][iy, ix] = c[slot]
        end
    end
    return fields
end

"""
    basin_map(fields) -> Matrix{Int}

The index of the dominant pigment at every vertex, 1–4. A tie goes to the
lowest slot, so the map is deterministic.
"""
function basin_map(fields::AbstractVector{<:AbstractMatrix})
    n = size(fields[1], 1)
    out = Matrix{Int}(undef, n, n)
    for iy in 1:n, ix in 1:n
        best = 1
        for slot in 2:4
            if fields[slot][iy, ix] > fields[best][iy, ix]
                best = slot
            end
        end
        out[iy, ix] = best
    end
    return out
end

"""
    residual_field(model, axis, level) -> Matrix{Float64}

The magnitude `‖x − M(U(x))‖` at every vertex of the slice: the error the
signed residual in a `Latent` carries there.
"""
function residual_field(model::PigmentModel, axis::Symbol, level::Real)
    n = PaintMix.grid_n(model)
    k = slice_index(n, level)
    fixed = slice_level_value(n, level)
    out = Matrix{Float64}(undef, n, n)
    for iy in 1:n, ix in 1:n
        c = _inverse_concentrations(model, _slice_vertex(axis, ix - 1, iy - 1, k)...)
        m = PaintMix.forward_rgb(model, c[1], c[2], c[3])
        x, y, z = _slice_rgb(axis, fixed, (ix - 1) / (n - 1), (iy - 1) / (n - 1))
        out[iy, ix] = sqrt((x - m[1])^2 + (y - m[2])^2 + (z - m[3])^2)
    end
    return out
end

"""
    ForwardPlane

The forward table on a fixed concentration plane: the fourth concentration
`c4` is fixed, `c1` runs along `x`, `c2` along `y`, and `c3` is implied.
`valid` marks the simplex triangle `c1 + c2 <= 1 - c4`; `rgb` is `NaN` outside.
"""
struct ForwardPlane
    c4::Float64
    rgb::Matrix{SVector{3, Float64}}
    valid::Matrix{Bool}
end

"""
    forward_plane(model, c4; n = grid_n) -> ForwardPlane

Evaluate `M(c)` over the concentration plane with the fourth concentration
fixed at `c4`.
"""
function forward_plane(model::PigmentModel, c4::Real; n::Integer = PaintMix.grid_n(model))
    m = Int(n)
    rgb = Matrix{SVector{3, Float64}}(undef, m, m)
    valid = Matrix{Bool}(undef, m, m)
    for iy in 1:m, ix in 1:m
        c1 = (ix - 1) / (m - 1)
        c2 = (iy - 1) / (m - 1)
        c3 = 1.0 - Float64(c4) - c1 - c2
        if c3 < 0
            valid[iy, ix] = false
            rgb[iy, ix] = SVector(NaN, NaN, NaN)
        else
            valid[iy, ix] = true
            f = PaintMix.forward_rgb(model, c1, c2, c3)
            rgb[iy, ix] = SVector(Float64(f[1]), Float64(f[2]), Float64(f[3]))
        end
    end
    return ForwardPlane(Float64(c4), rgb, valid)
end

"""
    spectral_error_field(model, spectral, axis, level; n = TABLE_SPECTRAL_N)

The worst channel difference between the forward table and the spectral model
at the concentrations the runtime inverse recovers, sampled `n x n` on the
slice. This is the forward-table approximation error at the points the model
actually visits, separate from the inverse-table error the residual carries.
"""
function spectral_error_field(
        model::PigmentModel, spectral::SpectralModel, axis::Symbol, level::Real;
        n::Integer = TABLE_SPECTRAL_N,
    )
    m = Int(n)
    fixed = slice_level_value(PaintMix.grid_n(model), level)
    out = Matrix{Float64}(undef, m, m)
    for iy in 1:m, ix in 1:m
        x = (ix - 1) / (m - 1)
        y = (iy - 1) / (m - 1)
        rgb = SVector{3, Float32}(_slice_rgb(axis, fixed, x, y))
        c = PaintMix.concentrations(PaintMix.encode(model, rgb))
        f = PaintMix.forward_rgb(model, c[1], c[2], c[3])
        s = mix_rgb_simplex(
            spectral, Float64(c[1]), Float64(c[2]), Float64(c[3])
        )
        out[iy, ix] = max(
            abs(Float64(f[1]) - Float64(s[1])),
            abs(Float64(f[2]) - Float64(s[2])),
            abs(Float64(f[3]) - Float64(s[3])),
        )
    end
    return out
end

# --- sidecar validation accessors -----------------------------------------

"""
    validation_report(d, key) -> Union{Nothing, Any}

One report from the sidecar's `validation` section, or `nothing` when the
sidecar or the report is absent.
"""
function validation_report(d::ExplorerData, key::AbstractString)
    validation = sidecar_section(d, "validation")
    validation isa AbstractDict || return nothing
    return Base.get(validation, key, nothing)
end

"""
    acceptance_gate_rows(d) -> Vector{Tuple{String, Bool}}

The acceptance gates as `(name, pass)` pairs, sorted by name. Empty without a
sidecar.
"""
function acceptance_gate_rows(d::ExplorerData)
    gates = sidecar_section(d, "acceptance_gates")
    gates isa AbstractDict || return Tuple{String, Bool}[]
    return [(String(k), Bool(v)) for (k, v) in sort(collect(gates); by = first)]
end

"""
    promoted_gate_names(d) -> Vector{String}

The gates the promoted payload overrode, from
`sidecar["promoted_with_failed_gates"]`. Empty when none were overridden.
"""
function promoted_gate_names(d::ExplorerData)
    failed = sidecar_section(d, "promoted_with_failed_gates")
    failed isa AbstractVector || return String[]
    return String[String(x) for x in failed]
end

"""
    number(x) -> String

A short number for a server-rendered table: four significant digits, or a dash
for a missing value.
"""
function number(x)
    x === nothing && return "—"
    return @sprintf("%.4g", Float64(x))
end
