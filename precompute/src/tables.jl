# Plan step 5: float table generation, simplex quantization, cube padding,
# and slab checkpointing.
#
# Two rules carry most of the subtlety:
#
#   * The inverse table stores three of four concentrations. They must be
#     quantized *jointly* to integers summing to 255. Independent rounding can
#     make the sum exceed 255, which would imply a negative fourth
#     concentration.
#   * The forward table's domain of interest is the concentration simplex, but
#     trilinear interpolation reads whole cubes. Off-simplex vertices are
#     padded by Euclidean projection onto the simplex before evaluation. The
#     paper does not specify this rule; `validate.jl` quantifies the error it
#     introduces near every face, edge, and vertex.

"""
    project_to_simplex(c) -> NTuple{4}

Euclidean projection of a four-vector onto `{c >= 0, sum(c) == 1}`.

The standard sort-and-threshold algorithm: with `v` sorted descending and
`css` its cumulative sums, the projection is `max(c - theta, 0)` where
`theta = (css[rho] - 1) / rho` for the largest `rho` that keeps
`v[rho] - theta > 0`.
"""
function project_to_simplex(c::NTuple{4,T}) where {T<:AbstractFloat}
    v = (c[1], c[2], c[3], c[4])
    # Sort descending; four elements, no allocation.
    w = [v[1], v[2], v[3], v[4]]
    sort!(w; rev = true)
    rho = 0
    theta = zero(T)
    css = zero(T)
    @inbounds for j in 1:4
        css += w[j]
        t = (css - one(T)) / j
        if w[j] - t > 0
            rho = j
            theta = t
        end
    end
    rho == 0 && return (T(0.25), T(0.25), T(0.25), T(0.25))
    return ntuple(i -> max(c[i] - theta, zero(T)), Val(4))
end

"""
    quantize_simplex(c) -> NTuple{3,UInt8}

Quantize all four concentrations jointly to non-negative integers summing to
255, then return the first three. The fourth is implied by the runtime and
equals `255 - (b1 + b2 + b3) >= 0`.

Uses largest-remainder rounding with deterministic tie-breaking: ties go to
the lower concentration index. Any leftover or excess after the floor is
distributed one unit at a time, largest remainder first.
"""
function quantize_simplex(c::NTuple{4,T}) where {T<:AbstractFloat}
    scaled = (c[1] * 255, c[2] * 255, c[3] * 255, c[4] * 255)
    base = ntuple(i -> floor(Int, scaled[i]), Val(4))
    rem = ntuple(i -> scaled[i] - base[i], Val(4))
    diff = 255 - (base[1] + base[2] + base[3] + base[4])
    order = _remainder_order(rem)
    if diff >= 0
        for i in 1:min(diff, 4)
            idx = order[i]
            base = _add_at(base, idx, 1)
        end
    else
        # sum(c) exceeded one slightly: remove from the smallest remainders.
        for i in 1:min(-diff, 4)
            idx = order[5 - i]
            if base[idx] > 0
                base = _add_at(base, idx, -1)
            end
        end
    end
    # Clamp any residual negative to zero and repair the sum, which can only
    # be needed if `c` was not exactly on the simplex.
    b = (max(base[1], 0), max(base[2], 0), max(base[3], 0), max(base[4], 0))
    s = b[1] + b[2] + b[3] + b[4]
    if s != 255
        b = _repair_byte_sum(b, s)
    end
    return (UInt8(b[1]), UInt8(b[2]), UInt8(b[3]))
end

function _remainder_order(rem::NTuple{4,T}) where {T}
    idx = (1, 2, 3, 4)
    # Insertion sort, descending remainder, ascending index on ties.
    for i in 2:4
        j = i
        while j > 1 && (rem[idx[j]] > rem[idx[j - 1]] ||
                   (rem[idx[j]] == rem[idx[j - 1]] && idx[j] < idx[j - 1]))
            idx = _swap_at(idx, j, j - 1)
            j -= 1
        end
    end
    return idx
end

@inline function _swap_at(t::NTuple{N,T}, i::Int, j::Int) where {N,T}
    return ntuple(k -> k == i ? t[j] : (k == j ? t[i] : t[k]), Val(N))
end

@inline function _add_at(t::NTuple{N,T}, i::Int, d::T) where {N,T}
    return ntuple(k -> k == i ? t[k] + d : t[k], Val(N))
end

function _repair_byte_sum(b::NTuple{4,Int}, s::Int)
    if s < 255
        return _add_at(b, 1, 255 - s)
    end
    excess = s - 255
    out = b
    for i in 1:4
        excess == 0 && break
        take = min(excess, out[i])
        out = _add_at(out, i, -take)
        excess -= take
    end
    return out
end

"""
    FloatTables{T}

Unquantized candidates for both tables, in the payload's channel-fastest
order: `3 * (i + n * (j + n * k)) + ch` with zero-based indices.

`forward` maps `(c1, c2, c3)` to linear sRGB; `inverse` maps `(r, g, b)` to
the first three concentrations. Both hold `3 n^3` values.
"""
struct FloatTables{T<:AbstractFloat}
    n::Int
    forward::Vector{T}
    inverse::Vector{T}
end

FloatTables(n::Integer, ::Type{T}) where {T} =
    FloatTables{T}(Int(n), zeros(T, 3 * Int(n)^3), zeros(T, 3 * Int(n)^3))

@inline _table_offset(n::Int, i::Int, j::Int, k::Int) = 3 * (i + n * (j + n * k))

"""
    forward_vertex!(dest, model, n, i, j, k)

Evaluate `mix(c)` at one forward-table vertex, padding with Euclidean
simplex projection when `(c1, c2, c3)` lies outside the simplex.
"""
@inline function forward_vertex!(
        dest::AbstractVector{T}, model::SpectralModel, n::Int, i::Int, j::Int, k::Int
    ) where {T}
    d = n - 1
    c1 = T(i) / d
    c2 = T(j) / d
    c3 = T(k) / d
    c4 = one(T) - c1 - c2 - c3
    if c4 < zero(T)
        c = project_to_simplex((c1, c2, c3, c4))
        c1, c2, c3, c4 = c[1], c[2], c[3], c[4]
    end
    rgb = mix_rgb(model, (c1, c2, c3, c4))
    o = _table_offset(n, i, j, k)
    @inbounds begin
        dest[o + 1] = rgb[1]
        dest[o + 2] = rgb[2]
        dest[o + 3] = rgb[3]
    end
    return nothing
end

"""
    generate_forward(model, n; threads) -> Vector{T}

Fill the forward table. Independent slabs are parallelized across `k`.
"""
function generate_forward(
        model::SpectralModel{T}, n::Integer; threads::Integer = Threads.nthreads()
    ) where {T}
    n = Int(n)
    n >= 2 || throw(ArgumentError("table size must be >= 2, got $n"))
    out = Vector{T}(undef, 3 * n^3)
    krange = 0:(n - 1)
    if threads <= 1 || n < 8
        for k in krange, j in 0:(n - 1), i in 0:(n - 1)
            forward_vertex!(out, model, n, i, j, k)
        end
    else
        Threads.@threads for k in krange
            for j in 0:(n - 1), i in 0:(n - 1)
                forward_vertex!(out, model, n, i, j, k)
            end
        end
    end
    return out
end

"""
    inverse_slab!(dest, model, n, k, scratch, settings, seedrow, prevrow)

Fill the slab `k` of the inverse table, traversing `j` then `i` in a fixed
order and seeding each solve from the two already-solved in-slab neighbours
`(i-1, j)` and `(i, j-1)`. When neither exists the uniform seed is used.

`seedrow` and `prevrow` are caller-owned `n`-element buffers of concentration
tuples. Keeping the traversal inside a slab makes the result independent of
how slabs are distributed across threads.
"""
function inverse_slab!(
        dest::AbstractVector{T}, model::SpectralModel{T}, n::Int, k::Int,
        scratch::SolverScratch, settings::UnmixSettings{T},
        seedrow::Vector{NTuple{4,T}}, prevrow::Vector{NTuple{4,T}},
    ) where {T}
    d = n - 1
    z = T(k) / d
    uniform = (T(0.25), T(0.25), T(0.25), T(0.25))
    cur = seedrow
    prev = prevrow
    @inbounds for j in 0:(n - 1)
        y = T(j) / d
        for i in 0:(n - 1)
            x = T(i) / d
            rgb = (x, y, z)
            # `cur[i]` is (i-1, j) from earlier in this row; `prev[i+1]` is
            # (i, j-1) from the previous row. Keeping the two buffers distinct
            # is what makes the traversal a proper 4-neighbour continuation.
            c = if j > 0 && i > 0
                unmix_bulk!(scratch, model, rgb, (cur[i], prev[i + 1]), settings).c
            elseif j > 0
                unmix_bulk!(scratch, model, rgb, (prev[i + 1],), settings).c
            elseif i > 0
                unmix_bulk!(scratch, model, rgb, (cur[i],), settings).c
            else
                unmix_bulk!(scratch, model, rgb, (uniform,), settings).c
            end
            cur[i + 1] = c
            o = _table_offset(n, i, j, k)
            dest[o + 1] = c[1]
            dest[o + 2] = c[2]
            dest[o + 3] = c[3]
        end
        cur, prev = prev, cur
    end
    return nothing
end

"""
    reference_slab!(dest, model, n, k, scratch, settings)

Fill slab `k` of the inverse table with the independent active-face
enumeration, one solve per vertex and no dependence on neighbouring cells.

This is the accuracy path. It costs roughly 15 times a single Newton solve
per vertex, but it is deterministic, parallel across slabs, and immune to the
local minima a continuation seed can fall into.
"""
function reference_slab!(
        dest::AbstractVector{T}, model, n::Int, k::Int,
        scratch::SolverScratch, settings::UnmixSettings{T},
    ) where {T}
    d = n - 1
    z = T(k) / d
    @inbounds for j in 0:(n - 1)
        y = T(j) / d
        for i in 0:(n - 1)
            rgb = (T(i) / d, y, z)
            r = unmix_reference(model, rgb; settings = settings, scratch = scratch)
            o = _table_offset(n, i, j, k)
            dest[o + 1] = r.c[1]
            dest[o + 2] = r.c[2]
            dest[o + 3] = r.c[3]
        end
    end
    return nothing
end

"""
    generate_inverse(model, n; threads, settings, on_slab, resume) -> Vector{T}

Fill the inverse table. Slabs are independent; each is passed to `on_slab`
when given, which is how checkpointing is wired in without making table
generation depend on the filesystem.

`resume` is an optional function `k -> Union{Nothing,Vector{T}}` returning an
already-completed slab. `on_slab(k, slab)` is called after a slab is filled.
"""
function generate_inverse(
        model, n::Integer;
        threads::Integer = Threads.nthreads(),
        settings::Union{Nothing,UnmixSettings} = nothing,
        solver::Symbol = :reference,
        resume = nothing, on_slab = nothing,
    )
    T = _scalar_type(model)
    st = settings === nothing ? UnmixSettings{T}(100, T(1.0e-10), T(1.0e-6), 4) : settings
    n = Int(n)
    n >= 2 || throw(ArgumentError("table size must be >= 2, got $n"))
    solver in (:reference, :bulk) || throw(ArgumentError(
        "solver must be :reference or :bulk, got $(repr(solver))"
    ))
    out = Vector{T}(undef, 3 * n^3)
    krange = collect(0:(n - 1))
    if threads <= 1 || n < 8
        seedrow = Vector{NTuple{4,T}}(undef, n)
        prevrow = Vector{NTuple{4,T}}(undef, n)
        scratch = SolverScratch()
        for k in krange
            _slab_or_generate!(
                out, model, n, k, scratch, st, seedrow, prevrow, resume, on_slab, solver
            )
        end
    else
        Threads.@threads for k in krange
            seedrow = Vector{NTuple{4,T}}(undef, n)
            prevrow = Vector{NTuple{4,T}}(undef, n)
            scratch = SolverScratch()
            _slab_or_generate!(
                out, model, n, k, scratch, st, seedrow, prevrow, resume, on_slab, solver
            )
        end
    end
    return out
end

function _slab_or_generate!(
        out, model, n, k, scratch, settings, seedrow, prevrow, resume, on_slab, solver,
    )
    if resume !== nothing
        loaded = resume(k)
        if loaded !== nothing
            o = 3 * k * n * n
            copyto!(out, o + 1, loaded, 1, 3 * n * n)
            return nothing
        end
    end
    if solver === :reference
        reference_slab!(out, model, n, k, scratch, settings)
    else
        inverse_slab!(out, model, n, k, scratch, settings, seedrow, prevrow)
    end
    if on_slab !== nothing
        o = 3 * k * n * n
        on_slab(k, view(out, o + 1:o + 3 * n * n))
    end
    return nothing
end

# --- slab checkpointing ----------------------------------------------------

"""
    CheckpointStore

Deterministic slab checkpoints with atomic completion markers. A store is
tied to one job hash: resuming a directory whose `job.json` records a
different configuration, input set, grid size, or surrogate parameters is an
error rather than a silent mix.
"""
struct CheckpointStore
    dir::String
    job_hash::String
    n::Int
    kind::String

    function CheckpointStore(
            dir::AbstractString, job_hash::AbstractString, n::Integer, kind::AbstractString,
        )
        d = abspath(dir)
        mkpath(d)
        store = new(String(d), String(job_hash), Int(n), String(kind))
        _write_job_manifest(store)
        return store
    end
end

function _job_manifest_path(store::CheckpointStore)
    return joinpath(store.dir, "job.json")
end

function _write_job_manifest(store::CheckpointStore)
    path = _job_manifest_path(store)
    if isfile(path)
        existing = _read_manifest(path)
        get(existing, "job_hash", "") == store.job_hash || throw(InputError(
            "checkpoint directory $(store.dir) belongs to job " *
                "$(get(existing, "job_hash", "?")) but this run has $(store.job_hash); " *
                "remove it or choose another --out"
        ))
        return nothing
    end
    open(path, "w") do io
        println(io, "job_hash = ", repr(store.job_hash))
        println(io, "grid_n = ", store.n)
        println(io, "kind = ", repr(store.kind))
        println(io, "schema_version = 1")
    end
    return nothing
end

function _read_manifest(path)
    d = Dict{String,Any}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line, "="; limit = 2)
        length(parts) == 2 || continue
        d[strip(parts[1])] = strip(strip(parts[2]), ['"'])
    end
    return d
end

_slab_path(store::CheckpointStore, k::Int) =
    joinpath(store.dir, "$(store.kind)_k$(lpad(k, 5, '0')).bin")
_slab_marker(store::CheckpointStore, k::Int) = _slab_path(store, k) * ".done"

"""
    slab_complete(store, k) -> Bool

Whether slab `k` has an atomic completion marker and a payload of the right
size.
"""
function slab_complete(store::CheckpointStore, k::Integer)
    path = _slab_path(store, Int(k))
    isfile(_slab_marker(store, Int(k))) || return false
    return isfile(path) && filesize(path) == 24 * store.n * store.n
end

"""
    write_slab(store, k, values)

Write slab `k` to a temporary file, then move it into place and create the
completion marker. A crash leaves either a completed slab or no marker.
"""
function write_slab(store::CheckpointStore, k::Integer, values)
    path = _slab_path(store, Int(k))
    tmp = path * ".tmp"
    open(tmp, "w") do io
        for v in values
            write(io, Float64(v))
        end
    end
    mv(tmp, path; force = true)
    open(_slab_marker(store, Int(k)), "w") do io
        print(io, store.job_hash)
    end
    return nothing
end

"""
    read_slab(store, k) -> Vector{Float64}

Load a completed slab, or `nothing` when it is absent.
"""
function read_slab(store::CheckpointStore, k::Integer)
    k = Int(k)
    slab_complete(store, k) || return nothing
    path = _slab_path(store, k)
    values = Vector{Float64}(undef, 3 * store.n * store.n)
    open(path, "r") do io
        for i in eachindex(values)
            values[i] = read(io, Float64)
        end
    end
    return values
end

"""
    slab_resume_function(store) -> Function

`resume(k)` for [`generate_inverse`](@ref), reading completed slabs from the
store.
"""
slab_resume_function(store::CheckpointStore) = k -> read_slab(store, k)

"""
    slab_callback_function(store) -> Function

`on_slab(k, values)` for [`generate_inverse`](@ref), writing each slab after
it is generated.
"""
slab_callback_function(store::CheckpointStore) = (k, values) -> write_slab(store, k, values)

"""
    completed_slabs(store) -> Vector{Int}

The indices of slabs already present, for progress reporting.
"""
completed_slabs(store::CheckpointStore) =
    Int[k for k in 0:(store.n - 1) if slab_complete(store, k)]

# --- float forward table as a fast evaluator -------------------------------
#
# The inverse table has to be filled at 16.7 million RGB vertices. Evaluating
# the spectral model costs about a microsecond per call; trilinear
# interpolation of the float forward table costs tens of nanoseconds and is
# the same function the runtime decodes with. The inverse is therefore solved
# against the float forward table, and the spectral model is used only to
# build that table. `validate.jl` measures the resulting LUT-versus-spectral
# mixing error, which is what the substitution could affect.

"""
    ForwardFloatLUT{T}

The forward table in floating point, in the payload's channel-fastest order.
Serves as the evaluator for the inverse solver through `_mix` and
`_jacobian4!`.
"""
struct ForwardFloatLUT{T<:AbstractFloat}
    n::Int
    data::Vector{T}
end

ForwardFloatLUT(n::Integer, ::Type{T}) where {T} =
    ForwardFloatLUT{T}(Int(n), Vector{T}(undef, 3 * Int(n)^3))

"""
    _scalar_type(ev) -> Type

The floating-point scalar type of a forward evaluator. The inverse solver
needs it to size its output before any solve has run.
"""
_scalar_type(::SpectralModel{T}) where {T} = T
_scalar_type(::ForwardFloatLUT{T}) where {T} = T

@inline function _axis_f(f::T, n::Int) where {T<:AbstractFloat}
    i = unsafe_trunc(Int, f)
    i > n - 2 && (i = n - 2)
    return i, f - T(i)
end

@inline _sub3(a::NTuple{3,T}, b::NTuple{3,T}) where {T} = (a[1] - b[1], a[2] - b[2], a[3] - b[3])
@inline _mul3(a::NTuple{3,T}, s::T) where {T} = (a[1] * s, a[2] * s, a[3] * s)

@inline function _mix(f::ForwardFloatLUT{T}, c::NTuple{4,S}) where {T,S<:Real}
    n = f.n
    d = f.data
    if n == 1
        return (d[1], d[2], d[3])
    end
    gx = clamp(S(c[1]), zero(S), one(S)) * (n - 1)
    gy = clamp(S(c[2]), zero(S), one(S)) * (n - 1)
    gz = clamp(S(c[3]), zero(S), one(S)) * (n - 1)
    i, fx = _axis_f(T(gx), n)
    j, fy = _axis_f(T(gy), n)
    k, fz = _axis_f(T(gz), n)
    b000 = 3 * (i + n * (j + n * k)) + 1
    b100 = b000 + 3
    b010 = b000 + 3n
    b110 = b010 + 3
    b001 = b000 + 3n * n
    b101 = b001 + 3
    b011 = b001 + 3n
    b111 = b011 + 3
    @inbounds begin
        c000 = (T(d[b000]), T(d[b000 + 1]), T(d[b000 + 2]))
        c100 = (T(d[b100]), T(d[b100 + 1]), T(d[b100 + 2]))
        c010 = (T(d[b010]), T(d[b010 + 1]), T(d[b010 + 2]))
        c110 = (T(d[b110]), T(d[b110 + 1]), T(d[b110 + 2]))
        c001 = (T(d[b001]), T(d[b001 + 1]), T(d[b001 + 2]))
        c101 = (T(d[b101]), T(d[b101 + 1]), T(d[b101 + 2]))
        c011 = (T(d[b011]), T(d[b011 + 1]), T(d[b011 + 2]))
        c111 = (T(d[b111]), T(d[b111 + 1]), T(d[b111 + 2]))
    end
    c00 = _lerp3f(c000, c100, fx)
    c10 = _lerp3f(c010, c110, fx)
    c01 = _lerp3f(c001, c101, fx)
    c11 = _lerp3f(c011, c111, fx)
    c0 = _lerp3f(c00, c10, fy)
    c1 = _lerp3f(c01, c11, fy)
    return _lerp3f(c0, c1, fz)
end

@inline function _lerp3f(a::NTuple{3,T}, b::NTuple{3,T}, w::T) where {T}
    s = one(T) - w
    return (s * a[1] + w * b[1], s * a[2] + w * b[2], s * a[3] + w * b[3])
end

@inline function _jacobian4!(J::AbstractMatrix, f::ForwardFloatLUT{T}, c::NTuple{4,S}) where {T,S<:Real}
    n = f.n
    if n == 1
        @inbounds for i in 1:4
            J[1, i] = zero(T)
            J[2, i] = zero(T)
            J[3, i] = zero(T)
        end
        return J
    end
    d = f.data
    gx = clamp(S(c[1]), zero(S), one(S)) * (n - 1)
    gy = clamp(S(c[2]), zero(S), one(S)) * (n - 1)
    gz = clamp(S(c[3]), zero(S), one(S)) * (n - 1)
    i, fx = _axis_f(T(gx), n)
    j, fy = _axis_f(T(gy), n)
    k, fz = _axis_f(T(gz), n)
    b000 = 3 * (i + n * (j + n * k)) + 1
    b100 = b000 + 3
    b010 = b000 + 3n
    b110 = b010 + 3
    b001 = b000 + 3n * n
    b101 = b001 + 3
    b011 = b001 + 3n
    b111 = b011 + 3
    @inbounds begin
        c000 = (T(d[b000]), T(d[b000 + 1]), T(d[b000 + 2]))
        c100 = (T(d[b100]), T(d[b100 + 1]), T(d[b100 + 2]))
        c010 = (T(d[b010]), T(d[b010 + 1]), T(d[b010 + 2]))
        c110 = (T(d[b110]), T(d[b110 + 1]), T(d[b110 + 2]))
        c001 = (T(d[b001]), T(d[b001 + 1]), T(d[b001 + 2]))
        c101 = (T(d[b101]), T(d[b101 + 1]), T(d[b101 + 2]))
        c011 = (T(d[b011]), T(d[b011 + 1]), T(d[b011 + 2]))
        c111 = (T(d[b111]), T(d[b111 + 1]), T(d[b111 + 2]))
    end
    c00 = _lerp3f(c000, c100, fx)
    c10 = _lerp3f(c010, c110, fx)
    c01 = _lerp3f(c001, c101, fx)
    c11 = _lerp3f(c011, c111, fx)
    c0 = _lerp3f(c00, c10, fy)
    c1 = _lerp3f(c01, c11, fy)
    z = one(T)
    x = one(T)
    # d/dx
    dx00 = _sub3(c100, c000)
    dx10 = _sub3(c110, c010)
    dx01 = _sub3(c101, c001)
    dx11 = _sub3(c111, c011)
    dx = _mul3(_lerp3f(_lerp3f(dx00, dx10, fy), _lerp3f(dx01, dx11, fy), fz), T(n - 1))
    # d/dy
    dy0 = _sub3(c10, c00)
    dy1 = _sub3(c11, c01)
    dy = _mul3(_lerp3f(dy0, dy1, fz), T(n - 1))
    # d/dz
    dz = _mul3(_sub3(c1, c0), T(n - 1))
    @inbounds begin
        J[1, 1], J[2, 1], J[3, 1] = dx
        J[1, 2], J[2, 2], J[3, 2] = dy
        J[1, 3], J[2, 3], J[3, 3] = dz
        # The table does not depend on the fourth concentration, but the
        # softmax directions sum to zero, so a zero column is the consistent
        # choice (see `_jacobian4!` in kubelka_munk.jl).
        J[1, 4] = zero(T)
        J[2, 4] = zero(T)
        J[3, 4] = zero(T)
    end
    return J
end

"""
    forward_float_lut(n, data) -> ForwardFloatLUT

Wrap an already-generated float forward table as an evaluator.
"""
forward_float_lut(n::Integer, data::AbstractVector{T}) where {T<:AbstractFloat} =
    ForwardFloatLUT{T}(Int(n), convert(Vector{T}, data))

# --- coarse-to-fine inverse ------------------------------------------------
#
# The reference solver costs minutes per slab at 256^3. Solving the same
# problem on a coarse RGB grid first and using the result as a seed gives the
# bulk solver a starting point in the right basin, so a handful of Newton
# steps per fine vertex suffice. The fine grid is then independent per vertex
# except for the seed, which keeps slab parallelism deterministic.

"""
    coarse_seed(coarse, cn, r, g, b) -> NTuple{4,T}

Trilinearly interpolate a coarse concentration field at a fine RGB position.
The stored coordinates are the first three concentrations; the fourth is
reconstructed.
"""
@inline function coarse_seed(
        coarse::AbstractVector{T}, cn::Int, r::T, g::T, b::T
    ) where {T<:AbstractFloat}
    d = cn - 1
    if cn == 1
        o = 1
        c1, c2, c3 = coarse[o], coarse[o + 1], coarse[o + 2]
        return (c1, c2, c3, one(T) - c1 - c2 - c3)
    end
    x = clamp(r, zero(T), one(T)) * d
    y = clamp(g, zero(T), one(T)) * d
    z = clamp(b, zero(T), one(T)) * d
    i, fx = _axis_f(x, cn)
    j, fy = _axis_f(y, cn)
    k, fz = _axis_f(z, cn)
    b000 = 3 * (i + cn * (j + cn * k)) + 1
    b100 = b000 + 3
    b010 = b000 + 3cn
    b110 = b010 + 3
    b001 = b000 + 3cn * cn
    b101 = b001 + 3
    b011 = b001 + 3cn
    b111 = b011 + 3
    @inbounds begin
        c1 = _tri(
            coarse[b000], coarse[b100], coarse[b010], coarse[b110],
            coarse[b001], coarse[b101], coarse[b011], coarse[b111], fx, fy, fz,
        )
        c2 = _tri(
            coarse[b000 + 1], coarse[b100 + 1], coarse[b010 + 1], coarse[b110 + 1],
            coarse[b001 + 1], coarse[b101 + 1], coarse[b011 + 1], coarse[b111 + 1], fx, fy, fz,
        )
        c3 = _tri(
            coarse[b000 + 2], coarse[b100 + 2], coarse[b010 + 2], coarse[b110 + 2],
            coarse[b001 + 2], coarse[b101 + 2], coarse[b011 + 2], coarse[b111 + 2], fx, fy, fz,
        )
    end
    return (c1, c2, c3, one(T) - c1 - c2 - c3)
end

@inline function _tri(
        c000, c100, c010, c110, c001, c101, c011, c111, fx, fy, fz,
    )
    a = c000 + fx * (c100 - c000)
    b = c010 + fx * (c110 - c010)
    c = c001 + fx * (c101 - c001)
    d = c011 + fx * (c111 - c011)
    e = a + fy * (b - a)
    f = c + fy * (d - c)
    return e + fz * (f - e)
end

"""
    coarse_to_fine_slab!(dest, model, n, k, scratch, settings, coarse, cn)

Fill fine slab `k` from a coarse concentration field plus a short bulk solve
per vertex.
"""
function coarse_to_fine_slab!(
        dest::AbstractVector{T}, model, n::Int, k::Int, scratch::SolverScratch,
        settings::UnmixSettings{T}, coarse::AbstractVector{T}, cn::Int,
    ) where {T}
    d = n - 1
    z = T(k) / d
    uniform = (T(0.25), T(0.25), T(0.25), T(0.25))
    @inbounds for j in 0:(n - 1)
        y = T(j) / d
        for i in 0:(n - 1)
            x = T(i) / d
            rgb = (x, y, z)
            seed = coarse_seed(coarse, cn, x, y, z)
            r = unmix_bulk!(scratch, model, rgb, (seed, uniform), settings)
            o = _table_offset(n, i, j, k)
            dest[o + 1] = r.c[1]
            dest[o + 2] = r.c[2]
            dest[o + 3] = r.c[3]
        end
    end
    return nothing
end

"""
    generate_inverse_coarse_to_fine(model, n; coarse_n, threads, settings, coarse_settings, resume, on_slab)

Solve the inverse on a `coarse_n^3` grid with the reference solver, then use
the trilinearly interpolated field as the seed for a short bulk solve at
`n^3`.

The coarse solve carries the global structure; the fine solve only polishes,
so it costs a few evaluations per vertex instead of an enumeration. The
coarse field is deterministic and the fine slabs are independent, so the
result does not depend on how slabs are distributed across threads.
"""
function generate_inverse_coarse_to_fine(
        model, n::Integer; coarse_n::Integer = 64,
        threads::Integer = Threads.nthreads(),
        settings::Union{Nothing,UnmixSettings} = nothing,
        coarse_settings::Union{Nothing,UnmixSettings} = nothing,
        resume = nothing, on_slab = nothing,
    )
    T = _scalar_type(model)
    n = Int(n)
    coarse_n = Int(coarse_n)
    2 <= coarse_n < n || throw(ArgumentError(
        "coarse_n must be at least 2 and smaller than n, got $coarse_n and $n"
    ))
    fine_settings = settings === nothing ?
        UnmixSettings{T}(15, T(1.0e-10), T(1.0e-6), 1) : settings
    cs = coarse_settings === nothing ?
        UnmixSettings{T}(50, T(1.0e-10), T(1.0e-6), 4) : coarse_settings
    coarse = generate_inverse(model, coarse_n; threads = threads, settings = cs, solver = :reference)
    out = Vector{T}(undef, 3 * n^3)
    krange = collect(0:(n - 1))
    if threads <= 1 || n < 8
        scratch = SolverScratch()
        for k in krange
            _coarse_slab_or_generate!(
                out, model, n, k, scratch, fine_settings, coarse, coarse_n, resume, on_slab
            )
        end
    else
        Threads.@threads for k in krange
            scratch = SolverScratch()
            _coarse_slab_or_generate!(
                out, model, n, k, scratch, fine_settings, coarse, coarse_n, resume, on_slab
            )
        end
    end
    return out
end

function _coarse_slab_or_generate!(
        out, model, n, k, scratch, settings, coarse, cn, resume, on_slab,
    )
    if resume !== nothing
        loaded = resume(k)
        if loaded !== nothing
            o = 3 * k * n * n
            copyto!(out, o + 1, loaded, 1, 3 * n * n)
            return nothing
        end
    end
    coarse_to_fine_slab!(out, model, n, k, scratch, settings, coarse, cn)
    if on_slab !== nothing
        o = 3 * k * n * n
        on_slab(k, view(out, o + 1:o + 3 * n * n))
    end
    return nothing
end
