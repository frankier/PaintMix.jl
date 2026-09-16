# Concrete C ABI entry points.
#
# Everything here is deliberately simple: fixed-width scalars, borrowed
# caller-owned buffers, and `JLWStatus` returns. The kernels validate their
# arguments first and then call the same `PaintMix` functions the Julia API
# calls, so a C client and a Julia client cannot drift apart.
#
# See `build/README.md` for the contract this file implements.

module PaintMixLib

using JLWInterop
using PaintMix

const ABI_VERSION = Int32(1)
const CHANNELS = 3
const LATENT_SCALARS = 7

# The payload is read while the image is built and becomes part of it.
#
# Its bytes were simplex-validated by `compile.jl` in an optimized process,
# so this load skips that check on purpose: the image-building and
# ABI-probing subprocesses run without optimization (`--compile=min`), where
# re-verifying 100 MB costs minutes of interpreted byte loops. Structural
# checks that are O(1) — magic, version, dimensions, offsets, storage
# conventions — still run.
const PAYLOAD_PATH = joinpath(@__DIR__, "data", "payload.pmx")
const MODEL = PaintMix.model_from_bytes(read(PAYLOAD_PATH); validate = false)

"Compiled-in ABI version. Bump when a signature changes."
Base.@ccallable function paintmix_abi_version()::Int32
    return ABI_VERSION
end

"""
    ModelInfo

Provenance of the compiled-in table pair, plus a status field so a query can
fail in band like every other entrypoint.
"""
struct ModelInfo
    status::JLWStatus
    abi_version::Int32
    format_version::Int32
    grid_n::Int32
    channel_count::Int32
    flags::UInt32
    model_id_lo::UInt64
    model_id_hi::UInt64
end

"""
    paintmix_model_info() -> ModelInfo

Grid size, generation flags, and the 16-byte model
identifier (split into two little-endian `UInt64` halves) of the embedded
tables.
"""
Base.@ccallable function paintmix_model_info()::ModelInfo
    id = MODEL.id
    lo = UInt64(0)
    hi = UInt64(0)
    for i in 1:8
        lo |= UInt64(id[i]) << (8 * (i - 1))
        hi |= UInt64(id[8 + i]) << (8 * (i - 1))
    end
    return ModelInfo(
        jlw_ok(), ABI_VERSION, Int32(MODEL.format_version), Int32(PaintMix.grid_n(MODEL)),
        Int32(CHANNELS), MODEL.flags, lo, hi,
    )
end

# --- argument checks -------------------------------------------------------

# Pointer comparison: `Ptr{T} === C_NULL` is false for any `T` other than
# `Nothing`, because `===` also compares the pointee type. `==` compares the
# address, which is what the contract needs.
@inline function _check_vec(v::CVector{:borrowed,T}, len::Int)::Int32 where {T}
    v.data == C_NULL && return PM_ERR_NULL
    v.dims[1] >= len || return PM_ERR_LENGTH
    return PM_OK
end

@inline function _status(code::Int32)::JLWStatus
    return code == PM_OK ? jlw_ok() : jlw_error(code, error_message(code))
end

@inline function _load3(v::CVector{:borrowed,T}, i::Int) where {T}
    return (
        unsafe_load(v.data, 3i - 2), unsafe_load(v.data, 3i - 1), unsafe_load(v.data, 3i),
    )
end

@inline function _store3!(v::CVector{:borrowed,T}, i::Int, rgb::NTuple{3,T}) where {T}
    unsafe_store!(v.data, rgb[1], 3i - 2)
    unsafe_store!(v.data, rgb[2], 3i - 1)
    unsafe_store!(v.data, rgb[3], 3i)
    return nothing
end

@inline function _finite3(c::NTuple{3,T}) where {T}
    return isfinite(c[1]) && isfinite(c[2]) && isfinite(c[3])
end

# --- encode and decode -----------------------------------------------------

function _encode!(
        rgb::CVector{:borrowed,T}, latent::CVector{:borrowed,T}
    ) where {T<:AbstractFloat}
    st = _check_vec(rgb, CHANNELS)
    st == PM_OK || return st
    st = _check_vec(latent, LATENT_SCALARS)
    st == PM_OK || return st
    c = (unsafe_load(rgb.data, 1), unsafe_load(rgb.data, 2), unsafe_load(rgb.data, 3))
    _finite3(c) || return PM_ERR_NONFINITE
    z = PaintMix.encode(MODEL, c)
    unsafe_store!(latent.data, z.c[1], 1)
    unsafe_store!(latent.data, z.c[2], 2)
    unsafe_store!(latent.data, z.c[3], 3)
    unsafe_store!(latent.data, z.c[4], 4)
    unsafe_store!(latent.data, z.r[1], 5)
    unsafe_store!(latent.data, z.r[2], 6)
    unsafe_store!(latent.data, z.r[3], 7)
    return PM_OK
end

function _decode!(
        latent::CVector{:borrowed,T}, rgb::CVector{:borrowed,T}
    ) where {T<:AbstractFloat}
    st = _check_vec(latent, LATENT_SCALARS)
    st == PM_OK || return st
    st = _check_vec(rgb, CHANNELS)
    st == PM_OK || return st
    for i in 1:LATENT_SCALARS
        isfinite(unsafe_load(latent.data, i)) || return PM_ERR_NONFINITE
    end
    z = Latent(
        (
            unsafe_load(latent.data, 1), unsafe_load(latent.data, 2),
            unsafe_load(latent.data, 3), unsafe_load(latent.data, 4),
        ),
        (
            unsafe_load(latent.data, 5), unsafe_load(latent.data, 6),
            unsafe_load(latent.data, 7),
        ),
    )
    _store3!(rgb, 1, PaintMix.decode(MODEL, z))
    return PM_OK
end

# --- mixing ----------------------------------------------------------------

function _mix!(
        a::CVector{:borrowed,T}, b::CVector{:borrowed,T}, t::T,
        out::CVector{:borrowed,T}
    ) where {T<:AbstractFloat}
    st = _check_vec(a, CHANNELS)
    st == PM_OK || return st
    st = _check_vec(b, CHANNELS)
    st == PM_OK || return st
    st = _check_vec(out, CHANNELS)
    st == PM_OK || return st
    ca = (unsafe_load(a.data, 1), unsafe_load(a.data, 2), unsafe_load(a.data, 3))
    cb = (unsafe_load(b.data, 1), unsafe_load(b.data, 2), unsafe_load(b.data, 3))
    (_finite3(ca) && _finite3(cb) && isfinite(t)) || return PM_ERR_NONFINITE
    _store3!(out, 1, mix(MODEL, ca, cb, t))
    return PM_OK
end

function _bulk_mix!(
        a::CVector{:borrowed,T}, b::CVector{:borrowed,T}, t::CVector{:borrowed,T},
        out::CVector{:borrowed,T}, count::Int64
    ) where {T<:AbstractFloat}
    count < 0 && return PM_ERR_LENGTH
    count == 0 && return PM_OK
    n = Int(count)
    st = _check_vec(a, 3n)
    st == PM_OK || return st
    st = _check_vec(b, 3n)
    st == PM_OK || return st
    st = _check_vec(t, n)
    st == PM_OK || return st
    st = _check_vec(out, 3n)
    st == PM_OK || return st
    return PaintMix.bulk_mix_kernel!(
        _flat(out, 3n), MODEL, _flat(a, 3n), _flat(b, 3n), _flat(t, n), n
    )
end

function _weighted_mix!(
        colors::CVector{:borrowed,T}, weights::CVector{:borrowed,T},
        out::CVector{:borrowed,T}, count::Int64
    ) where {T<:AbstractFloat}
    count < 0 && return PM_ERR_LENGTH
    count == 0 && return PM_ERR_TOTAL
    n = Int(count)
    st = _check_vec(colors, 3n)
    st == PM_OK || return st
    st = _check_vec(weights, n)
    st == PM_OK || return st
    st = _check_vec(out, CHANNELS)
    st == PM_OK || return st
    return PaintMix.weighted_mix_kernel!(
        _flat(out, CHANNELS), MODEL, _flat(colors, 3n), _flat(weights, n), n
    )
end

# A `CVector` is already a dense linear view; this exists to give the kernels
# a plain `CArray` with an explicit length so their own length checks are the
# only ones in play.
@inline function _flat(v::CVector{:borrowed,T}, len::Int) where {T}
    return CVector{:borrowed,T}((len,), v.data)
end

# --- exported entrypoints --------------------------------------------------

# Encode a linear-light sRGB color into seven scalars: c1, c2, c3, c4, r, g, b.
Base.@ccallable function paintmix_encode_f64(
        rgb::CVector{:borrowed,Float64}, latent::CVector{:borrowed,Float64}
    )::JLWStatus
    return _status(_encode!(rgb, latent))
end

Base.@ccallable function paintmix_encode_f32(
        rgb::CVector{:borrowed,Float32}, latent::CVector{:borrowed,Float32}
    )::JLWStatus
    return _status(_encode!(rgb, latent))
end

Base.@ccallable function paintmix_decode_f64(
        latent::CVector{:borrowed,Float64}, rgb::CVector{:borrowed,Float64}
    )::JLWStatus
    return _status(_decode!(latent, rgb))
end

Base.@ccallable function paintmix_decode_f32(
        latent::CVector{:borrowed,Float32}, rgb::CVector{:borrowed,Float32}
    )::JLWStatus
    return _status(_decode!(latent, rgb))
end

Base.@ccallable function paintmix_mix_f64(
        a::CVector{:borrowed,Float64}, b::CVector{:borrowed,Float64}, t::Float64,
        out::CVector{:borrowed,Float64}
    )::JLWStatus
    return _status(_mix!(a, b, t, out))
end

Base.@ccallable function paintmix_mix_f32(
        a::CVector{:borrowed,Float32}, b::CVector{:borrowed,Float32}, t::Float32,
        out::CVector{:borrowed,Float32}
    )::JLWStatus
    return _status(_mix!(a, b, t, out))
end

Base.@ccallable function paintmix_bulk_mix_f64(
        a::CVector{:borrowed,Float64}, b::CVector{:borrowed,Float64},
        t::CVector{:borrowed,Float64}, out::CVector{:borrowed,Float64}, count::Int64
    )::JLWStatus
    return _status(_bulk_mix!(a, b, t, out, count))
end

Base.@ccallable function paintmix_bulk_mix_f32(
        a::CVector{:borrowed,Float32}, b::CVector{:borrowed,Float32},
        t::CVector{:borrowed,Float32}, out::CVector{:borrowed,Float32}, count::Int64
    )::JLWStatus
    return _status(_bulk_mix!(a, b, t, out, count))
end

Base.@ccallable function paintmix_weighted_mix_f64(
        colors::CVector{:borrowed,Float64}, weights::CVector{:borrowed,Float64},
        out::CVector{:borrowed,Float64}, count::Int64
    )::JLWStatus
    return _status(_weighted_mix!(colors, weights, out, count))
end

Base.@ccallable function paintmix_weighted_mix_f32(
        colors::CVector{:borrowed,Float32}, weights::CVector{:borrowed,Float32},
        out::CVector{:borrowed,Float32}, count::Int64
    )::JLWStatus
    return _status(_weighted_mix!(colors, weights, out, count))
end

end # module PaintMixLib
