#!/usr/bin/env julia
# Generate a synthetic `.pmx` payload for build and relocation checks.
#
# These payloads are *not* the release model: they exist so the compiled ABI
# can be exercised, and so the "does a full-size table survive trimming and
# relocation" question can be answered before any spectral work happens.
#
#   julia --project=build scripts/make_dummy_payload.jl tiny   out.pmx
#   julia --project=build scripts/make_dummy_payload.jl full   out.pmx
#   julia --project=build scripts/make_dummy_payload.jl grid 64 out.pmx

using PaintMix
using Random: Xoshiro

const PROFILES = Dict(
    "tiny" => 2,
    "small" => 32,
    "full" => 256,
)

# `Xoshiro` fills the 100 MB of a 256^3 payload fast enough, and gives the
# same bytes on every run.

"""
    dummy_model(n; seed = 0x7061696e74) -> PigmentModel

A model with deterministic pseudo-random forward values and an inverse table
whose every vertex satisfies the simplex constraint by construction.
"""
function dummy_model(n::Integer; seed::UInt64 = 0x0000007061696e74)
    m = Int(n)
    payload = 3 * m^3
    rng = Xoshiro(seed)
    fwd = Vector{UInt8}(undef, payload)
    inv = Vector{UInt8}(undef, payload)
    for i in 1:payload
        fwd[i] = rand(rng, UInt8)
    end
    for v in 0:(m^3 - 1)
        a = Int(rand(rng, UInt8))
        b = Int(rand(rng, UInt8)) % (256 - a)
        c = Int(rand(rng, UInt8)) % (256 - a - b)
        base = 3v + 1
        inv[base] = UInt8(a)
        inv[base + 1] = UInt8(b)
        inv[base + 2] = UInt8(c)
    end
    id = ntuple(_ -> rand(rng, UInt8), Val(16))
    flags = PaintMix.FLAG_FORWARD_SIMPLEX_PROJECTED |
        PaintMix.FLAG_INVERSE_LARGEST_REMAINDER
    return PigmentModel(id, ByteLUT(n, inv), ByteLUT(n, fwd), FORMAT_VERSION, flags)
end

function main(args)
    isempty(args) && error("usage: make_dummy_payload.jl <tiny|small|full|grid N> <out.pmx>")
    profile = args[1]
    n = if haskey(PROFILES, profile)
        PROFILES[profile]
    elseif profile == "grid"
        length(args) >= 2 || error("grid profile needs an edge length")
        parse(Int, args[2])
    else
        error("unknown profile \"$profile\"; use tiny, small, full, or grid N")
    end
    out = haskey(PROFILES, profile) ? args[2] : args[3]
    model = dummy_model(n)
    write_model(out, model)
    println(
        "wrote $out: n = $n, $(3 * n^3) bytes per table, ",
        "$(length(model_to_bytes(model))) bytes total, id = $(model_id(model))",
    )
    return nothing
end

main(ARGS)
