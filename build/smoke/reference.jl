#!/usr/bin/env julia
# Emit reference values for the compiled-ABI smoke tests.
#
#   julia --project=build smoke/reference.jl [payload] [out]
#
# The values are computed by the ordinary Julia package, from the same payload
# the shared library embeds. `client.c` and `client.py` read this file and
# check the compiled entrypoints against it, which is what makes "the C and
# Python clients agree with Julia" a measurable claim rather than a claim
# about the same code being called twice.
#
# Token stream (all whitespace separated):
#
#   PAINTMIX-REF 1
#   <grid_n> <abi_version> <id_lo> <id_hi>
#   mix <count>
#     <dtype> <a_r g b> <b_r g b> <t> <out_r g b>
#   encode <count>
#     <dtype> <rgb r g b> <c1 c2 c3 c4> <r r g b> <out_r g b>
#   wmix <count>
#     <dtype> <n> <colors 3n> <weights n> <out_r g b>
#
# `dtype` is 0 for Float64 and 1 for Float32.

using PaintMix
using Random: Xoshiro

const BUILD_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_PAYLOAD = joinpath(BUILD_DIR, "data", "payload.pmx")

const MIX_COUNT = 200
const ENCODE_COUNT = 100
const WMIX_COUNT = 50

function emit!(io, values...)
    for v in values
        print(io, v, ' ')
    end
    println(io)
    return nothing
end

function emit3!(io, v)
    print(io, v[1], ' ', v[2], ' ', v[3], ' ')
    return nothing
end

function main(args)
    payload = length(args) >= 1 ? args[1] : DEFAULT_PAYLOAD
    out = length(args) >= 2 ? args[2] : joinpath(BUILD_DIR, "out", "reference.txt")
    model = read_model(payload)
    n = grid_n(model)
    id = model.id
    # Print the two identifier halves as signed Int64 so a plain `%lld` /
    # Python `int()` round trips them; comparing the unsigned bits is what
    # matters, and both clients reinterpret back to UInt64.
    lo = reinterpret(Int64, sum(UInt64(id[i]) << (8 * (i - 1)) for i in 1:8))
    hi = reinterpret(Int64, sum(UInt64(id[8 + i]) << (8 * (i - 1)) for i in 1:8))
    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "PAINTMIX-REF 1")
        emit!(io, n, 1, lo, hi)
        rng = Xoshiro(20240716)

        println(io, "mix ", MIX_COUNT)
        for i in 1:MIX_COUNT
            if isodd(i)
                a = (rand(rng), rand(rng), rand(rng))
                b = (rand(rng), rand(rng), rand(rng))
                t = rand(rng)
                print(io, 0, ' ')
                emit3!(io, a)
                emit3!(io, b)
                print(io, t, ' ')
                emit3!(io, mix(model, a, b, t))
                println(io)
            else
                a = (rand(rng, Float32), rand(rng, Float32), rand(rng, Float32))
                b = (rand(rng, Float32), rand(rng, Float32), rand(rng, Float32))
                t = rand(rng, Float32)
                print(io, 1, ' ')
                emit3!(io, a)
                emit3!(io, b)
                print(io, t, ' ')
                emit3!(io, mix(model, a, b, t))
                println(io)
            end
        end

        println(io, "encode ", ENCODE_COUNT)
        for i in 1:ENCODE_COUNT
            if isodd(i)
                x = (rand(rng), rand(rng), rand(rng))
                z = encode(model, x)
                print(io, 0, ' ')
            else
                x = (rand(rng, Float32), rand(rng, Float32), rand(rng, Float32))
                z = encode(model, x)
                print(io, 1, ' ')
            end
            emit3!(io, x)
            print(io, z.c[1], ' ', z.c[2], ' ', z.c[3], ' ', z.c[4], ' ')
            emit3!(io, z.r)
            emit3!(io, decode(model, z))
            println(io)
        end

        println(io, "wmix ", WMIX_COUNT)
        for i in 1:WMIX_COUNT
            count = rand(rng, 1:8)
            dtype32 = iseven(i)
            colors = if dtype32
                [
                    (rand(rng, Float32), rand(rng, Float32), rand(rng, Float32))
                        for _ in 1:count
                ]
            else
                [(rand(rng), rand(rng), rand(rng)) for _ in 1:count]
            end
            weights = dtype32 ? rand(rng, Float32, count) : rand(rng, count)
            print(io, dtype32 ? 1 : 0, ' ', count, ' ')
            for c in colors
                emit3!(io, c)
            end
            for w in weights
                print(io, w, ' ')
            end
            emit3!(io, weighted_mix(model, colors, weights))
            println(io)
        end
    end
    println("wrote $out for model $(model_id(model)) (n = $n)")
    return nothing
end

main(ARGS)
