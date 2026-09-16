#!/usr/bin/env julia
# Runtime benchmarks.
#
#   julia --project=benchmarks run.jl [--payload PATH] [--n 256] [--samples 100000]
#
# With no payload the benchmarks use a synthetic table of the requested size,
# so the numbers describe the kernel and the cache behavior rather than a
# particular model. Pass `--payload data/default/default.pmx` to measure the
# release payload, including its load time.
#
# A hand-rolled timer is used instead of BenchmarkTools on purpose: the model
# is a large value, and `@benchmark`'s `$` interpolation splices it into the
# expression as a literal (hundreds of thousands of bytes per sample), while
# leaving it uninterpolated makes the macro resolve it as a `Main` global.
# A closure over typed locals measures the same thing without either problem.

using PaintMix
using Printf: @printf
using Random: Xoshiro

function synthetic_model(n::Integer)
    rng = Xoshiro(0x62656e63686d6172)
    m = Int(n)
    payload = 3 * m^3
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
    return PigmentModel(
        ntuple(_ -> rand(rng, UInt8), Val(16)),
        ByteLUT(m, inv), ByteLUT(m, fwd), FORMAT_VERSION,
        FLAG_FORWARD_SIMPLEX_PROJECTED | FLAG_INVERSE_LARGEST_REMAINDER,
    )
end

function parse_args(args)
    opts = Dict{String,Any}("payload" => nothing, "n" => 256, "samples" => 100_000)
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--payload", "--n", "--samples")
            i < length(args) || error("$a needs a value")
            key = a[3:end]
            opts[key] = key == "payload" ? args[i + 1] : parse(Int, args[i + 1])
            i += 1
        else
            error("unknown argument $a")
        end
        i += 1
    end
    return opts
end

"""
    timeit(f; target = 0.05) -> (seconds_per_call, iterations)

Time `f` with a warmup call, growing the iteration count until the batch takes
at least `target` seconds.
"""
function timeit(f; target::Float64 = 0.05)
    f()
    n = 1
    while true
        t0 = time_ns()
        for _ in 1:n
            f()
        end
        elapsed = (time_ns() - t0) / 1e9
        if elapsed >= target || n >= 10^9
            return elapsed / n, n
        end
        n *= 4
    end
end

# `elements` is how many colors one call processes, so a batch kernel reports
# a per-element cost alongside the per-call one.
function report(label, f; elements::Int = 1)
    per_call, n = timeit(f)
    allocs = @allocated f()
    per_element = per_call / elements
    @printf(
        "%-38s %10.4f us/call  %8.1f ns/color  %8.2f Mcolor/s  %5d allocs  (%d iters)\n",
        label, per_call * 1e6, per_element * 1e9, elements / per_call * 1e-6, allocs, n,
    )
    return nothing
end

function run_benchmarks(model::PigmentModel, samples::Int)
    n = grid_n(model)
    rng = Xoshiro(7)
    a = (rand(rng), rand(rng), rand(rng))
    b = (rand(rng), rand(rng), rand(rng))
    dest = zeros(3)
    z = encode(model, a)
    as = rand(rng, 3samples)
    bs = rand(rng, 3samples)
    ts = rand(rng, samples)
    out = zeros(3samples)
    weights = rand(rng, samples)
    wout = zeros(3)

    report("encode (Float64)", () -> encode(model, a))
    report("decode (Float64)", () -> decode(model, z))
    report("mix (Float64)", () -> mix(model, a, b, 0.5))
    report("mix! (Float64, preallocated)", () -> mix!(dest, model, a, b, 0.5))
    report(
        "bulk_mix! ($(samples) pairs)", () -> bulk_mix!(out, model, as, bs, ts);
        elements = samples,
    )
    report(
        "weighted_mix! ($(samples) colors)",
        () -> weighted_mix!(wout, model, as, weights);
        elements = samples,
    )
    return nothing
end

function main(args)
    opts = parse_args(args)
    load_time = 0.0
    model = if opts["payload"] === nothing
        synthetic_model(opts["n"])
    else
        load_time = @elapsed m = read_model(opts["payload"])
        m
    end
    @printf("model: n = %d, id = %s\n", grid_n(model), model_id(model))
    opts["payload"] === nothing ||
        @printf("%-40s %10.4f s\n", "payload load and validation", load_time)
    return run_benchmarks(model, Int(opts["samples"]))
end

main(ARGS)
