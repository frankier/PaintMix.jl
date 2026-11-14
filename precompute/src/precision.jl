"""Quantize one floating inverse with a joint four-component integer sum."""
function quantize_inverse(values, n, ::Type{U}) where {U <: Union{UInt8, UInt16}}
    length(values) == 3n^3 || throw(DimensionMismatch("inverse grid length mismatch"))
    out = Vector{U}(undef, length(values))
    for v in 0:(n^3 - 1)
        o = 3v
        c1, c2, c3 = values[o + 1], values[o + 2], values[o + 3]
        c4 = 1.0 - c1 - c2 - c3
        if c4 < 0
            repaired = project_to_simplex(SVector(c1, c2, c3, 0.0))
            c1, c2, c3 = repaired[1], repaired[2], repaired[3]
        end
        c = SVector(c1, c2, c3, max(1.0 - c1 - c2 - c3, 0.0))
        q = quantize_simplex(c, U)
        out[o + 1], out[o + 2], out[o + 3] = q
    end
    return out
end

"""Generate only valid concentration vertices in packed cumulative order."""
function generate_simplex_forward(model, n; threads = Threads.nthreads())
    PaintMix._experimental_grid(n)
    data = Vector{UInt8}(undef, 3 * PaintMix._simplex_vertices(n))
    body = function (_, k)
        for j in 0:k, i in 0:j
            c = SVector(i, j - i, k - j, n - 1 - k) / (n - 1)
            rgb = mix_rgb(model, c)
            o = 3 * PaintMix._simplex_index(i, j, k)
            for ch in 1:3
                data[o + ch] = round(UInt8, clamp(rgb[ch], 0.0, 1.0) * 255)
            end
        end
        return
    end
    _run_slabs!(() -> nothing, body, n, threads)
    return PaintMix.SimplexLUT(n, data)
end

"""Build an opt-in experimental model, identifying its actual bytes and conventions."""
function experimental_model(inverse::Union{PaintMix.ByteLUT, PaintMix.WordLUT}, forward)
    spec = Dict(
        "format" => "pmxe-1", "grid_n" => inverse.n,
        "inverse_bits" => PaintMix._inverse_bits(inverse), "forward_mode" => PaintMix._forward_mode(forward),
        "inverse_sha256" => bytes2hex(sha256(reinterpret(UInt8, htol.(inverse.data)))),
        "forward_sha256" => bytes2hex(sha256(forward.data)),
    )
    return PaintMix.ExperimentalModel(_model_id_from(canonical_text(spec)), inverse, forward)
end

"""Compare full mixing and separate forward/concentration errors on a common corpus."""
function precision_report(model, spectral, cfg; samples = validation_corpus())
    forward, concentration = Float64[], Float64[]
    U = eltype(model.inverse.data)
    scale = Float64(typemax(U))
    for c in samples.concentrations
        exact = mix_rgb(spectral, c)
        push!(forward, maximum(abs, PaintMix.forward_rgb(model, c) - exact))
        q = quantize_simplex(c, U)
        cc = SVector(Float64(q[1]) / scale, Float64(q[2]) / scale, Float64(q[3]) / scale, (scale - sum(Int, q)) / scale)
        push!(concentration, maximum(abs, mix_rgb(spectral, cc) - exact))
    end
    return Dict(
        "quality" => quality_report(model, spectral, cfg; samples),
        "forward_channel" => _error_stats(forward),
        "concentration_channel" => _error_stats(concentration),
        "roundtrip" => roundtrip_report(model, cfg),
        "quantization" => quantization_report(model, cfg),
        "runtime" => _mixing_probe(model, samples),
    )
end

function _mixing_checksum(model, samples)
    total = SVector(0.0, 0.0, 0.0)
    for (ia, ib, t) in samples.pairs
        total += PaintMix.mix(model, samples.colors[ia], samples.colors[ib], t)
    end
    return total
end
function _mixing_probe(model, samples)
    checksum = _mixing_checksum(model, samples)
    times = [@elapsed(_mixing_checksum(model, samples)) for _ in 1:3]
    allocations = @allocated _mixing_checksum(model, samples)
    return Dict("warm_scalar_ns" => minimum(times) * 1.0e9 / length(samples.pairs), "allocated_bytes" => allocations, "pairs" => length(samples.pairs), "checksum" => collect(checksum))
end
