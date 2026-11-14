"""Seven coordinates per wavelength, fixing mean scattering to one.

The positivity floor is applied to measured inputs before normalization.
Fitted K uses positive softplus values; S uses three ratios to a fourth
reference pigment. A common K/S scale has no effect on reflectance.
"""
struct GaugeParameters
    epsilon::Float64
end
parameterization(::Val{:gauge}, epsilon) = GaugeParameters(epsilon)
_parameter_bounds(::GaugeParameters, x) = (;)
_gradient_floor(::GaugeParameters) = -Inf

function initial_parameters(p::GaugeParameters, spectra)
    K, S = parameter_matrices(SoftplusParameters(p.epsilon), initial_theta(spectra, p.epsilon), length(spectra.wavelength))
    W = size(K, 2)
    x = Vector{eltype(K)}(undef, 7W)
    for j in 1:W
        scale = sum(view(S, :, j)) / 4
        for i in 1:4
            y = K[i, j] / scale
            x[(i - 1) * W + j] = y + log(-expm1(-y))
        end
        for i in 1:3
            x[4W + (i - 1) * W + j] = log(S[i, j]) - log(S[4, j])
        end
    end
    return x
end

function _parameter_matrices_wl(::GaugeParameters, x::AbstractVector{T}, W) where {T}
    length(x) == 7W || throw(DimensionMismatch("expected 7W gauge coordinates"))
    K, S = Matrix{T}(undef, W, 4), Matrix{T}(undef, W, 4)
    tiny = oftype(x[1], floatmin(Float64))
    for j in 1:W
        logits = SVector(x[4W + j], x[5W + j], x[6W + j], zero(T))
        e = exp.(logits .- maximum(logits)) .+ tiny
        scattering = 4 .* e ./ sum(e)
        for i in 1:4
            K[j, i] = softplus(x[(i - 1) * W + j]) + tiny
            S[j, i] = scattering[i]
        end
    end
    return K, S
end
parameter_matrices(p::GaugeParameters, x, W) = _transpose_coefficients(_parameter_matrices_wl(p, x, W)...)
