"""Experimental dense inverse table with `UInt16` concentrations scaled by 65535."""
struct WordLUT
    n::Int
    data::Vector{UInt16}
    function WordLUT(n::Integer, data::Vector{UInt16})
        _experimental_grid(n)
        length(data) == 3n^3 || throw(ArgumentError("wrong inverse table length"))
        return new(Int(n), data)
    end
end

"""Packed forward RGB lattice indexed by ordered cumulative concentrations.

Vertices `(i,j,k)` satisfy `0 <= i <= j <= k < n`. Each stores three bytes.
Their concentrations are `(i,j-i,k-j,n-1-k)/(n-1)`.
"""
struct SimplexLUT
    n::Int
    data::Vector{UInt8}
    function SimplexLUT(n::Integer, data::Vector{UInt8})
        _experimental_grid(n)
        length(data) == 3 * _simplex_vertices(n) || throw(ArgumentError("wrong simplex table length"))
        return new(Int(n), data)
    end
end
_experimental_grid(n) = 2 <= n <= 1024 ? nothing : throw(ArgumentError("experimental grid must be between 2 and 1024"))
_simplex_vertices(n) = n * (n + 1) * (n + 2) ÷ 6
@inline _simplex_index(i, j, k) = k * (k + 1) * (k + 2) ÷ 6 + j * (j + 1) ÷ 2 + i

"""Opt-in `.pmxe` model with independent inverse precision and forward layout.

Supports the ordinary Julia encoding and mixing API. Production defaults and
compiled ABI payloads continue to use `PigmentModel` and `.pmx`.
"""
struct ExperimentalModel{I <: Union{ByteLUT, WordLUT}, F <: Union{ByteLUT, SimplexLUT}} <: AbstractPigmentModel
    id::NTuple{16, UInt8}
    inverse::I
    forward::F
    function ExperimentalModel(id::NTuple{16, UInt8}, inverse::I, forward::F) where {I <: Union{ByteLUT, WordLUT}, F <: Union{ByteLUT, SimplexLUT}}
        model = new{I, F}(id, inverse, forward)
        validate_model(model)
        return model
    end
end
function Base.show(io::IO, model::ExperimentalModel)
    return print(io, "ExperimentalModel(id = ", model_id(model), ", n = ", grid_n(model), ", inverse = ", eltype(model.inverse.data), ", forward = ", nameof(typeof(model.forward)), ")")
end
