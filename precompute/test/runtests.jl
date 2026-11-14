# Precompute tests.
#
# The fast tests run on synthetic inputs and cover the numerical contracts
# the pipeline depends on: the spectral reference, the objective derivatives,
# the simplex tools, the inverse solvers, the table kernels, and
# checkpointing. Tests that need the measurement workbooks skip when the
# files are not present.

using Test
using ForwardDiff
using LeastSquaresOptim
using PaintMix
using PaintMixPrecompute
using StaticArrays: SVector

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include("fixtures.jl")

# The evaluator interface is exercised through these aliases so the tests do
# not reach into the module's internals.
const eval_mix = PaintMixPrecompute._mix
const eval_jacobian4! = PaintMixPrecompute._jacobian4!

@testset "PaintMixPrecompute" begin
    include("test_config.jl")
    include("test_inputs.jl")
    include("test_spectral.jl")
    include("test_surrogates.jl")
    include("test_tables.jl")
    include("test_inverse.jl")
    include("test_validation.jl")
    include("test_experiments.jl")
    include("test_precision_cli.jl")
end
