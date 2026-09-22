# Structural quality checks: concrete struct fields and the conventions that
# Aqua enforces. These read only the package's own definitions, so they need
# no lookup tables.

using Aqua: Aqua
using CheckConcreteStructs: all_concrete

@testset "struct fields are concretely typed" begin
    # The module form also covers types added later. Warnings name the
    # offending field when the assertion fails.
    @test all_concrete(PaintMix)
end

@testset "Aqua" begin
    Aqua.test_all(PaintMix)
end
