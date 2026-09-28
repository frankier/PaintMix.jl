# Data-assembly tests. These exercise the degradation paths directly, without
# the 96 MiB release payload: a synthetic `.pmx` is written to a temp file and
# pinned with `payload_path`.

using Test
using PaintMix
using Explorer

@testset "data assembly" begin
    path = tempname() * ".pmx"
    PaintMix.write_model(path, Explorer.synthetic_model())
    try
        data = Explorer.load_explorer_data(; with_sidecar = false, payload_path = path)
        @test data.model.id == ntuple(_ -> 0x00, Val(16))
        @test data.sidecar === nothing
        @test data.source == "explicit ($path)"
        @test any(n -> occursin("--no-sidecar", n), data.notes)
        @test data.fitted === nothing
    finally
        rm(path; force = true)
    end
end
