@testset "provenance text is deterministic and convention-sensitive" begin
    cfg = load_config()
    inputs = Dict{String, String}("config" => "c", "spectral_k_s" => "k")
    a = provenance_text(cfg, "c", inputs)
    b = provenance_text(cfg, "c", inputs)
    @test a == b
    cfg2 = deepcopy(cfg)
    cfg2["quantization"]["rule"] = "other"
    @test provenance_text(cfg2, "c", inputs) != a
end


@testset "acceptance gates read the reports" begin
    cfg = load_config()
    reports = Dict{String, Any}(
        "quality" => Dict{String, Any}(
            "channel_error" => Dict{String, Any}("p99" => 0.001, "max" => 0.002),
            "oklab_error" => Dict{String, Any}("p99" => 0.001, "max" => 0.002),
        ),
        "padding" => Dict{String, Any}(
            "by_depth" => Dict{String, Any}("0.0" => Dict{String, Any}("p99" => 0.001)),
        ),
        "roundtrip" => Dict{String, Any}(
            "float32" => Dict{String, Any}("max" => 1.0e-8),
            "float64" => Dict{String, Any}("max" => 1.0e-14),
        ),
        "quantization" => Dict{String, Any}("invalid_vertices" => 0),
    )
    gates = acceptance_gates(cfg, reports)
    @test all(values(gates))
    reports["quantization"]["invalid_vertices"] = 1
    @test !acceptance_gates(cfg, reports)["quantization_simplex"]
end
