@testset "default configuration is valid" begin
    cfg = load_config()
    @test pigment_codes(cfg) == ("PB15:4", "PR122", "PY74", "PW6")
    grid = wavelength_grid(cfg)
    @test length(grid) == 38
    @test first(grid) == 380
    @test last(grid) == 750
    @test step(grid) == 10
    @test cfg["color"]["xyz_to_rgb"][2][1] == -0.9689
    @test length(config_hash()) == 64
end


@testset "configuration errors are caught" begin
    good = load_config()
    @test_throws ConfigError validate_config(merge(good, Dict("schema_version" => 2)))

    bad = deepcopy(good)
    bad["spectra"]["wavelength_step_nm"] = 7
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["color"]["storage_transfer"] = "srgb-encoded"
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["saunderson"]["k2"] = 1.5
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    push!(bad["pigments"], deepcopy(bad["pigments"][1]))
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["grid"]["storage"] = "f32"
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["surrogate"]["epsilon"] = 0.0
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["unmix"]["objective"] = "oklab-least-squares"
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["color"]["xyz_to_rgb"] = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    @test_throws ConfigError validate_config(bad)

    bad = deepcopy(good)
    bad["inputs"]["observer_file"] = 42
    @test_throws ConfigError validate_config(bad)

    for (section, key) in (("grid", "interpolation"), ("grid", "index_order"), ("quantization", "rule"), ("quantization", "padding"))
        bad = deepcopy(good)
        bad[section][key] = "unsupported"
        @test_throws ConfigError validate_config(bad)
    end

    @test_throws ConfigError load_config(joinpath(@__DIR__, "no-such-config.toml"))
end


@testset "selected inputs match their recorded checksums" begin
    checksum_file = joinpath(ROOT, "precompute", "inputs", "checksums.sha256")
    @test isfile(checksum_file)
    checked, failures = check_checksums(ROOT, checksum_file)
    mismatches = [f for f in failures if startswith(f, "checksum mismatch")]
    @test isempty(mismatches)
    if checked < 4
        @info "some input files not present; checksum verification is partial" checked
    end
end
