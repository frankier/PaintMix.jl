# Precompute tests.
#
# The fast tests run on synthetic inputs and cover the numerical contracts
# the pipeline depends on: the spectral reference, the objective derivatives,
# the simplex tools, the inverse solvers, the table kernels, and
# checkpointing. Tests that need the measurement workbooks skip when the
# files are not present.

using Test
using ForwardDiff
using PaintMix
using PaintMixPrecompute

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include("fixtures.jl")

@testset "PaintMixPrecompute" begin
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

    @testset "synthetic inputs load" begin
        cfg, db, spectra, quad, model = synthetic_model()
        @test size(spectra.K) == (4, length(wavelength_grid(cfg)))
        @test all(>(0), spectra.K)
        @test all(>(0), spectra.S)
        @test length(quad.weights) == length(quad.wavelength)
        @test quad.norm > 0
    end

    @testset "Kubelka-Munk limits" begin
        @test km_reflectance(0.0) == 1.0
        @test km_reflectance(1.0e12) < 1.0e-6
        @test km_reflectance(1.0) ≈ 1 / (2 + sqrt(3))
        # Equation (2) and its cancellation-free form agree.
        q = 0.3
        @test km_reflectance(q) ≈ 1 + q - sqrt(q * q + 2q)
        @test saunderson(1.0, 0.03, 0.65) ≈ (0.97 * 0.35) / 0.35
        @test saunderson(0.0, 0.03, 0.65) == 0.0
    end

    @testset "Kubelka-Munk is invariant to a common K and S scale" begin
        _, _, spectra, quad, model = synthetic_model()
        scaled = deepcopy(spectra)
        scaled.K .= 37.0 .* spectra.K
        scaled.S .= 37.0 .* spectra.S
        other = spectral_model(scaled, quad)
        for c in ((0.25, 0.25, 0.25, 0.25), (0.7, 0.1, 0.1, 0.1), (0.0, 0.5, 0.5, 0.0))
            @test maximum(abs.(collect(mix_rgb(model, c)) .- collect(mix_rgb(other, c)))) < 1.0e-12
        end
    end

    @testset "spectral Jacobian matches finite differences" begin
        _, _, _, _, model = synthetic_model()
        for c in ((0.25, 0.25, 0.25, 0.25), (0.7, 0.1, 0.1, 0.1), (0.0, 0.4, 0.3, 0.3))
            J = zeros(3, 4)
            mix_rgb_jacobian!(J, model, c)
            h = 1.0e-6
            for i in 1:4
                cp = ntuple(j -> j == i ? c[j] + h : c[j], 4)
                cm = ntuple(j -> j == i ? c[j] - h : c[j], 4)
                rp = mix_rgb(model, cp)
                rm = mix_rgb(model, cm)
                for r in 1:3
                    @test J[r, i] ≈ (rp[r] - rm[r]) / (2h) rtol = 1.0e-5 atol = 1.0e-7
                end
            end
        end
    end

    @testset "Oklab conventions" begin
        black = linear_srgb_to_oklab((0.0, 0.0, 0.0))
        white = linear_srgb_to_oklab((1.0, 1.0, 1.0))
        @test maximum(abs, black) < 1.0e-12
        @test white[1] ≈ 1.0 atol = 1.0e-6
        @test abs(white[2]) < 1.0e-6
        @test abs(white[3]) < 1.0e-6
        # Signed cube root keeps out-of-gamut mixtures finite.
        o = linear_srgb_to_oklab((1.2, -0.1, 0.5))
        @test all(isfinite, o)
        @test oklab_distance_squared(white, white) == 0.0
    end

    @testset "cube penalty equals the signed distance squared" begin
        for p in ((0.5, 0.5, 0.5), (1.2, 0.5, 0.5), (-0.1, -0.2, 0.5), (1.1, 1.2, 1.3))
            d = cube_signed_distance(p)
            want = d > 0 ? d^2 : 0.0
            @test cube_outside_penalty(p) ≈ want atol = 1.0e-12
        end
    end

    @testset "surface quadrature" begin
        sq = surface_quadrature(5)
        @test length(sq.points) == length(sq.weights) == length(sq.face)
        @test sum(sq.weights) ≈ 1.0
        for c in sq.points
            @test sum(c) ≈ 1.0 atol = 1.0e-12
            @test all(>=(0), c)
            @test count(==(0.0), c) >= 1
        end
    end

    @testset "softplus round trip and alpha schedule" begin
        for x in (-8.0, -1.0, 0.0, 1.0, 6.0)
            @test inv_softplus(softplus(x)) ≈ x atol = 1.0e-9
        end
        cfg = load_config()
        α = alpha_schedule(cfg)
        @test issorted(α; rev = true)
        @test first(α) == cfg["surrogate"]["alpha_initial"]
        @test last(α) >= cfg["surrogate"]["alpha_final"]
    end

    @testset "projection onto the simplex" begin
        for c in ((0.25, 0.25, 0.25, 0.25), (1.0, 1.0, 1.0, 1.0), (-1.0, 0.5, 0.5, 0.5),
                  (2.0, -1.0, 0.0, 0.0))
            p = project_to_simplex(c)
            @test sum(p) ≈ 1.0 atol = 1.0e-12
            @test all(>=(0), p)
            # Projection is idempotent on the simplex.
            @test maximum(abs.(collect(project_to_simplex(p)) .- collect(p))) < 1.0e-12
        end
        # A point already on the simplex is unchanged.
        c = (0.1, 0.2, 0.3, 0.4)
        @test maximum(abs.(collect(project_to_simplex(c)) .- collect(c))) < 1.0e-12
    end

    @testset "joint simplex quantization" begin
        for c in ((0.25, 0.25, 0.25, 0.25), (1.0, 0.0, 0.0, 0.0), (0.0, 0.0, 0.0, 1.0),
                  (1 / 3, 1 / 3, 1 / 3, 0.0), (0.001, 0.001, 0.001, 0.997))
            b = quantize_simplex(c)
            @test length(b) == 3
            @test all(x -> 0 <= Int(x) <= 255, b)
            @test Int(b[1]) + Int(b[2]) + Int(b[3]) <= 255
            # Deterministic.
            @test quantize_simplex(c) == b
        end
        # Exact concentrations quantize exactly.
        @test quantize_simplex((0.5, 0.25, 0.25, 0.0)) == (UInt8(127), UInt8(64), UInt8(64))
        # Ties break towards the lower index.
        b = quantize_simplex((0.5, 0.5, 0.0, 0.0))
        @test Int(b[1]) + Int(b[2]) == 255
        @test b[1] >= b[2]
    end

    @testset "affine table interpolation and Jacobian" begin
        lut, f = affine_lut(9)
        for c in ((0.3, 0.4, 0.5, 0.0), (0.0, 0.0, 0.0, 1.0), (1.0, 0.0, 0.0, 0.0))
            got = eval_mix(lut, c)
            want = f(c)
            @test maximum(abs, (got[1] - want[1], got[2] - want[2], got[3] - want[3])) < 1.0e-12
        end
        J = zeros(3, 4)
        c = (0.3, 0.4, 0.5, 0.0)
        eval_jacobian4!(J, lut, c)
        h = 1.0e-7
        for i in 1:3
            cp = ntuple(j -> j == i ? c[j] + h : c[j], 4)
            cm = ntuple(j -> j == i ? c[j] - h : c[j], 4)
            rp = eval_mix(lut, cp)
            rm = eval_mix(lut, cm)
            for r in 1:3
                @test J[r, i] ≈ (rp[r] - rm[r]) / (2h) atol = 1.0e-6
            end
        end
        @test all(iszero, J[:, 4])
    end

    @testset "forward table reproduces the spectral model at vertices" begin
        _, _, _, _, model = synthetic_model()
        lut, _ = affine_lut(5)
        n = 5
        fwd = generate_forward(model, n)
        grid = forward_float_lut(n, fwd)
        d = n - 1
        for k in 0:d, j in 0:d, i in 0:d
            c1, c2, c3 = i / d, j / d, k / d
            c4 = 1 - c1 - c2 - c3
            if c4 < 0
                p = project_to_simplex((c1, c2, c3, c4))
                c1, c2, c3, c4 = p[1], p[2], p[3], p[4]
            end
            want = mix_rgb(model, (c1, c2, c3, c4))
            got = eval_mix(grid, (i / d, j / d, k / d, 1 - i / d - j / d - k / d))
            @test maximum(abs, (got[1] - want[1], got[2] - want[2], got[3] - want[3])) < 1.0e-12
        end
    end

    @testset "inverse solver recovers concentrations" begin
        cfg, _, _, _, model = synthetic_model()
        settings = unmix_settings(cfg)
        scratch = SolverScratch()
        for c in ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0),
                  (0.0, 0.0, 0.0, 1.0), (0.25, 0.25, 0.25, 0.25), (0.5, 0.0, 0.5, 0.0),
                  (0.6, 0.2, 0.1, 0.1))
            rgb = mix_rgb(model, c)
            r = unmix_reference(model, rgb; settings = settings, scratch = scratch)
            @test r.sse < 1.0e-16
            @test maximum(abs.(collect(r.c) .- collect(c))) < 1.0e-6
        end
        # Bulk solver agrees with the reference on interior and boundary points.
        for c in ((0.25, 0.25, 0.25, 0.25), (0.5, 0.0, 0.5, 0.0), (0.8, 0.1, 0.05, 0.05))
            rgb = mix_rgb(model, c)
            ref = unmix_reference(model, rgb; settings = settings, scratch = scratch)
            bulk = unmix_bulk!(scratch, model, rgb, (c,), settings)
            @test bulk.sse <= ref.sse + 1.0e-10
        end
    end

    @testset "table generation, quantization, and payload round trip" begin
        cfg, _, _, _, model = synthetic_model()
        n = 9
        forward = generate_forward(model, n)
        inverse = generate_inverse(model, n; threads = 1, solver = :reference)
        ft = FloatTables(n, forward, inverse)
        inv, fwd = quantize_tables(ft, cfg)
        @test inv.n == n && fwd.n == n
        # Every inverse vertex is a valid simplex point.
        for v in 0:(n^3 - 1)
            s = Int(inv.data[3v + 1]) + Int(inv.data[3v + 2]) + Int(inv.data[3v + 3])
            @test s <= 255
        end
        provenance = Dict{String,Any}(
            "config_hash" => repeat("0", 64),
            "inputs" => Dict{String,String}("config" => repeat("0", 64)),
            "grid" => Dict{String,Any}("n" => n),
        )
        m = build_model(ft, cfg, provenance)
        bytes = PaintMix.model_to_bytes(m)
        back = PaintMix.model_from_bytes(bytes)
        @test PaintMix.model_id(back) == PaintMix.model_id(m)
        @test back.inverse.data == m.inverse.data
        @test back.forward.data == m.forward.data

        cont = continuity_report(m, cfg; lines = 16, points = 17)
        @test cont["concentration_jump"]["max"] >= 0
        @test cont["decode_error"]["max"] <= 2.0e-6
        beh = behavior_report(m, model, cfg)
        @test beh["white_tint_monotone"]
        @test beh["same_color_max_error"] <= 2.0e-6
        @test beh["reversal_symmetry_error"] <= 2.0e-6
    end

    @testset "coarse-to-fine inverse matches the reference" begin
        cfg, _, _, _, model = synthetic_model()
        n = 9
        coarse_n = 6
        settings = unmix_settings(cfg)
        fine = generate_inverse_coarse_to_fine(
            model, n; coarse_n = coarse_n, threads = 1,
            settings = UnmixSettings{Float64}(15, 1.0e-10, 1.0e-6, 1),
            coarse_settings = settings,
        )
        # Every vertex is a valid simplex point and decodes to the requested
        # color at least as well as the coarse seeding allows.
        sc = SolverScratch()
        d = n - 1
        for k in 0:d, j in 0:d, i in 0:d
            o = 3 * (i + n * (j + n * k))
            c1, c2, c3 = fine[o + 1], fine[o + 2], fine[o + 3]
            @test c1 >= -1.0e-9 && c2 >= -1.0e-9 && c3 >= -1.0e-9
            @test c1 + c2 + c3 <= 1 + 1.0e-9
            rgb = (i / d, j / d, k / d)
            res = unmix_reference(model, rgb; settings = settings, scratch = sc)
            # The coarse seed is a starting point, not the answer; the fine
            # solve should land in the reference basin even when it does not
            # reach the same local minimum.
            @test sum(abs2, mix_rgb(model, (c1, c2, c3, 1 - c1 - c2 - c3)) .- rgb) <=
                res.sse + 1.0e-3
        end
    end

    @testset "slab checkpointing" begin
        dir = mktempdir()
        store = CheckpointStore(joinpath(dir, "checkpoints"), "abc123", 3, "inverse")
        @test completed_slabs(store) == Int[]
        @test !slab_complete(store, 0)
        write_slab(store, 1, fill(0.5, 3 * 3 * 3))
        @test slab_complete(store, 1)
        @test completed_slabs(store) == [1]
        @test read_slab(store, 1) == fill(0.5, 27)
        @test read_slab(store, 0) === nothing
        # A different job hash in the same directory is rejected.
        @test_throws InputError CheckpointStore(joinpath(dir, "checkpoints"), "other", 3, "inverse")
    end

    @testset "resume reproduces an uninterrupted run" begin
        cfg, _, _, _, model = synthetic_model()
        n = 7
        full = generate_inverse(model, n; threads = 1, solver = :reference)
        dir = mktempdir()
        store = CheckpointStore(joinpath(dir, "cp"), "job", n, "inverse")
        # Complete one slab, then resume.
        write_slab(store, 0, view(full, 1:3n^2))
        resumed = generate_inverse(
            model, n; threads = 1, solver = :reference,
            resume = slab_resume_function(store), on_slab = slab_callback_function(store),
        )
        @test resumed == full
        @test completed_slabs(store) == collect(0:(n - 1))
    end

    @testset "provenance text is deterministic and convention-sensitive" begin
        cfg = load_config()
        inputs = Dict{String,String}("config" => "c", "spectral_k_s" => "k")
        a = provenance_text(cfg, "c", inputs, nothing)
        b = provenance_text(cfg, "c", inputs, nothing)
        @test a == b
        cfg2 = deepcopy(cfg)
        cfg2["quantization"]["rule"] = "other"
        @test provenance_text(cfg2, "c", inputs, nothing) != a
    end

    @testset "acceptance gates read the reports" begin
        cfg = load_config()
        reports = Dict{String,Any}(
            "quality" => Dict{String,Any}(
                "channel_error" => Dict{String,Any}("p99" => 0.001, "max" => 0.002),
                "oklab_error" => Dict{String,Any}("p99" => 0.001, "max" => 0.002),
            ),
            "padding" => Dict{String,Any}(
                "by_depth" => Dict{String,Any}("0.0" => Dict{String,Any}("p99" => 0.001)),
            ),
            "roundtrip" => Dict{String,Any}(
                "float32" => Dict{String,Any}("max" => 1.0e-8),
                "float64" => Dict{String,Any}("max" => 1.0e-14),
            ),
            "quantization" => Dict{String,Any}("invalid_vertices" => 0),
        )
        gates = acceptance_gates(cfg, reports)
        @test all(values(gates))
        reports["quantization"]["invalid_vertices"] = 1
        @test !acceptance_gates(cfg, reports)["quantization_simplex"]
    end
end
