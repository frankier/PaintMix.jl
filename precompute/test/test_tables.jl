@testset "projection onto the simplex" begin
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(1.0, 1.0, 1.0, 1.0),
            SVector(-1.0, 0.5, 0.5, 0.5), SVector(2.0, -1.0, 0.0, 0.0),
        )
        p = project_to_simplex(c)
        @test sum(p) ≈ 1.0 atol = 1.0e-12
        @test all(>=(0), p)
        # Projection is idempotent on the simplex.
        @test maximum(abs.(collect(project_to_simplex(p)) .- collect(p))) < 1.0e-12
    end
    # A point already on the simplex is unchanged.
    c = SVector(0.1, 0.2, 0.3, 0.4)
    @test maximum(abs.(collect(project_to_simplex(c)) .- collect(c))) < 1.0e-12
end


@testset "joint simplex quantization" begin
    for c in (
            SVector(0.25, 0.25, 0.25, 0.25), SVector(1.0, 0.0, 0.0, 0.0),
            SVector(0.0, 0.0, 0.0, 1.0), SVector(1 / 3, 1 / 3, 1 / 3, 0.0),
            SVector(0.001, 0.001, 0.001, 0.997),
        )
        b = quantize_simplex(c)
        @test length(b) == 3
        @test all(x -> 0 <= Int(x) <= 255, b)
        @test Int(b[1]) + Int(b[2]) + Int(b[3]) <= 255
        # Deterministic.
        @test quantize_simplex(c) == b
    end
    # Exact concentrations quantize exactly.
    @test quantize_simplex(SVector(0.5, 0.25, 0.25, 0.0)) == SVector(0x7f, 0x40, 0x40)
    # Ties break towards the lower index.
    b = quantize_simplex(SVector(0.5, 0.5, 0.0, 0.0))
    @test Int(b[1]) + Int(b[2]) == 255
    @test b[1] >= b[2]
end


@testset "affine table interpolation and Jacobian" begin
    lut, f = affine_lut(9)
    for c in (
            SVector(0.3, 0.4, 0.5, 0.0), SVector(0.0, 0.0, 0.0, 1.0),
            SVector(1.0, 0.0, 0.0, 0.0),
        )
        got = eval_mix(lut, c)
        want = f(c)
        @test maximum(
            abs, (got[1] - want[1], got[2] - want[2], got[3] - want[3])
        ) < 1.0e-12
    end
    J = zeros(3, 4)
    c = SVector(0.3, 0.4, 0.5, 0.0)
    eval_jacobian4!(J, lut, c)
    h = 1.0e-7
    for i in 1:3
        cp = SVector{4}(ntuple(j -> j == i ? c[j] + h : c[j], 4))
        cm = SVector{4}(ntuple(j -> j == i ? c[j] - h : c[j], 4))
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
            p = project_to_simplex(SVector(c1, c2, c3, c4))
            c1, c2, c3, c4 = p[1], p[2], p[3], p[4]
        end
        want = mix_rgb(model, SVector(c1, c2, c3, c4))
        got = eval_mix(grid, SVector(i / d, j / d, k / d, 1 - i / d - j / d - k / d))
        @test maximum(
            abs, (got[1] - want[1], got[2] - want[2], got[3] - want[3])
        ) < 1.0e-12
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
    provenance = Dict{String, Any}(
        "config_hash" => repeat("0", 64),
        "inputs" => Dict{String, String}("config" => repeat("0", 64)),
        "grid" => Dict{String, Any}("n" => n),
    )
    m = build_model(ft, cfg, provenance)
    bytes = PaintMix.model_to_bytes(m)
    back = PaintMix.PigmentModel(bytes)
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
    @test_throws InputError CheckpointStore(
        joinpath(dir, "checkpoints"), "other", 3, "inverse"
    )
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

    # The bulk solver over the threaded driver resumes to the same table.
    # `n8 > 8` forces the threaded branch even on a single-threaded run.
    n8 = 8
    full_bulk = generate_inverse(model, n8; threads = 1, solver = :bulk)
    dir2 = mktempdir()
    store2 = CheckpointStore(joinpath(dir2, "cp"), "job", n8, "inverse")
    write_slab(store2, 2, view(full_bulk, (3 * 2 * n8^2 + 1):(3 * 3 * n8^2)))
    resumed_bulk = generate_inverse(
        model, n8; threads = 2, solver = :bulk,
        resume = slab_resume_function(store2), on_slab = slab_callback_function(store2),
    )
    @test resumed_bulk == full_bulk
    @test completed_slabs(store2) == collect(0:(n8 - 1))
end
