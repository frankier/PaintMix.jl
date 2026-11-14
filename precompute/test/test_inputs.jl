@testset "synthetic inputs load" begin
    cfg, db, spectra, quad, model = synthetic_model()
    @test size(spectra.K) == (4, length(wavelength_grid(cfg)))
    @test all(>(0), spectra.K)
    @test all(>(0), spectra.S)
    @test length(quad.weights) == length(quad.wavelength)
    @test quad.norm > 0
end


@testset "input database round trip" begin
    cfg = synthetic_config()
    db = synthetic_database(cfg)
    path = joinpath(mktempdir(), "roundtrip.duckdb")
    save_database(db, path)
    back = open_database(path)
    by_role = f -> f.role
    # `save_database` records the observer source alongside the other
    # source files, and `open_database` sorts each table.
    expected_files = sort(vcat(db.source_files, [db.observer_source]); by = by_role)
    @test sort(back.source_files; by = by_role) == expected_files
    @test back.pigments == sort(db.pigments; by = p -> p.slot)
    @test sort(back.spectra; by = s -> (s.code, s.quantity, s.wavelength_nm)) ==
        sort(db.spectra; by = s -> (s.code, s.quantity, s.wavelength_nm))
    @test back.saunderson == db.saunderson
    @test back.observer == sort(db.observer; by = o -> o.wavelength_nm)
    @test back.observer_source == db.observer_source
    @test back.build_info == db.build_info
    @test_throws InputError open_database(joinpath(mktempdir(), "absent.duckdb"))
end
