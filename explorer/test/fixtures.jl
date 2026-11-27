module PlotFixtures

using Explorer
using PaintMix
using PaintMixPrecompute

# Reuse the physical synthetic inputs rather than requiring the workbooks.
include("../../precompute/test/fixtures.jl")

function explorer_data()
    cfg, db, spectra, quad, spectral = synthetic_model(; samples = 16)
    model = Explorer.synthetic_model(; n = 5)
    # An identity forward table cannot distinguish paint from RGB blending.
    forward = generate_forward(spectral, 5; threads = 1)
    bytes = UInt8.(clamp.(round.(255 .* forward), 0, 255))
    model = PaintMix.PigmentModel(
        model.id, model.inverse, PaintMix.ByteLUT(5, bytes), model.format_version, model.flags,
    )
    curves = Explorer.derive_curves(spectra, quad, spectral, spectral, model, db)
    return Explorer.ExplorerData(;
        cfg, db, spectra, quad, spectral, fitted = spectral, model, derived = curves,
    )
end

function traced_data(dir)
    cfg, db, spectra, quad, _ = synthetic_model()
    cfg["grid"]["dev_n"] = 3
    cfg["surrogate"]["surface_divisions"] = 2
    cfg["surrogate"]["max_iterations"] = 2
    cfg["surrogate"]["alpha_halvings"] = 1
    cfg["surrogate"]["push_tolerance"] = 0.0
    cfg["unmix"]["max_iterations"] = 5
    settings = TableSettings(cfg; strategy = :bulk, threads = 1)
    run_pipeline(
        cfg; inputs = (db, spectra, quad), settings, out_dir = dir,
        stage = :fit, trace_fit = true, log = nothing
    )
    result = run_pipeline(
        cfg; inputs = (db, spectra, quad), settings, out_dir = dir,
        stage = :tables, trace_inverse = true, log = nothing
    )
    sidecar = PaintMixPrecompute.TOML.parsefile(result["sidecar_path"])
    model = result["model"]
    notes = String[]
    traces = Explorer.load_traces(sidecar, dir, model, notes; payload_path = result["payload_path"])
    return Explorer.ExplorerData(; cfg, model, sidecar, notes, traces)
end

end
