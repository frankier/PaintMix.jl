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

end
