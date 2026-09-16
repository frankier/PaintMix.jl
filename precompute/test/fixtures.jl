# Synthetic inputs for the fast precompute tests.
#
# The tests must run without the measurement workbooks and without the
# DuckDB import, so they build a tiny but structurally valid `InputDatabase`
# on a shortened wavelength grid. Values are chosen to be physically sane
# (positive K and S, white-ish titanium, absorbing chromatic pigments) rather
# than to reproduce the paper.

using PaintMixPrecompute

"""
    synthetic_config(; samples = 6) -> Dict

A copy of the default configuration on a shorter, evenly spaced wavelength
grid. Everything else — pigment order, matrix, Saunderson constants,
quantization rules — is the shipped configuration, so the code paths under
test are the real ones.
"""
function synthetic_config(; samples::Integer = 6)
    cfg = deepcopy(load_config())
    lo, hi = 400.0, 700.0
    step = (hi - lo) / (samples - 1)
    cfg["spectra"]["wavelength_min_nm"] = lo
    cfg["spectra"]["wavelength_max_nm"] = hi
    cfg["spectra"]["wavelength_step_nm"] = step
    return cfg
end

"""
    synthetic_database(cfg; rng_seed = 1) -> InputDatabase

Build a structurally valid input database for `cfg`: four pigments, both
quantities on the configuration grid, a D65-like observer, and the
configuration's Saunderson constants.
"""
function synthetic_database(cfg::AbstractDict; rng_seed::Integer = 1)
    grid = collect(Float64, wavelength_grid(cfg))
    codes = pigment_codes(cfg)
    spectra = SpectrumRecord[]
    # Pigment 4 is a near-ideal white: low absorption, high scattering. The
    # others absorb strongly in one third of the visible range, which gives a
    # wide, non-degenerate gamut.
    for (i, code) in enumerate(codes)
        for (j, wl) in enumerate(grid)
            frac = (j - 1) / max(length(grid) - 1, 1)
            k = i == 4 ? 0.01 + 0.005 * frac : 0.05 + 2.5 * exp(-((frac - (i - 1) / 3)^2) / 0.05)
            s = i == 4 ? 1.0 + 0.1 * frac : 0.2 + 0.1 * frac
            push!(spectra, (
                code = code, quantity = "K", wavelength_nm = wl, value = k,
                source_file = "synthetic", sheet = "synthetic", cell_range = "synthetic",
            ))
            push!(spectra, (
                code = code, quantity = "S", wavelength_nm = wl, value = s,
                source_file = "synthetic", sheet = "synthetic", cell_range = "synthetic",
            ))
        end
    end
    observer = ObserverRecord[]
    for (j, wl) in enumerate(grid)
        x = 1.0 - abs((wl - 450) / 250)
        y = 1.0 - abs((wl - 550) / 250)
        z = 1.0 - abs((wl - 500) / 250)
        push!(observer, (
            wavelength_nm = wl, x_bar = max(x, 0.01), y_bar = max(y, 0.01),
            z_bar = max(z, 0.01), d65 = 1.0,
        ))
    end
    pigments = [
        (
            slot = i, code = codes[i], name = "synthetic $(codes[i])", ci = "CI $(codes[i])",
            column = cfg["pigments"][i]["column"],
        ) for i in 1:4
    ]
    saunderson = (
        k1 = Float64(cfg["saunderson"]["k1"]), k2 = Float64(cfg["saunderson"]["k2"]),
        source_file = "synthetic", sheet = "synthetic", note = "synthetic",
    )
    db = InputDatabase(
        [("synthetic", "synthetic.xlsx", repeat("0", 64))], pigments, spectra, saunderson,
        observer, ("observer", "synthetic.csv", repeat("1", 64)), Dict{String,String}("tool" => "test"),
    )
    return validate_database(db)
end

"""
    synthetic_model(; samples = 6) -> (cfg, db, spectra, quad, model)

Everything a solver test needs, built from [`synthetic_database`](@ref).
"""
function synthetic_model(; samples::Integer = 6)
    cfg = synthetic_config(; samples = samples)
    db = synthetic_database(cfg)
    spectra = load_spectra(cfg, db)
    quad = build_quadrature(cfg, db)
    return cfg, db, spectra, quad, spectral_model(spectra, quad)
end

"""
    affine_lut(n) -> (ForwardFloatLUT, f)

A float forward table holding an exactly trilinear function
`f(c) = a + B*c`, used to check the interpolation and Jacobian kernels
against the analytic answer.
"""
function affine_lut(n::Integer = 9)
    d = n - 1
    a = (0.1, 0.2, 0.3)
    b = ((0.5, 0.1, 0.0), (0.05, 0.4, 0.02), (0.0, 0.1, 0.6))
    data = Vector{Float64}(undef, 3 * n^3)
    for k in 0:(n - 1), j in 0:(n - 1), i in 0:(n - 1)
        c = (i / d, j / d, k / d)
        v = (
            a[1] + b[1][1] * c[1] + b[1][2] * c[2] + b[1][3] * c[3],
            a[2] + b[2][1] * c[1] + b[2][2] * c[2] + b[2][3] * c[3],
            a[3] + b[3][1] * c[1] + b[3][2] * c[2] + b[3][3] * c[3],
        )
        o = 3 * (i + n * (j + n * k))
        data[o + 1] = v[1]
        data[o + 2] = v[2]
        data[o + 3] = v[3]
    end
    f = c -> (
        a[1] + b[1][1] * c[1] + b[1][2] * c[2] + b[1][3] * c[3],
        a[2] + b[2][1] * c[1] + b[2][2] * c[2] + b[2][3] * c[3],
        a[3] + b[3][1] * c[1] + b[3][2] * c[2] + b[3][3] * c[3],
    )
    return forward_float_lut(n, data), f
end
