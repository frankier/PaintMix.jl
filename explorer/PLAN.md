# explorer — interactive viewer plan

An interactive, server-rendered viewer for the PaintMix input database, the
precompute intermediates, and the runtime payload. It answers "what is in the
data, and what does each processing step do to it?" by exposing the measured
spectra, the fitted surrogate, the two lookup tables, the color spaces, and a
small painting demo on one local web app.

This is a **development tool**, not part of the runtime package or its
dependency graph. It is the only place in the repository allowed to depend on
plotting, a web server, and a browser toolkit.

## Decisions

| Question | Choice |
| --- | --- |
| Deliverable | This plan only, in this pass. No code yet. |
| Data layer | Depend on `../precompute` (`PaintMixPrecompute`) and reuse `open_database`, `PigmentSpectra`, `Quadrature`, `SpectralModel`, `mix_rgb`, `km_reflectance`, `saunderson`, `theta_parameters`, `linear_srgb_to_oklab`. Do not duplicate colorimetry. |
| Payload | Load the promoted release model `data/default/default.pmx` (256³, 96 MiB). Degrade to a precompute dev payload, then to a synthetic model, rather than failing. |
| Painting demo | A server-side WGLMakie brush canvas **and** a two/three-pigment mixer ramp. |
| Host stack | The `Oxygen` **release** plus upstream `Bonnie.jl` pinned by URL. No local fork, no HTTP-2 branch. |
| Browser tests | [playwright-julia](https://github.com/frankier/playwright-julia), opt-in via `EXPLORER_E2E=1`, mirroring Bonnie's suite. |

## Verified numbers

Computed from the real inputs and payload while writing this plan
(`julia --project=precompute`), so the plan can be checked against the
viewer:

| Pigment | Raw K/S chromaticity (xy) | Spectral pure `c = eᵢ` (linear sRGB) | Runtime forward table pure vertex | Runtime chromaticity (xy) |
| --- | --- | --- | --- | --- |
| PB15:4 | (0.1972, 0.1186) | (0.0115, 0.0059, 0.0706) | (0.0118, 0.0039, 0.0706) | (0.1954, 0.1070) |
| PR122 | (0.5292, 0.2731) | (0.2048, 0.0013, 0.0313) | (0.2039, 0.0000, 0.0275) | (0.5416, 0.2757) |
| PY74 | (0.4747, 0.4929) | (1.0389, 0.6332, **−0.0542**) | (1.0000, 0.6510, 0.0039) | (0.4533, 0.4761) |
| PW6 | (0.3124, 0.3301) | (0.9450, 0.9578, 0.9483) | (0.9529, 0.9529, 0.9059) | (0.3160, 0.3344) |

sRGB primaries from the pinned `xyz_to_rgb` matrix: R (0.6401, 0.3300),
G (0.3000, 0.6000), B (0.1500, 0.0600), D65 (0.3127, 0.3290). Nominal paper
colors (linear and encoded): blue `(0.02, 0.09, 0.42)` → `#2755ad`,
magenta `(0.55, 0.02, 0.12)` → `#c42761`, yellow `(0.71, 0.62, 0.02)` →
`#dbce27`, white `(1, 1, 1)`.

Three facts the viewer should make visible:

1. PY74 is outside the sRGB cube in linear light (blue channel is negative).
   The residual in `Latent` is what lets the model reproduce it; the 8-bit
   forward table clips it. Plot both.
2. The runtime pure-pigment vertices are the **surrogate** pigments, so they
   differ from the raw K/S chromaticity. That difference is the whole point
   of the fit; show the three layers side by side.
3. The 256³ inverse table has a documented quality tail near pure pigments.
   Show it directly (basin map, residual field, LUT-vs-spectral error), not
   just as the aggregate table in the sidecar.

## Architecture

Standalone project, **not** a member of the root workspace. The tracked root
`Manifest.toml` stays free of Makie, Oxygen, and Bonito. `explorer/` depends
one way: `explorer -> precompute -> PaintMix`, plus a URL pin on
`Bonnie.jl` and the registered `Oxygen` release.

```
explorer/
  PLAN.md               this file
  Project.toml          deps + [sources] paths
  run.jl                CLI entry point: load data, build routes, serve
  README.md             added with the code (run instructions)
  src/
    Explorer.jl         module, exports, includes
    data.jl             ExplorerData assembly, caching, degradation
    colorimetry.jl      locus, chromaticity, gamut, derived spectra
    palette.jl          swatches, ramps, nominal colors, 24-pigment context
    templates.jl        OteraEngine Template loading and render helpers
    routes.jl           Oxygen routes and JSON endpoints
    apps/
      spectra.jl        spectra figure
      cie.jl            CIE 1931 figure
      palette.jl        sRGB / palette figure
      paint.jl          painting figure
      tables.jl         table-slice and validation figures
      fit.jl            surrogate-fit figure
  templates/
    base.html
    index.html  spectra.html  cie.html  palette.html
    paint.html  tables.html   fit.html  provenance.html
  static/explorer.css
  test/
    runtests.jl  test_colorimetry.jl  test_server.jl
```

`Project.toml` (deps by UUID, sources by path):

```toml
[deps]
PaintMix            = "baf22a60-1c55-4971-a207-55391bb51cf2"
PaintMixPrecompute  = "eeffd428-b4c8-40bd-84d3-0fe91856487d"
Bonnie              = "6f62f684-a5f8-44b1-bc4e-a9317abb7de4"
Bonito              = "824d6782-a2ef-11e9-3a09-e5662e0c26f8"
Oxygen              = "df9a0d86-3283-4920-82dc-4555fc0d1d8b"
OteraEngine         = "b2d7f28f-acd6-4007-8b26-bc27716e5513"
WGLMakie            = "276b4fcb-3e11-5398-bf8b-a0c2d153d008"
HTTP                = "cd3eb016-35fb-5094-929b-558a96fad6f3"
StaticArrays        = "90137ffa-7385-5640-81b9-e52037218182"
TOML                = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
Printf              = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[sources]
PaintMix           = { path = ".." }
# PaintMixPrecompute itself sources PaintMix; declare it here too so the
# resolver sees the path regardless of which project is active.
PaintMixPrecompute = { path = "../precompute" }
# Upstream Bonnie, which targets a released Oxygen. Do not use a sibling
# checkout and do not add an Oxygen [sources] entry: Oxygen resolves from the
# registry. Pin a commit SHA once Phase 0 confirms the release resolves.
Bonnie             = { url = "https://github.com/frankier/Bonnie.jl.git" }
# Test-only; see [extras]/[targets] below. Also unregistered.
Playwright         = { url = "https://github.com/frankier/playwright-julia" }
```

Test-only dependencies stay out of the app environment:

```toml
[extras]
Test       = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
Playwright = "d777c63b-3b21-49c6-865f-99d0e068a423"

[targets]
test = ["Test", "Playwright"]
```

Add `/explorer/Manifest.toml` to the root `.gitignore` (it is a generated
dev environment, like `/test/Manifest.toml`). Do **not** add `explorer` to
the root `[workspace] projects`.

Request flow:

```
browser ──GET /cie──► BONNIE.middleware ──► Oxygen handler
                                            │  builds WGLMakie App
                                            │  head_content()  (root session bootstrap)
                                            │  app_html(app)   (fragment, subsession)
                                            ▼
                                       OteraEngine template ──► HTML
browser ──WS /bonito/ws/<id>──► Bonnie registry ──► Bonito session ──► figure
```

`Bonnie.setup!(Val(:oxygen))` registers `/bonito/assets/{key}` and
`/bonito/ws/{session_id}` on the Oxygen module and returns
`(; middleware, context, router, prefix, app)`. Install the middleware with
`serve(; middleware = [handle.middleware])`. Call `WGLMakie.activate!()`
before building any figure.

## Data layer

One `ExplorerData` value is assembled at startup and cached for the process.

```julia
struct ExplorerData
    cfg::Dict{String, Any}            # precompute/config/default.toml
    db::InputDatabase                 # precompute/output/paintmix.duckdb
    spectra::PigmentSpectra{Float64}  # measured K and S, 4 x 38
    quad::Quadrature{Float64}         # weights, CIE 1931, D65, xyz_to_rgb
    spectral::SpectralModel{Float64}  # measured, wavelength-major
    fitted::Union{Nothing, SpectralModel}   # from sidecar theta
    model::PaintMix.PigmentModel      # release payload
    source::String                    # which payload path was used
    sidecar::Union{Nothing, Dict}     # release TOML, or nothing
    derived::DerivedCurves            # precomputed per-pigment curves/points
    notes::Vector{String}             # degradation and data warnings
end
```

Sources, in order:

* `precompute/config/default.toml` via `load_config()`. The config's paths
  are relative to the repository root, so the explorer resolves them the way
  `precompute/scripts/generate.jl` does: `ROOT = normpath(joinpath(@__DIR__, "..", ".."))`,
  then `joinpath(ROOT, cfg["inputs"]["database"])`. Never rely on the
  process's working directory.
* `precompute/output/paintmix.duckdb` via `open_database()`; the observer and
  Saunderson rows this holds are the colorimetric ground truth.
* `PigmentSpectra(cfg, db)` and `Quadrature(cfg, db)`, then
  `SpectralModel(spectra, quad)`.
* Payload, first found wins:
  `data/default/default.pmx` (release) →
  `precompute/output/release/default.pmx` →
  `precompute/output/dev/default.pmx` →
  the build's dummy payload →
  synthetic in-memory tables. Record the choice in `source` and show it in
  the provenance banner. Load with `PaintMix.read_model`; an untracked
  96 MiB payload is normal, so absence is not an error.
* `precompute/output/release/<model-id>.toml` if present. `precompute/output/`
  is gitignored, so the sidecar is optional. When it is absent: no fitted
  curves, no fit history, no validation tables; the viewer shows a "not
  available" note instead of failing.

Derived per-pigment curves, all computed once in `colorimetry.jl`:

| Curve | Formula |
| --- | --- |
| Absorption `K` | measured, `PigmentSpectra.K` |
| Scattering `S` | measured, `PigmentSpectra.S` |
| Ratio `K/S` | `K ./ S` |
| Reflectance `R∞` | `km_reflectance.(K ./ S)` |
| Corrected `R′` | `saunderson.(R∞, k1, k2)` |
| Tristimulus integrand | `w .* d65 .* R′ .* cmf` for each of `x̄, ȳ, z̄` |
| Fit curves | `K, S = theta_parameters(sidecar["surrogate"]["theta"], 38, 1e-6)`, then `SpectralModel(spectral; K, S)` |

Chromaticities in `cie.jl`:

* **Spectral locus.** For each grid wavelength, `XYZ ∝ (x̄, ȳ, z̄)`, so
  `xy = (x̄, ȳ) / (x̄ + ȳ + z̄)`. Close the locus with the line of purples.
* **sRGB gamut.** Columns of `inv(quad.xyz_to_rgb)`, normalized per column.
* **Pigment points.** Integrate `R′·D65` against the observer with the
  trapezoidal weights to get `XYZ`, then `xy`. This equals the spectral pure
  `mix_rgb(m, eᵢ)` chromaticity, which is a useful self-check. Also compute
  the runtime forward-table `PaintMix.forward_rgb(model, eᵢ)` chromaticity.
* **Mixture gamut.** Convex hull of the chromaticities of a dense simplex
  sample (`SurfaceQuadrature(d)` points), clipped to the sRGB triangle for
  display. Optionally also the hull of the runtime model.

Two conversion helpers, used everywhere, so linear and encoded are never
mixed up:

```julia
encoded(c) = PaintMix.srgb8_from_linear(SVector{3, Float32}(c...))
linear(b)  = PaintMix.linear_from_srgb8(b)

# Chromaticity of a linear-light sRGB triple: map back to XYZ, then
# sum-normalize. Uses the pinned matrix, so it is the inverse of the same
# transform `mix_rgb` applies.
function xy_of_linear(c)
    xyz = inv(quad.xyz_to_rgb) * collect(Float64, c)
    s = sum(xyz)
    return (xyz[1] / s, xyz[2] / s)
end
```

## Web layer

### Routes (`Oxygen`)

| Route | Contents |
| --- | --- |
| `GET /` | overview: provenance banner, data-card, links |
| `GET /spectra` | spectra figure + readout table |
| `GET /cie` | CIE 1931 figure + gamut/legend |
| `GET /palette` | sRGB/palette figure + ramps + swatch table |
| `GET /paint` | painting figure + mixer ramp |
| `GET /tables` | table slices + validation reports |
| `GET /fit` | surrogate-fit figure + history |
| `GET /provenance` | full metadata tables |
| `GET /healthz` | JSON: payload id, source, sidecar present, notes |
| `GET /fig/<name>` | standalone `figure_page` of one app (debug/iframe) |

Templates are constructed once at startup and called inside the handler:

```julia
const TPL = Dict(n => Template(joinpath(TEMPLATES, "$n.html")) for n in PAGE_NAMES)

@get "/cie" function (req::HTTP.Request)
    d = EXPLORER_DATA[]
    app = cie_app(d)                 # a Bonito.App
    head = head_content()            # emits root-session bootstrap, once per page
    body = app_html(app)             # fragment; becomes a subsession
    Oxygen.html(TPL["cie"](; init = Dict{Symbol, Any}(
        :title => "CIE 1931", :head => head, :body => body,
        :model_id => PaintMix.model_id(d.model), :source => d.source,
    )))
end
```

`head_content()` must be called **before** `app_html()` so the fragment does
not carry the bootstrap twice. Both run inside the request scope that
`BONNIE.middleware` installs. Pass every fragment through the `safe` filter —
`head_content`, `app_html`, and any pre-rendered HTML table must not be
escaped by OteraEngine's default autoescape:

```html
<head>{{ head|safe }}</head>
...
<main>{{ body|safe }}</main>
```

Several figures may share one page: each `app_html` call adds a subsession
of the page's root session over the single websocket. Prefer one `App` per
figure and let the template lay them out with CSS.

`base.html` holds the nav and the shared banner; page templates use
`{% extends %}` and `{% block %}`. Readout tables are rendered server-side
into HTML strings (with `to_html`) and inserted with `|safe`, so no
client-side templating is needed for non-figure content.

### Session lifecycle

`setup!(Val(:oxygen); session_ttl = ..., reconnect_window = ...)` defaults are
fine for a dev tool. Close the server with Ctrl-C; sessions are swept by the
registry. Add a `--prefix` flag only if the viewer is put behind a reverse
proxy; do not hard-code one.

## Pages and figures

### `/spectra` — measured, fitted, and derived spectra

One figure, four axes or one shared axis plus a mode toggle:

```
┌───────────────────────────────────────────────┬──────────────────┐
│ K (log) / S (log) vs λ    ─ measured          │ toggle: pigment  │
│ K/S (log)                 ─ fitted (dashed)   │ toggle: quantity │
│ R∞, R′ (0..1)             ▏λ crosshair        │ λ = 550 nm ▲▼    │
│                                               │ readout table    │
└───────────────────────────────────────────────┴──────────────────┘
```

* Pigment toggles (4), quantity checkboxes (`K`, `S`, `K/S`, `R∞`, `R′`),
  and a "measured vs fitted" toggle when the sidecar is present.
* A crosshair at a selected wavelength (slider or mouse hover); the readout
  table shows the value at that wavelength for each enabled pigment and the
  source cell reference (`SpectrumRecord.cell_range`, e.g. `Q6`).
* `K` and `S` range over six orders of magnitude, so plot them on a log axis
  by default, with a linear toggle.
* The Saunderson constants `k1`, `k2` and the grid (380–750 nm, 10 nm) are
  shown in a caption, and the `kins` contradiction (`inputs/README.md`) is
  surfaced as a warning tooltip.

### `/cie` — where the pigments and the gamut sit

The chromaticity diagram. This is the page the request calls out, and it is
the natural home for the "spectra around the outside" idea.

```
                        λ=450   λ=500
                    ╭──────────────────╮
       λ=400 ──────╯       ·  R (0.64,0.33)
              \         ╱ sRGB triangle
               \  PW6 ●───────● G (0.30,0.60)
        PR122 ●  \    ╱   D65 ✛
                 \ ● PY74
        PB15:4 ●─╲
                    ╰──────╯   B (0.15,0.06)
                        λ=650   λ=700
```

* **Wavelength crown.** The spectral locus is drawn as a line, closed by the
  line of purples. At every grid sample, draw a short radial tick outward
  from the locus, filled with that wavelength's encoded sRGB color, and label
  every other sample (380–750 nm). Because the locus is a curve of
  monochromatic chromaticities, this turns the outside of the horseshoe into
  a labelled spectrum ring — the "cool" version of the request. Out-of-gamut
  monochromatic colors are desaturated for the tick fill and marked; the
  chromaticity point is exact regardless.
* **sRGB gamut.** The triangle R–G–B closed over the D65 white point, drawn
  under the locus. Note in the caption that chromaticity folds out lightness,
  so the sRGB *chromaticity* gamut is exactly this triangle; the cube itself
  lives on `/palette`.
* **Pigment points**, each with three markers and a hover readout:
  raw K/S chromaticity (filled), fitted surrogate (open ring), runtime
  forward table (cross). Connect each triple with a faint segment so the
  effect of the fit and the quantization is visible as a displacement.
* **Nominal paper colors** as small outlined squares at their own
  chromaticity, with a table on the side: pigment, nominal linear RGB,
  encoded hex, spectral pure RGB, runtime forward RGB, per-channel residual
  at `c = eᵢ`.
* **Mixture gamut** toggle: the convex hull of a dense simplex sample, both
  spectral and runtime. This shows how much of the triangle the four pigments
  can actually reach.
* **Click to probe.** Clicking the diagram sets a target xy, converts to a
  displayable sRGB color, runs `PaintMix.encode`, and shows the recovered
  four concentrations, the reconstructed color, and the residual. This is
  the inverse mapping made tangible.

### `/palette` — sRGB space, primaries, and mixtures

* The sRGB cube (WGLMakie 3D): edges, the pure-pigment runtime vertices, the
  nominal colors, and the D65 white point.
* A **mixing curve** pane: for a selected pair (or a weighted set), draw
  `mix(a, b, t)` for `t ∈ [0, 1]` as a path in linear RGB, with the naive
  `(1-t)a + t b` straight line for contrast. It is the clearest single
  demonstration that paint mixing is not RGB interpolation.
* A 4×4 **pairwise ramp matrix**: each cell is a small encoded-sRGB strip of
  the mixture of two pigments.
* The **mixer ramp** used by the paint page as well: pigment selection,
  per-pigment weight sliders, a fraction slider, a large swatch, and the
  encoded hex plus the decoded concentrations.
* Optional **database context**: the workbook's `k and s data` sheet names 24
  pigments (C3:Z3); the config selects four. Read the `Details` sheet and the
  full K/S block with DuckDB (the same `read_xlsx` path as
  `load_spreadsheet_inputs`) and show the other twenty as a faint nominal
  swatch row, clearly labelled "database context, not part of the model". This
  is the literal reading of "nominal values of the other pigments" and costs
  one extra spreadsheet read.

### `/paint` — painting demo

A brush canvas plus the mixer ramp.

* **Canvas.** A coarse RGB bitmap (`--canvas-size`, default 128×128) held as
  `Matrix{PaintMix.Latent{Float32}}` and a derived `Observable{Matrix{RGBf}}`
  rendered with `image!`. Mouse drag paints.
* **Dab.** For each dab, blend the canvas latent toward the brush latent:
  `z ← Latent((1-α)z.c + α z_brush.c, (1-α)z.r + α z_brush.r)`, then
  `decode`. This is exactly `mix`, and it is associative, so overlapping
  strokes compose correctly without re-encoding. The brush latent is either a
  pure pigment `eᵢ`, a `weighted_mix` of the selected pigments, or the
  `encode` of a picked color.
* **Controls.** Pigment weight sliders, brush radius, flow `α`, a "paint
  model / linear RGB" toggle that switches the blend to naive channel
  interpolation for side-by-side comparison, clear, and a color picker whose
  `encode` result is shown as concentrations and residual.
* **Mixer ramp** beside the canvas: the same two/three-pigment ramp as
  `/palette`, so a single pigment or mixture can be loaded as the brush.
* **Performance.** The bitmap is the only large buffer that crosses the
  websocket. At 128² `RGBf` one frame is ~192 KiB; throttle `canvas_obs`
  updates to ~15–20 Hz and update only the dirty rectangle. Measure on the
  target machine and expose the size as a flag. The dab-accumulation variant
  (draw each dab as a `scatter!` marker with the mixed color) is the fallback
  if streaming is too slow; it trades image quality for a fixed per-dab cost.

### `/tables` — the intermediate and output data

The full pair is resident, so 2D slices are cheap. A slice axis selector and
a fixed channel (e.g. `b = 0.5`) drive a column of small heatmaps:

* **Concentration field** `U(rgb)`: one map per pigment.
* **Basin map**: the dominant pigment or active set at each vertex. This is
  the picture behind `continuity_report`'s concentration jumps.
* **Residual field** `‖x − M(U(x))‖`: shows where the signed residual is
  carrying the error, exactly the region the quality report flags.
* **Forward table** `M(c)` on a fixed concentration plane.
* **LUT vs spectral error** on the slice, computed at a reduced resolution
  (e.g. 96²) because it evaluates the spectral model per sample.
* **Validation reports** from the sidecar rendered as tables and small
  plots: `Epush`/`Epull`/objective vs `alpha` from the fit history,
  padding error by depth, history distribution of the quality report, and the
  `acceptance_gates` pass/fail table with the overridden gates highlighted.

### `/fit` — the surrogate fit

* Measured vs fitted `K` and `S` per pigment, and the resulting reflectance
  and chromaticity displacement.
* `Epush`, `Epull`, objective, and seconds vs `alpha` from
  `sidecar["surrogate"]["history"]`, marking the step where `Epush` first
  meets `push_tolerance` and the step the run stopped on.
* Fit diagnostics from `sidecar["surrogate"]["diagnostics"]`, including the
  dense-cube violation and Oklab deviation quoted in the root README.

### `/provenance` and `/`

* Source files with role, path, and SHA-256 prefix.
* Config hash, model id, grid, flags, payload bytes, format version.
* The acceptance-gate table with failures called out, and the note that the
  promoted payload overrode `lut_vs_spectral_channel_max`,
  `lut_vs_spectral_channel_p99`, and `lut_vs_spectral_oklab_p99`.
* Environment versions from the sidecar and the live process.
* Data caveats: PY74 versus the plan's PY73, the contradictory `kins`, the
  unpromoted cross-check workbook, and whichever payload actually loaded.

## Testing

Keep the fast lane browser-free; the browser lane is opt-in and heavy.

* `test_colorimetry.jl` — unit tests against the numbers above: sRGB
  primaries from `inv(xyz_to_rgb)`, D65 xy, raw pigment chromaticity, and the
  identity `xy(mix_rgb(m, eᵢ)) == xy` from the `R′·D65` integral. These pin
  the two independent routes to the same chromaticity.
* `test_server.jl` — start the server on an ephemeral port, `GET /` and
  `GET /healthz`, assert 200 and that the HTML contains the Bonito bootstrap
  and the model id. Model the websocket canary on Bonnie's
  `test/test_canary.jl` so a broken session is caught without a browser.
* A degradation test: start with `--no-sidecar` and check the pages still
  return 200 with the "not available" note.
* An escape test: assert the page contains the Bonito bootstrap as live
  markup (a `setup_connection` call and the `/bonito/...` script URL), not
  `&lt;script&gt;`. A forgotten `|safe` is the most likely silent rendering
  regression.

### Browser end-to-end (`EXPLORER_E2E=1`)

An HTTP-status test proves almost nothing about the figure pages: WGLMakie
renders through WebGL and Bonito opens a websocket. Use
[playwright-julia](https://github.com/frankier/playwright-julia) in a
separate, opt-in lane, structured like Bonnie's suite.

* `test/test_e2e.jl`, skipped unless `EXPLORER_E2E=1`. Playwright.jl is
  unregistered, so it is a `[sources]` URL pin and a test-only dependency;
  only `Pkg.test()` can load it.
* One-time install (driver + browsers):
  `git clone https://github.com/frankier/playwright-julia /tmp/playwright-julia && julia /tmp/playwright-julia/bin/install.jl`.
  `bin/install.jl` takes an engine name to install only one.
* Coverage:
  * every route returns 200 and writes a screenshot;
  * `/` and `/provenance` show the model id and the gate table;
  * `/healthz` matches the JSON contract;
  * `/spectra`: toggle a pigment and a quantity, assert the figure redraws;
  * `/cie`: click the diagram, assert the probe readout appears;
  * `/paint`: drag across the canvas and assert it changes, then flip the
    paint-vs-RGB toggle and assert the result differs.
* WGLMakie can trip browser "slow script" protection; Bonnie carries a
  separate Firefox lane for exactly that. Default to Chromium and expose
  `EXPLORER_E2E_BROWSER` to swap engines.
* Write a screenshot, the console log, and the page errors per page to
  `EXPLORER_E2E_ARTIFACTS`. Set a global action timeout
  (`EXPLORER_E2E_TIMEOUT`, default 20000 ms): figure boot waits on WebGL and
  the first websocket round trip, so 15 s is tight.
* Serve a small deterministic payload (a tiny synthetic `.pmx`) for the E2E
  lane, so it is fast and does not depend on the 96 MiB release artifact.

## Phases

| Phase | Contents |
| --- | --- |
| 0 | Project, `.gitignore`, standalone-environment decision, `ExplorerData` with an injectable constructor, server skeleton, `base.html`, `/healthz`, `/provenance`, and the split test lanes (fast + `EXPLORER_E2E` scaffolding). Verify the registered Oxygen release resolves with the upstream Bonnie pin before anything else. |
| 1 | `/spectra`: derived curves, crosshair, readout, measured-vs-fitted toggle. |
| 2 | `/cie`: locus, wavelength crown, sRGB triangle, pigment triples, click-to-probe. |
| 3 | `/palette`: cube, mixing curves, ramp matrix, mixer ramp, nominal table. |
| 4 | `/paint`: canvas spike (throughput measured) then brushes, mixer ramp, RGB-vs-paint toggle. |
| 5 | `/tables` and `/fit`: slices, basin/residual maps, validation and history plots. |
| 6 | Playwright E2E coverage, optional 24-pigment context overlay, README, polish. |

Phases 1–2 cover the explicit request; 3–4 the palette and painting demo;
5–6 the "intermediate steps" exploration the request invites.

## Big remaining problems

These can sink the viewer or the schedule. They are decisions to make, in
rough order of weight.

1. **WGLMakie websocket throughput for the paint canvas is unmeasured.**
   Each update is a binary frame of `size² × 3` Float32s; at 256² that is
   ~786 KiB per frame, which will not sustain an interactive brush. The plan
   defaults to 128² and throttling, but that is a guess. *Need:* spike and
   measure in Phase 4 before building the rest of `/paint`, and be ready to
   fall back to dab accumulation.
2. **Playwright E2E of WebGL figures is fragile and slow.** The driver and
   browsers install out of band; WGLMakie can trip browser "slow script"
   protection; clicking a WebGL canvas is timing-sensitive; a screenshot of a
   live figure is nondeterministic. This lane will need generous timeouts,
   failure artifacts, and possibly retries, and it will dominate CI time.
   *Need:* decide whether E2E is a required gate or an advisory lane, and
   keep the fast tests independent of it.
3. **Dropping the fork leaves an unregistered, unpinned dependency chain.**
   `Oxygen` is now a normal registered dependency, but `Bonnie.jl` is pinned
   by URL and is pre-release, so Pkg gives no SemVer protection: an upstream
   change or a new Oxygen release can break the build silently.
   `PaintMixPrecompute` and `Playwright` are also unregistered. *Need:* verify
   at Phase 0 that the registered Oxygen release satisfies Bonnie's HTTP 2.6
   compat, and pin commit SHAs for Bonnie and Playwright instead of tracking
   default branches.
4. **The `ExplorerData` singleton blocks test isolation.** A process-wide
   cache is right for serving and wrong for tests that must swap payloads or
   exercise the missing-sidecar path. *Need:* an injectable constructor
   (`ExplorerData(; model, sidecar, ...)`) and a `--no-sidecar` path the
   tests drive, so degradation is tested rather than assumed.
5. **Template escaping is a silent-failure surface.** OteraEngine
   autoescapes, so one forgotten `|safe` renders the Bonito bootstrap as
   visible text and no figure ever appears. *Need:* pass every embedded
   fragment through one helper returning `OteraEngine.SafeString`, and assert
   in `test_server.jl` that the bootstrap is live markup.
6. **The 24-pigment overlay and the cross-check workbook are unverified.**
   The `Details` C.I.-name-to-column mapping is exercised only for the four
   selected pigments, and the reflectance cross-check workbook has never been
   imported. Either can silently mislabel swatches. *Need:* confirm the
   mapping for all 24 columns before showing any, or drop the overlay.
