# explorer — interactive viewer

A local web app that shows the PaintMix input database, the precompute
intermediates, and the runtime payload. It answers "what is in the data, and
what does each processing step do to it?" on one page.

This is a development tool. It is not part of the runtime package and not part
of the root workspace. The tracked root `Manifest.toml` stays free of Makie,
Oxygen, and Bonito. The dependency direction is one way:
`explorer -> precompute -> PaintMix`.

`PLAN.md` holds the full design and the phase list. This README covers how to
run what exists now.

## Status

Phases 0–3 are implemented: the data layer, the server skeleton, the
`/spectra` figure, the `/cie` chromaticity diagram, and the `/palette` cube
and mixer.

| Route | Contents |
| --- | --- |
| `GET /` | overview: provenance banner, data card, links, live Bonito probe |
| `GET /spectra` | measured, fitted, and derived K, S, K/S, R∞, R′ curves with a wavelength crosshair and readout |
| `GET /cie` | CIE 1931 locus with a wavelength crown, sRGB triangle, the three pigment chromaticity layers, nominal colors, mixture gamut, and click-to-probe |
| `GET /palette` | sRGB cube with the pure-pigment vertices and the mixing path of a selected pair, the pairwise mixture ramp matrix, the mixer ramp, and the nominal table |
| `GET /fig/spectra` | the spectra figure alone (debug/iframe) |
| `GET /fig/cie` | the CIE figure alone (debug/iframe) |
| `GET /fig/palette` | the palette cube alone (debug/iframe) |
| `GET /provenance` | configuration, inputs and hashes, environment, payload header, acceptance gates, caveats |
| `GET /healthz` | JSON: model id, grid, source, sidecar state, notes |
| `GET /static/...` | CSS |

`/spectra` is one WGLMakie figure built inside a Bonito `App`. Pigment and
quantity checkboxes, a measured/fitted toggle (shown only when a sidecar
supplied surrogate parameters), and a log-scale toggle drive it. A wavelength
slider moves the crosshair; a server-rendered readout table follows it and
carries the source spreadsheet cell of each K and S sample. The caption lists
the grid and the Saunderson constants and surfaces the `kins` contradiction.
The figure degrades to a note when the input database is absent.

`/cie` draws the spectral locus closed by the line of purples, with a radial
crown of wavelength ticks and labels; the sRGB chromaticity triangle and D65;
and each pigment three times — the raw K/S point, the fitted surrogate ring,
and the runtime forward-table cross — joined so the fit and quantization
displacement is visible. Toggles hide pigments and layers, and switch the
spectral and runtime mixture-gamut hulls. Clicking the diagram runs the
inverse lookup at that chromaticity and shows the recovered concentrations,
the reconstructed color, the residual, and the nearest pigment. The caption
lists the primaries, D65, the pigment layers, and the nominal paper colors.
It degrades to a note when the input database is absent.

`/palette` draws the sRGB cube with its edges, the four pure-pigment runtime
vertices, the nominal paper colors, and the D65 white point. A dropdown picks
a pigment pair and draws the paint mixing path through the cube, with the
naive linear-RGB line dashed for contrast. A 4×4 matrix of server-rendered
strips shows every pairwise mixture ramp. The mixer has one weight slider per
pigment, a paper-to-mixture fraction, a large swatch, and the encoded hex and
recovered concentrations. The nominal table lists each slot's design color
against its runtime vertex and the per-channel residual. The cube and mixer
need only the runtime payload, so the page still works without the input
database; the chromaticity column then shows a dash. The workbook's other
twenty pigments are deliberately not shown, because the
C.I.-name-to-column mapping is unverified (PLAN.md, problem 6).

The remaining figure pages (`/paint`, `/tables`, `/fit`) come in later phases.

## Run

```sh
julia --project=explorer explorer/run.jl
julia --project=explorer explorer/run.jl --port 8081 --no-sidecar
julia --project=explorer explorer/run.jl --payload precompute/output/dev/default.pmx
```

Then open <http://127.0.0.1:8080/>. Stop it with Ctrl-C.

Flags:

| Flag | Meaning |
| --- | --- |
| `--host HOST` | bind address (default `127.0.0.1`) |
| `--port PORT` | TCP port (default `8080`) |
| `--payload PATH` | pin a specific `.pmx` instead of searching |
| `--no-sidecar` | skip the provenance sidecar (exercises the degraded path) |

The first run resolves and precompiles the environment. That includes Makie and
WGLMakie, which the later figure phases use; it is the slow part.

## Data sources

Sources are read from the repository root, resolved from the source file, never
from the working directory.

| Source | Path |
| --- | --- |
| Configuration | `precompute/config/default.toml` |
| Input database | `precompute/output/paintmix.duckdb` |
| Payload, first found | `data/default/default.pmx`, `precompute/output/release/default.pmx`, `precompute/output/dev/default.pmx`, `build/data/payload.pmx`, else synthetic |
| Sidecar | `<payload dir>/<model-id>.toml`, else `precompute/output/release/<model-id>.toml` |

The 96 MiB release payload and the measurement database are untracked. When a
source is absent the viewer degrades instead of failing: it loads the next
payload, or synthetic tables, and records a note. The provenance banner and
`/provenance` show the choice and every note.

`ExplorerData` also has a keyword constructor, so tests inject a synthetic
model without touching disk:

```julia
using Explorer
data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
viewer = Explorer.serve_explorer(data; port = 8080)
```

## Tests

The fast lane is browser-free:

```sh
julia --project=explorer -e 'using Pkg; Pkg.test()'
```

`test_colorimetry.jl` pins the sRGB primaries and D65 from the configuration
matrix, and checks that `xy(mix_rgb(m, eᵢ))` agrees with a direct `R′·D65`
integral. The pigment tests skip when the input database is absent.
`test_spectra.jl` exercises the quantity accessors, the nearest-wavelength
lookup, and the readout on a hand-built `DerivedCurves`, then checks the
caption against the real database when it is present. `test_cie.jl` checks the
hull and clipper on hand-built points, then the locus, primaries, pigment
layers, mixture gamut, probe, and caption against the real database.
`test_palette.jl` checks the pure vertices, the mixing curve endpoints and
naive line, the ramp matrix, the mixer, and the nominal rows on the synthetic
model, then the spectral pure RGB and the chromaticity column against the real
database. `test_server.jl` starts a server on an ephemeral port and checks the
routes, the JSON contract, the degraded sidecar, spectra, CIE, and palette
pages, and a Bonito websocket handshake modeled on Bonnie's canary.

The browser lane is opt-in and heavy. Install the driver and browsers once:

```sh
git clone https://github.com/frankier/playwright-julia /tmp/playwright-julia
julia /tmp/playwright-julia/bin/install.jl
```

Then:

```sh
EXPLORER_E2E=1 julia --project=explorer -e 'using Pkg; Pkg.test()'
```

`EXPLORER_E2E_BROWSER` swaps the engine (default `chromium`),
`EXPLORER_E2E_TIMEOUT` sets the action timeout in milliseconds (default
`20000`), and `EXPLORER_E2E_ARTIFACTS` sets where screenshots are written.

## Dependency notes

`Bonnie.jl` and `playwright-julia` are unregistered and pinned by commit SHA.
`Oxygen` is the registered release (`1.11`), which resolves against Bonnie's
HTTP 2.6 requirement. Do not add an `Oxygen` `[sources]` entry.

The routes live in this module's own Oxygen context (the `@oxidize` macro), so
the viewer does not share Oxygen's global route table with anything else in the
process.

Templates use OteraEngine, whose filter syntax is `|>` (not Jinja's `|`).
Autoescape is on, so every embedded fragment passes through the `safe` filter
as live markup. `test_server.jl` asserts this: a forgotten filter renders the
Bonito bootstrap as visible text and no figure ever appears.
