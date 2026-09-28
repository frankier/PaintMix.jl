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

Phases 0 and 1 are implemented: the data layer, the server skeleton, and the
`/spectra` figure.

| Route | Contents |
| --- | --- |
| `GET /` | overview: provenance banner, data card, links, live Bonito probe |
| `GET /spectra` | measured, fitted, and derived K, S, K/S, R∞, R′ curves with a wavelength crosshair and readout |
| `GET /fig/spectra` | the spectra figure alone (debug/iframe) |
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

The remaining figure pages (`/cie`, `/palette`, `/paint`, `/tables`, `/fit`)
come in later phases.

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
caption against the real database when it is present. `test_server.jl` starts a
server on an ephemeral port and checks the routes, the JSON contract, the
degraded sidecar and spectra pages, and a Bonito websocket handshake modeled on
Bonnie's canary.

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
