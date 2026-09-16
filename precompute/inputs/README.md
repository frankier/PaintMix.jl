# PaintMix precompute inputs

Everything the model depends on, and where it comes from. This file is the
record of the *data source selection*; the machine-readable form is
`precompute/config/default.toml` and `precompute/inputs/checksums.sha256`.

## Selected sources

| Role | File | SHA-256 |
| --- | --- | --- |
| Absorption and scattering spectra, Saunderson constants | `data/Final_artist_database.xlsx` | `2ab4ef751c781a58b6ada5f424eac7f7e0ca191b7644c22ad53524ee64c7d847` |
| Drawdown reflectances (cross-check only) | `data/Reflectance Data for Golden HB 10 mil Drawdowns over White.xlsx` | `584a38368c4af637a1253b6465b9f71493e38c65340092a0cfe9f73b3ed227cf` |
| Algorithm reference | `mixbox.pdf` | see `checksums.sha256` |
| CIE 1931 2-degree observer and D65 | `precompute/inputs/cie_1931_2deg_d65_10nm.csv` | `8add4c8c68ee83491e0fd41e2c0a21debb456e17706a9217d9319868fb4abc61` |

`Final_artist_database.xlsx` is the primary source. Its `k and s data` sheet
holds absorption `K` and scattering `S` directly on the wavelength grid, plus
the Saunderson constants, which is exactly what equations (1)-(7) need.
`Mix data` holds reflectance spectra for the same paints and is the natural
cross-check. The `Reflectance Data …` workbook is a second, independent
measurement (400-700 nm, `%` reflectance plus `K/S`) and is not used to build
the default model; it is kept as evidence for future validation.

Nothing here may be replaced by RGB swatches or the upstream `mixbox`
artifacts. The upstream lookup tables are not established as redistributable,
and the surrogate fit makes the published tables an input in their own right.

## Sheet layout of the primary source

`k and s data`:

| Cell range | Contents |
| --- | --- |
| `B1:D2` | Saunderson constants: `k1 = 0.03`, `k2 = 0.65`, `kins = 1.0` |
| `C3:Z3` | 24 pigment names, one per column |
| `A5` | label `Absorption (k)` |
| `B6:B43` | wavelengths 380, 390, … 750 nm |
| `C6:Z43` | absorption `K` per pigment per wavelength |
| `A44` | label `Scattering (s)` |
| `B45:B82` | the same wavelength grid |
| `C45:Z82` | scattering `S` per pigment per wavelength |

The `Details` sheet maps the same pigments to C.I. names. Row 31-32 of that
sheet repeats the Saunderson constants as `k1 = 0.03`, `k2 = 0.65`,
`kins = 0.0`, contradicting `kins = 1.0` in `k and s data`. The plan's
section 4.2 uses the paper's convention, which has no added specular term, so
neither `kins` value is used; record the discrepancy rather than silently
picking one.

## Pigment order and the one substitution

The paper names Phthalo Blue (PB15:4), Quinacridone Magenta (PR122), Hansa
Yellow, and Titanium White (PW6), and the plan fixes that order. The primary
source provides:

| Slot | Planned pigment | Column | Name in source | C.I. |
| --- | --- | --- | --- | --- |
| 1 | Phthalo Blue | `Q` | Phathalo blue green shade tints | PB 15:4 |
| 2 | Quinacridone Magenta | `K` | k quinacridone magenta | PR 122 |
| 3 | Hansa Yellow | `D` | Arylide (Hansa) yellow opaque | **PY 74** |
| 4 | Titanium White | `Z` | Titanium White | PW 6 |

**Open item:** the plan names PY73, the source has PY74. PY73 and PY74 are
both arylide/Hansa yellows from the same manufacturer with near-identical
hue, but they are not the same pigment and their spectra are not identical.
The default config selects PY74 with the code `PY74` so that no report can
accidentally claim PY73 fidelity. If PY73 spectra are obtained later, add
them to the input set, switch the config entry, and re-run the pipeline; the
table format records the pigment codes, so a different yellow is a different
model.

## The CIE observer and D65 input

The source workbook carries absorption and scattering data but no colorimetry,
so the standard CIE tables are an explicit input rather than something read
from the measurement file. `cie_1931_2deg_d65_10nm.csv` holds `x_bar`,
`y_bar`, `z_bar` (CIE 1931 2-degree standard observer) and the D65 relative
spectral power distribution at 380, 390, ... 750 nm.

The values were extracted from the `colour-science` datasets, which reproduce
the published CIE tables; at 10 nm they agree with the standard tabulation
(for example `x_bar(380) = 0.001368`, `y_bar(555) = 1.0`). D65 is stored as
the standard relative SPD (unscaled). The shipped configuration pins the
observer name and the illuminant, and `build_quadrature` rejects an observer
grid that does not match `spectra` exactly, so resampling the file would be a
new model, not a runtime option.

## Wavelength grid

380-750 nm at 10 nm, 38 samples, taken from the source as stored. The
spectral model must not resample: quadrature weights and the colorimetric
transform are derived from this exact array. A different grid is a different
configuration, not a runtime option.

## Color conventions

  * Core coordinate system: linear-light sRGB, sRGB primaries, D65 white
    point, matching equation (7).
  * Observer: CIE 1931 2-degree, as tabulated with the source. The `Mix data`
    sheet's `L*a*b*` columns are for D50 and are not used.
  * `XYZ -> RGB` matrix and the `1 / Y_D65` normalization come from paper
    equation (7) and are pinned in the config, not recomputed at run time.
  * The model never stores encoded sRGB. The runtime package's
    `linear_from_srgb8` / `srgb8_from_linear` are the only transfer-function
    adapters, and they live outside the mixing path.

## Import step

`precompute/scripts/import_inputs.jl` loads these ranges into a DuckDB
database under `precompute/output/` together with source metadata (file,
sheet, cell range, checksum, units). It is implemented as part of plan step 3;
until then the config and this file are the authoritative record of the
selection.
