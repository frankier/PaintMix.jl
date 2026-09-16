#!/usr/bin/env julia
# Import the selected spreadsheet ranges into a DuckDB database.
#
#   julia --project=precompute scripts/import_inputs.jl [--config PATH] [--db PATH] [--force]
#
# Plan step 3. The import is specified in precompute/inputs/README.md:
#
#   * `k and s data` B6:B43 / C6:Z43  absorption K, 380-750 nm at 10 nm
#   * `k and s data` B45:B82 / C45:Z82 scattering S on the same grid
#   * `k and s data` B2/C2            Saunderson k1 and k2
#   * `Details`      E3:G25           source paint name and C.I. mapping
#   * `precompute/inputs/cie_1931_2deg_d65_10nm.csv`  observer and D65
#
# Each spectrum records its source file, sheet, and cell, so a value in the
# database can always be traced back to the measurement. Requires the `duckdb`
# command line client with the `excel` extension.

using PaintMixPrecompute
using Printf: @printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

function parse_args(args)
    opts = Dict{String,Any}(
        "config" => PaintMixPrecompute.DEFAULT_CONFIG_PATH,
        "db" => nothing,
        "force" => false,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("--config", "--db")
            i < length(args) || error("$a needs a value")
            opts[a[3:end]] = args[i + 1]
            i += 1
        elseif a == "--force"
            opts["force"] = true
        elseif a == "--help" || a == "-h"
            print("""
            usage: import_inputs.jl [options]

              --config PATH   configuration file (default: config/default.toml)
              --db PATH       output database (default: from the configuration)
              --force         overwrite an existing database
            """)
            exit(0)
        else
            error("unknown argument $a")
        end
        i += 1
    end
    return opts
end

function main(args)
    opts = parse_args(args)
    cfg_path = opts["config"]
    cfg = load_config(cfg_path)
    db_path = something(
        opts["db"], joinpath(ROOT, get(cfg["inputs"], "database", "precompute/output/paintmix.duckdb"))
    )
    if isfile(db_path) && !opts["force"]
        @printf("database already exists: %s (use --force to overwrite)\n", db_path)
        return nothing
    end
    @printf("reading inputs from %s\n", joinpath(ROOT, cfg["inputs"]["primary"]))
    db = load_spreadsheet_inputs(cfg, cfg_path; root = ROOT)
    @printf(
        "read %d spectra for %d pigments, %d observer samples\n",
        length(db.spectra), length(db.pigments), length(db.observer)
    )
    save_database(db, db_path)
    @printf("wrote %s\n", db_path)
    for f in db.source_files
        @printf("  %-24s %s\n", f.role, f.sha256[1:16] * "…")
    end
    return nothing
end

main(ARGS)
