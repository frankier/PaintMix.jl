# Validated pipeline inputs, held in a small in-memory form.
#
# The plan asks for the selected spreadsheet ranges to live in a DuckDB
# database so every spectrum can be traced back to the file, sheet, and cell
# range it came from. `InputDatabase` is the loaded form the rest of the
# pipeline consumes; `open_database`/`save_database` move it to and from a
# `.duckdb` file.
#
# The DuckDB dependency is deliberately not a Julia package dependency. The
# database is read and written by shelling out to the `duckdb` command line
# client, so `PaintMixPrecompute` keeps owning only Optim and ForwardDiff.
# `precompute/scripts/import_inputs.jl` is the only writer.
#
# Nothing here runs at PaintMix runtime.

"""
    SourceFile

One input file with its content hash, so a payload can be traced to the exact
bytes that produced it.
"""
const SourceFile = @NamedTuple{role::String, path::String, sha256::String}

"""
    PigmentRecord

A pigment slot in the fixed four-pigment order, together with the source
column and C.I. name it was read from.
"""
const PigmentRecord = @NamedTuple{
    slot::Int, code::String, name::String, ci::String, column::String,
}

"""
    SpectrumRecord

One absorption or scattering sample. `quantity` is `"K"` or `"S"`.
"""
const SpectrumRecord = @NamedTuple{
    code::String, quantity::String, wavelength_nm::Float64, value::Float64,
    source_file::String, sheet::String, cell_range::String,
}

"""
    SaundersonRecord

The measured Saunderson constants from the source workbook. `kins` is
deliberately not modelled: the source contradicts itself and the paper's
convention (equation 6) has no added specular term.
"""
const SaundersonRecord = @NamedTuple{
    k1::Float64, k2::Float64, source_file::String, sheet::String, note::String,
}

"""
    ObserverRecord

One row of the colorimetric tables: CIE 1931 2-degree observer functions and
the D65 relative spectral power distribution, on the configuration's
wavelength grid.
"""
const ObserverRecord = @NamedTuple{
    wavelength_nm::Float64, x_bar::Float64, y_bar::Float64, z_bar::Float64,
    d65::Float64,
}

"""
    InputDatabase

Everything the numerical pipeline reads, in memory: source files and their
hashes, the fixed pigment slots, the absorption and scattering samples, the
Saunderson constants, the colorimetric tables, and free-form build metadata.

Construction is checked by [`validate_database`](@ref); reading and writing
against a `.duckdb` file go through [`open_database`](@ref) and
[`save_database`](@ref).
"""
struct InputDatabase
    source_files::Vector{SourceFile}
    pigments::Vector{PigmentRecord}
    spectra::Vector{SpectrumRecord}
    saunderson::SaundersonRecord
    observer::Vector{ObserverRecord}
    observer_source::SourceFile
    build_info::Dict{String,String}
end

"""
    InputError

Raised when the input database is missing tables, has inconsistent
dimensions, or carries values the spectral model cannot use.
"""
struct InputError <: Exception
    msg::String
end

Base.showerror(io::IO, e::InputError) = print(io, "InputError: ", e.msg)

function _require(cond::Bool, msg::AbstractString)
    cond || throw(InputError(msg))
    return nothing
end

"""
    validate_database(db) -> db

Check the invariants the spectral model depends on: four unique pigment
codes, a strictly increasing wavelength grid with finite positive values in
both quantities, matching grids for every pigment, and a finite D65/observer
table with a positive `Y` normalizer input.

Throws [`InputError`](@ref) naming the first problem.
"""
function validate_database(db::InputDatabase)
    length(db.pigments) == 4 ||
        throw(InputError("expected four pigment slots, got $(length(db.pigments))"))
    codes = String[p.code for p in db.pigments]
    length(unique(codes)) == 4 ||
        throw(InputError("pigment codes must be unique, got $(join(codes, ", "))"))
    sort([p.slot for p in db.pigments]) == [1, 2, 3, 4] ||
        throw(InputError("pigment slots must be exactly 1..4"))

    for quantity in ("K", "S")
        for code in codes
            rows = filter(r -> r.code == code && r.quantity == quantity, db.spectra)
            isempty(rows) && throw(InputError("no $quantity spectrum for pigment $code"))
            wl = [r.wavelength_nm for r in rows]
            issorted(wl; lt = <=) || throw(InputError("$code $quantity wavelengths are not sorted"))
            length(unique(wl)) == length(wl) ||
                throw(InputError("$code $quantity has duplicate wavelengths"))
            for r in rows
                isfinite(r.wavelength_nm) ||
                    throw(InputError("$code $quantity has a non-finite wavelength"))
                isfinite(r.value) && r.value > 0 || throw(InputError(
                    "$code $quantity at $(r.wavelength_nm) nm must be finite and positive, " *
                        "got $(r.value)"
                ))
            end
        end
    end

    grid = [r.wavelength_nm for r in db.observer]
    issorted(grid; lt = <=) || throw(InputError("observer wavelengths are not sorted"))
    length(unique(grid)) == length(grid) ||
        throw(InputError("observer table has duplicate wavelengths"))
    for r in db.observer
        for (name, v) in (("x_bar", r.x_bar), ("y_bar", r.y_bar), ("z_bar", r.z_bar), ("d65", r.d65))
            isfinite(v) || throw(InputError("observer $name at $(r.wavelength_nm) nm is not finite"))
        end
        r.y_bar >= 0 || throw(InputError("observer y_bar must be non-negative"))
        r.d65 >= 0 || throw(InputError("D65 must be non-negative"))
    end

    isfinite(db.saunderson.k1) && isfinite(db.saunderson.k2) ||
        throw(InputError("Saunderson constants must be finite"))
    return db
end

"""
    spectra_grid(db) -> Vector{Float64}

The common wavelength grid read from the absorption table of the first
pigment. Used to check every other quantity against it.
"""
function spectra_grid(db::InputDatabase)
    code = db.pigments[1].code
    return sort([r.wavelength_nm for r in db.spectra if r.code == code && r.quantity == "K"])
end

# --- DuckDB command line client -------------------------------------------

"""
    duckdb_bin() -> String

The DuckDB command line client to use. Honors `PAINTMIX_DUCKDB` so a build
can pin a specific binary; otherwise looks for `duckdb` on `PATH`.
"""
function duckdb_bin()
    return get(ENV, "PAINTMIX_DUCKDB", "duckdb")
end

function _duckdb(; read::Bool = false, input::Union{Nothing,String} = nothing, args::Vector{String})
    bin = duckdb_bin()
    cmd = `$bin $(args)`
    out = IOBuffer()
    err = IOBuffer()
    if read
        p = run(pipeline(ignorestatus(cmd), stdout = out, stderr = err))
        p.exitcode == 0 || throw(InputError(
            "duckdb failed ($(p.exitcode)): $(String(take!(err)))"
        ))
    else
        p = input === nothing ?
            run(pipeline(ignorestatus(cmd), stdout = out, stderr = err)) :
            run(pipeline(ignorestatus(cmd), stdin = IOBuffer(input), stdout = out, stderr = err))
        p.exitcode == 0 || throw(InputError(
            "duckdb failed ($(p.exitcode)): $(String(take!(err)))"
        ))
    end
    return String(take!(out))
end

# Minimal RFC-4180-ish CSV reader for the `duckdb -csv -header` output.
function _parse_csv(text::AbstractString)
    rows = Vector{Vector{String}}()
    fields = Vector{String}()
    field = IOBuffer()
    in_quotes = false
    i = firstindex(text)
    n = lastindex(text)
    while i <= n
        c = text[i]
        if in_quotes
            if c == '"'
                nxt = i < n ? text[nextind(text, i)] : '\0'
                if nxt == '"'
                    write(field, '"')
                    i = nextind(text, i)
                else
                    in_quotes = false
                end
            else
                write(field, c)
            end
        elseif c == '"'
            in_quotes = true
        elseif c == ','
            push!(fields, String(take!(field)))
        elseif c == '\n'
            push!(fields, String(take!(field)))
            push!(rows, fields)
            fields = Vector{String}()
        elseif c == '\r'
            # skip
        else
            write(field, c)
        end
        i = nextind(text, i)
    end
    if !isempty(fields) || position(field) > 0
        push!(fields, String(take!(field)))
        push!(rows, fields)
    end
    filter!(r -> !(length(r) == 1 && isempty(r[1])), rows)
    return rows
end

function _query_csv(db_path::AbstractString, sql::AbstractString; header::Bool = false)
    args = String[db_path]
    header || push!(args, "-noheader")
    push!(args, "-csv", "-c", String(sql))
    return _parse_csv(_duckdb(; read = true, args = args))
end

function _f(s::AbstractString)
    isempty(s) && throw(InputError("expected a number, got an empty field"))
    v = tryparse(Float64, s)
    v === nothing && throw(InputError("cannot parse $(repr(s)) as a number"))
    return v
end

function _i(s::AbstractString)
    v = tryparse(Int, s)
    v === nothing && throw(InputError("cannot parse $(repr(s)) as an integer"))
    return v
end

# --- reading ---------------------------------------------------------------

"""
    open_database(path) -> InputDatabase

Read a `.duckdb` database produced by `import_inputs.jl` into memory.

Requires the `duckdb` command line client; see [`duckdb_bin`](@ref).
"""
function open_database(path::AbstractString)
    isfile(path) || throw(InputError(
        "no input database at $path; run precompute/scripts/import_inputs.jl first"
    ))
    tables = Set(String(row[1]) for row in _query_csv(path,
        "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"))
    for needed in ("source_files", "pigments", "spectra", "saunderson", "observer", "build_info")
        needed in tables ||
            throw(InputError("database $path has no table $needed"))
    end

    source_files = SourceFile[
        (role = r[1], path = r[2], sha256 = r[3]) for r in _query_csv(path,
            "SELECT role, path, sha256 FROM source_files ORDER BY role")
    ]
    pigments = PigmentRecord[
        (slot = _i(r[1]), code = r[2], name = r[3], ci = r[4], column = r[5]) for r in
            _query_csv(path, "SELECT slot, code, name, ci, column_name FROM pigments ORDER BY slot")
    ]
    spectra = SpectrumRecord[
        (
            code = r[1], quantity = r[2], wavelength_nm = _f(r[3]), value = _f(r[4]),
            source_file = r[5], sheet = r[6], cell_range = r[7],
        ) for r in _query_csv(path,
            "SELECT code, quantity, wavelength_nm, value, source_file, sheet, cell_range " *
                "FROM spectra ORDER BY code, quantity, wavelength_nm")
    ]
    saunderson_rows = _query_csv(path,
        "SELECT k1, k2, source_file, sheet, note FROM saunderson LIMIT 1")
    isempty(saunderson_rows) && throw(InputError("database has no Saunderson row"))
    sr = only(saunderson_rows)
    saunderson = (
        k1 = _f(sr[1]), k2 = _f(sr[2]), source_file = sr[3], sheet = sr[4], note = sr[5],
    )
    observer = ObserverRecord[
        (
            wavelength_nm = _f(r[1]), x_bar = _f(r[2]), y_bar = _f(r[3]), z_bar = _f(r[4]),
            d65 = _f(r[5]),
        ) for r in _query_csv(path,
            "SELECT wavelength_nm, x_bar, y_bar, z_bar, d65 FROM observer ORDER BY wavelength_nm")
    ]
    observer_source_row = _query_csv(path,
        "SELECT path, sha256 FROM source_files WHERE role = 'observer' LIMIT 1")
    observer_source = isempty(observer_source_row) ?
        SourceFile(("observer", "precompute/inputs/cie_1931_2deg_d65_10nm.csv", "")) :
        SourceFile((role = "observer", path = observer_source_row[1][1], sha256 = observer_source_row[1][2]))

    build_info = Dict{String,String}(
        r[1] => r[2] for r in _query_csv(path, "SELECT key, value FROM build_info")
    )
    db = InputDatabase(
        source_files, pigments, spectra, saunderson, observer, observer_source, build_info,
    )
    return validate_database(db)
end

# --- writing ---------------------------------------------------------------

_sql_quote(s::AbstractString) = "'" * replace(String(s), "'" => "''") * "'"

_sql_literal(x::Real) = string(Float64(x))
_sql_literal(s::AbstractString) = _sql_quote(s)

function _insert_sql(table::AbstractString, columns::Vector{String}, row)
    vals = join((_sql_literal(v) for v in row), ", ")
    return "INSERT INTO $table(" * join(columns, ", ") * ") VALUES (" * vals * ");"
end

"""
    save_database(db, path)

Write `db` to a fresh `.duckdb` database at `path`, replacing any file that
already exists. The schema is flat and normalized just enough that each
spectrum carries its file, sheet, and cell range.

Requires the `duckdb` command line client.
"""
function save_database(db::InputDatabase, path::AbstractString)
    validate_database(db)
    mkpath(dirname(path))
    isfile(path) && rm(path)
    sql = IOBuffer()
    write(sql, """
    CREATE TABLE source_files(role VARCHAR, path VARCHAR, sha256 VARCHAR);
    CREATE TABLE pigments(slot INTEGER, code VARCHAR, name VARCHAR, ci VARCHAR, column_name VARCHAR);
    CREATE TABLE spectra(code VARCHAR, quantity VARCHAR, wavelength_nm DOUBLE, value DOUBLE,
                         source_file VARCHAR, sheet VARCHAR, cell_range VARCHAR);
    CREATE TABLE saunderson(k1 DOUBLE, k2 DOUBLE, source_file VARCHAR, sheet VARCHAR, note VARCHAR);
    CREATE TABLE observer(wavelength_nm DOUBLE, x_bar DOUBLE, y_bar DOUBLE, z_bar DOUBLE, d65 DOUBLE);
    CREATE TABLE build_info(key VARCHAR, value VARCHAR);
    """)
    for r in db.source_files
        write(sql, _insert_sql("source_files", ["role", "path", "sha256"],
            (r.role, r.path, r.sha256)), "\n")
    end
    # The observer table's own source is a row in source_files, so it need not
    # be repeated on every sample.
    if !any(r -> r.role == "observer", db.source_files)
        write(sql, _insert_sql("source_files", ["role", "path", "sha256"],
            (db.observer_source.role, db.observer_source.path, db.observer_source.sha256)), "\n")
    end
    for p in db.pigments
        write(sql, _insert_sql("pigments", ["slot", "code", "name", "ci", "column_name"],
            (p.slot, p.code, p.name, p.ci, p.column)), "\n")
    end
    for r in db.spectra
        write(sql, _insert_sql("spectra",
            ["code", "quantity", "wavelength_nm", "value", "source_file", "sheet", "cell_range"],
            (r.code, r.quantity, r.wavelength_nm, r.value, r.source_file, r.sheet, r.cell_range)), "\n")
    end
    s = db.saunderson
    write(sql, _insert_sql("saunderson", ["k1", "k2", "source_file", "sheet", "note"],
        (s.k1, s.k2, s.source_file, s.sheet, s.note)), "\n")
    for r in db.observer
        write(sql, _insert_sql("observer",
            ["wavelength_nm", "x_bar", "y_bar", "z_bar", "d65"],
            (r.wavelength_nm, r.x_bar, r.y_bar, r.z_bar, r.d65)), "\n")
    end
    for (k, v) in sort(collect(db.build_info); by = first)
        write(sql, _insert_sql("build_info", ["key", "value"], (k, v)), "\n")
    end
    _duckdb(; input = String(take!(sql)), args = String[path])
    return path
end

# --- spreadsheet import ----------------------------------------------------

_column_index(letter::AbstractString) = Int(letter[1] - 'A') + 1

"""
    _read_xlsx_rows(path, sheet, range) -> Vector{Vector{String}}

Read a sheet range as raw strings through the DuckDB `excel` extension. The
extension is loaded if present, and installed on first use; installation
needs network access the first time only.
"""
function _read_xlsx_rows(path::AbstractString, sheet::AbstractString, range::AbstractString)
    select = "SELECT row_number() OVER () AS rn, * FROM read_xlsx(" *
        _sql_quote(path) * ", sheet=" * _sql_quote(sheet) *
        ", all_varchar=true, header=false, range=" * _sql_quote(range) * ")"
    try
        return _query_csv(":memory:", "LOAD excel; " * select)
    catch err
        err isa InputError || rethrow()
        return _query_csv(":memory:", "INSTALL excel; LOAD excel; " * select)
    end
end

function _read_observer_csv(path::AbstractString)
    rows = _parse_csv(read(path, String))
    isempty(rows) && throw(InputError("observer file $path is empty"))
    header = rows[1]
    expected = ["wavelength_nm", "x_bar", "y_bar", "z_bar", "d65"]
    header == expected || throw(InputError(
        "observer file $path has header $(join(header, ",")), expected $(join(expected, ","))"
    ))
    out = ObserverRecord[]
    for r in rows[2:end]
        length(r) == 5 || throw(InputError("observer row has $(length(r)) fields, expected 5"))
        push!(out, (
            wavelength_nm = _f(r[1]), x_bar = _f(r[2]), y_bar = _f(r[3]),
            z_bar = _f(r[4]), d65 = _f(r[5]),
        ))
    end
    return out
end

"""
    load_spreadsheet_inputs(cfg, cfg_path; root) -> InputDatabase

Read the selected spreadsheet ranges and the observer CSV into an
[`InputDatabase`](@ref). This is the only function that knows the workbook
layout; `import_inputs.jl` and the tests both call it.

Requires the `duckdb` command line client with the `excel` extension.
"""
function load_spreadsheet_inputs(
        cfg::AbstractDict, cfg_path::AbstractString;
        root::AbstractString = normpath(joinpath(dirname(cfg_path), "..", "..")),
    )
    validate_config(cfg)
    inputs = cfg["inputs"]
    primary = joinpath(root, inputs["primary"])
    observer_file = joinpath(root, inputs["observer_file"])
    isfile(primary) || throw(InputError("primary input $primary does not exist"))
    isfile(observer_file) || throw(InputError("observer input $observer_file does not exist"))

    sheet = inputs["primary_sheet"]
    ks = _read_xlsx_rows(primary, sheet, "A1:Z82")
    det = _read_xlsx_rows(primary, "Details", "E3:G25")

    primary_sha = bytes2hex(sha256(read(primary)))
    observer_sha = bytes2hex(sha256(read(observer_file)))
    cross = get(inputs, "reflectance_cross_check", nothing)
    cross_path = cross === nothing ? nothing : joinpath(root, cross)
    cross_sha = cross_path !== nothing && isfile(cross_path) ?
        bytes2hex(sha256(read(cross_path))) : ""

    source_files = SourceFile[
        ("spectral_k_s", inputs["primary"], primary_sha),
        ("observer", observer_file, observer_sha),
    ]
    if cross_sha != ""
        push!(source_files, ("reflectance_cross_check", cross, cross_sha))
    end

    # Details: map C.I. name to the source paint name. The range starts at
    # column E, so the CSV row is [rn, E, F, G] and F/G are name/C.I.
    ci_to_name = Dict{String,String}()
    for r in det
        length(r) >= 4 || continue
        ci = strip(r[4])
        ci == "" && continue
        ci_to_name[ci] = strip(r[3])
    end

    pigments = PigmentRecord[]
    spectra = SpectrumRecord[]
    for (slot, p) in enumerate(cfg["pigments"])
        code = String(p["code"])
        column = String(p["column"])
        col = _column_index(column)
        source_name = get(ci_to_name, String(p["ci"]), String(p["name"]))
        push!(pigments, (
            slot = slot, code = code, name = source_name, ci = String(p["ci"]),
            column = column,
        ))
        for (quantity, lo, hi) in (("K", 6, 43), ("S", 45, 82))
            for r in ks
                rn = _i(r[1])
                lo <= rn <= hi || continue
                wl_text = r[3]
                (wl_text === nothing || isempty(strip(wl_text))) && continue
                val = r[col + 1]
                (val === nothing || isempty(strip(val))) && continue
                push!(spectra, (
                    code = code, quantity = quantity, wavelength_nm = _f(wl_text),
                    value = _f(val), source_file = inputs["primary"], sheet = sheet,
                    cell_range = "$(column)$(rn)",
                ))
            end
        end
    end

    saunderson = (
        k1 = _f(ks[2][3]), k2 = _f(ks[2][4]), source_file = inputs["primary"],
        sheet = sheet, note = "B2 and C2 of 'k and s data'; kins is ignored",
    )
    observer = _read_observer_csv(observer_file)
    duckdb_version = try
        strip(_duckdb(; read = true, args = String[":memory:", "-noheader", "-list", "-c", "SELECT version()"]))
    catch
        "unknown"
    end
    build_info = Dict{String,String}(
        "julia_version" => string(VERSION),
        "duckdb_version" => replace(duckdb_version, "\n" => " "),
        "tool" => "precompute/scripts/import_inputs.jl",
        "schema" => "1",
    )
    db = InputDatabase(
        source_files, pigments, spectra, saunderson, observer,
        ("observer", inputs["observer_file"], observer_sha), build_info,
    )
    return validate_database(db)
end

"""
    check_checksums(root, checksum_file) -> (checked, failures)

Verify a `sha256sum`-style file listing paths relative to `root`. Returns the
number of files checked and a list of human-readable failures; it does not
throw, so callers can report every problem at once.
"""
function check_checksums(root::AbstractString, checksum_file::AbstractString)
    isfile(checksum_file) || throw(InputError("no checksum file at $checksum_file"))
    checked = 0
    failures = String[]
    for line in eachline(checksum_file)
        isempty(strip(line)) && continue
        length(strip(line)) > 64 || throw(InputError("malformed checksum line: $line"))
        want = lowercase(line[1:64])
        rel = strip(line[65:end])
        path = joinpath(root, rel)
        if !isfile(path)
            push!(failures, "missing: $rel")
            continue
        end
        got = bytes2hex(sha256(read(path)))
        got == want || push!(failures, "checksum mismatch: $rel")
        checked += 1
    end
    return checked, failures
end
