# Validated pipeline inputs, held in a small in-memory form.
#
# The plan asks for the selected spreadsheet ranges to live in a DuckDB
# database so every spectrum can be traced back to the file, sheet, and cell
# range it came from. `InputDatabase` is the loaded form the rest of the
# pipeline consumes; `open_database`/`save_database` move it to and from a
# `.duckdb` file.
#
# Database access goes through DuckDB.jl and DataFrames.jl. DuckDB's own
# spreadsheet and CSV readers turn the inputs into DataFrames, and the same
# connection writes the normalized `.duckdb` file, so no external client or
# hand-written parser is involved. The `excel` extension is installed on
# first use, which needs network access once.
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
    build_info::Dict{String, String}
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
                isfinite(r.value) && r.value > 0 || throw(
                    InputError(
                        "$code $quantity at $(r.wavelength_nm) nm must be finite and positive, " *
                            "got $(r.value)"
                    )
                )
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

# --- DuckDB and DataFrame helpers -----------------------------------------

"""
    _connect(path = ":memory:") -> DuckDB.DB

Open a DuckDB database, or a scratch in-memory one for reading a spreadsheet
or CSV file.
"""
function _connect(path::AbstractString = ":memory:")
    return try
        DuckDB.DB(path)
    catch err
        throw(InputError("cannot open DuckDB database $(path): $(sprint(showerror, err))"))
    end
end

_sql_str(s::AbstractString) = "'" * replace(String(s), "'" => "''") * "'"

_query(con, sql::AbstractString) = DataFrame(DBInterface.execute(con, sql))

"""
    _records(::Type{T}, df) -> Vector{T}

Turn each row of `df` into a `T`. The column names must match the field names
of `T`; use `rename!` first when the database column is spelled differently.
"""
_records(::Type{T}, df::DataFrame) where {T} = T[NamedTuple(row) for row in eachrow(df)]

"""
    _cell(value) -> String

A trimmed spreadsheet cell as a string. Empty cells arrive as `missing`.
"""
_cell(::Missing) = ""
_cell(value::AbstractString) = strip(value)

function _float(s::AbstractString)
    v = tryparse(Float64, s)
    v === nothing && throw(InputError("cannot parse $(repr(s)) as a number"))
    return v
end

# --- writing ---------------------------------------------------------------

function _write_table(con, name::AbstractString, df::DataFrame)
    view = "src_" * name
    DuckDB.register_data_frame(con, df, view)
    try
        DBInterface.execute(con, "CREATE TABLE $name AS SELECT * FROM \"$view\"")
    finally
        DuckDB.unregister_data_frame(con, view)
    end
    return nothing
end

function _source_file_records(db::InputDatabase)
    files = copy(db.source_files)
    any(f -> f.role == "observer", files) || push!(files, db.observer_source)
    return files
end

function _build_info_frame(db::InputDatabase)
    ks = sort!(collect(Base.keys(db.build_info)))
    return DataFrame(key = ks, value = String[db.build_info[k] for k in ks])
end

"""
    save_database(db, path)

Write `db` to a fresh `.duckdb` database at `path`, replacing any file that
already exists. The schema is the flat table-per-record-type form
[`open_database`](@ref) reads back.

The observer table's own source is a row in `source_files`, so it need not be
repeated on every sample.
"""
function save_database(db::InputDatabase, path::AbstractString)
    validate_database(db)
    mkpath(dirname(path))
    isfile(path) && rm(path)
    con = _connect(path)
    try
        _write_table(con, "source_files", DataFrame(_source_file_records(db)))
        _write_table(
            con, "pigments",
            rename!(DataFrame(db.pigments), :column => :column_name)
        )
        _write_table(con, "spectra", DataFrame(db.spectra))
        _write_table(con, "saunderson", DataFrame([db.saunderson]))
        _write_table(con, "observer", DataFrame(db.observer))
        _write_table(con, "build_info", _build_info_frame(db))
    finally
        DBInterface.close!(con)
    end
    return path
end

# --- reading ---------------------------------------------------------------

const _TABLES = ("source_files", "pigments", "spectra", "saunderson", "observer", "build_info")

"""
    open_database(path) -> InputDatabase

Read a `.duckdb` database produced by `import_inputs.jl` into memory.
"""
function open_database(path::AbstractString)
    isfile(path) || throw(
        InputError(
            "no input database at $path; run precompute/scripts/import_inputs.jl first"
        )
    )
    con = _connect(path)
    try
        present = Set(
            String.(
                _query(
                    con,
                    "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
                ).table_name
            )
        )
        for needed in _TABLES
            needed in present || throw(InputError("database $path has no table $needed"))
        end

        source_files = _records(
            SourceFile, _query(
                con,
                "SELECT role, path, sha256 FROM source_files ORDER BY role"
            )
        )
        pigments = _records(
            PigmentRecord, rename!(
                _query(
                    con,
                    "SELECT slot, code, name, ci, column_name FROM pigments ORDER BY slot"
                ),
                :column_name => :column
            )
        )
        spectra = _records(
            SpectrumRecord, _query(
                con,
                "SELECT code, quantity, wavelength_nm, value, source_file, sheet, cell_range " *
                    "FROM spectra ORDER BY code, quantity, wavelength_nm"
            )
        )

        saunderson_rows = _query(
            con,
            "SELECT k1, k2, source_file, sheet, note FROM saunderson LIMIT 1"
        )
        nrow(saunderson_rows) == 1 || throw(InputError("database has no Saunderson row"))
        saunderson = only(_records(SaundersonRecord, saunderson_rows))

        observer = _records(
            ObserverRecord, _query(
                con,
                "SELECT wavelength_nm, x_bar, y_bar, z_bar, d65 FROM observer ORDER BY wavelength_nm"
            )
        )

        observer_rows = _query(
            con,
            "SELECT path, sha256 FROM source_files WHERE role = 'observer' LIMIT 1"
        )
        observer_source = nrow(observer_rows) == 0 ? SourceFile(
                (
                    role = "observer",
                    path = "precompute/inputs/cie_1931_2deg_d65_10nm.csv",
                    sha256 = "",
                )
            ) : SourceFile(
                (
                    role = "observer",
                    path = String(observer_rows.path[1]),
                    sha256 = String(observer_rows.sha256[1]),
                )
            )

        info = _query(con, "SELECT key, value FROM build_info")
        build_info = Dict{String, String}(
            String(info.key[i]) => String(info.value[i]) for i in 1:nrow(info)
        )

        db = InputDatabase(
            source_files, pigments, spectra, saunderson, observer, observer_source,
            build_info,
        )
        return validate_database(db)
    finally
        DBInterface.close!(con)
    end
end

# --- spreadsheet import ----------------------------------------------------

"""
    _load_excel!(con)

Make DuckDB's `excel` extension available. The extension is installed on
first use, which needs network access once.
"""
function _load_excel!(con)
    try
        DBInterface.execute(con, "LOAD excel")
    catch
        DBInterface.execute(con, "INSTALL excel; LOAD excel")
    end
    return con
end

"""
    _read_xlsx(con, path, sheet, range) -> DataFrame

Read a sheet range as raw strings, with a leading `rn` column holding the
sheet row number. Every cell is a `String` or `missing`, so the caller does
the parsing.
"""
function _read_xlsx(con, path::AbstractString, sheet::AbstractString, range::AbstractString)
    return _query(
        con,
        "SELECT row_number() OVER () AS rn, * FROM read_xlsx(" *
            _sql_str(path) * ", sheet=" * _sql_str(sheet) *
            ", all_varchar=true, header=false, range=" * _sql_str(range) * ")"
    )
end

"""
    _read_observer(con, path) -> Vector{ObserverRecord}

Read the observer CSV, checking that its header is exactly the five expected
columns.
"""
function _read_observer(con, path::AbstractString)
    df = _query(con, "SELECT * FROM read_csv(" * _sql_str(path) * ", header = true)")
    expected = ["wavelength_nm", "x_bar", "y_bar", "z_bar", "d65"]
    names(df) == expected || throw(
        InputError(
            "observer file $path has header $(join(names(df), ",")), expected $(join(expected, ","))"
        )
    )
    return _records(ObserverRecord, df)
end

"""
    load_spreadsheet_inputs(cfg, cfg_path; root) -> InputDatabase

Read the selected spreadsheet ranges and the observer CSV into an
[`InputDatabase`](@ref). This is the only function that knows the workbook
layout; `import_inputs.jl` and the tests both call it.
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

    primary_sha = bytes2hex(sha256(read(primary)))
    observer_sha = bytes2hex(sha256(read(observer_file)))
    cross = get(inputs, "reflectance_cross_check", nothing)
    cross_path = cross === nothing ? nothing : joinpath(root, cross)
    cross_sha = cross_path !== nothing && isfile(cross_path) ?
        bytes2hex(sha256(read(cross_path))) : ""

    source_files = SourceFile[
        ("spectral_k_s", inputs["primary"], primary_sha),
        ("observer", inputs["observer_file"], observer_sha),
    ]
    cross_sha == "" || push!(source_files, ("reflectance_cross_check", cross, cross_sha))

    sheet = inputs["primary_sheet"]
    con = _connect()
    local ks, details, observer, duckdb_version
    try
        _load_excel!(con)
        ks = _read_xlsx(con, primary, sheet, "A1:Z82")
        details = _read_xlsx(con, primary, "Details", "E3:G25")
        observer = _read_observer(con, observer_file)
        duckdb_version = try
            String(only(_query(con, "SELECT version() AS version").version))
        catch
            "unknown"
        end
    finally
        DBInterface.close!(con)
    end

    # Details: map the C.I. name in column G to the source paint name in
    # column F. The range starts at column E, so the frame has columns E, F,
    # and G.
    ci_to_name = Dict{String, String}()
    for row in eachrow(details)
        ci = _cell(row.G)
        isempty(ci) && continue
        ci_to_name[ci] = _cell(row.F)
    end

    pigments = PigmentRecord[]
    for (slot, p) in enumerate(cfg["pigments"])
        code = String(p["code"])
        column = String(p["column"])
        push!(
            pigments, (
                slot = slot, code = code,
                name = get(ci_to_name, String(p["ci"]), String(p["name"])),
                ci = String(p["ci"]), column = column,
            )
        )
    end
    column_to_code = Dict(p.column => p.code for p in pigments)

    # Melt the pigment columns of each quantity block into one long form, then
    # drop the empty cells. `rn` is the sheet row number, so it is also the
    # cell reference.
    pigment_columns = String[p["column"] for p in cfg["pigments"]]
    spectra = SpectrumRecord[]
    for (quantity, rows) in (("K", 6:43), ("S", 45:82))
        block = select(ks[rows, :], :rn, :B, pigment_columns...)
        long = stack(
            block, pigment_columns;
            variable_name = :column, value_name = :value
        )
        for row in eachrow(long)
            wavelength = _cell(row.B)
            value = _cell(row.value)
            (isempty(wavelength) || isempty(value)) && continue
            column = String(row.column)
            push!(
                spectra, (
                    code = column_to_code[column], quantity = quantity,
                    wavelength_nm = _float(wavelength), value = _float(value),
                    source_file = inputs["primary"], sheet = sheet,
                    cell_range = "$(column)$(row.rn)",
                )
            )
        end
    end
    sort!(spectra; by = r -> (r.code, r.quantity, r.wavelength_nm))

    saunderson = (
        k1 = _float(_cell(ks.B[2])), k2 = _float(_cell(ks.C[2])),
        source_file = inputs["primary"], sheet = sheet,
        note = "B2 and C2 of 'k and s data'; kins is ignored",
    )
    build_info = Dict{String, String}(
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
