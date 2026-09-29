# OteraEngine template loading and the server-side HTML helpers.
#
# Templates are constructed once at startup. Fragments produced by Bonnie
# (`head_content`, `app_html`) and every pre-rendered HTML table are passed
# through `safe_html`, because OteraEngine autoescapes by default and a
# forgotten `|> safe` renders the Bonito bootstrap as visible text with no
# error. The page templates use OteraEngine's `|>` filter syntax, not Jinja's
# `|`.

const PAGE_NAMES = (
    "index", "provenance", "spectra", "cie", "palette", "paint", "tables", "fit",
)

const TPL = Dict{String, Template}()

"""
    load_templates!() -> TPL

Parse every page template once. `base.html` is only extended, never rendered
on its own, so it is not in [`PAGE_NAMES`](@ref).
"""
function load_templates!()
    isempty(TPL) || return TPL
    for name in PAGE_NAMES
        TPL[name] = Template(joinpath(TEMPLATES, "$name.html"))
    end
    return TPL
end

"""
    safe_html(s) -> OteraEngine.SafeString

Mark a server-rendered fragment as trusted markup for the template engine.
"""
safe_html(s::AbstractString) = OteraEngine.safe(String(s))

function _esc(s::AbstractString)
    return replace(s, '&' => "&amp;", '<' => "&lt;", '>' => "&gt;", '"' => "&quot;")
end

"""
    render_page(name; title, body, head = "", extra = Dict()) -> String

Render one page template. `head` and `body` are trusted fragments. `extra`
supplies page-specific variables; extra keys a template does not use are
ignored.
"""
function render_page(
        name::AbstractString;
        title::AbstractString, body::AbstractString, head::AbstractString = "",
        extra::AbstractDict = Dict{Symbol, Any}(),
    )
    load_templates!()
    init = Dict{Symbol, Any}(
        :title => title,
        :head => safe_html(head),
        :body => safe_html(body),
        :banner => safe_html(provenance_banner(explorer_data())),
    )
    for (k, v) in extra
        init[Symbol(k)] = v
    end
    return TPL[name](; init = init)
end

# --- server-rendered fragments --------------------------------------------

"""
    html_table(rows; headers = nothing) -> String

A minimal two-column table for flat metadata, or a header row plus rows when
`headers` is given. Every cell is escaped.
"""
function html_table(rows::AbstractVector; headers::Union{Nothing, AbstractVector} = nothing)
    io = IOBuffer()
    print(io, "<table class=\"data\">")
    if headers !== nothing
        print(io, "<thead><tr>")
        for h in headers
            print(io, "<th>", _esc(string(h)), "</th>")
        end
        print(io, "</tr></thead>")
    end
    print(io, "<tbody>")
    for row in rows
        print(io, "<tr>")
        for cell in row
            print(io, "<td>", _esc(string(cell)), "</td>")
        end
        print(io, "</tr>")
    end
    print(io, "</tbody></table>")
    return String(take!(io))
end

function flat_rows(dict::AbstractDict)
    return [(k, string(v)) for (k, v) in sort(collect(dict); by = first)]
end

# A sidecar section, or `nothing` when the sidecar is absent.
sidecar_section(d::ExplorerData, key::AbstractString) =
    d.sidecar === nothing ? nothing : Base.get(d.sidecar, key, nothing)

"""
    provenance_banner(d) -> String

The one-line banner every page carries: which payload loaded, its grid, and
whether the sidecar is present.
"""
function provenance_banner(d::ExplorerData)
    id = PaintMix.model_id(d.model)
    sidecar = d.sidecar === nothing ? "sidecar absent" : "sidecar present"
    return """
    <div class="banner">
      <span class="badge">$(_esc(d.source))</span>
      <span>model <code>$(id)</code></span>
      <span>grid $(PaintMix.grid_n(d.model))<sup>3</sup></span>
      <span>$(sidecar)</span>
    </div>
    """
end

function state_word(ok::Bool)
    return ok ? "available" : "not available"
end

"""
    index_body(d) -> String

The overview data card and the page list.
"""
function index_body(d::ExplorerData)
    id = PaintMix.model_id(d.model)
    rows = [
        ("Model id", id),
        ("Grid", "$(PaintMix.grid_n(d.model))^3 ($(PaintMix.grid_n(d.model)^3) vertices)"),
        ("Source", d.source),
        ("Format version", string(d.model.format_version)),
        ("Flags", "0x" * string(d.model.flags, base = 16, pad = 8)),
        ("Sidecar", state_word(d.sidecar !== nothing)),
        ("Input database", state_word(d.db !== nothing)),
        ("Spectral reference", state_word(d.spectral !== nothing)),
        ("Fitted surrogate", state_word(d.fitted !== nothing)),
    ]
    return """
    <section class="cards">
      <article class="card">
        <h2>Loaded data</h2>
        $(html_table(rows))
      </article>
      <article class="card">
        <h2>Explore</h2>
        <ul class="links">
          <li><a href="/provenance">Provenance</a> — sources, hashes, gates, caveats</li>
          <li><a href="/spectra">Spectra</a> — measured, fitted, and derived curves</li>
          <li><a href="/cie">CIE 1931</a> — locus, gamut, pigment points, click-to-probe</li>
          <li><a href="/palette">Palette</a> — sRGB cube, mixing curves, ramps, mixer</li>
          <li><a href="/paint">Paint</a> — brush canvas and mixer ramp</li>
          <li><a href="/tables">Tables</a> — lookup-table slices and validation</li>
          <li><a href="/fit">Fit</a> — surrogate fit and continuation history</li>
          <li><a href="/healthz">/healthz</a> — JSON status</li>
        </ul>
      </article>
    </section>
    """
end

"""
    provenance_body(d) -> String

Full metadata: configuration, inputs and hashes, environment, payload header,
acceptance gates, and the degradation notes.
"""
function provenance_body(d::ExplorerData)
    io = IOBuffer()
    print(io, "<section class=\"cards\">")
    print(io, "<article class=\"card\"><h2>Payload</h2>")
    print(
        io, html_table(
            [
                ("Source", d.source),
                ("Model id", PaintMix.model_id(d.model)),
                ("Grid edge", string(PaintMix.grid_n(d.model))),
                ("Format version", string(d.model.format_version)),
                ("Flags", "0x" * string(d.model.flags, base = 16, pad = 8)),
            ]
        )
    )
    print(io, "</article>")

    print(io, "<article class=\"card\"><h2>Configuration</h2>")
    if isempty(d.cfg)
        print(io, "<p class=\"muted\">No precompute configuration was found.</p>")
    else
        color = Base.get(d.cfg, "color", Dict())
        grid = Base.get(d.cfg, "grid", Dict())
        spectra = Base.get(d.cfg, "spectra", Dict())
        rows = [
            ("schema_version", Base.get(d.cfg, "schema_version", "?")),
            ("color_space", Base.get(color, "color_space", "?")),
            ("storage_transfer", Base.get(color, "storage_transfer", "?")),
            ("release_n", Base.get(grid, "release_n", "?")),
            ("dev_n", Base.get(grid, "dev_n", "?")),
            ("integration", Base.get(spectra, "integration", "?")),
        ]
        print(io, html_table(rows))
    end
    print(io, "</article>")

    print(io, "<article class=\"card\"><h2>Environment (sidecar)</h2>")
    env = sidecar_section(d, "environment")
    if env isa AbstractDict
        print(io, html_table(flat_rows(env)))
    else
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    end
    print(io, "</article></section>")

    print(io, "<section class=\"cards\">")
    print(io, "<article class=\"card\"><h2>Inputs and hashes</h2>")
    inputs = sidecar_section(d, "inputs")
    if inputs isa AbstractDict
        rows = [
            (k, first(string(v), 16) * "…") for (k, v) in sort(collect(inputs); by = first)
        ]
        print(io, html_table(rows))
    else
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    end
    print(io, "</article>")

    print(io, "<article class=\"card\"><h2>Acceptance gates</h2>")
    gates = sidecar_section(d, "acceptance_gates")
    if gates isa AbstractDict
        rows = [(k, v ? "pass" : "FAIL") for (k, v) in sort(collect(gates); by = first)]
        print(io, html_table(rows))
    else
        print(io, "<p class=\"muted\">Not available: no sidecar.</p>")
    end
    failed = sidecar_section(d, "promoted_with_failed_gates")
    if failed isa AbstractVector && !isempty(failed)
        print(
            io, "<p class=\"warn\">Promoted with overridden gates: ",
            _esc(join(failed, ", ")), "</p>"
        )
    end
    print(io, "</article></section>")

    print(io, "<section class=\"cards\"><article class=\"card\"><h2>Notes and caveats</h2>")
    if isempty(d.notes)
        print(io, "<p class=\"muted\">None.</p>")
    else
        print(io, "<ul class=\"notes\">")
        for note in d.notes
            print(io, "<li>", _esc(note), "</li>")
        end
        print(io, "</ul>")
    end
    print(io, "</article></section>")
    return String(take!(io))
end
