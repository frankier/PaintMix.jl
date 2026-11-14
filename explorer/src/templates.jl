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
        :intro => safe_html(page_intro(name)),
        :evidence => safe_html(evidence_panel(name)),
        :next => safe_html(next_step(name)),
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
    <details class="banner">
      <summary>Loaded: $(_esc(d.source)) · $(PaintMix.grid_n(d.model))³ grid · $(sidecar)</summary>
      <p>Model <code>$(id)</code> · <a href="/provenance">Sources, acceptance gates &amp; caveats →</a></p>
    </details>
    """
end

function state_word(ok::Bool)
    return ok ? "available" : "not available"
end

"""
    index_body(d) -> String

The pipeline landing page, with data availability separated from the story.
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
    steps = IOBuffer()
    for (i, stage) in enumerate(PIPELINE)
        available = i == 1 ? d.spectral !== nothing : i == 2 ? d.fitted !== nothing : true
        status = available ? (i >= 3 ? "Loaded tables" : "Data available") : "Data unavailable"
        print(
            steps, """
            <li class="pipeline-stage">
              <span class="step-number">0$i</span>
              <h3><a href="/$(stage.route)">$(stage.title) <span aria-hidden="true">↗</span></a></h3>
              <p>$(stage.description)</p>
              <p class="stage-output">$(stage.output)</p>
              <span class="stage-status">$status</span>
              $(i == 4 ? "<a class=\"text-link\" href=\"/paint\">Open the paint canvas →</a>" : "")
            </li>
            """
        )
    end
    synthetic = d.source == "synthetic" ? "<p class=\"data-notice\">You are viewing synthetic demonstration tables. Measured pigment behaviour requires the input database and a generated payload.</p>" : ""
    curve = mix_curve(d.model, nominal_rgb(1), nominal_rgb(3); n = 32)
    paint_ramp = ramp_gradient([hex_from_linear(curve.paint[:, j]) for j in 1:32])
    rgb_ramp = ramp_gradient([hex_from_linear(curve.naive[:, j]) for j in 1:32])
    return """
    <section class="hero">
      <div>
      <p class="eyebrow">Inside PaintMix</p>
      <h1>From pigment spectra<br>to paint on screen.</h1>
      <p class="lede">Follow how measured pigments become a fast RGB mixing model.</p>
      <a class="button" href="/spectra">Start with the pigments →</a>
      <a class="text-link" href="/paint">Try painting first</a>
      </div>
      <figure class="hero-comparison">
        <figcaption>Blue + yellow · loaded model</figcaption>
        <span>PaintMix</span><div class="hero-ramp" role="img" aria-label="Blue to yellow mixed with the loaded PaintMix tables" style="background:$paint_ramp"></div>
        <span>Linear RGB</span><div class="hero-ramp" role="img" aria-label="The same blue to yellow mixed in linear RGB" style="background:$rgb_ramp"></div>
        <a class="text-link" href="/palette">Explore the difference →</a>
      </figure>
    </section>
    $synthetic
    <section aria-labelledby="pipeline-title">
      <div class="section-heading"><h2 id="pipeline-title">The pipeline</h2><span>Four stages, from measurement to mixing</span></div>
      <div class="pipeline-phases"><span>Precompute · build once</span><span>Runtime · every mix</span></div>
      <ol class="pipeline">$(String(take!(steps)))</ol>
    </section>
    <aside class="cross-stage">
      <div><p class="eyebrow">Compare across stages</p>
        <h2><a href="/cie">Where do the colours move?</a></h2>
        <p>The CIE 1931 view connects measured pigments, fitted surrogates, and runtime colours.</p>
      </div><a class="text-link" href="/cie">Open the gamut comparison →</a>
    </aside>
    <details class="diagnostics"><summary>Loaded data &amp; validation</summary>
      $(html_table(rows))
      <p><a href="/provenance">Inspect provenance, acceptance gates &amp; caveats →</a></p>
      $(provenance_gates(d))
    </details>
    <p class="paper-credit">Based on Sochorová &amp; Jamriška (2021),
      <a href="$(PAPER_URL)">Practical Pigment Mixing for Digital Painting</a>.
      Expand “paper &amp; code” beside each stage for the relevant passage and implementation.</p>
    """
end

function provenance_gates(d::ExplorerData)
    gates = sidecar_section(d, "acceptance_gates")
    gates isa AbstractDict || return "<p class=\"muted\">Acceptance gates not available: no sidecar.</p>"
    return html_table([(k, v ? "pass" : "FAIL") for (k, v) in sort(collect(gates); by = first)])
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
