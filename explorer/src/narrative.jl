# The same stage map drives the landing page and each figure's context.
const PAPER_URL = "https://dcgi.fel.cvut.cz/wp-content/wpallimport-dist/publications/pdf/publications-2021-sochorova-tog-pigments-paper.pdf"
const CODE_URL = "https://github.com/frankier/PaintMix.jl/blob/main/"

const PIPELINE = (
    (
        title = "Measure the pigments", route = "spectra", label = "Spectra",
        question = "Where does a paint's colour come from?",
        description = "Measured absorption and scattering predict the colour of a pigment mixture.",
        output = "K, S → spectral reference", reference = "§2 · Eqs. (1–7)", page = 3,
        passage = "absorption and scattering coefficients",
        path = "precompute/src/kubelka_munk.jl",
        first = "        k = c1 * K[l, 1]", last = "        rc = saunderson",
    ),
    (
        title = "Fit a displayable palette", route = "fit", label = "Fit",
        question = "How can pigment mixtures fit inside sRGB?",
        description = "Adjust the spectra to reduce out-of-gamut colours while preserving the original mixing behaviour.",
        output = "Spectral reference → surrogate pigments", reference = "§3.1 · Eqs. (15–17)", page = 6,
        passage = "inside the RGB cube",
        path = "precompute/src/surrogates.jl",
        first = "    for (step, alpha) in enumerate(alpha_schedule(st))", last = "        push!(history, entry)",
    ),
    (
        title = "Bake the lookup tables", route = "tables", label = "Tables",
        question = "How does the expensive model become fast?",
        description = "Solve RGB-to-concentration inverses, sample concentrations back to RGB, then quantize and validate both tables.",
        output = "Surrogate pigments → inverse + forward tables", reference = "§3.2 · Eq. (9)", page = 7,
        passage = "precompute the expensive operations",
        path = "precompute/src/solvers.jl",
        first = "            rgb = SVector{3, T}(i / (n - 1)", last = "            _record!(summary, r)",
    ),
    (
        title = "Mix and paint", route = "palette", label = "Palette",
        question = "How do arbitrary RGB colours mix like paint?",
        description = "Encode colours as concentrations plus a residual, blend those latents, then decode to RGB.",
        output = "RGB → latents → blend → RGB", reference = "§3 · Eqs. (12–14), Fig. 7", page = 5,
        passage = "additive residual term",
        path = "src/mixing.jl",
        first = "    x, y, z = rgb", last = "    return Latent(SVector(c[1]",
    ),
)

function stage_index(name)
    name == "paint" && return 4
    name == "cie" && return 2
    return findfirst(s -> s.route == name, PIPELINE)
end

function pipeline_nav(name)
    io = IOBuffer()
    print(io, "<nav class=\"pipeline-nav\" aria-label=\"Pipeline stages\">")
    for (i, stage) in enumerate(PIPELINE)
        current = name != "cie" && stage_index(name) == i ? " aria-current=\"step\"" : ""
        print(io, "<a href=\"/$(stage.route)\"$current><span>0$i</span> $(stage.label)</a>")
    end
    print(io, "</nav>")
    return String(take!(io))
end

function page_intro(name)
    name == "provenance" && return """
    <div class="page-intro">
      <p class="eyebrow">Supporting evidence</p>
      <h1>Know what you are looking at.</h1>
      <p class="lede">Sources, model identity, and acceptance gates for the loaded data.</p>
      <a href="/">Back to the pipeline →</a>
    </div>
    """
    i = stage_index(name)
    i === nothing && return ""
    stage = PIPELINE[i]
    heading = name == "cie" ? "Compare the colour gamuts" :
        name == "paint" ? "Try the mixing model" : stage.question
    description = name == "cie" ?
        "Compare measured pigments, fitted surrogates, and runtime tables in one chromaticity diagram. Follow the colours across stages 01–04." :
        name == "paint" ? "Paint with the table-based model and compare it with linear RGB blending." : stage.description
    return """
    $(pipeline_nav(name))
    <div class="page-intro">
      <p class="eyebrow">$(name == "cie" ? "Across the pipeline" : "Stage 0$i · $(i == 4 ? "Runtime" : "Precompute")")</p>
      <h1>$heading</h1><p class="lede">$description</p>
    </div>
    """
end

# Read small, allow-listed excerpts from this checkout, so displayed code
# cannot silently diverge from the evaluator being explored.
function source_excerpt(stage)
    lines = readlines(joinpath(ROOT, stage.path))
    a = findfirst(line -> startswith(line, stage.first), lines)
    a === nothing && return nothing
    b = findnext(line -> startswith(line, stage.last), lines, a)
    (b === nothing || b - a > 16) && return nothing
    return (; code = join(lines[a:b], '\n'), first = a, last = b)
end

function evidence_panel(name)
    i = stage_index(name)
    i === nothing && return ""
    stage = PIPELINE[i]
    excerpt = source_excerpt(stage)
    code = excerpt === nothing ? "<p>Source excerpt unavailable in this checkout.</p>" : """
        <pre><code>$(_esc(excerpt.code))</code></pre>
        <a href="$(CODE_URL)$(stage.path)#L$(excerpt.first)-L$(excerpt.last)">$(stage.path):$(excerpt.first)–$(excerpt.last) ↗</a>
        """
    caveat = i == 2 ? "<p class=\"muted\">This implementation defaults to L-BFGS in softplus coordinates; the paper uses L-BFGS-B. Saved history contains continuation summaries, not a replay of pigment motion.</p>" :
        i == 3 ? "<p class=\"muted\">Simplex padding and joint concentration rounding are implementation choices. The validation reports below measure their effects.</p>" : ""
    return """
    <details class="explanation">
      <summary>Why this works · paper &amp; code</summary>
      <div class="evidence-grid">
        <div><p class="eyebrow">From the paper</p>
          <blockquote>“$(stage.passage)”</blockquote>
          <a href="$(PAPER_URL)#page=$(stage.page)">$(stage.reference) · PDF p. $(stage.page) ↗</a>
          $caveat
        </div>
        <div><p class="eyebrow">In this implementation</p>$code</div>
      </div>
    </details>
    """
end

function next_step(name)
    i = stage_index(name)
    i === nothing && return ""
    route, label = name == "palette" ? ("paint", "Try painting") :
        i == 4 ? ("", "Back to the pipeline") :
        (PIPELINE[i + 1].route, PIPELINE[i + 1].title)
    related = i == 2 && name != "cie" ? "<a href=\"/cie\">Compare gamuts in CIE 1931 →</a>" : ""
    return "<nav class=\"next-step\" aria-label=\"Continue exploring\">$related<a href=\"/$route\">$label →</a></nav>"
end

function diagnostic_details(body; label = "Readouts, conventions & validation")
    return "<details class=\"diagnostics\"><summary>$(_esc(label))</summary>$body</details>"
end
