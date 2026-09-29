# Browser end-to-end lane. Opt-in: set EXPLORER_E2E=1.
#
# Playwright.jl is a test-only, unregistered dependency. Install the driver
# and browsers once, out of band:
#
#   git clone https://github.com/frankier/playwright-julia /tmp/playwright-julia
#   julia /tmp/playwright-julia/bin/install.jl
#
# WGLMakie figures boot through WebGL and the first websocket round trip, so
# the default action timeout is generous. `EXPLORER_E2E_BROWSER` swaps the
# engine; Chromium is the default. WGLMakie can trip Firefox's slow-script
# limit, so that pref is disabled. `EXPLORER_E2E_ARTIFACTS` names where the
# per-page screenshot, console log, and page errors are written.
#
# The lane serves a tiny synthetic payload, never the 96 MiB release artifact,
# so it is fast and deterministic. The untracked input database is still read
# when present, because `/spectra` and `/cie` need the spectral reference; the
# tests that drive those figures skip on a bare checkout. A small synthetic
# sidecar makes the gate, validation, and fit-history rendering deterministic.

using Test

if get(ENV, "EXPLORER_E2E", "") != "1"
    @testset "e2e (skipped)" begin
        @test_skip "set EXPLORER_E2E=1 to run the browser lane"
    end
else
    using HTTP
    using PaintMix
    using Playwright
    using Explorer
    include("helpers.jl")

    const E2E_ARTIFACTS = get(
        ENV, "EXPLORER_E2E_ARTIFACTS", joinpath(@__DIR__, "artifacts")
    )
    const E2E_TIMEOUT = parse(Int, get(ENV, "EXPLORER_E2E_TIMEOUT", "20000"))
    const E2E_BROWSER = Symbol(get(ENV, "EXPLORER_E2E_BROWSER", "chromium"))

    # What containerised CI running as root needs, plus Firefox's slow-script
    # limit off: WGLMakie's first frame can outlast its 10 s default and leave
    # the Bonito spinner up forever.
    const E2E_LAUNCH = (;
        chromium_sandbox = false,
        args = ["--disable-dev-shm-usage"],
        firefox_user_prefs = Dict("dom.max_script_run_time" => 0),
    )

    # A small deterministic sidecar: enough gates, validation history, and
    # diagnostics for `/provenance`, `/tables`, and `/fit` to render their
    # populated state without the untracked release sidecar.
    function synthetic_sidecar()
        return Dict{String, Any}(
            "acceptance_gates" => Dict{String, Any}(
                "roundtrip_f32" => true, "lut_channel_p99" => false,
            ),
            "promoted_with_failed_gates" => ["lut_channel_p99"],
            "surrogate" => Dict{String, Any}(
                "history" => [
                    Dict{String, Any}(
                        "step" => s, "alpha" => 1.0e5 / 2.0^(s - 1),
                        "Epush" => 10.0^(-s), "Epull" => 0.1 / s,
                        "objective" => 10.0^(-s) + 0.1 / s,
                        "seconds" => Float64(s), "iterations" => 5 + s,
                        "converged" => true,
                    ) for s in 1:3
                ],
                "diagnostics" => Dict{String, Any}(
                    "Epush_fit" => 1.0e-9, "max_cube_violation" => 0.0,
                    "max_oklab_dilation" => 0.0, "max_oklab_deviation" => 0.02,
                ),
            ),
        )
    end

    # A tiny `.pmx` pinned by path, with the real database when it exists.
    function e2e_data()
        path = joinpath(mktempdir(), "synthetic.pmx")
        PaintMix.write_model(path, Explorer.synthetic_model())
        base = Explorer.load_explorer_data(; payload_path = path)
        return Explorer.ExplorerData(;
            cfg = base.cfg, db = base.db, spectra = base.spectra, quad = base.quad,
            spectral = base.spectral, fitted = base.fitted, model = base.model,
            source = base.source, sidecar = synthetic_sidecar(),
            derived = base.derived, notes = base.notes,
        )
    end

    # The Bonito websocket is the readiness signal: driving a widget before it
    # connects is a silent no-op.
    bonito_ready(target) =
        wait_for_function(target, "() => window.Bonito?.can_send_to_julia() === true")

    checkboxes(target) = locator(target, "input[type=checkbox]"; strict = false)

    # A page error is a real failure; console errors are logged as warnings
    # because WebGL and context warnings are environment-dependent. Both are
    # in the artifact dump either way.
    function check_clean(page, url)
        errors = page_errors(page)
        isempty(errors) ||
            @warn "uncaught page errors" url errors = [e.message for e in errors]
        @test isempty(errors)
        console = filter(m -> m.type == "error", console_messages(page))
        isempty(console) ||
            @warn "console errors" url errors = [m.text for m in console]
        return
    end

    function with_browser(f)
        return playwright() do pw
            browser = launch(getfield(pw, E2E_BROWSER); headless = true, E2E_LAUNCH...)
            try
                f(browser)
            finally
                close!(browser)
            end
        end
    end

    # One page per route, in its own context, with a screenshot, console log,
    # and page-error dump written unconditionally.
    function with_route_page(f, browser, base, route, name; artifacts = true)
        ctx = new_context(browser)
        set_default_timeout!(ctx, E2E_TIMEOUT)
        set_default_navigation_timeout!(ctx, 60_000)
        url = base * route
        dir = artifacts ? joinpath(E2E_ARTIFACTS, name) : nothing
        try
            result = with_page(
                ctx, url; artifacts = dir, artifacts_on = :always
            ) do page
                r = f(page)
                check_clean(page, url)
                return r
            end
            dir === nothing || @test isfile(joinpath(dir, "screenshot.png"))
            return result
        finally
            close!(ctx)
        end
    end

    # Playwright.jl exposes no high-level mouse API yet, so the drag drives the
    # generated channel layer directly: real, trusted input events. The rev is
    # pinned, so these names are fixed. The first move must land before the
    # press, because the paint app reads `Makie.mouseposition` on mousedown.
    mouse_move!(page, x, y) = Playwright._page_mouse_move(page; x = x, y = y)
    mouse_down!(page) = Playwright._page_mouse_down(page; button = "left")
    mouse_up!(page) = Playwright._page_mouse_up(page; button = "left")

    function canvas_rect(page)
        return evaluate(
            page, """
            () => {
              const cs = [...document.querySelectorAll('canvas')]
                .filter(c => c.offsetParent !== null);
              const c = cs.sort((x, y) =>
                y.width * y.height - x.width * x.height)[0];
              if (!c) return null;
              const r = c.getBoundingClientRect();
              return {x: r.left, y: r.top, w: r.width, h: r.height};
            }
            """
        )
    end

    function drag_canvas!(page, x0, y0, x1, y1; steps = 6)
        # Clicking a control below the fold scrolls the page, so the canvas may
        # be off-viewport. Bring it back before reading its rect.
        evaluate(
            page, """
            () => {
              const cs = [...document.querySelectorAll('canvas')]
                .filter(c => c.offsetParent !== null);
              const c = cs.sort((x, y) =>
                y.width * y.height - x.width * x.height)[0];
              c.scrollIntoView({block: 'center'});
            }
            """
        )
        sleep(0.2)
        rect = canvas_rect(page)
        rect === nothing && return false
        px(f) = rect["x"] + f * rect["w"]
        py(f) = rect["y"] + f * rect["h"]
        mouse_move!(page, px(x0), py(y0))
        mouse_down!(page)
        for i in 1:steps
            t = i / steps
            mouse_move!(page, px(x0 + (x1 - x0) * t), py(y0 + (y1 - y0) * t))
            sleep(0.05)
        end
        mouse_up!(page)
        return true
    end

    # The canvas content as a data URL. Comparing this isolates the WebGL
    # canvas from the rest of the page, whose screenshots are not byte-stable.
    canvas_data(page) = evaluate(page, "() => document.querySelector('canvas').toDataURL()")

    # Wait until two consecutive canvas reads agree, so the first WGLMakie
    # paint has landed before the before/after compare.
    function settled_canvas(page)
        previous = canvas_data(page)
        for _ in 1:40
            sleep(0.25)
            current = canvas_data(page)
            current == previous && return current
            previous = current
        end
        return previous
    end

    # A single drag is one frame update, and WGLMakie occasionally drops the
    # first update after a plot is created. Retry the same stroke until the
    # canvas moves, so the lane is not flaky about that race.
    function drag_until_change!(page, before, x0, y0, x1, y1)
        for _ in 1:3
            drag_canvas!(page, x0, y0, x1, y1) || return false
            changed = retry_until(; timeout = 5_000, on_timeout = :false) do
                canvas_data(page) != before
            end
            changed && return true
        end
        return false
    end

    @testset "e2e ($E2E_BROWSER)" begin
        mkpath(E2E_ARTIFACTS)
        data = e2e_data()
        has_spectra = data.quad !== nothing
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        base = "http://127.0.0.1:$port"
        try
            wait_for(
                () -> HTTP.get(base; status_exception = false, retry = false).status == 200;
                timeout = 30,
            )

            @testset "routes return 200 and screenshot" begin
                pages = [
                    ("/", "index", "Loaded data"),
                    (
                        "/spectra", "spectra",
                        has_spectra ? "Spectra caption" :
                            "Spectral reference not available",
                    ),
                    (
                        "/cie", "cie",
                        has_spectra ? "Click to probe" :
                            "Spectral reference not available",
                    ),
                    ("/palette", "palette", "Pairwise mixture ramps"),
                    ("/paint", "paint", "How the dab works"),
                    ("/tables", "tables", "Acceptance gates"),
                    ("/fit", "fit", "Continuation history"),
                    ("/provenance", "provenance", "Acceptance gates"),
                ]
                with_browser() do browser
                    for (route, name, text) in pages
                        with_route_page(browser, base, route, name) do page
                            expect(page; to_have_title = r"^PaintMix explorer")
                            expect(locator(page, "body"); to_contain_text = text)
                        end
                    end
                end
            end

            @testset "overview and provenance show the model and gates" begin
                with_browser() do browser
                    with_route_page(browser, base, "/", "overview") do page
                        expect(
                            locator(page, "body");
                            to_contain_text = PaintMix.model_id(data.model),
                        )
                    end
                    with_route_page(browser, base, "/provenance", "gates") do page
                        body = locator(page, "body")
                        expect(body; to_contain_text = PaintMix.model_id(data.model))
                        expect(body; to_contain_text = "lut_channel_p99")
                        expect(body; to_contain_text = "FAIL")
                        expect(
                            body; to_contain_text = "Promoted with overridden gates"
                        )
                    end
                end
            end

            if has_spectra
                @testset "/spectra toggles a pigment and a quantity" begin
                    with_browser() do browser
                        with_route_page(browser, base, "/spectra", "spectra-toggle") do page
                            bonito_ready(page)
                            boxes = checkboxes(page)
                            rows = locator(
                                page, "table.readout tbody tr"; strict = false
                            )
                            heads = locator(
                                page, "table.readout thead th"; strict = false
                            )
                            expect(boxes; to_have_count = 10)
                            expect(rows; to_have_count = 4)
                            expect(heads; to_have_count = 8)
                            # Pigment 1 off: one fewer readout row.
                            click!(nth(boxes, 1))
                            expect(rows; to_have_count = 3)
                            # Quantity K off: one fewer readout column.
                            click!(nth(boxes, 5))
                            expect(heads; to_have_count = 7)
                        end
                    end
                end

                @testset "/cie click-to-probe fills the readout" begin
                    with_browser() do browser
                        with_route_page(browser, base, "/cie", "cie-probe") do page
                            bonito_ready(page)
                            body = locator(page, "body")
                            expect(body; to_contain_text = "No probe yet")
                            click!(nth(locator(page, "canvas"; strict = false), 1))
                            expect(body; to_contain_text = "target xy")
                            expect(body; to_contain_text = "concentrations")
                        end
                    end
                end
            end

            @testset "/paint drags and differs under the RGB blend" begin
                with_browser() do browser
                    with_route_page(browser, base, "/paint", "paint-drag") do page
                        bonito_ready(page)
                        before = settled_canvas(page)

                        # Paint the same stroke in each mode, clearing between,
                        # and compare the two canvases. That isolates the blend
                        # from the checkbox and from the stroke position.
                        @test drag_until_change!(page, before, 0.3, 0.3, 0.7, 0.7)
                        paint_stroke = canvas_data(page)

                        cleared = false
                        for _ in 1:3
                            click!(locator(page, "button"))
                            cleared = retry_until(;
                                timeout = 5_000, on_timeout = :false
                            ) do
                                canvas_data(page) == before
                            end
                            cleared && break
                        end
                        @test cleared

                        click!(nth(checkboxes(page), 1))
                        sleep(0.5)
                        @test drag_until_change!(page, before, 0.3, 0.3, 0.7, 0.7)
                        @test canvas_data(page) != paint_stroke
                    end
                end
            end
        finally
            Explorer.close_explorer(viewer)
        end
    end
end
