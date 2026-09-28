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
# engine; Chromium is the default.

using Test

if get(ENV, "EXPLORER_E2E", "") != "1"
    @testset "e2e (skipped)" begin
        @test_skip "set EXPLORER_E2E=1 to run the browser lane"
    end
else
    using HTTP
    using Playwright
    using Explorer
    include("helpers.jl")

    const E2E_ARTIFACTS = get(
        ENV, "EXPLORER_E2E_ARTIFACTS", joinpath(@__DIR__, "artifacts")
    )
    const E2E_TIMEOUT = parse(Int, get(ENV, "EXPLORER_E2E_TIMEOUT", "20000"))
    const E2E_BROWSER = get(ENV, "EXPLORER_E2E_BROWSER", "chromium")

    @testset "e2e" begin
        mkpath(E2E_ARTIFACTS)
        data = Explorer.ExplorerData(; model = Explorer.synthetic_model())
        port = free_port()
        viewer = Explorer.serve_explorer(data; port = port, async = true)
        url = "http://127.0.0.1:$port"
        try
            wait_for(
                () -> HTTP.get(url; status_exception = false, retry = false).status == 200;
                timeout = 30,
            )
            playwright() do pw
                browser = launch(engine(pw, E2E_BROWSER); headless = true)
                try
                    page = new_page(browser)
                    set_default_timeout!(page, E2E_TIMEOUT)
                    set_default_navigation_timeout!(page, E2E_TIMEOUT)

                    goto!(page, "$url/")
                    expect(page; to_have_title = "PaintMix explorer")
                    expect(
                        locator(page, "body");
                        to_have_text = r"00000000000000000000000000000000",
                    )
                    screenshot(page; path = joinpath(E2E_ARTIFACTS, "index.png"))

                    goto!(page, "$url/provenance")
                    expect(
                        locator(page, "body"); to_have_text = r"Acceptance gates"
                    )
                    screenshot(page; path = joinpath(E2E_ARTIFACTS, "provenance.png"))

                    goto!(page, "$url/cie")
                    expect(locator(page, "body"); to_have_text = r"CIE 1931")
                    screenshot(page; path = joinpath(E2E_ARTIFACTS, "cie.png"))

                    goto!(page, "$url/palette")
                    expect(
                        locator(page, "body"); to_have_text = r"Pairwise mixture ramps"
                    )
                    screenshot(page; path = joinpath(E2E_ARTIFACTS, "palette.png"))
                finally
                    close!(browser)
                end
            end
        finally
            Explorer.close_explorer(viewer)
        end
    end
end
