const AXIS_NAVIGATION = (:dragpan, :scrollzoom, :rectanglezoom, :limitreset)

function navigation!(ax, enabled::Bool)
    for name in AXIS_NAVIGATION
        if enabled
            Makie.activate_interaction!(ax, name)
        else
            Makie.deactivate_interaction!(ax, name)
        end
    end
    return ax
end

function tool_icon(path)
    return Bonito.HTML(
        """
        <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.7"
             stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
          <path d="$path"/>
        </svg>
        """
    )
end

function cie_tools(session, ax, mode)
    navigation!(ax, mode[] == "navigate")
    Bonito.on(session, mode) do value
        navigation!(ax, value == "navigate")
    end
    tools = (
        ("navigate", "Pan / zoom", "M12 3v18M3 12h18M9 6l3-3 3 3M9 18l3 3 3-3M6 9l-3 3 3 3M18 9l3 3-3 3"),
        ("pick", "Pick point", "M12 3v4m0 10v4M3 12h4m10 0h4M17 12a5 5 0 1 1-10 0a5 5 0 1 1 10 0"),
    )
    choices = map(tools) do (value, label, path)
        checked = Makie.lift(==(value), mode)
        input = Bonito.DOM.input(
            ; type = "radio", name = "cie-interaction", value = value, checked = checked,
            onchange = js"event => { if (event.target.checked) $(mode).notify($(value)); }",
        )
        return Bonito.DOM.label(
            input, Bonito.DOM.span(tool_icon(path), label; class = "plot-tool-face");
            class = "plot-tool",
        )
    end
    reset = Bonito.Observable(false)
    Bonito.on(session, reset) do _
        Makie.limits!(ax, 0.0, 0.8, 0.0, 0.9)
    end
    return Bonito.DOM.div(
        Bonito.DOM.fieldset(
            Bonito.DOM.legend("Diagram tools"), choices...; class = "plot-mode-group",
        ),
        Bonito.DOM.button(
            tool_icon("M3 10a9 9 0 1 1 2 8M3 4v6h6"), "Reset";
            type = "button", title = "Reset diagram view", class = "plot-reset",
            onclick = js"event => $(reset).notify(true);",
        ); class = "plot-toolbar",
    )
end
