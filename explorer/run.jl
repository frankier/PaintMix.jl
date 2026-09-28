#!/usr/bin/env julia
# Entry point for the PaintMix explorer.
#
#   julia --project=explorer explorer/run.jl
#   julia --project=explorer explorer/run.jl --port 8081 --no-sidecar
#
# Loads the data, registers the routes, and serves until interrupted.

using Explorer
using PaintMix
using Printf: @printf

const USAGE = """
Usage: run.jl [--host HOST] [--port PORT] [--payload PATH] [--no-sidecar]

  --host HOST     bind address (default 127.0.0.1)
  --port PORT     TCP port (default 8080)
  --payload PATH  pin a specific .pmx payload instead of searching
  --no-sidecar    skip the provenance sidecar (exercises the degraded path)
  -h, --help      show this message
"""

function parse_args(args)
    opts = Dict{String, Any}(
        "host" => "127.0.0.1", "port" => 8080, "payload" => nothing, "sidecar" => true,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--host"
            i += 1
            opts["host"] = args[i]
        elseif arg == "--port"
            i += 1
            opts["port"] = parse(Int, args[i])
        elseif arg == "--payload"
            i += 1
            opts["payload"] = args[i]
        elseif arg == "--no-sidecar"
            opts["sidecar"] = false
        elseif arg in ("-h", "--help")
            print(USAGE)
            exit(0)
        else
            error("unknown argument: $arg\n$USAGE")
        end
        i += 1
    end
    return opts
end

function main(args = ARGS)
    opts = parse_args(args)
    data = load_explorer_data(;
        with_sidecar = opts["sidecar"], payload_path = opts["payload"]
    )
    for note in data.notes
        @warn note
    end
    viewer = serve_explorer(data; host = opts["host"], port = opts["port"], async = true)
    @printf(
        "PaintMix explorer on http://%s:%d/ (model %s, source %s)\n",
        opts["host"], opts["port"], model_id(data.model), data.source
    )
    wait(viewer.server)
    return viewer
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
