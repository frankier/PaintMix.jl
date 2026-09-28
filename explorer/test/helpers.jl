# Shared fast-lane test helpers. The websocket client mirror is modeled on
# Bonnie's `test/test_canary.jl`, so a broken Bonito embedding surface fails
# here before it fails in a browser.

using Sockets: Sockets
using HTTP: HTTP
using HTTP.WebSockets: WebSockets

function free_port()
    server = Sockets.listen(Sockets.InetAddr(Sockets.ip"127.0.0.1", 0))
    _, port = Sockets.getsockname(server)
    close(server)
    return Int(port)
end

function wait_for(f; timeout = 10.0)
    deadline = time() + timeout
    while !f() && time() < deadline
        sleep(0.05)
    end
    return f()
end

# The browser side of Bonito's protocol: plain msgpack, gzipped only if the
# session negotiated compression.
function client_message(session, msg::AbstractDict)
    bytes = Bonito.MsgPack.pack(msg)
    session.compression_enabled && (bytes = Bonito.transcode(Bonito.GzipCompressor, bytes))
    return bytes
end

"Session id of the page's root session, parsed from the Bonito bootstrap."
function root_session_id(body::AbstractString)
    m = match(r"init_session\(\"([0-9a-f-]+)\"", body)
    m === nothing && error("no Bonito init_session call in the rendered page")
    return m.captures[1]
end
