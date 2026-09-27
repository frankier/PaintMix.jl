#!/usr/bin/env julia
# Build the shared library and its C header, Python bindings, and R package.
#
#   julia --project=build-env -e 'using Pkg; Pkg.instantiate()'
#   julia --project=build -e 'using Pkg; Pkg.instantiate()'
#   julia --project=build compile.jl [options]
#
# Options:
#   --payload PATH        payload to embed (default build/data/payload.pmx)
#   --profile NAME        regenerate the payload first: tiny, small, full,
#                         or "grid N"
#   --out DIR             artifact directory (default build/out)
#   --libname NAME        library name (default paintmix)
#   --rname NAME          R package name (default the library name)
#   --bundle              also produce a self-contained runtime bundle
#   --no-python           skip the Python package target
#   --no-r                skip the R package target
#   --no-c                skip the C header target
#   --keep-going          do not delete an existing out/ directory
#
# `juliac` relocates the project into a temporary copy, so this script copies
# the selected payload to `build/data/payload.pmx` first: that path is what
# `library.jl` reads while the image is built.

push!(LOAD_PATH, joinpath(@__DIR__, "build-env"))
using ArgParse
using JuliaLibWrapping
using JuliaC
using Libdl: Libdl
using PaintMix
using TOML: TOML

const BUILD_DIR = @__DIR__
const DEFAULT_PAYLOAD = joinpath(BUILD_DIR, "data", "payload.pmx")

function parse_args(args)
    s = ArgParseSettings(
        prog = "compile.jl",
        description = "Build the shared library and its C header, Python bindings, and R package.",
    )
    @add_arg_table! s begin
        "--payload"
        help = "payload to embed"
        arg_type = String
        default = DEFAULT_PAYLOAD
        "--profile"
        help = "regenerate the payload first: tiny, small, full, or \"grid N\""
        arg_type = String
        "--out"
        help = "artifact directory"
        arg_type = String
        default = joinpath(BUILD_DIR, "out")
        "--libname"
        help = "library name"
        arg_type = String
        default = "paintmix"
        "--rname"
        help = "R package name (default: the library name)"
        arg_type = String
        "--bundle"
        help = "also produce a self-contained runtime bundle"
        action = :store_true
        "--no-c"
        help = "skip the C header target"
        action = :store_false
        dest_name = "c"
        "--no-python"
        help = "skip the Python package target"
        action = :store_false
        dest_name = "python"
        "--no-r"
        help = "skip the R package target"
        action = :store_false
        dest_name = "r"
    end
    return ArgParse.parse_args(args, s)
end

function make_payload(profile::AbstractString)
    jl = joinpath(BUILD_DIR, "scripts", "make_dummy_payload.jl")
    project = BUILD_DIR
    julia = Base.julia_cmd()
    cmd = `$julia --project=$project --startup-file=no $jl $profile $DEFAULT_PAYLOAD`
    println("Generating payload: $cmd")
    run(cmd)
    return DEFAULT_PAYLOAD
end

"""
    compile(opts)

Build `paintmix.so`, `paintmix.h`, the `paintmix_py` ctypes package, and the
`paintmix` R package.
"""
function compile(opts)
    profile = opts["profile"]
    payload = profile === nothing ? String(opts["payload"]) : make_payload(profile)
    isfile(payload) || error(
        "no payload at $payload; pass --profile tiny|small|full to generate one"
    )
    mkpath(dirname(DEFAULT_PAYLOAD))
    if payload != DEFAULT_PAYLOAD
        cp(payload, DEFAULT_PAYLOAD; force = true)
    end
    # Validate here, in an optimized process, so the embedded model can skip
    # the byte loops in the unoptimized image-building interpreter. This is
    # the check that makes `validate = false` in library.jl honest.
    model = read_model(DEFAULT_PAYLOAD)
    println(
        "Validated payload: $(filesize(DEFAULT_PAYLOAD)) bytes, n = $(grid_n(model)), " *
            "id = $(model_id(model))",
    )

    out = String(opts["out"])
    libname = String(opts["libname"])
    targets = AbstractTarget[]
    if opts["c"]
        push!(targets, CTarget(out, libname))
    end
    if opts["python"]
        bundle = opts["bundle"]
        push!(
            targets,
            PythonTarget(
                out, libname * "_py", libname;
                bundle_subdir = bundle ? "bundle" : nothing,
            ),
        )
    end
    rname = opts["rname"] === nothing ? libname : String(opts["rname"])
    if opts["r"]
        bundle = opts["bundle"]
        push!(
            targets,
            RTarget(
                out, rname, libname;
                bundle_subdir = bundle ? "bundle" : nothing,
            ),
        )
    end
    isempty(targets) && error("no targets selected")

    result = build_library(
        joinpath(BUILD_DIR, "library.jl"), targets;
        project = BUILD_DIR,
        libname,
        libdir = out,
        trim = :safe,
        bundle = opts["bundle"],
        verbose = true,
    )
    if opts["python"]
        _finish_python_target(out, libname, opts["bundle"])
    end
    if opts["r"]
        _finish_r_target(out, rname, libname, opts["bundle"])
    end
    println("library:  ", result.library)
    println("abi:      ", result.abi_path)
    if result.bundle_dir !== nothing
        println("bundle:   ", result.bundle_dir)
    end
    for t in result.target_outputs
        println("target:   ", t.dir)
    end
    if !opts["bundle"]
        println(
            "note: this build is not bundled, so a consumer needs Julia's lib " *
                "directory on the loader path; pass --bundle for a self-contained " *
                "package.",
        )
    end
    return result
end

"""
    _finish_python_target(out, libname, bundle)

Install the maintained façade over the generated starter and, for a
non-bundled build, place the shared library inside the package so the
generated loader finds it without an environment variable. The generated
`_lowlevel.py` is never edited.
"""
function _finish_python_target(out::AbstractString, libname::AbstractString, bundle::Bool)
    pkg = joinpath(out, libname * "_py")
    facade = joinpath(BUILD_DIR, "python", "_facade.py")
    if isfile(facade)
        cp(facade, joinpath(pkg, "_facade.py"); force = true)
        println("façade:   ", joinpath(pkg, "_facade.py"))
    end
    if !bundle
        cp(
            joinpath(out, libname * "." * Libdl.dlext),
            joinpath(pkg, libname * "." * Libdl.dlext);
            force = true,
        )
    end
    return nothing
end

"""
    _finish_r_target(out, rname, libname, bundle)

Install the maintained façade over the generated starter and, for a
non-bundled build, place the shared library in `inst/` so the generated
loader finds it once `R CMD INSTALL` has moved that tree to the package root.
The generated `R/lowlevel.R` is never edited.
"""
function _finish_r_target(
        out::AbstractString, rname::AbstractString, libname::AbstractString, bundle::Bool
    )
    pkg = joinpath(out, rname)
    facade = joinpath(BUILD_DIR, "r", "_facade.R")
    if isfile(facade)
        cp(facade, joinpath(pkg, "R", "facade.R"); force = true)
        println("façade:   ", joinpath(pkg, "R", "facade.R"))
    end
    if !bundle
        inst = joinpath(pkg, "inst")
        mkpath(inst)
        cp(
            joinpath(out, libname * "." * Libdl.dlext),
            joinpath(inst, libname * "." * Libdl.dlext);
            force = true,
        )
    end
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    compile(parse_args(ARGS))
end

pop!(LOAD_PATH)
