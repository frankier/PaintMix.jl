#!/usr/bin/env bash
# Run the compiled-ABI smoke tests: C, Python, and R clients against Julia
# reference values derived from the same embedded payload.
#
#   build/smoke/run.sh [out-dir] [payload]
#
# Default out-dir is build/out, default payload is build/data/payload.pmx.
# A non-bundled library still needs Julia's lib directory on the loader path;
# the script finds it from `julia` unless JULIA_LIB_DIR is already set. A
# bundled build (compile.jl --bundle) does not need that.

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
build_dir=$(dirname "$here")

out=${1:-"$build_dir/out"}
payload=${2:-"$build_dir/data/payload.pmx"}

if [ -z "${JULIA_LIB_DIR:-}" ]; then
    bindir=$(julia -e 'print(Sys.BINDIR)')
    JULIA_LIB_DIR="$bindir/../lib"
fi
export LD_LIBRARY_PATH="$out:$JULIA_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PAINTMIX_PY_LIBRARY="${PAINTMIX_PY_LIBRARY:-$out/paintmix.so}"

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

echo "== reference values from $payload"
julia --project="$build_dir" "$here/reference.jl" "$payload" "$out/reference.txt"

echo "== C client"
# The library is named paintmix.so, not libpaintmix.so, so link it by path.
gcc -std=c11 -Wall -Wextra -I"$out" "$here/client.c" \
    "$out/paintmix.so" -Wl,-rpath,"$out" -lm -o "$work/client"
status=0
"$work/client" "$out/reference.txt" || status=1

echo "== Python client"
python3 "$here/client.py" "$out" "$out/reference.txt" || status=1

# The R package needs rdyncall and the R toolchain. Skip with a message when
# they are not installed, so the C and Python checks still run on a host
# without R.
if command -v Rscript >/dev/null 2>&1 && \
    Rscript --vanilla -e 'quit(status = !requireNamespace("rdyncall", quietly = TRUE))' >/dev/null 2>&1; then
    echo "== R client"
    rlib="$work/rlib"
    mkdir -p "$rlib"
    R CMD INSTALL --no-multiarch --no-docs -l "$rlib" "$out/paintmix" >/dev/null || status=1
    R_LIBS="$rlib" Rscript "$here/client.R" "$out" "$out/reference.txt" || status=1
else
    echo "== R client (skipped: Rscript or rdyncall is not available)"
fi

if [ "$status" -ne 0 ]; then
    echo "smoke tests FAILED" >&2
    exit "$status"
fi
echo "smoke tests passed"
