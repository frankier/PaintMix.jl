#!/usr/bin/env bash
# Run the compiled-ABI smoke tests: C and Python clients against Julia
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

if [ "$status" -ne 0 ]; then
    echo "smoke tests FAILED" >&2
    exit "$status"
fi
echo "smoke tests passed"
