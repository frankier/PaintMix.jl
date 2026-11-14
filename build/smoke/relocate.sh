#!/usr/bin/env bash
# Recheck the non-bundled library and bindings from a separate directory.
set -euo pipefail
: "${TMPDIR:?Set TMPDIR to disk-backed scratch storage}"
: "${JULIA_LIB_DIR:?Set JULIA_LIB_DIR to the tested runtime library directory}"
[ "$#" = 1 ] || { echo "usage: relocate.sh OUT" >&2; exit 2; }
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
out=$(cd "$1" && pwd)
[ -f "$out/reference.txt" ] || { echo "run smoke/run.sh first" >&2; exit 2; }
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
dest="$work/relocated"
mkdir -p "$dest"
cp "$out/paintmix.so" "$out/paintmix.h" "$out/reference.txt" "$dest/"
cp -a "$out/paintmix_py" "$out/paintmix" "$dest/"
gcc -std=c11 -Wall -Wextra -I"$dest" "$here/client.c" \
    "$dest/paintmix.so" -Wl,-rpath,"$dest" -lm -o "$work/client"
python=$(command -v "${PYTHON:-python3}")
if command -v Rscript >/dev/null 2>&1 && \
    Rscript --vanilla -e 'quit(status = !requireNamespace("rdyncall", quietly = TRUE))' >/dev/null 2>&1; then
    rscript=$(command -v Rscript)
    mkdir -p "$work/rlib"
    R CMD INSTALL --no-multiarch --no-docs -l "$work/rlib" "$dest/paintmix" >/dev/null
else
    rscript=""
    [ "${PAINTMIX_REQUIRE_R:-0}" != 1 ] || { echo "R smoke coverage is required" >&2; exit 1; }
fi
export LD_LIBRARY_PATH="$dest:$JULIA_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PAINTMIX_PY_LIBRARY="$dest/paintmix.so"
export PAINTMIX_R_LIBRARY="$dest/paintmix.so"
export PATH="${PAINTMIX_CLEAN_PATH:-$(getconf PATH)}"
if command -v julia >/dev/null 2>&1; then
    echo "set PAINTMIX_CLEAN_PATH to a command path without Julia" >&2
    exit 1
fi
"$work/client" "$dest/reference.txt"
"$python" "$here/client.py" "$dest" "$dest/reference.txt"
if [ -n "$rscript" ]; then
    R_LIBS="$work/rlib${R_LIBS:+:$R_LIBS}" "$rscript" "$here/client.R" "$dest" "$dest/reference.txt"
fi
echo "relocated smoke tests passed without Julia on PATH"
