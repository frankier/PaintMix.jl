#!/usr/bin/env python3
"""Python smoke test for the compiled PaintMix ABI.

    python3 client.py <out-dir> <reference.txt>

Imports the generated ``paintmix_py`` package from ``<out-dir>`` and checks
every entrypoint against values computed by Julia from the same embedded
payload, plus the documented error contract. Run it with the Julia library
directory on the loader path unless the library was built with ``--bundle``.
"""

import ctypes
import math
import os
import sys

import numpy as np

TOL64 = 1e-12
TOL32 = 2e-6

failures = 0
checks = 0


def fail(what: str) -> None:
    global failures
    print(f"FAIL: {what}", file=sys.stderr)
    failures += 1


def ok() -> None:
    global checks
    checks += 1


def close(what: str, got, want, tol: float) -> None:
    global failures, checks
    got = float(got)
    want = float(want)
    if not math.isfinite(got) or abs(got - want) > tol + tol * abs(want):
        fail(f"{what}: got {got!r} want {want!r}")
        return
    checks += 1


def expect_status(what: str, code: int, expected: int) -> None:
    global failures, checks
    if code != expected:
        fail(f"{what}: status {code}, expected {expected}")
        return
    checks += 1


def main(argv) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    out_dir, ref_path = argv[1], argv[2]
    sys.path.insert(0, os.path.abspath(out_dir))

    import paintmix_py as pm
    from paintmix_py import JLWError, _lowlevel
    from paintmix_py._lowlevel import CVector_borrowed_Float64

    with open(ref_path, "r", encoding="utf-8") as f:
        tokens = f.read().split()
    pos = 0

    def word():
        nonlocal pos
        v = tokens[pos]
        pos += 1
        return v

    def number():
        return float(word())

    def integer():
        return int(word())

    assert word() == "PAINTMIX-REF"
    assert integer() == 1
    grid_n = integer()
    abi = integer()
    id_lo = integer()
    id_hi = integer()

    info = pm.model_info()
    if info.grid_n != grid_n:
        fail("model_info grid_n")
    else:
        ok()
    if info.abi_version != abi:
        fail("model_info abi_version")
    else:
        ok()
    got_id = (
        int.from_bytes(info.model_id[:8], "little", signed=True),
        int.from_bytes(info.model_id[8:], "little", signed=True),
    )
    if got_id != (id_lo, id_hi):
        fail("model_info model_id")
    else:
        ok()
    if pm.abi_version() != abi:
        fail("abi_version")
    else:
        ok()

    assert word() == "mix"
    count = integer()
    for _ in range(count):
        dtype = integer()
        a = np.array([number() for _ in range(3)])
        b = np.array([number() for _ in range(3)])
        t = number()
        want = [number() for _ in range(3)]
        npdtype = np.float32 if dtype else np.float64
        got = pm.mix(a.astype(npdtype), b.astype(npdtype), npdtype(t))
        tol = TOL32 if dtype else TOL64
        for k in range(3):
            close(f"mix[{'f32' if dtype else 'f64'}][{k}]", got[k], want[k], tol)

    assert word() == "encode"
    count = integer()
    for _ in range(count):
        dtype = integer()
        rgb = np.array([number() for _ in range(3)])
        c = [number() for _ in range(4)]
        r = [number() for _ in range(3)]
        want = [number() for _ in range(3)]
        npdtype = np.float32 if dtype else np.float64
        tol = TOL32 if dtype else TOL64
        latent = pm.encode(rgb.astype(npdtype), dtype=npdtype)
        for k in range(4):
            close(f"encode.c[{k}]", latent[k], c[k], tol)
        for k in range(3):
            close(f"encode.r[{k}]", latent[4 + k], r[k], tol)
        got = pm.decode(latent, dtype=npdtype)
        for k in range(3):
            close(f"decode[{k}]", got[k], want[k], tol)

    assert word() == "wmix"
    count = integer()
    for _ in range(count):
        dtype = integer()
        n = integer()
        colors = np.array([number() for _ in range(3 * n)]).reshape(n, 3)
        weights = np.array([number() for _ in range(n)])
        want = [number() for _ in range(3)]
        npdtype = np.float32 if dtype else np.float64
        got = pm.weighted_mix(colors.astype(npdtype), weights.astype(npdtype), dtype=npdtype)
        tol = TOL32 if dtype else TOL64
        for k in range(3):
            close(f"weighted_mix[{k}]", got[k], want[k], tol)

    # bulk_mix must agree with the scalar entrypoint exactly: same kernel.
    n = 4
    a = np.array([[0.1, 0.2, 0.3], [0.9, 0.1, 0.4], [0.05, 0.05, 0.05], [0.8, 0.3, 0.2]])
    b = np.array([[0.7, 0.1, 0.2], [0.2, 0.2, 0.2], [0.4, 0.4, 0.9], [0.1, 0.9, 0.1]])
    ts = np.array([0.0, 0.5, 1.0, 0.25])
    bulk = pm.bulk_mix(a, b, ts)
    for i in range(n):
        scalar = pm.mix(a[i], b[i], ts[i])
        for k in range(3):
            close(f"bulk_mix[{i}][{k}] agrees with mix", bulk[i][k], scalar[k], 0.0)
    for k in range(3):
        close("bulk endpoint t=0", bulk[0][k], a[0][k], 0.0)
        close("bulk endpoint t=1", bulk[2][k], b[2][k], 0.0)

    # Error contract.
    def vec64(values):
        return CVector_borrowed_Float64.from_numpy(values)

    a3 = np.array([0.1, 0.2, 0.3])
    b3 = np.array([0.6, 0.5, 0.4])
    out3 = np.zeros(3)
    encoded = pm.encode(a3)

    try:
        pm.mix(a3, b3, float("nan"))
        fail("NaN fraction did not raise")
    except JLWError as e:
        expect_status("NaN fraction", e.code, pm.PM_ERR_NONFINITE)

    try:
        pm.weighted_mix(a3.reshape(1, 3), np.array([-1.0]))
        fail("negative weight did not raise")
    except JLWError as e:
        expect_status("negative weight", e.code, pm.PM_ERR_WEIGHT)

    try:
        pm.weighted_mix(a3.reshape(1, 3), np.array([0.0]))
        fail("zero total did not raise")
    except JLWError as e:
        expect_status("zero total", e.code, pm.PM_ERR_TOTAL)

    # A short output buffer must be reported as PM_ERR_LENGTH, not corrupt
    # memory. The generated low-level wrapper turns any non-zero status into
    # `JLWError`, which carries the code and message.
    short = np.zeros(2)
    try:
        _lowlevel.paintmix_mix_f64(vec64(a3), vec64(b3), 0.5, vec64(short))
        fail("short output buffer did not raise")
    except JLWError as e:
        expect_status("short output buffer", e.code, pm.PM_ERR_LENGTH)
        if "shorter than element count" not in e.message:
            fail(f"short buffer message: {e.message!r}")
        else:
            ok()

    # count == 0 succeeds and writes nothing.
    sentinel = np.array([7.0, 7.0, 7.0])
    empty = np.zeros(0)
    try:
        _lowlevel.paintmix_bulk_mix_f64(
            vec64(empty), vec64(empty), vec64(empty), vec64(sentinel), 0
        )
        ok()
    except JLWError as e:
        fail(f"zero-length bulk mix raised {e.code}")
    if not np.all(sentinel == 7.0):
        fail("zero-length bulk mix wrote output")
    else:
        ok()

    # Null pointer with a positive count is PM_ERR_NULL.
    null = CVector_borrowed_Float64((3,), ctypes.POINTER(ctypes.c_double)())
    try:
        _lowlevel.paintmix_mix_f64(null, vec64(b3), 0.5, vec64(out3))
        fail("null input did not raise")
    except JLWError as e:
        expect_status("null input", e.code, pm.PM_ERR_NULL)

    # Decoding the encoded value reproduces the input to float precision.
    for k in range(3):
        close(f"round trip[{k}]", pm.decode(encoded)[k], a3[k], TOL64)

    # The sRGB byte adapters are host-side helpers; check the documented
    # rounding rule on a value that needs it.
    assert (pm.linear_to_rgb8([1.0, 0.0, 0.0]) == np.array([255, 0, 0], dtype=np.uint8)).all()

    print(f"Python client: {checks} checks, {failures} failures")
    return 0 if failures == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
