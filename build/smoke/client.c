/*
 * C smoke test for the compiled PaintMix ABI.
 *
 *   gcc -std=c11 -I<out> client.c -L<out> -lpaintmix -Wl,-rpath,<out> \
 *       -o client && ./client <out>/reference.txt
 *
 * Checks the compiled entrypoints against values computed by Julia
 * (smoke/reference.jl, which reads the same embedded payload) and against the
 * documented error contract. Exits non-zero on the first failure summary.
 *
 * This file also serves as the worked example of how to call the library:
 * every buffer is caller-owned, every call returns a status, and the only
 * memory the library touches is the caller's.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "paintmix.h"

static int failures = 0;
static int checks = 0;

static void fail(const char *what) {
    fprintf(stderr, "FAIL: %s\n", what);
    failures++;
}

static void ok(void) {
    checks++;
}

/* Reference values differ from the compiled library only by rounding; both
 * run the same table lookups. */
static const double TOL64 = 1e-12;
static const double TOL32 = 2e-6;

static void close_enough(const char *what, double got, double want, double tol) {
    if (!(fabs(got - want) <= tol + tol * fabs(want))) {
        fprintf(stderr, "FAIL: %s: got %.17g want %.17g\n", what, got, want);
        failures++;
        return;
    }
    ok();
}

static int token_string(FILE *f, char *buf, size_t n) {
    (void)n;
    return fscanf(f, "%63s", buf) == 1;
}

static int token_int(FILE *f, long long *v) {
    return fscanf(f, "%lld", v) == 1;
}

static int token_double(FILE *f, double *v) {
    return fscanf(f, "%lf", v) == 1;
}

static CVector_borrowed_Float64 vec64(double *data, int64_t n) {
    CVector_borrowed_Float64 v;
    v.dims[0] = n;
    v.data = data;
    return v;
}

static CVector_borrowed_Float32 vec32(float *data, int64_t n) {
    CVector_borrowed_Float32 v;
    v.dims[0] = n;
    v.data = data;
    return v;
}

static int status_is_ok(JLWStatus s) {
    return s.code == 0;
}

static void report_status(const char *what, JLWStatus s) {
    if (s.code == 0) {
        ok();
        return;
    }
    fprintf(stderr, "FAIL: %s: unexpected status %d: %.256s\n", what, s.code, s.message);
    failures++;
}

/* ---- error contract ---------------------------------------------------- */

static void check_error_contract(void) {
    double a[3] = {0.1, 0.2, 0.3};
    double b[3] = {0.6, 0.5, 0.4};
    double out[3] = {0.0, 0.0, 0.0};
    double big_a[24];
    double big_b[24];
    double ts[8];
    double weights[8];
    float fa[3] = {0.1f, 0.2f, 0.3f};
    float fout[3] = {0.0f, 0.0f, 0.0f};
    float flatent[7] = {0, 0, 0, 0, 0, 0, 0};
    JLWStatus s;

    /* A null data pointer with a positive count is PM_ERR_NULL (1). */
    CVector_borrowed_Float64 null_vec;
    null_vec.dims[0] = 3;
    null_vec.data = NULL;
    s = paintmix_mix_f64(null_vec, vec64(b, 3), 0.5, vec64(out, 3));
    if (s.code != 1) {
        fprintf(stderr, "FAIL: null input gave status %d, expected 1\n", s.code);
        failures++;
    } else {
        ok();
    }

    /* A short buffer is PM_ERR_LENGTH (2). */
    s = paintmix_mix_f64(vec64(a, 2), vec64(b, 3), 0.5, vec64(out, 3));
    if (s.code != 2) {
        fprintf(stderr, "FAIL: short input gave status %d, expected 2\n", s.code);
        failures++;
    } else {
        ok();
    }
    s = paintmix_mix_f64(vec64(a, 3), vec64(b, 3), 0.5, vec64(out, 2));
    if (s.code != 2) {
        fprintf(stderr, "FAIL: short output gave status %d, expected 2\n", s.code);
        failures++;
    } else {
        ok();
    }

    /* Non-finite input is PM_ERR_NONFINITE (3), and the output is untouched. */
    out[0] = 42.0;
    s = paintmix_mix_f64(vec64(a, 3), vec64(b, 3), NAN, vec64(out, 3));
    if (s.code != 3 || out[0] != 42.0) {
        fprintf(stderr, "FAIL: NaN fraction gave status %d, output %g\n", s.code, out[0]);
        failures++;
    } else {
        ok();
    }

    /* A short latent buffer for decode is PM_ERR_LENGTH. */
    s = paintmix_decode_f64(vec64(a, 3), vec64(out, 3));
    if (s.code != 2) {
        fprintf(stderr, "FAIL: short latent gave status %d, expected 2\n", s.code);
        failures++;
    } else {
        ok();
    }

    /* Negative weight is PM_ERR_WEIGHT (4); a zero total is PM_ERR_TOTAL (5). */
    double colors6[6] = {0.1, 0.2, 0.3, 0.6, 0.5, 0.4};
    weights[0] = 1.0;
    weights[1] = -1.0;
    s = paintmix_weighted_mix_f64(vec64(colors6, 6), vec64(weights, 2), vec64(out, 3), 2);
    if (s.code != 4) {
        fprintf(stderr, "FAIL: negative weight gave status %d, expected 4\n", s.code);
        failures++;
    } else {
        ok();
    }
    weights[0] = 0.0;
    weights[1] = 0.0;
    s = paintmix_weighted_mix_f64(vec64(colors6, 6), vec64(weights, 2), vec64(out, 3), 2);
    if (s.code != 5) {
        fprintf(stderr, "FAIL: zero total gave status %d, expected 5\n", s.code);
        failures++;
    } else {
        ok();
    }

    /* count == 0 is a success and writes nothing. */
    out[0] = 7.0;
    s = paintmix_bulk_mix_f64(vec64(big_a, 0), vec64(big_b, 0), vec64(ts, 0),
                              vec64(out, 3), 0);
    if (s.code != 0 || out[0] != 7.0) {
        fprintf(stderr, "FAIL: zero-length bulk mix: status %d output %g\n", s.code, out[0]);
        failures++;
    } else {
        ok();
    }

    /* The f32 path exists and shares the contract. */
    s = paintmix_mix_f32(vec32(fa, 2), vec32(fa, 3), 0.5f, vec32(fout, 3));
    if (s.code != 2) {
        fprintf(stderr, "FAIL: f32 short input gave status %d, expected 2\n", s.code);
        failures++;
    } else {
        ok();
    }
    s = paintmix_encode_f32(vec32(fa, 3), vec32(flatent, 7));
    report_status("encode_f32", s);
}

/* ---- main -------------------------------------------------------------- */

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <reference.txt>\n", argv[0]);
        return 2;
    }
    FILE *f = fopen(argv[1], "r");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        return 2;
    }

    char magic[64];
    long long version = 0;
    if (!token_string(f, magic, sizeof magic) || strcmp(magic, "PAINTMIX-REF") != 0 ||
        !token_int(f, &version) || version != 1) {
        fprintf(stderr, "unrecognized reference file\n");
        fclose(f);
        return 2;
    }

    long long grid_n = 0, abi = 0, id_lo = 0, id_hi = 0;
    if (token_int(f, &grid_n) != 1 || token_int(f, &abi) != 1 ||
        token_int(f, &id_lo) != 1 || token_int(f, &id_hi) != 1) {
        fprintf(stderr, "truncated reference header\n");
        fclose(f);
        return 2;
    }

    ModelInfo info = paintmix_model_info();
    report_status("paintmix_model_info", info.status);
    if (info.status.code == 0) {
        if (info.grid_n != (int32_t)grid_n) fail("model_info grid_n");
        else ok();
        if (info.abi_version != (int32_t)abi) fail("model_info abi_version");
        else ok();
        if ((long long)info.model_id_lo != id_lo || (long long)info.model_id_hi != id_hi) {
            fail("model_info model_id");
        } else {
            ok();
        }
    }
    if ((long long)paintmix_abi_version() != abi) fail("paintmix_abi_version");
    else ok();

    /* mix */
    char section[64];
    long long count = 0;
    if (!token_string(f, section, sizeof section) || strcmp(section, "mix") != 0 ||
        !token_int(f, &count)) {
        fprintf(stderr, "expected mix section\n");
        fclose(f);
        return 2;
    }
    for (long long i = 0; i < count; i++) {
        long long dtype = 0;
        double a[3], b[3], t, want[3];
        if (token_int(f, &dtype) != 1) break;
        for (int k = 0; k < 3; k++) token_double(f, &a[k]);
        for (int k = 0; k < 3; k++) token_double(f, &b[k]);
        token_double(f, &t);
        for (int k = 0; k < 3; k++) token_double(f, &want[k]);
        if (dtype == 0) {
            double out[3] = {0, 0, 0};
            JLWStatus s = paintmix_mix_f64(vec64(a, 3), vec64(b, 3), t, vec64(out, 3));
            if (!status_is_ok(s)) {
                fail("mix_f64 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("mix_f64", out[k], want[k], TOL64);
        } else {
            float fa[3], fb[3], fout[3] = {0, 0, 0};
            for (int k = 0; k < 3; k++) {
                fa[k] = (float)a[k];
                fb[k] = (float)b[k];
            }
            JLWStatus s = paintmix_mix_f32(vec32(fa, 3), vec32(fb, 3), (float)t, vec32(fout, 3));
            if (!status_is_ok(s)) {
                fail("mix_f32 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("mix_f32", fout[k], want[k], TOL32);
        }
    }

    /* encode */
    if (!token_string(f, section, sizeof section) || strcmp(section, "encode") != 0 ||
        !token_int(f, &count)) {
        fprintf(stderr, "expected encode section\n");
        fclose(f);
        return 2;
    }
    for (long long i = 0; i < count; i++) {
        long long dtype = 0;
        double rgb[3], c[4], r[3], want[3];
        if (token_int(f, &dtype) != 1) break;
        for (int k = 0; k < 3; k++) token_double(f, &rgb[k]);
        for (int k = 0; k < 4; k++) token_double(f, &c[k]);
        for (int k = 0; k < 3; k++) token_double(f, &r[k]);
        for (int k = 0; k < 3; k++) token_double(f, &want[k]);
        if (dtype == 0) {
            double latent[7] = {0, 0, 0, 0, 0, 0, 0};
            double out[3] = {0, 0, 0};
            JLWStatus s = paintmix_encode_f64(vec64(rgb, 3), vec64(latent, 7));
            if (!status_is_ok(s)) {
                fail("encode_f64 status");
                continue;
            }
            for (int k = 0; k < 4; k++) close_enough("encode_f64.c", latent[k], c[k], TOL64);
            for (int k = 0; k < 3; k++) close_enough("encode_f64.r", latent[4 + k], r[k], TOL64);
            s = paintmix_decode_f64(vec64(latent, 7), vec64(out, 3));
            if (!status_is_ok(s)) {
                fail("decode_f64 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("decode_f64", out[k], want[k], TOL64);
        } else {
            float frgb[3], latent[7] = {0, 0, 0, 0, 0, 0, 0}, fout[3] = {0, 0, 0};
            for (int k = 0; k < 3; k++) frgb[k] = (float)rgb[k];
            JLWStatus s = paintmix_encode_f32(vec32(frgb, 3), vec32(latent, 7));
            if (!status_is_ok(s)) {
                fail("encode_f32 status");
                continue;
            }
            for (int k = 0; k < 4; k++) close_enough("encode_f32.c", latent[k], c[k], TOL32);
            for (int k = 0; k < 3; k++) close_enough("encode_f32.r", latent[4 + k], r[k], TOL32);
            s = paintmix_decode_f32(vec32(latent, 7), vec32(fout, 3));
            if (!status_is_ok(s)) {
                fail("decode_f32 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("decode_f32", fout[k], want[k], TOL32);
        }
    }

    /* weighted mix */
    double colors[24];
    double weights[8];
    if (!token_string(f, section, sizeof section) || strcmp(section, "wmix") != 0 ||
        !token_int(f, &count)) {
        fprintf(stderr, "expected wmix section\n");
        fclose(f);
        return 2;
    }
    for (long long i = 0; i < count; i++) {
        long long dtype = 0, n = 0;
        double want[3];
        if (token_int(f, &dtype) != 1 || token_int(f, &n) != 1) break;
        if (n > 8) {
            fprintf(stderr, "reference file has n = %lld, test supports 8\n", n);
            failures++;
            break;
        }
        for (long long k = 0; k < 3 * n; k++) token_double(f, &colors[k]);
        for (long long k = 0; k < n; k++) token_double(f, &weights[k]);
        for (int k = 0; k < 3; k++) token_double(f, &want[k]);
        if (dtype == 0) {
            double out[3] = {0, 0, 0};
            JLWStatus s =
                paintmix_weighted_mix_f64(vec64(colors, 3 * n), vec64(weights, n), vec64(out, 3), n);
            if (!status_is_ok(s)) {
                fail("weighted_mix_f64 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("weighted_mix_f64", out[k], want[k], TOL64);
        } else {
            float fc[24], fw[8], fout[3] = {0, 0, 0};
            for (long long k = 0; k < 3 * n; k++) fc[k] = (float)colors[k];
            for (long long k = 0; k < n; k++) fw[k] = (float)weights[k];
            JLWStatus s =
                paintmix_weighted_mix_f32(vec32(fc, 3 * n), vec32(fw, n), vec32(fout, 3), n);
            if (!status_is_ok(s)) {
                fail("weighted_mix_f32 status");
                continue;
            }
            for (int k = 0; k < 3; k++) close_enough("weighted_mix_f32", fout[k], want[k], TOL32);
        }
    }

    /* bulk mix: check a handful of records against the scalar entrypoint. */
    {
        double a[9] = {0.1, 0.2, 0.3, 0.9, 0.1, 0.4, 0.05, 0.05, 0.05};
        double b[9] = {0.7, 0.1, 0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.9};
        double t[3] = {0.0, 0.5, 1.0};
        double out[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        double scalar[3] = {0, 0, 0};
        JLWStatus s = paintmix_bulk_mix_f64(vec64(a, 9), vec64(b, 9), vec64(t, 3), vec64(out, 9), 3);
        report_status("bulk_mix_f64", s);
        for (int i = 0; i < 3; i++) {
            s = paintmix_mix_f64(vec64(&a[3 * i], 3), vec64(&b[3 * i], 3), t[i], vec64(scalar, 3));
            if (!status_is_ok(s)) {
                fail("bulk scalar reference status");
                continue;
            }
            for (int k = 0; k < 3; k++) {
                close_enough("bulk_mix_f64 agrees with mix_f64", out[3 * i + k], scalar[k], 0.0);
            }
        }
        /* Endpoints are exact by contract. */
        for (int k = 0; k < 3; k++) {
            close_enough("bulk endpoint t=0", out[k], a[k], 0.0);
            close_enough("bulk endpoint t=1", out[6 + k], b[6 + k], 0.0);
        }
    }

    fclose(f);
    check_error_contract();

    printf("C client: %d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
