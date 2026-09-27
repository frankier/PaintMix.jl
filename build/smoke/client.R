#!/usr/bin/env Rscript
# R smoke test for the compiled PaintMix ABI.
#
#     Rscript client.R <out-dir> <reference.txt>
#
# Loads the generated `paintmix` package from the installed R library and
# checks every entrypoint against values computed by Julia from the same
# embedded payload, plus the documented error contract. Run it with the Julia
# library directory on the loader path unless the library was built with
# `--bundle`.
#
# The package must already be installed, e.g.
# `R CMD INSTALL -l <libdir> <out-dir>/paintmix` and `R_LIBS=<libdir>`.

TOL64 <- 1e-12
TOL32 <- 2e-6

failures <- 0L
checks <- 0L

fail <- function(what) {
  cat(sprintf("FAIL: %s\n", what), file = stderr())
  failures <<- failures + 1L
}

ok <- function() checks <<- checks + 1L

close <- function(what, got, want, tol) {
  got <- as.numeric(got)
  want <- as.numeric(want)
  if (length(got) != 1L || !is.finite(got) || abs(got - want) > tol + tol * abs(want)) {
    fail(sprintf("%s: got %s want %s", what, format(got), format(want)))
    return(invisible(NULL))
  }
  ok()
}

expect_code <- function(what, e, expected) {
  if (is.null(e) || !identical(as.integer(e$code), as.integer(expected))) {
    fail(sprintf("%s: expected code %d", what, expected))
    return(invisible(NULL))
  }
  ok()
}

# A decimal integer string, interpreted as an unsigned 64-bit value, to its
# eight little-endian bytes. Repeated division by 256 avoids the 64-bit
# arithmetic R does not have.
pmx_dec_le <- function(s) {
  digits <- as.integer(strsplit(s, "")[[1L]])
  out <- integer(8L)
  for (i in 1:8) {
    rem <- 0L
    for (j in seq_along(digits)) {
      cur <- rem * 10L + digits[j]
      digits[j] <- cur %/% 256L
      rem <- cur %% 256L
    }
    out[i] <- rem
  }
  as.raw(out)
}

# The reference prints the two identifier halves as signed Int64 decimals.
pmx_i64_le <- function(s) {
  negative <- startsWith(s, "-")
  if (negative) {
    s <- substring(s, 2L)
  }
  bytes <- pmx_dec_le(s)
  if (!negative) {
    return(bytes)
  }
  # Two's complement: invert and add one, least significant byte first.
  bytes <- as.raw(bitwXor(as.integer(bytes), 255L))
  carry <- 1L
  for (i in 1:8) {
    v <- as.integer(bytes[[i]]) + carry
    bytes[[i]] <- as.raw(v %% 256L)
    carry <- v %/% 256L
    if (carry == 0L) break
  }
  bytes
}

main <- function(argv) {
  if (length(argv) < 2L) {
    cat("usage: client.R <out-dir> <reference.txt>\n", file = stderr())
    return(2L)
  }
  out_dir <- argv[[1L]]
  ref_path <- argv[[2L]]
  # Prefer the library the build just wrote over the copy `R CMD INSTALL`
  # placed in the package. The variable is read by `.onLoad`, so it must be
  # set before the package is attached.
  lib <- file.path(out_dir, paste0("paintmix", .Platform$dynlib.ext))
  if (file.exists(lib)) {
    Sys.setenv(PAINTMIX_R_LIBRARY = lib)
  }
  suppressPackageStartupMessages(library(paintmix))

  tokens <- scan(ref_path, what = character(), quiet = TRUE)
  pos <- 1L
  word <- function() {
    v <- tokens[[pos]]
    pos <<- pos + 1L
    v
  }
  num <- function() as.numeric(word())
  int <- function() as.integer(word())
  vec <- function(n) vapply(seq_len(n), function(i) num(), numeric(1))

  stopifnot(word() == "PAINTMIX-REF")
  stopifnot(int() == 1L)
  grid_n <- int()
  abi <- int()
  id_lo <- word()
  id_hi <- word()

  info <- model_info()
  if (info$grid_n != grid_n) fail("model_info grid_n") else ok()
  if (info$abi_version != abi) fail("model_info abi_version") else ok()
  want_id <- c(pmx_i64_le(id_lo), pmx_i64_le(id_hi))
  if (!identical(info$model_id, want_id)) fail("model_info model_id") else ok()
  if (abi_version() != abi) fail("abi_version") else ok()

  stopifnot(word() == "mix")
  count <- int()
  for (i in seq_len(count)) {
    dtype <- int()
    a <- vec(3L)
    b <- vec(3L)
    t <- num()
    want <- vec(3L)
    precision <- if (dtype == 1L) "single" else "double"
    tol <- if (dtype == 1L) TOL32 else TOL64
    got <- mix(a, b, t, precision)
    for (k in 1:3) close(sprintf("mix[%s][%d]", precision, k), got[[k]], want[[k]], tol)
  }

  stopifnot(word() == "encode")
  count <- int()
  for (i in seq_len(count)) {
    dtype <- int()
    rgb <- vec(3L)
    c4 <- vec(4L)
    r3 <- vec(3L)
    want <- vec(3L)
    precision <- if (dtype == 1L) "single" else "double"
    tol <- if (dtype == 1L) TOL32 else TOL64
    latent <- encode(rgb, precision)
    for (k in 1:4) close(sprintf("encode.c[%d]", k), latent[[k]], c4[[k]], tol)
    for (k in 1:3) close(sprintf("encode.r[%d]", k), latent[[4L + k]], r3[[k]], tol)
    got <- decode(latent, precision)
    for (k in 1:3) close(sprintf("decode[%d]", k), got[[k]], want[[k]], tol)
  }

  stopifnot(word() == "wmix")
  count <- int()
  for (i in seq_len(count)) {
    dtype <- int()
    n <- int()
    colors <- matrix(vec(3L * n), nrow = n, ncol = 3L, byrow = TRUE)
    weights <- vec(n)
    want <- vec(3L)
    precision <- if (dtype == 1L) "single" else "double"
    tol <- if (dtype == 1L) TOL32 else TOL64
    got <- weighted_mix(colors, weights, precision)
    for (k in 1:3) close(sprintf("weighted_mix[%d]", k), got[[k]], want[[k]], tol)
  }

  # bulk_mix must agree with the scalar entrypoint exactly: same kernel.
  a <- rbind(
    c(0.1, 0.2, 0.3), c(0.9, 0.1, 0.4),
    c(0.05, 0.05, 0.05), c(0.8, 0.3, 0.2)
  )
  b <- rbind(
    c(0.7, 0.1, 0.2), c(0.2, 0.2, 0.2),
    c(0.4, 0.4, 0.9), c(0.1, 0.9, 0.1)
  )
  ts <- c(0.0, 0.5, 1.0, 0.25)
  bulk <- bulk_mix(a, b, ts)
  for (i in 1:4) {
    scalar <- mix(a[i, ], b[i, ], ts[[i]])
    for (k in 1:3) {
      close(sprintf("bulk_mix[%d][%d] agrees with mix", i, k), bulk[i, k], scalar[[k]], 0.0)
    }
  }
  for (k in 1:3) {
    close("bulk endpoint t=0", bulk[1L, k], a[1L, k], 0.0)
    close("bulk endpoint t=1", bulk[3L, k], b[3L, k], 0.0)
  }

  a3 <- c(0.1, 0.2, 0.3)
  b3 <- c(0.6, 0.5, 0.4)

  e <- tryCatch({mix(a3, b3, NaN); NULL}, error = function(e) e)
  expect_code("NaN fraction", e, PM_ERR_NONFINITE)

  e <- tryCatch({weighted_mix(rbind(a3), -1.0); NULL}, error = function(e) e)
  expect_code("negative weight", e, PM_ERR_WEIGHT)

  e <- tryCatch({weighted_mix(rbind(a3), 0.0); NULL}, error = function(e) e)
  expect_code("zero total", e, PM_ERR_TOTAL)

  arg <- paintmix:::.jlr_CVector_borrowed_Float64_arg
  mix64 <- paintmix:::.jlr_paintmix_mix_f64

  # A short output buffer must be reported as PM_ERR_LENGTH, not corrupt
  # memory. The generated low-level binding raises a condition for any
  # non-zero status.
  e <- tryCatch({mix64(arg(a3), arg(b3), 0.5, arg(numeric(2L))); NULL}, error = function(e) e)
  expect_code("short output buffer", e, PM_ERR_LENGTH)
  if (!is.null(e) && grepl("shorter than element count", conditionMessage(e), fixed = TRUE)) {
    ok()
  } else {
    fail("short buffer message")
  }

  # count == 0 succeeds and writes nothing.
  bm64 <- paintmix:::.jlr_paintmix_bulk_mix_f64
  sentinel <- c(7.0, 7.0, 7.0)
  empty <- numeric(0L)
  e <- tryCatch({
    bm64(arg(empty), arg(empty), arg(empty), arg(sentinel), 0)
    NULL
  }, error = function(e) e)
  if (is.null(e)) ok() else fail(sprintf("zero-length bulk mix raised %s", e$code))
  if (!identical(sentinel, c(7.0, 7.0, 7.0))) {
    fail("zero-length bulk mix wrote output")
  } else {
    ok()
  }

  # A null pointer with a positive count is PM_ERR_NULL.
  null_carrier <- rdyncall::cdata("CVector_borrowed_Float64")
  null_carrier$dims <- 3L
  null_carrier$data <- new("externalptr")
  e <- tryCatch({mix64(null_carrier, arg(b3), 0.5, arg(numeric(3L))); NULL}, error = function(e) e)
  expect_code("null input", e, PM_ERR_NULL)

  # Decoding the encoded value reproduces the input to double precision.
  encoded <- encode(a3)
  for (k in 1:3) close("round trip", decode(encoded)[[k]], a3[[k]], TOL64)

  # The sRGB byte adapters are host-side helpers; check the documented
  # rounding rule on a value that needs it.
  if (identical(linear_to_rgb8(c(1.0, 0.0, 0.0)), c(255L, 0L, 0L))) {
    ok()
  } else {
    fail("linear_to_rgb8")
  }

  cat(sprintf("R client: %d checks, %d failures\n", checks, failures))
  if (failures == 0L) 0L else 1L
}

status <- main(commandArgs(trailingOnly = TRUE))
quit(status = status)
