# Public façade for the paintmix R package.
#
# `compile.jl` installs this over the starter the generator writes, on every
# build, so `R/lowlevel.R` is never edited and this file is the maintained
# API. It mirrors the Python façade in `build/python/_facade.py`: small
# numeric-returning helpers over the flat, caller-allocated ABI, plus the
# model-info query that the generator leaves for a human because it returns a
# struct.
#
# R has no float32 vector type. `precision = "single"` routes a call through
# the f32 entrypoints using a `floatraw` buffer; the default `"double"` uses
# the f64 entrypoints and R's native `numeric`. Nothing is computed in the
# other precision unless it is asked for.
#
# See `build/README.md` for the buffer contract and `build/library.jl` for
# the entrypoints.

# Status codes from `PaintMix`. A failed call raises a `jlw_*` condition that
# carries the code; the class is the generator's generic mapping, so match on
# `e$code` to tell the PaintMix failures apart.
PM_OK <- 0L
PM_ERR_NULL <- 1L
PM_ERR_LENGTH <- 2L
PM_ERR_NONFINITE <- 3L
PM_ERR_WEIGHT <- 4L
PM_ERR_TOTAL <- 5L

# A borrowed buffer plus what is needed to read it back. For `"double"` the
# buffer is the caller's vector, which the library writes in place. For
# `"single"` it is a `floatraw` copy, because R cannot hold float32 directly.
.pmx_buffer <- function(values, precision) {
  values <- as.double(values)
  if (precision == "single") {
    buffer <- as.floatraw(values)
    carrier <- cdata("CVector_borrowed_Float32")
    carrier$dims <- as.integer(length(values))
    carrier$data <- .jlr_extptr(buffer)
    list(carrier = carrier, buffer = buffer, n = length(values), single = TRUE)
  } else {
    list(
      carrier = .jlr_CVector_borrowed_Float64_arg(values),
      buffer = values,
      n = length(values),
      single = FALSE
    )
  }
}

.pmx_read <- function(b) {
  if (!b$single) {
    return(b$buffer)
  }
  vapply(
    seq_len(b$n),
    function(i) unpack(b$buffer, (i - 1L) * 4L, "f"),
    numeric(1)
  )
}

# An `n x 3` matrix, in which each row is a color, or a flat channel-fastest
# vector: `r0, g0, b0, r1, g1, b1, ...`. R flattens column-major, so a matrix
# is transposed before flattening.
.pmx_flat <- function(x, len) {
  if (is.matrix(x)) {
    if (ncol(x) != 3L) {
      stop("a matrix argument must have 3 columns", call. = FALSE)
    }
    x <- as.double(t(x))
  } else {
    x <- as.double(x)
  }
  if (length(x) != len) {
    stop(sprintf("expected %d values, got %d", len, length(x)), call. = FALSE)
  }
  x
}

#' Compiled-in ABI version
#'
#' @return An integer scalar.
#' @export
abi_version <- function() as.integer(.jlr_paintmix_abi_version())

#' Provenance of the model compiled into the library
#'
#' @return A list with `abi_version`, `format_version`, `grid_n`,
#'   `channel_count`, `flags`, `model_id` (16 raw bytes), `model_id_hex`, and
#'   the two decoded flag booleans.
#' @export
model_info <- function() {
  info <- .jlr_paintmix_model_info()
  .jlr_check_status(info$status)
  # The identifier is a pair of little-endian UInt64 fields. R's `double`
  # cannot hold them, so read the bytes straight out of the returned struct;
  # the layout check done at load time guarantees the offset.
  fields <- attr(info, "typeinfo")$fields
  lo <- fields$offset[fields$name == "model_id_lo"] + 1L
  model_id <- info[lo:(lo + 15L)]
  flags <- as.integer(info$flags)
  list(
    abi_version = as.integer(info$abi_version),
    format_version = as.integer(info$format_version),
    grid_n = as.integer(info$grid_n),
    channel_count = as.integer(info$channel_count),
    flags = flags,
    model_id = model_id,
    model_id_hex = paste(sprintf("%02x", as.integer(model_id)), collapse = ""),
    forward_table_simplex_projected = bitwAnd(flags, 1L) != 0L,
    inverse_table_largest_remainder = bitwAnd(flags, 2L) != 0L
  )
}

#' Encode a linear-light RGB triple
#'
#' @param rgb A numeric vector of length 3.
#' @param precision `"double"` (default) or `"single"`.
#' @return Seven scalars: `c1, c2, c3, c4, r, g, b`.
#' @export
encode <- function(rgb, precision = c("double", "single")) {
  precision <- match.arg(precision)
  src <- .pmx_buffer(rgb, precision)
  dst <- .pmx_buffer(double(7), precision)
  if (precision == "single") {
    .jlr_paintmix_encode_f32(src$carrier, dst$carrier)
  } else {
    .jlr_paintmix_encode_f64(src$carrier, dst$carrier)
  }
  .pmx_read(dst)
}

#' Decode seven scalars into a linear-light RGB triple
#'
#' @param latent A numeric vector of length 7: `c1, c2, c3, c4, r, g, b`.
#' @param precision `"double"` (default) or `"single"`.
#' @return A numeric vector of length 3.
#' @export
decode <- function(latent, precision = c("double", "single")) {
  precision <- match.arg(precision)
  src <- .pmx_buffer(latent, precision)
  dst <- .pmx_buffer(double(3), precision)
  if (precision == "single") {
    .jlr_paintmix_decode_f32(src$carrier, dst$carrier)
  } else {
    .jlr_paintmix_decode_f64(src$carrier, dst$carrier)
  }
  .pmx_read(dst)
}

#' Mix two linear-light RGB colors
#'
#' @param a,b Numeric vectors of length 3.
#' @param t The share of `b`, between 0 and 1.
#' @param precision `"double"` (default) or `"single"`.
#' @return A numeric vector of length 3.
#' @export
mix <- function(a, b, t, precision = c("double", "single")) {
  precision <- match.arg(precision)
  ca <- .pmx_buffer(a, precision)
  cb <- .pmx_buffer(b, precision)
  out <- .pmx_buffer(double(3), precision)
  if (precision == "single") {
    .jlr_paintmix_mix_f32(ca$carrier, cb$carrier, as.double(t), out$carrier)
  } else {
    .jlr_paintmix_mix_f64(ca$carrier, cb$carrier, as.double(t), out$carrier)
  }
  .pmx_read(out)
}

#' Vectorized mixing
#'
#' @param a,b `n x 3` matrices or flat channel-fastest vectors of length `3n`.
#' @param t A numeric vector of length `n`.
#' @param precision `"double"` (default) or `"single"`.
#' @return An `n x 3` matrix.
#' @export
bulk_mix <- function(a, b, t, precision = c("double", "single")) {
  precision <- match.arg(precision)
  t <- as.double(t)
  n <- length(t)
  ca <- .pmx_buffer(.pmx_flat(a, 3L * n), precision)
  cb <- .pmx_buffer(.pmx_flat(b, 3L * n), precision)
  ct <- .pmx_buffer(t, precision)
  out <- .pmx_buffer(double(3L * n), precision)
  if (precision == "single") {
    .jlr_paintmix_bulk_mix_f32(
      ca$carrier, cb$carrier, ct$carrier, out$carrier, n
    )
  } else {
    .jlr_paintmix_bulk_mix_f64(
      ca$carrier, cb$carrier, ct$carrier, out$carrier, n
    )
  }
  matrix(.pmx_read(out), nrow = n, ncol = 3L, byrow = TRUE)
}

#' Weighted average of colors in latent space
#'
#' @param colors An `n x 3` matrix or a flat channel-fastest vector.
#' @param weights A numeric vector of length `n`.
#' @param precision `"double"` (default) or `"single"`.
#' @return A numeric vector of length 3.
#' @export
weighted_mix <- function(colors, weights, precision = c("double", "single")) {
  precision <- match.arg(precision)
  weights <- as.double(weights)
  n <- length(weights)
  cc <- .pmx_buffer(.pmx_flat(colors, 3L * n), precision)
  cw <- .pmx_buffer(weights, precision)
  out <- .pmx_buffer(double(3), precision)
  if (precision == "single") {
    .jlr_paintmix_weighted_mix_f32(cc$carrier, cw$carrier, out$carrier, n)
  } else {
    .jlr_paintmix_weighted_mix_f64(cc$carrier, cw$carrier, out$carrier, n)
  }
  .pmx_read(out)
}

#' Encoded-sRGB bytes to linear light
#'
#' Host-side helper for image plumbing, evaluated in double precision. It
#' mirrors `PaintMix.linear_from_srgb8` but is not a bit-for-bit port.
#'
#' @param rgb8 Integer or numeric bytes in `[0, 255]`.
#' @return A numeric vector.
#' @export
rgb8_to_linear <- function(rgb8) {
  c <- as.double(rgb8) / 255
  ifelse(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055)^2.4)
}

#' Linear light to encoded-sRGB bytes
#'
#' Clips, applies the transfer function, rounds, and clamps. Host-side helper,
#' as for [rgb8_to_linear()].
#'
#' @param rgb A numeric vector in `[0, 1]`.
#' @return An integer vector.
#' @export
linear_to_rgb8 <- function(rgb) {
  c <- pmin(pmax(as.double(rgb), 0), 1)
  e <- ifelse(c <= 0.0031308, 12.92 * c, 1.055 * c^(1 / 2.4) - 0.055)
  as.integer(pmin(pmax(round(255 * e), 0), 255))
}
