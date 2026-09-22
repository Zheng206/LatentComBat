test_that("empty_matrix has requested dimensions", {
  z <- LatentComBat:::empty_matrix(7L)
  expect_equal(dim(z), c(7L, 0L))
})

test_that("reduce_design_qr removes zero and collinear columns", {
  x <- seq_len(8)
  Z <- cbind(intercept = 1, x = x, twice_x = 2 * x, zero = 0)
  out <- LatentComBat:::reduce_design_qr(Z)

  expect_equal(out$rank, 2L)
  expect_equal(nrow(out$Z), 8L)
  expect_equal(ncol(out$Z), 2L)
  expect_equal(crossprod(out$Q), diag(2), tolerance = 1e-10)
})

test_that("reduce_design_qr handles a zero-column design", {
  Z <- matrix(numeric(0), nrow = 5L, ncol = 0L)
  out <- LatentComBat:::reduce_design_qr(Z)
  expect_equal(out$rank, 0L)
  expect_equal(dim(out$Z), c(5L, 0L))
  expect_equal(dim(out$Q), c(5L, 0L))
})

test_that("build_Z_preserve includes or excludes batch as requested", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 3L)

  with_batch <- LatentComBat:::build_Z_preserve(
    model = stats::lm,
    formula = fx$formula,
    covar = fx$covar,
    bat = fx$bat,
    preserve_batch = TRUE,
    return_Q = TRUE
  )
  without_batch <- LatentComBat:::build_Z_preserve(
    model = stats::lm,
    formula = fx$formula,
    covar = fx$covar,
    bat = fx$bat,
    preserve_batch = FALSE,
    return_Q = TRUE
  )

  expect_gt(ncol(with_batch$Z_raw), ncol(without_batch$Z_raw))
  expect_equal(crossprod(with_batch$Q), diag(with_batch$rank), tolerance = 1e-9)
})

test_that("build_Z_preserve validates protected covariates", {
  fx <- make_lc_fixture(n_per_batch = 6L, p = 3L)

  expect_error(
    LatentComBat:::build_Z_preserve(
      stats::lm, fx$formula, fx$covar, fx$bat,
      protected_covar = data.frame(z = rep(1, nrow(fx$data) - 1L))
    ),
    "same number of rows"
  )

  bad <- data.frame(z = rep(1, nrow(fx$data)))
  bad$z[1] <- NA
  expect_error(
    LatentComBat:::build_Z_preserve(
      stats::lm, fx$formula, fx$covar, fx$bat,
      protected_covar = bad
    ),
    "missing values"
  )
})

test_that("orthogonalize_latent removes protected-space projection", {
  set.seed(3)
  Z <- cbind(1, seq(-1, 1, length.out = 20))
  Q <- qr.Q(qr(Z))[, 1:2, drop = FALSE]
  H <- matrix(rnorm(60), 20, 3)

  H_orth <- LatentComBat:::orthogonalize_latent(H, Q)
  expect_equal(crossprod(Q, H_orth), matrix(0, 2, 3), tolerance = 1e-10)
  expect_equal(dim(H_orth), dim(H))
})
