test_that("variance_step_quasi decomposes pooled and batch log variance", {
  set.seed(1)
  bat <- factor(rep(c("A", "B"), each = 20))
  e <- cbind(
    c(rnorm(20, sd = 1), rnorm(20, sd = 2)),
    c(rnorm(20, sd = 0.7), rnorm(20, sd = 1.5))
  )

  out <- LatentComBat:::variance_step_quasi(e, bat)
  expect_equal(dim(out$ell), dim(e))
  expect_equal(dim(out$theta_bg), c(2L, 2L))
  expect_equal(out$log_delta_bg, 0.5 * out$theta_bg)
  expect_true(is.finite(out$eps) && out$eps > 0)

  w <- as.numeric(table(bat)) / length(bat)
  expect_equal(as.numeric(crossprod(w, out$theta_bg)), c(0, 0), tolerance = 1e-10)

  b_int <- as.integer(bat)
  expected_ell <- matrix(rep(out$c_g, each = nrow(e)), nrow(e), ncol(e)) +
    out$theta_bg[b_int, , drop = FALSE]
  expect_equal(out$ell, expected_ell)
})

test_that("variance_step_quasi handles zero residuals with a positive offset", {
  e <- matrix(0, 8, 3)
  bat <- factor(rep(c("A", "B"), each = 4))
  out <- LatentComBat:::variance_step_quasi(e, bat)

  expect_equal(out$eps, 1e-8)
  expect_true(all(is.finite(out$ell)))
})

test_that("em_var_stabilize_only follows the damped fixed-mean update", {
  fx <- make_lc_fixture(n_per_batch = 10L, p = 4L)
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)

  damp <- 0.4
  n_iter <- 3L
  out <- LatentComBat:::em_var_stabilize_only(
    data = fx$data,
    bat = fx$bat,
    fitted_model = fit,
    n_iter = n_iter,
    damp = damp
  )

  expect_equal(dim(out$ell), dim(fx$data))
  expect_equal(dim(out$mu_hat), dim(fx$data))
  multiplier <- 1 - (1 - damp)^n_iter
  expect_equal(out$ell, multiplier * out$vfit$ell, tolerance = 1e-10)
})

test_that("biweight_midvar is less sensitive to a gross outlier", {
  set.seed(2)
  x <- c(rnorm(100), 30)
  expect_lt(biweight_midvar(x), stats::var(x))
})
