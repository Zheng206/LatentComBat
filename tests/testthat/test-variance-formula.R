test_that("variance formula is applied in the log-variance fit", {
  fx <- make_lc_fixture(n_per_batch = 24L, p = 4L)
  br <- batch_matrix(fx$bat)
  design <- model.matrix(~ splines::ns(x, df = 3) + sex, data = fx$covar)
  design <- cbind(model.matrix(~ fx$bat - 1), design[, -1, drop = FALSE])

  fake_select <- function(Z_resid, k_max, B_perm, quant, seed,
                          X_batch = NULL) {
    expect_equal(dim(X_batch), dim(design))
    list(k = 1L, d_obs = 1, thr = 0)
  }
  testthat::local_mocked_bindings(
    select_k_var_parallel = fake_select, .package = "LatentComBat"
  )

  out <- suppressMessages(LatentComBat:::.stabilize_latent_variance(
    data = fx$data, bat = fx$bat, covar = fx$covar,
    batch_result = br, model = stats::lm, formula = fx$formula,
    variance_formula = ~ splines::ns(x, df = 3) + sex
  ))

  expect_equal(out$K_var, 1L)
  expect_equal(dim(out$ell_bio), dim(fx$data))
  expect_equal(nrow(out$coef_bio), ncol(design) - 2L)
  expect_lt(max(abs(crossprod(design, out$ell_latent))), 1e-6)
  expect_lte(max(abs(out$ell_latent)), 2 + 1e-8)
  expect_true(all(is.finite(out$data)))
})

test_that("variance formula is rejected when the variance stage is off", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 3L)
  expect_error(com_harm(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula,
    var_stable = FALSE, variance_formula = ~ x
  ), "var_stable")
})
