test_that("latent mean stage skips removal when no SVs are detected", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 4L)
  br <- batch_matrix(fx$bat)
  fitted <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)

  fake_sva <- function(data_nb, Z, B, alpha, preserve_Z) {
    list(H = matrix(numeric(0), nrow(data_nb), 0L), K = 0L)
  }

  testthat::local_mocked_bindings(
    sva_per_measurement_np = fake_sva,
    .package = "LatentComBat"
  )

  expect_message(
    out <- LatentComBat:::.com_harm_sva_mean(
      data = fx$data,
      bat = fx$bat,
      covar = fx$covar,
      fitted_model = fitted,
      model = stats::lm,
      formula = fx$formula,
      B = 3L
    ),
    "found 0 hidden factors"
  )

  expect_equal(out$data, fx$data)
  expect_null(out$latent$H_orth)
  expect_null(out$latent$fitted)
  expect_true(out$protection$preserve_batch)
})

test_that("latent mean stage removes only the orthogonalized latent component", {
  fx <- make_lc_fixture(n_per_batch = 10L, p = 4L)
  br <- batch_matrix(fx$bat)
  fitted <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  n <- nrow(fx$data)

  H <- cbind(
    seq(-1, 1, length.out = n),
    sin(seq(0, 2 * pi, length.out = n))
  )

  fake_sva <- function(data_nb, Z, B, alpha, preserve_Z) list(H = H, K = 2L)

  testthat::local_mocked_bindings(
    sva_per_measurement_np = fake_sva,
    .package = "LatentComBat"
  )

  out <- suppressMessages(LatentComBat:::.com_harm_sva_mean(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    fitted_model = fitted,
    model = stats::lm,
    formula = fx$formula,
    B = 3L
  ))

  expect_equal(dim(out$latent$H), c(n, 2L))
  expect_equal(dim(out$latent$fitted), dim(fx$data))
  expect_equal(out$data, fx$data - out$latent$fitted, tolerance = 1e-10)
  expect_lt(max(abs(crossprod(out$protection$Z_raw, out$latent$H_orth))), 1e-8)
})
