test_that("batch_matrix constructs design, indices, sizes, and reference", {
  bat <- factor(c("B", "A", "B", "A", "A"), levels = c("unused", "A", "B"))
  out <- batch_matrix(bat, ref.batch = "B")

  expect_identical(levels(out$batch_vector), c("A", "B"))
  expect_equal(dim(out$batch_matrix), c(5L, 2L))
  expect_equal(out$batch_index$A, c(2L, 4L, 5L))
  expect_equal(out$batch_index$B, c(1L, 3L))
  expect_equal(unname(out$n_batches), c(3L, 2L))
  expect_identical(out$ref, out$batch_vector == "B")
})

test_that("batch_matrix rejects an unknown reference batch", {
  expect_error(
    batch_matrix(factor(c("A", "B")), ref.batch = "C"),
    "Reference batch"
  )
})

test_that("batch_matrix accepts character batch labels", {
  character_batch <- batch_matrix(c("B", "A", "B"))
  factor_batch <- batch_matrix(factor(c("B", "A", "B")))

  expect_identical(character_batch$batch_vector, factor_batch$batch_vector)
  expect_equal(character_batch$batch_matrix, factor_batch$batch_matrix)
  expect_identical(character_batch$batch_index, factor_batch$batch_index)
})

test_that("model_fitting fits one model per feature", {
  fx <- make_lc_fixture()
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(
    data = fx$data,
    batch = br$batch_matrix,
    covar = fx$covar,
    model = stats::lm,
    formula = fx$formula
  )

  expect_length(fit$fits, ncol(fx$data))
  expect_identical(names(fit$fits), colnames(fx$data))
  expect_true(all(vapply(fit$fits, inherits, logical(1), what = "lm")))
  expect_equal(dim(fit$data), dim(fx$data))
})

test_that("model_fitting validates formula and weights", {
  fx <- make_lc_fixture()
  br <- batch_matrix(fx$bat)

  expect_error(
    model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, formula = NULL),
    "provide a formula"
  )
  expect_error(
    model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm,
                  formula = fx$formula, weights = rep(1, nrow(fx$data) - 1L)),
    "Weight vector"
  )
  expect_error(
    model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm,
                  formula = fx$formula,
                  weights = matrix(1, nrow(fx$data), ncol(fx$data) - 1L)),
    "same number of columns"
  )
})

test_that("model_fitting supports vector and feature-specific weights", {
  fx <- make_lc_fixture()
  br <- batch_matrix(fx$bat)
  w <- seq(1, 2, length.out = nrow(fx$data))

  fit_vec <- model_fitting(
    fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula,
    weights = w
  )
  fit_mat <- model_fitting(
    fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula,
    weights = matrix(rep(w, ncol(fx$data)), nrow = nrow(fx$data))
  )

  expect_length(fit_vec$fits, ncol(fx$data))
  expect_length(fit_mat$fits, ncol(fx$data))
})

test_that("model_fitting handles a null covariate model", {
  fx <- make_lc_fixture()
  br <- batch_matrix(fx$bat)

  expect_warning(
    fit <- model_fitting(fx$data, br$batch_matrix, NULL, stats::lm),
    "No covariates"
  )
  expect_length(fit$fits, ncol(fx$data))
})

test_that("standardize_data returns conformable finite matrices", {
  fx <- make_lc_fixture()
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  st <- standardize_data(fit, br)

  expect_finite_matrix(st$data_stand, nrow(fx$data), ncol(fx$data))
  expect_finite_matrix(st$stand_mean, nrow(fx$data), ncol(fx$data))
  expect_finite_matrix(st$sd_mat, nrow(fx$data), ncol(fx$data))
  expect_true(all(st$sd_mat > 0))
})

test_that("reference-batch standardization uses reference residual variance", {
  fx <- make_lc_fixture(n_per_batch = 15L)
  br <- batch_matrix(fx$bat, ref.batch = "A")
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  st <- standardize_data(fit, br)

  expect_equal(dim(st$data_stand), dim(fx$data))
  expect_true(all(is.finite(st$data_stand)))
})
