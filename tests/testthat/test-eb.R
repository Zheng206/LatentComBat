test_that("combat_type creates an S3 dispatch tag", {
  x <- combat_type("univariate")
  expect_s3_class(x, "univariate")
  expect_identical(unclass(x), "univariate")
})

test_that("univariate gamma and delta estimates have batch-by-feature shape", {
  fx <- make_lc_fixture(n_per_batch = 12L, p = 5L)
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  st <- standardize_data(fit, br)
  type <- combat_type("univariate")

  g <- gamm_hat_gen(type, st, br)
  d <- delta_hat_gen(type, st, br)

  expect_equal(dim(g), c(nlevels(fx$bat), ncol(fx$data)))
  expect_equal(dim(d), c(nlevels(fx$bat), ncol(fx$data)))
  expect_identical(rownames(g), levels(fx$bat))
  expect_true(all(is.finite(g)))
  expect_true(all(is.finite(d)))
  expect_true(all(d >= 0))
})

test_that("mom_calculation returns finite EB hyperparameters on variable inputs", {
  type <- combat_type("univariate")
  out <- mom_calculation(type,
                         gamma_hat = c(-0.4, 0.1, 0.5, 0.8),
                         delta_hat = c(0.5, 0.8, 1.2, 1.7))
  expect_named(out, c("g_bar", "g_var", "d_bar", "d_var", "d_a", "d_b"))
  expect_true(all(vapply(out, is.finite, logical(1))))
})

test_that("one EB iteration returns conformable finite updates", {
  type <- combat_type("univariate")
  set.seed(4)
  bdat <- matrix(rnorm(40), 10, 4)
  g <- colMeans(bdat)
  d <- apply(bdat, 2, stats::var)
  mom <- mom_calculation(type, g + c(-.2, .1, .3, .5), d + c(.1, .2, .4, .6))

  out <- eb_one_iteration(type, bdat, g, g, d, mom)
  expect_length(out$g_new, 4L)
  expect_length(out$d_new, 4L)
  expect_true(all(is.finite(out$g_new)))
  expect_true(all(is.finite(out$d_new)))
  expect_true(is.finite(out$change))
})

test_that("EB disabled returns initial estimates exactly", {
  fx <- make_lc_fixture(n_per_batch = 12L, p = 5L)
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  st <- standardize_data(fit, br)
  type <- combat_type("univariate")

  out <- eb_algorithm(type, st, br, eb = FALSE)
  expect_equal(out$gamma_star, out$gamma_hat)
  expect_equal(out$delta_star, out$delta_hat)
  expect_null(out$mom)
})

test_that("EB enabled returns finite shrunken batch parameters", {
  fx <- make_lc_fixture(n_per_batch = 18L, p = 8L, seed = 44)
  br <- batch_matrix(fx$bat)
  fit <- model_fitting(fx$data, br$batch_matrix, fx$covar, stats::lm, fx$formula)
  st <- standardize_data(fit, br)
  type <- combat_type("univariate")

  out <- eb_algorithm(type, st, br, eb = TRUE)
  expect_equal(dim(out$gamma_star), c(2L, 8L))
  expect_equal(dim(out$delta_star), c(2L, 8L))
  expect_true(all(is.finite(out$gamma_star)))
  expect_true(all(is.finite(out$delta_star)))
  expect_true(all(out$delta_star > 0))
  expect_length(out$mom, 2L)
})
