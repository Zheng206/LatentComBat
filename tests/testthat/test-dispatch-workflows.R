test_that("com_harm validates method", {
  fx <- make_lc_fixture(n_per_batch = 5L, p = 3L)
  expect_error(
    com_harm(
      bat = fx$bat, data = fx$data, covar = fx$covar,
      formula = fx$formula, method = "not-a-method"
    ),
    "arg"
  )
})

test_that("legacy bayes_sva conflict warns when method is explicit", {
  fx <- make_lc_fixture(n_per_batch = 5L, p = 3L)
  expect_warning(
    com_harm(
      bat = fx$bat, data = fx$data, covar = fx$covar,
      formula = fx$formula, method = "sequential", bayes_sva = TRUE,
      eb = FALSE
    ),
    "overrides legacy"
  )
})

test_that("basic sequential ComBat returns a valid fit object", {
  fx <- make_lc_fixture(n_per_batch = 15L, p = 6L)
  out <- com_harm(
    bat = fx$bat,
    data = fx$data,
    covar = fx$covar,
    model = stats::lm,
    formula = fx$formula,
    method = "sequential",
    eb = FALSE,
    sva = FALSE,
    var_stable = FALSE
  )

  expect_s3_class(out, "latentcombat_sequential")
  expect_s3_class(out, "latentcombat_fit")
  expect_identical(out$method, "sequential")
  expect_equal(dim(out$harm_data), dim(fx$data))
  expect_null(out$latent)
  expect_null(out$variance)
  expect_true("eb_result" %in% names(out))
})

test_that("sequential reference batch is unchanged", {
  fx <- make_lc_fixture(n_per_batch = 15L, p = 5L)
  out <- com_harm(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula,
    ref.batch = "A", eb = FALSE
  )
  expect_equal(
    out$harm_data[fx$bat == "A", , drop = FALSE],
    fx$data[fx$bat == "A", , drop = FALSE],
    tolerance = 1e-10
  )
})

test_that("var_stable branch is behaviorally independent of sva", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 4L)
  fake_variance <- function(data, bat, covar, batch_result, model, formula, ...) {
    list(data = as.matrix(data) + 10, marker = "variance-ran")
  }
  fake_combat <- function(data, ...) {
    list(
      harm_data = as.matrix(data),
      resid = as.matrix(data),
      data_stand_result = list(data_stand = as.matrix(data)),
      fit = list(), engine = "eb"
    )
  }

  testthat::local_mocked_bindings(
    .stabilize_latent_variance = fake_variance,
    .com_harm_combat = fake_combat,
    .package = "LatentComBat"
  )

  out <- com_harm.sequential(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula,
    sva = FALSE, var_stable = TRUE
  )

  expect_equal(out$harm_data, fx$data + 10)
  expect_identical(out$variance$marker, "variance-ran")
  expect_null(out$latent)
})

test_that("when both branches are enabled SVA output feeds variance stabilization", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 4L)

  fake_sva <- function(data, ...) {
    list(
      data = as.matrix(data) + 1,
      latent = list(marker = "sva-ran"),
      protection = list(rank = 1L)
    )
  }
  fake_variance <- function(data, ...) {
    expect_equal(data, fx$data + 1)
    list(data = as.matrix(data) * 2, marker = "variance-ran")
  }
  fake_combat <- function(data, ...) {
    list(
      harm_data = as.matrix(data), resid = as.matrix(data),
      data_stand_result = list(), fit = list(), engine = "eb"
    )
  }

  testthat::local_mocked_bindings(
    .com_harm_sva_mean = fake_sva,
    .stabilize_latent_variance = fake_variance,
    .com_harm_combat = fake_combat,
    .package = "LatentComBat"
  )

  out <- com_harm.sequential(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula,
    sva = TRUE, var_stable = TRUE
  )

  expect_equal(out$harm_data, (fx$data + 1) * 2)
  expect_identical(out$latent$marker, "sva-ran")
  expect_identical(out$variance$marker, "variance-ran")
})

test_that("Bayesian workflow can be tested without running CmdStan", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 4L)

  fake_detect <- function(R, Z, B, alpha) {
    list(K = 2L, pvals = c(0.01, 0.02), svd = svd(R))
  }
  fake_stan <- function(type, stan_data, ...) {
    list(
      R_adj = stan_data$R,
      H_hat = matrix(0, stan_data$N, stan_data$K),
      latent_hat = matrix(0, stan_data$N, stan_data$P),
      fit = structure(list(), class = "fake_cmdstan")
    )
  }

  testthat::local_mocked_bindings(
    detect_K_one = fake_detect,
    stan_algorithm = fake_stan,
    .package = "LatentComBat"
  )

  expect_message(
    out <- com_harm.bayesian(
      bat = fx$bat, data = fx$data, covar = fx$covar,
      model = stats::lm, formula = fx$formula,
      K_max = 4L, K_buffer = 1L
    ),
    "Harmonization Complete"
  )

  expect_s3_class(out, "latentcombat_bayesian")
  expect_equal(out$latent$K_target, 2L)
  expect_equal(out$latent$K_bayes, 3L)
  expect_false(out$protection$preserve_batch)
  expect_equal(dim(out$harm_data), dim(fx$data))
})
