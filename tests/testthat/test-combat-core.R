test_that(".com_harm_combat runs EB pathway", {
  fx <- make_lc_fixture(
    n_per_batch = 12L,
    p = 5L,
    seed = 100
  )

  br <- batch_matrix(fx$bat)

  out <- LatentComBat:::.com_harm_combat(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    batch_result = br,
    model = stats::lm,
    formula = fx$formula,
    eb = FALSE,
    stan = FALSE,
    cov = FALSE
  )

  expect_identical(out$engine, "eb")
  expect_equal(dim(out$harm_data), dim(fx$data))
  expect_equal(dim(out$resid), dim(fx$data))
  expect_true(all(is.finite(out$harm_data)))
  expect_true(all(is.finite(out$resid)))
})


test_that(".com_harm_combat preserves reference batch", {
  fx <- make_lc_fixture(
    n_per_batch = 12L,
    p = 4L,
    seed = 101
  )

  br <- batch_matrix(
    fx$bat,
    ref.batch = "A"
  )

  out <- LatentComBat:::.com_harm_combat(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    batch_result = br,
    model = stats::lm,
    formula = fx$formula,
    eb = FALSE,
    stan = FALSE,
    ref.batch = "A"
  )

  expect_equal(
    out$harm_data[fx$bat == "A", , drop = FALSE],
    fx$data[fx$bat == "A", , drop = FALSE],
    tolerance = 1e-12
  )
})


test_that(".com_harm_combat runs CovBat branch", {
  fx <- make_lc_fixture(
    n_per_batch = 10L,
    p = 4L,
    seed = 102
  )

  br <- batch_matrix(fx$bat)

  fake_covbat <- function(
    R,
    site,
    center,
    scale_scores,
    var_thresh,
    min_rblock,
    max_rblock,
    ref.batch
  ) {
    expect_equal(site, fx$bat)
    expect_true(center)
    expect_true(scale_scores)
    expect_equal(var_thresh, 0.90)
    expect_equal(min_rblock, 1)
    expect_equal(max_rblock, 2)

    R
  }

  testthat::local_mocked_bindings(
    covbat = fake_covbat,
    .package = "LatentComBat"
  )

  out <- LatentComBat:::.com_harm_combat(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    batch_result = br,
    model = stats::lm,
    formula = fx$formula,
    eb = FALSE,
    stan = FALSE,
    cov = TRUE,
    var_thresh = 0.90,
    min_rblock = 1,
    max_rblock = 2
  )

  expect_equal(
    dim(out$harm_data),
    dim(fx$data)
  )
})


test_that(".com_harm_combat runs Stan pathway using mocks", {
  fx <- make_lc_fixture(
    n_per_batch = 8L,
    p = 4L,
    seed = 103
  )

  br <- batch_matrix(fx$bat)

  fake_stan_data_prep <- function(
    type,
    data_stand_result,
    bat
  ) {
    list(fake = TRUE)
  }

  fake_stan_algorithm <- function(
    type,
    stan_data,
    batch_names,
    adapt_delta,
    iter_sampling,
    iter_warmup,
    max_treedepth
  ) {
    p <- ncol(fx$data)
    B <- nlevels(fx$bat)

    list(
      gamma_star = matrix(
        0,
        nrow = B,
        ncol = p,
        dimnames = list(
          levels(fx$bat),
          colnames(fx$data)
        )
      ),
      delta_star = matrix(
        1,
        nrow = B,
        ncol = p,
        dimnames = list(
          levels(fx$bat),
          colnames(fx$data)
        )
      )
    )
  }

  testthat::local_mocked_bindings(
    stan_data_prep = fake_stan_data_prep,
    stan_algorithm = fake_stan_algorithm,
    .package = "LatentComBat"
  )

  out <- LatentComBat:::.com_harm_combat(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    batch_result = br,
    model = stats::lm,
    formula = fx$formula,
    stan = TRUE
  )

  expect_identical(out$engine, "stan")
  expect_equal(dim(out$harm_data), dim(fx$data))
  expect_true(all(is.finite(out$harm_data)))
})
