test_that("combine_harm_draw_inference uses within plus between variance", {
  x <- data.frame(
    feature = rep("g1", 3),
    term = rep("x", 3),
    estimate = c(1, 2, 3),
    se = c(1, 1, 1)
  )

  out <- combine_harm_draw_inference(x)

  expect_equal(out$estimate, 2)
  expect_equal(out$within_var, 1)
  expect_equal(out$between_var, 1)
  expect_equal(out$total_var, 2)
  expect_equal(out$se, sqrt(2))
  expect_equal(out$harmonization_fraction, 0.5)
})


test_that("get_harm_draws retrieves retained Bayesian draws", {
  draws <- array(rnorm(3 * 10 * 2), dim = c(3, 10, 2))
  fit <- list(
    posterior_harmonization = list(harm_data_draws = draws)
  )
  class(fit) <- c("latentcombat_bayesian", "latentcombat_fit", "list")

  expect_equal(get_harm_draws(fit), draws)
})


test_that("Bayesian harm_inference propagates posterior harmonization uncertainty", {
  set.seed(10)
  n <- 30
  x <- rnorm(n)
  Y0 <- cbind(g1 = 0.7 * x + rnorm(n), g2 = -0.4 * x + rnorm(n))
  draws <- array(NA_real_, dim = c(4, n, 2))
  for (s in 1:4) draws[s, , ] <- Y0 + matrix(rnorm(n * 2, sd = 0.1), n, 2)

  fit <- list(
    harm_data = apply(draws, c(2, 3), mean),
    posterior_harmonization = list(harm_data_draws = draws),
    inference_data = list(
      covar = data.frame(x = x),
      formula = y ~ x,
      model = stats::lm
    ),
    method = "bayesian"
  )
  class(fit) <- c("latentcombat_bayesian", "latentcombat_fit", "list")

  out <- harm_inference(fit, terms = "x")

  expect_s3_class(out, "latentcombat_inference")
  expect_equal(nrow(out), 2L)
  expect_true(all(out$n_draws == 4L))
  expect_true(all(out$total_var >= out$within_var))
  expect_true(all(out$inference_method == "posterior_propagation"))
})


test_that("Sequential harm_inference uses full-pipeline bootstrap estimates", {
  set.seed(11)
  n <- 40
  x <- rnorm(n)
  Y <- cbind(g1 = 0.5 * x + rnorm(n), g2 = -0.2 * x + rnorm(n))
  bat <- factor(rep(c("A", "B"), each = n / 2))

  fit <- list(
    harm_data = Y,
    inference_data = list(
      data = Y,
      bat = bat,
      covar = data.frame(x = x),
      model = stats::lm,
      formula = y ~ x
    ),
    method = "sequential"
  )
  class(fit) <- c("latentcombat_sequential", "latentcombat_fit", "list")

  fake_refit <- function(fit, idx) {
    list(harm_data = fit$inference_data$data[idx, , drop = FALSE])
  }

  testthat::local_mocked_bindings(
    .refit_sequential_bootstrap = fake_refit,
    .package = "LatentComBat"
  )

  out <- harm_inference(fit, terms = "x", B = 10, seed = 3, progress = FALSE)

  expect_s3_class(out, "latentcombat_inference")
  expect_equal(nrow(out), 2L)
  expect_true(all(out$n_boot == 10L))
  expect_true(all(is.finite(out$se)))
  expect_true(all(out$inference_method == "full_pipeline_bootstrap"))
})
