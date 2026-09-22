testthat::skip_if_not_installed("cmdstanr")
test_that("univariate Stan data are flattened consistently", {
  Y <- matrix(1:12, nrow = 4, ncol = 3)
  colnames(Y) <- paste0("g", 1:3)
  ds <- list(data_stand = Y)
  bat <- factor(c("A", "A", "B", "B"))
  type <- combat_type("univariate")

  out <- stan_data_prep(type, ds, bat)
  expect_equal(out$N, 12L)
  expect_equal(out$G, 3L)
  expect_equal(out$I, 2L)
  expect_equal(out$y, as.vector(Y))
  expect_equal(out$g, rep(1:3, each = 4))
  expect_equal(out$i, rep(c(1L, 1L, 2L, 2L), times = 3))
})

test_that("univariate Stan data validates batch-vector length", {

  data_stand_result <- list(
    data_stand = matrix(rnorm(20), nrow = 5, ncol = 4)
  )

  type <- combat_type("univariate")

  expect_error(
    stan_data_prep(
      type = type,
      data_stand_result = data_stand_result,
      bat = c("A", "B")
    ),
    "one entry per observation"
  )
})

test_that("univariate Stan data uses observation-level batch vector", {

  data_stand_result <- list(
    data_stand = matrix(1:12, nrow = 4, ncol = 3)
  )

  bat <- factor(c("A", "A", "B", "B"))

  out <- stan_data_prep(
    type = combat_type("univariate"),
    data_stand_result = data_stand_result,
    bat = bat
  )

  expect_equal(out$N, 12)
  expect_equal(out$G, 3)
  expect_equal(out$I, 2)

  expect_equal(
    out$i,
    c(
      1, 1, 2, 2,  # feature 1
      1, 1, 2, 2,  # feature 2
      1, 1, 2, 2   # feature 3
    )
  )

  expect_equal(
    out$g,
    c(
      1, 1, 1, 1,
      2, 2, 2, 2,
      3, 3, 3, 3
    )
  )
})

test_that("joint Stan data validate sigma_hat", {
  R <- matrix(rnorm(30), 10, 3)
  bat <- factor(rep(c("A", "B"), each = 5))
  type <- combat_type("sva_combat")

  expect_error(
    stan_data_prep(type, R = R, bat = bat, K_bayes = 2, K_target = 1),
    "sigma_hat is required"
  )
  expect_error(
    stan_data_prep(type, R = R, bat = bat, K_bayes = 2, K_target = 1,
                   sigma_hat = c(1, 1)),
    "length P"
  )
  expect_error(
    stan_data_prep(type, R = R, bat = bat, K_bayes = 2, K_target = 1,
                   sigma_hat = c(1, 0, 1)),
    "finite and > 0"
  )
})

test_that("cross-sectional joint Stan data have expected indices", {
  R <- matrix(rnorm(30), 10, 3)
  bat <- factor(rep(c("A", "B"), each = 5))
  type <- combat_type("sva_combat")
  out <- stan_data_prep(
    type, R = R, bat = bat, K_bayes = 2L, K_target = 1L,
    sigma_hat = rep(1, 3)
  )

  expect_equal(out$N, 10L)
  expect_equal(out$P, 3L)
  expect_equal(out$B, 2L)
  expect_equal(out$K, 2L)
  expect_equal(out$K_target, 1L)
  expect_false("S" %in% names(out))
  expect_equal(out$batch, as.integer(bat))
})

test_that("longitudinal joint Stan data construct subject-by-batch indices", {
  R <- matrix(rnorm(36), 12, 3)
  bat <- factor(rep(c("A", "B"), each = 6))
  subj <- factor(rep(paste0("s", 1:6), each = 2))
  type <- combat_type("sva_combat")
  out <- stan_data_prep(
    type, R = R, bat = bat, K_bayes = 2L, K_target = 1L,
    subj = subj, sigma_hat = rep(1, 3)
  )

  expect_equal(out$N, 12L)
  expect_equal(out$S, 6L)
  expect_equal(length(out$S_uniq_b), 2L)
  expect_true(all(out$N_bs >= 0))
  expect_equal(dim(out$idx_bs)[1], 2L)
})
