test_that("covariance recovery metrics are zero for identical matrices", {
  S <- matrix(c(2, .4, .4, 1), 2, 2)
  expect_equal(frobenius_norm(S, S), 0)
  expect_equal(mse_cov(S, S), 0)
  expect_equal(spectral_norm(S, S), 0)
  expect_equal(eigen_error(S, S), 0)

  all_metrics <- evaluate_cov_recovery(S, S)
  expect_true(all(unlist(all_metrics) == 0))
})

test_that("spectral_norm uses the largest absolute eigenvalue", {
  S <- diag(2)
  Shat <- diag(c(4, 1))
  # S - Shat has eigenvalues 0 and -3, so ||S-Shat||_2 = 3.
  expect_equal(spectral_norm(S, Shat), 3)
})

test_that("pick_r_from_pc respects variance threshold and bounds", {
  pca <- list(sdev = sqrt(c(5, 3, 1, 1)))
  expect_equal(pick_r_from_pc(pca, var_thresh = 0.70, min_r = 1), 2)
  expect_equal(pick_r_from_pc(pca, var_thresh = 0.10, min_r = 3), 3)
  expect_equal(pick_r_from_pc(pca, var_thresh = 0.99, min_r = 1, max_r = 2), 2)
})
