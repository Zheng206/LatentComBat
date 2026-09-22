test_that("compute_cov_matrix matches base cor and cov", {
  x <- matrix(c(1, 2, 3, 4, 2, 1, 4, 3, 5, 6, 7, 8), ncol = 3)

  expect_equal(
    LatentComBat:::compute_cov_matrix(x, "correlation"),
    stats::cor(x, use = "pairwise.complete.obs")
  )
  expect_equal(
    LatentComBat:::compute_cov_matrix(x, "covariance"),
    stats::cov(x, use = "pairwise.complete.obs")
  )
})

test_that("batch covariance distances return one row per batch pair", {
  set.seed(101)
  bat <- factor(rep(c("A", "B", "C"), each = 8))
  y <- matrix(rnorm(24 * 4), 24, 4)
  y[bat == "B", ] <- y[bat == "B", ] * c(1, 1.5, 0.8, 1.2)

  out <- batch_covariance_distance(y, bat, type = "covariance", metric = "frobenius")

  expect_equal(nrow(out), choose(3, 2))
  expect_named(out, c("batch1", "batch2", "distance"))
  expect_true(all(is.finite(out$distance)))
  expect_true(all(out$distance >= 0))
})

test_that("batch covariance distance validates batch input", {
  y <- matrix(rnorm(20), 10, 2)
  expect_error(batch_covariance_distance(y, rep("A", 10)), "At least two batches")
  expect_error(batch_covariance_distance(y, rep(c("A", "B"), each = 4)), "one entry per observation")
})

test_that("covariance summaries and reduction are internally consistent", {
  set.seed(102)
  bat <- factor(rep(c("A", "B"), each = 15))
  raw <- matrix(rnorm(30 * 4), 30, 4)
  raw[bat == "B", ] <- sweep(raw[bat == "B", , drop = FALSE], 2, c(1, 2, .5, 1.5), "*")
  harm <- raw
  harm[bat == "B", ] <- sweep(harm[bat == "B", , drop = FALSE], 2, c(1, 2, .5, 1.5), "/")

  sm <- summarize_batch_covariance(raw, bat, type = "covariance", metric = "mse", relative = FALSE)
  red <- covariance_batch_reduction(raw, harm, bat, type = "covariance", metric = "mse", relative = FALSE)

  expect_equal(sm$mean_distance, mean(sm$pairwise$distance))
  expect_equal(sm$median_distance, stats::median(sm$pairwise$distance))
  expect_named(red, c("before", "after", "relative_reduction"))
  expect_true(all(is.finite(red)))
})

test_that("structure perturbation is zero when data are unchanged", {
  set.seed(103)
  y <- matrix(rnorm(80), 20, 4)
  out <- evaluate_structure_perturbation(y, y, type = "correlation")

  expect_equal(out$raw, out$harmonized)
  expect_equal(out$difference, matrix(0, 4, 4), tolerance = 1e-12)
  expect_equal(out$Frobenius, 0, tolerance = 1e-12)
  expect_equal(out$MSE, 0, tolerance = 1e-12)
  expect_equal(out$Spectral, 0, tolerance = 1e-12)
  expect_equal(out$EigenError, 0, tolerance = 1e-12)
})

test_that("structure perturbation requires conformable matrices", {
  expect_error(
    evaluate_structure_perturbation(matrix(1, 4, 2), matrix(1, 5, 2)),
    "same dimensions"
  )
})

test_that("Box M returns a valid test and checks group sizes", {
  set.seed(104)
  y <- matrix(rnorm(60), 20, 3)
  g <- factor(rep(c("A", "B"), each = 10))
  out <- boxM_simple(y, g)

  expect_equal(out$method, "Box's M test")
  expect_true(is.finite(out$statistic))
  expect_true(is.finite(out$parameter))
  expect_true(out$p.value >= 0 && out$p.value <= 1)

  expect_error(boxM_simple(y, rep("A", 20)), "At least two groups")
  expect_error(boxM_simple(y[1:3, , drop = FALSE], c("A", "B", "B")), "at least two observations")
})

test_that("PCA batch preparation obeys requested and automatic dimensions", {
  set.seed(105)
  y <- matrix(rnorm(36 * 6), 36, 6)
  bat <- factor(rep(c("A", "B", "C"), each = 12))

  auto <- LatentComBat:::.prepare_batch_test_pca(y, bat, var_explained = 0.8)
  fixed <- LatentComBat:::.prepare_batch_test_pca(y, bat, n_pc = 2)

  expect_true(auto$n_pc >= 1L && auto$n_pc <= auto$max_pc)
  expect_equal(fixed$n_pc, 2L)
  expect_equal(dim(fixed$scores), c(36L, 2L))
  expect_equal(dim(fixed$loadings), c(6L, 2L))
  expect_error(LatentComBat:::.prepare_batch_test_pca(y, bat, n_pc = 0), "positive integer")
})

test_that("batch structure tests return expected classes and components", {
  set.seed(106)
  bat <- factor(rep(c("A", "B"), each = 20))
  raw <- matrix(rnorm(40 * 5), 40, 5)
  raw[bat == "B", 1:2] <- raw[bat == "B", 1:2] + 0.8
  harm <- raw
  harm[bat == "B", 1:2] <- harm[bat == "B", 1:2] - 0.8

  one <- test_batch_structure(raw, bat, n_pc = 2)
  cmp <- compare_batch_structure(raw, harm, bat, n_pc = 2)

  expect_s3_class(one, "batch_structure_test")
  expect_equal(one$n_pc, 2L)
  expect_named(one, c("n_pc", "variance_explained", "mean", "covariance", "pca"))

  expect_s3_class(cmp, "batch_structure_comparison")
  expect_equal(nrow(cmp$result), 4L)
  expect_setequal(cmp$result$stage, c("Before", "After"))
  expect_setequal(cmp$result$structure, c("Mean", "Covariance"))
  expect_equal(dim(cmp$scores_raw), dim(cmp$scores_harmonized))
})

test_that("covariance diagnostic plot functions return ggplots", {
  skip_if_not_installed("ggplot2")
  set.seed(107)
  bat <- factor(rep(c("A", "B"), each = 12))
  raw <- matrix(rnorm(24 * 4), 24, 4)
  harm <- raw + matrix(rnorm(24 * 4, sd = 0.05), 24, 4)

  expect_s3_class(plot_batch_distance(raw, harm, bat), "ggplot")
  expect_s3_class(plot_structure_preservation(raw, harm), "ggplot")
  expect_s3_class(plot_batch_covariance(raw, bat), "ggplot")
})
