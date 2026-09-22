test_that("diagnostic model generation and lm summary return expected residuals", {
  fx <- make_lc_fixture(n_per_batch = 15L, p = 5L)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )
  sm <- diag_model_summary(dm)

  expect_s3_class(dm, "lm_model")
  expect_equal(dim(sm$resid_add), dim(fx$data))
  expect_equal(dim(sm$resid_mul), dim(fx$data))
})

test_that("batch-only diagnostic summaries retain the observed data", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 3L)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = NULL,
    model = stats::lm
  )
  sm <- diag_model_summary(dm)

  expect_equal(sm$resid_add, fx$data)
  expect_equal(dim(sm$resid_mul), dim(fx$data))
})

test_that("PCA preparation forwards batch controls and plot ellipse option", {
  fx <- make_lc_fixture(n_per_batch = 10L, p = 4L)
  adjusted <- pca_prep(
    fx$bat, fx$data, fx$covar, stats::lm, fx$formula,
    ref.batch = "B"
  )
  unadjusted <- pca_prep(
    fx$bat, fx$data, fx$covar, stats::lm, fx$formula,
    bat_adjust = FALSE
  )

  expect_false(isTRUE(all.equal(adjusted$F_t, unadjusted$F_t)))

  plot_with_ellipse <- pca_plot(adjusted, ellipse = TRUE)
  plot_without_ellipse <- pca_plot(adjusted, ellipse = FALSE)
  is_ellipse <- function(layer) inherits(layer$stat, "StatEllipse")
  expect_equal(sum(vapply(plot_with_ellipse$layers, is_ellipse, logical(1))), 1L)
  expect_equal(sum(vapply(plot_without_ellipse$layers, is_ellipse, logical(1))), 0L)
})

test_that("feature-wise nonparametric variance diagnostics return valid tables", {
  set.seed(8)
  bat <- factor(rep(c("A", "B"), each = 20))
  R <- cbind(
    c(rnorm(20, sd = 1), rnorm(20, sd = 2)),
    c(rnorm(20), rnorm(20)),
    c(rnorm(20, sd = 0.5), rnorm(20, sd = 1.5))
  )
  colnames(R) <- paste0("f", 1:3)

  for (fun in list(kruskal_test, lv_test, bl_test, fk_test)) {
    out <- fun(R, bat)
    expect_equal(nrow(out$test_table), 3L)
    expect_match(out$perc.sig, "%$")
  }
})

test_that("anova_test detects the batch term in diagnostic models", {
  fx <- make_lc_fixture(n_per_batch = 20L, p = 5L, seed = 88)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )
  out <- anova_test(dm)

  expect_equal(nrow(out$test_table), ncol(fx$data))
  expect_true(all(as.numeric(out$test_table$p_value) >= 0))
  expect_match(out$perc.sig, "%$")
})

test_that("uni_test validates input dimensions", {
  fx <- make_lc_fixture(n_per_batch = 5L, p = 3L)
  expect_error(
    uni_test(fx$bat[-1], fx$data, fx$covar, stats::lm, fx$formula),
    "same number of observations"
  )
  expect_error(
    uni_test(fx$bat, as.list(fx$data), fx$covar, stats::lm, fx$formula),
    "data.frame or matrix"
  )
})

test_that("uni_test returns applicable lm diagnostics", {
  fx <- make_lc_fixture(n_per_batch = 18L, p = 4L)
  out <- uni_test(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c("ANOVA", "Kruskal", "Levene", "Bartlett", "FlignerKilleen") %in% names(out)))
  expect_false("KenwardRoger" %in% names(out))
  expect_equal(nrow(out), 1L)
})
