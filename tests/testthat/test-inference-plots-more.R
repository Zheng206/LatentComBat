make_plot_inference_fixture <- function() {
  x <- data.frame(
    feature = rep(paste0("g", 1:4), 2),
    term = rep(c("age", "diagnosis"), each = 4),
    estimate = c(0.25, 0.10, -0.08, 0.04, 0.40, 0.20, 0.05, -0.10),
    se = c(0.12, 0.10, 0.11, 0.09, 0.15, 0.14, 0.12, 0.13),
    conditional_se = c(0.08, 0.08, 0.09, 0.08, 0.10, 0.10, 0.10, 0.11),
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    bootstrap_bias = c(0.01, -0.02, 0.00, 0.01, 0.03, -0.01, 0.02, -0.02),
    harmonization_fraction = c(0.25, 0.15, 0.10, 0.05, 0.30, 0.20, 0.12, 0.08),
    inference_method = "full_pipeline_bootstrap",
    harmonization_method = "sequential",
    stringsAsFactors = FALSE
  )
  x$p_value <- 2 * stats::pnorm(-abs(x$estimate / x$se))
  zcrit <- stats::qnorm(0.975)
  x$conf_low <- x$estimate - zcrit * x$se
  x$conf_high <- x$estimate + zcrit * x$se
  class(x) <- c("latentcombat_inference", "data.frame")
  x
}

test_that("all inference plot types return ggplot objects", {
  skip_if_not_installed("ggplot2")
  x <- make_plot_inference_fixture()

  for (type in c(
    "estimate", "estimate_ci", "forest", "se", "se_ratio",
    "pvalue", "significance", "harmonization", "bootstrap_bias"
  )) {
    p <- plot(x, type = type)
    expect_s3_class(p, "ggplot")
  }
})

test_that("forest plot supports propagated and reference inference", {
  skip_if_not_installed("ggplot2")
  x <- make_plot_inference_fixture()

  expect_s3_class(
    plot(x, type = "forest", inference = "propagated", order_by = "abs_estimate"),
    "ggplot"
  )
  expect_s3_class(
    plot(x, type = "forest", inference = "reference", order_by = "pvalue", adjusted = FALSE),
    "ggplot"
  )
})

test_that("inference plot filtering accepts valid terms and features", {
  skip_if_not_installed("ggplot2")
  x <- make_plot_inference_fixture()

  p1 <- plot(x, type = "estimate", term = "age")
  p2 <- plot(x, type = "se", feature = c("g1", "g2"))

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
})

test_that("inference plotting validates bins and unknown filters", {
  skip_if_not_installed("ggplot2")
  x <- make_plot_inference_fixture()

  expect_error(plot(x, type = "harmonization", bins = 0), "positive integer")
  expect_error(plot(x, type = "estimate", term = "not_a_term"), "term|Term")
  expect_error(plot(x, type = "estimate", feature = "not_a_feature"), "feature|Feature")
})

test_that("plot helper labels and ordering are stable", {
  expect_equal(LatentComBat:::.title_case("plug-in"), "Plug-in")
  expect_equal(LatentComBat:::.forest_order_label("estimate"), "effect estimate")
  expect_equal(LatentComBat:::.forest_order_label("abs_estimate"), "absolute effect size")
})
