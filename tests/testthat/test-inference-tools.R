test_that("inference_table augments posterior-propagation results", {
  x <- data.frame(
    feature = paste0("g", 1:4),
    term = "diagnosis",
    estimate = c(0.30, 0.20, 0.05, 0.00),
    se = c(0.12, 0.12, 0.10, 0.10),
    p_value = 2 * pnorm(-abs(c(0.30, 0.20, 0.05, 0.00) / c(0.12, 0.12, 0.10, 0.10))),
    plugin_estimate = c(0.30, 0.20, 0.05, 0.00),
    plugin_se = c(0.08, 0.09, 0.09, 0.09),
    harmonization_fraction = c(0.25, 0.15, 0.05, 0.00),
    inference_method = "posterior_propagation",
    harmonization_method = "bayesian"
  )
  class(x) <- c("latentcombat_inference", "data.frame")

  tab <- inference_table(x, adjust_method = "BH")

  expect_equal(tab$reference_se, x$plugin_se)
  expect_equal(tab$se_ratio, x$se / x$plugin_se)
  expect_true(all(c("p_adjust", "reference_p_adjust") %in% names(tab)))
  expect_identical(attr(tab, "reference_label"), "plug-in")
})


test_that("inference_table uses conditional SE for sequential results", {
  x <- data.frame(
    feature = c("g1", "g2"),
    term = "diagnosis",
    estimate = c(0.2, 0.1),
    se = c(0.12, 0.11),
    conditional_se = c(0.08, 0.09),
    p_value = 2 * pnorm(-abs(c(0.2, 0.1) / c(0.12, 0.11))),
    bootstrap_bias = c(0.01, -0.02),
    inference_method = "full_pipeline_bootstrap",
    harmonization_method = "sequential"
  )
  class(x) <- c("latentcombat_inference", "data.frame")

  tab <- inference_table(x)

  expect_equal(tab$reference_estimate, x$estimate)
  expect_equal(tab$reference_se, x$conditional_se)
  expect_identical(attr(tab, "reference_label"), "conditional")
})


test_that("adjustment is performed separately by term", {
  x <- data.frame(
    feature = rep(c("g1", "g2"), 2),
    term = rep(c("age", "diagnosis"), each = 2),
    estimate = 1,
    se = 1,
    p_value = c(0.01, 0.20, 0.01, 0.20),
    plugin_estimate = 1,
    plugin_se = 1,
    inference_method = "posterior_propagation",
    harmonization_method = "bayesian"
  )
  class(x) <- c("latentcombat_inference", "data.frame")

  tab <- inference_table(x, adjust_method = "BH")

  expect_equal(tab$p_adjust[tab$term == "age"], c(0.02, 0.20))
  expect_equal(tab$p_adjust[tab$term == "diagnosis"], c(0.02, 0.20))
})


test_that("summary reports discovery transitions", {
  x <- data.frame(
    feature = c("g1", "g2"),
    term = "diagnosis",
    estimate = c(0.30, 0.02),
    se = c(0.20, 0.10),
    p_value = c(0.14, 0.84),
    plugin_estimate = c(0.30, 0.02),
    plugin_se = c(0.05, 0.10),
    inference_method = "posterior_propagation",
    harmonization_method = "bayesian"
  )
  class(x) <- c("latentcombat_inference", "data.frame")

  s <- summary(x, adjust_method = "none")

  expect_s3_class(s, "summary_latentcombat_inference")
  expect_equal(s$term_summary$lost, 1)
  expect_equal(s$term_summary$propagated_discoveries, 0)
})


test_that("method-specific plots fail informatively", {
  skip_if_not_installed("ggplot2")

  x <- data.frame(
    feature = "g1",
    term = "diagnosis",
    estimate = 0,
    se = 1,
    conditional_se = 1,
    p_value = 1,
    bootstrap_bias = 0,
    inference_method = "full_pipeline_bootstrap",
    harmonization_method = "sequential"
  )
  class(x) <- c("latentcombat_inference", "data.frame")

  expect_error(
    plot(x, type = "harmonization"),
    "No finite harmonization fractions"
  )
  expect_s3_class(plot(x, type = "bootstrap_bias"), "ggplot")
})
