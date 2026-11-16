test_that("technical partial R2 compares nested feature-wise models", {
  set.seed(12)
  n <- 40L
  technical <- data.frame(
    scanner = factor(sample(c("A", "B"), n, replace = TRUE)),
    field_strength = factor(sample(c("low", "high"), n, replace = TRUE))
  )
  covar <- data.frame(age = rnorm(n), sex = factor(sample(c("F", "M"), n, replace = TRUE)))
  bat <- factor(rep(c("site1", "site2", "site3", "site4"), each = 10L))
  data <- cbind(
    feature_1 = 2 * (technical$field_strength == "high") + rnorm(n),
    feature_2 = rnorm(n)
  )

  out <- technical_partial_r2(
    data = data,
    technical = technical,
    covar = covar,
    bat = bat,
    formula = y ~ age + sex
  )

  expect_equal(nrow(out), ncol(data) * ncol(technical))
  expect_true(all(c("Feature", "Technical", "partial_R2", "rank_deficient") %in% names(out)))
  expect_true(any(out$partial_R2[out$Technical == "field_strength"] > 0.5))
  expect_false(any(out$rank_deficient))
})

test_that("technical partial R2 flags collinear technical designs", {
  n <- 20L
  technical <- data.frame(scanner = factor(rep(c("A", "B"), each = 10L)))
  data <- matrix(rnorm(n * 2L), ncol = 2L)
  bat <- technical$scanner

  out <- technical_partial_r2(data, technical, bat = bat)

  expect_true(all(out$rank_deficient))
  expect_true(all(is.na(out$partial_R2)))
})

test_that("technical partial R2 supports GAMs and both plot types", {
  skip_if_not_installed("mgcv")
  set.seed(21)
  n <- 30L
  technical <- data.frame(scanner = factor(sample(c("A", "B"), n, replace = TRUE)))
  covar <- data.frame(age = runif(n, 18, 80))
  data <- cbind(feature_1 = sin(covar$age / 10) + rnorm(n, sd = 0.1))
  bat <- factor(rep(c("site1", "site2", "site3"), each = 10L))

  out <- technical_partial_r2(
    data, technical, covar, bat, model = mgcv::gam,
    formula = y ~ s(age, k = 5)
  )

  expect_true(is.numeric(out$partial_R2))
  expect_s3_class(plot_technical_r2(out, out, type = "aggregate"), "ggplot")
  expect_s3_class(plot_technical_r2(out, out, type = "feature"), "ggplot")
})
test_that("partial R2 matches direct OLS and reuses each full fit", {
  set.seed(52)
  n <- 60L
  covar <- data.frame(age = rnorm(n))
  technical <- data.frame(scanner = rnorm(n), field = rnorm(n))
  bat <- rep(1:3, each = 20)
  y <- covar$age + technical$scanner + rnorm(n)
  calls <- 0L
  counted_lm <- function(...) {
    calls <<- calls + 1L
    stats::lm(...)
  }
  out <- technical_partial_r2(cbind(a = y, b = y + 1), technical, covar, bat,
                              model = counted_lm, formula = y ~ age)
  frame <- data.frame(y, covar, technical, batch = factor(bat))
  rss_full <- sum(residuals(lm(y ~ age + batch + scanner + field, frame))^2)
  rss_reduced <- sum(residuals(lm(y ~ age + batch + field, frame))^2)
  expect_equal(out$partial_R2[1], 1 - rss_full / rss_reduced)
  expect_equal(calls, 6L)
  expect_true(all(out$rank_added == 1L))
})

test_that("nuisance redundancy does not invalidate an estimable technical effect", {
  set.seed(53)
  n <- 60L
  age <- rnorm(n)
  covar <- data.frame(age, duplicate_age = age)
  technical <- data.frame(scanner = rnorm(n))
  out <- technical_partial_r2(cbind(y = age + technical$scanner + rnorm(n)),
                              technical, covar, rep(1:3, each = 20),
                              formula = y ~ age + duplicate_age)
  expect_true(out$rank_deficient_full)
  expect_true(out$rank_deficient_reduced)
  expect_true(out$estimable)
  expect_equal(out$rank_added, 1L)
  expect_true(is.finite(out$partial_R2))
})

test_that("batch adjustment exposes confounding and additive batch invariance", {
  set.seed(54)
  n <- 60L
  bat <- factor(rep(1:3, each = 20))
  technical <- data.frame(scanner = rnorm(n))
  y <- technical$scanner + rnorm(n)
  shifted <- y + c(3, -4, 8)[as.integer(bat)]
  before <- technical_partial_r2(cbind(y), technical, bat = bat)
  after <- technical_partial_r2(cbind(shifted), technical, bat = bat)
  expect_equal(before$partial_R2, after$partial_R2)
  confounded <- data.frame(scanner = bat)
  conditional <- technical_partial_r2(cbind(shifted), confounded, bat = bat)
  marginal <- technical_partial_r2(cbind(shifted), confounded, adjust_batch = FALSE)
  expect_false(conditional$estimable)
  expect_equal(conditional$rank_added, 0L)
  expect_equal(conditional$status, "non_estimable")
  expect_true(is.na(conditional$partial_R2))
  expect_true(marginal$partial_R2 > 0.8)
  expect_false(attr(marginal, "adjust_batch"))
})

test_that("complete cases use only full-model variables and the same sample", {
  set.seed(55)
  n <- 40L
  technical <- data.frame(scanner = rnorm(n), field = rnorm(n))
  technical$field[1] <- NA_real_
  covar <- data.frame(age = rnorm(n), unused = rep(NA_real_, n))
  covar$age[2] <- NA_real_
  bat <- rep(1:2, each = 20)
  data <- cbind(a = rnorm(n), b = rnorm(n))
  data[3, "a"] <- NA_real_
  out <- technical_partial_r2(data, technical, covar, bat, formula = y ~ age)
  expect_equal(out$n_obs, c(37L, 37L, 38L, 38L))
  frame <- data.frame(y = data[, "a"], technical, age = covar$age, batch = factor(bat))
  frame <- frame[complete.cases(frame), ]
  full <- lm(y ~ age + batch + scanner + field, frame)
  reduced <- lm(y ~ age + batch + field, frame)
  expect_equal(out$partial_R2[1], 1 - sum(residuals(full)^2) / sum(residuals(reduced)^2))
  # Batch missingness is irrelevant when batch is not adjusted for.
  bat[] <- NA
  unadjusted <- technical_partial_r2(data, technical, covar, bat, formula = y ~ age,
                                    adjust_batch = FALSE)
  expect_equal(unadjusted$n_obs, out$n_obs)
  expect_true(all(is.finite(unadjusted$partial_R2)))
})

test_that("failed fits retain errors without implying rank deficiency", {
  technical <- data.frame(scanner = rep(1:2, 5))
  broken <- function(...) stop("Specific fitting problem")
  out <- technical_partial_r2(cbind(a = seq_len(10)), technical,
                              adjust_batch = FALSE, model = broken)
  expect_equal(out$status, "fit_error")
  expect_match(out$error, "Specific fitting problem")
  expect_true(is.na(out$rank_deficient))
  expect_true(is.na(out$estimable))
  expect_true(is.na(out$partial_R2))
  constant <- technical_partial_r2(cbind(a = seq_len(10)),
                                   data.frame(scanner = factor(rep("A", 10))),
                                   adjust_batch = FALSE)
  expect_match(constant$error, "contrasts")
  expect_equal(constant$status, "fit_error")
  missing <- technical_partial_r2(cbind(a = rep(NA_real_, 10)), technical,
                                  adjust_batch = FALSE)
  expect_equal(missing$n_obs, 0L)
  expect_match(missing$error, "No complete observations")
  dropped <- technical_partial_r2(cbind(a = seq_len(10)), technical,
                                  adjust_batch = FALSE, subset = 1:8)
  expect_match(dropped$error, "common analysis sample")
})

test_that("names are preserved and ambiguous metadata are rejected", {
  set.seed(56)
  technical <- data.frame("scanner model" = rnorm(20), check.names = FALSE)
  covar <- data.frame("age years" = rnorm(20), check.names = FALSE)
  out <- technical_partial_r2(cbind(a = rnorm(20)), technical, covar,
                              adjust_batch = FALSE, formula = y ~ `age years`)
  expect_equal(out$Technical, "scanner model")
  expect_true(is.finite(out$partial_R2))
  dotted <- technical_partial_r2(cbind(a = rnorm(20)), technical, covar,
                                 adjust_batch = FALSE, formula = y ~ .)
  expect_true(is.finite(dotted$partial_R2))
  expect_error(technical_partial_r2(cbind(a = 1:20), technical, technical,
                                  adjust_batch = FALSE), "unique across")
  for (name in c("y", ".primary_batch")) {
    expect_error(technical_partial_r2(cbind(a = 1:20), setNames(technical, name),
                                    adjust_batch = FALSE), "reserved")
  }
  expect_error(technical_partial_r2(cbind(a = 1:20), technical, covar,
                                  adjust_batch = FALSE, formula = y ~ missing), "columns of `covar`")
  expect_error(technical_partial_r2(cbind(a = 1:20), technical,
                                  adjust_batch = FALSE, formula = log(y) ~ 1), "response")
})

test_that("a technical effect orthogonal to the outcome has zero partial R2", {
  x <- rep(c(-1, 1), 10)
  y <- rep(c(-1, -1, 1, 1), 5)
  out <- technical_partial_r2(cbind(y), data.frame(scanner = x), adjust_batch = FALSE)
  expect_equal(out$partial_R2, 0, tolerance = 1e-12)
  zero <- technical_partial_r2(cbind(y = rep(0, 20)), data.frame(scanner = x),
                               adjust_batch = FALSE)
  expect_equal(zero$status, "zero_rss")
  expect_true(is.na(zero$partial_R2))
})

test_that("GAM RSS comparison matches independently fitted smooth models", {
  skip_if_not_installed("mgcv")
  set.seed(57)
  n <- 80L
  covar <- data.frame(age = runif(n))
  technical <- data.frame(scanner = rnorm(n))
  y <- sin(6 * covar$age) + technical$scanner + rnorm(n, sd = 0.2)
  out <- technical_partial_r2(cbind(y), technical, covar, adjust_batch = FALSE,
                              model = mgcv::gam, formula = y ~ s(age, k = 6), method = "REML")
  frame <- data.frame(y, covar, technical)
  full <- mgcv::gam(y ~ s(age, k = 6) + scanner, data = frame, method = "REML")
  reduced <- mgcv::gam(y ~ s(age, k = 6), data = frame, method = "REML")
  expect_true(is.finite(out$partial_R2))
  expect_equal(out$partial_R2, 1 - sum(residuals(full)^2) / sum(residuals(reduced)^2))
})

test_that("negative RSS comparisons and unavailable plot groups remain visible", {
  set.seed(58)
  # A deterministic fitter whose full-model residuals are deliberately worse
  # exercises the non-monotone RSS path without relying on a particular GAM version.
  nonmonotone <- function(formula, data, ...) {
    fit <- lm(formula, data = data, ...)
    fit$residuals <- rep(if ("scanner" %in% all.vars(formula)) 2 else 1, nrow(data))
    fit
  }
  out <- technical_partial_r2(cbind(a = rnorm(20)), data.frame(scanner = rnorm(20)),
                              adjust_batch = FALSE, model = nonmonotone)
  expect_equal(out$partial_R2, -3)
  expect_true(out$negative)
  expect_equal(out$status, "negative")
  result <- data.frame(Feature = c("a", "b", "a", "b"),
                       Technical = c("scanner", "scanner", "field", "field"),
                       partial_R2 = c(-0.2, NA, NA, NA))
  heatmap <- plot_technical_r2(result, result)
  built_heatmap <- ggplot2::ggplot_build(heatmap)
  expect_false(built_heatmap$data[[1]]$fill[1] == "grey85")
  expect_equal(heatmap$scales$scales[[1]]$limits, c(-1, 1))
  expect_equal(heatmap$data$median_delta_r2, c(0, NA))
  expect_error(plot_technical_r2(result, result, labels = c("Same", "Same")), "distinct")
})

test_that("summarize_r2 describes finite feature estimates without hiding failures", {
  x <- data.frame(
    Technical = factor(c(rep("scanner model", 7), "field", "field", "single")),
    partial_R2 = c(-0.2, 0, 0.2, 0.4, NA, Inf, NaN, NA, -Inf, 0.3)
  )
  out <- summarize_r2(x, method = "CALM")
  expect_named(out, c("Technical", "median_R2", "q25_R2", "q75_R2", "n_valid", "n_total", "Method"))
  expect_equal(out$Technical, c("scanner model", "field", "single"))
  expect_equal(out$Method, rep("CALM", 3))
  expect_equal(out$n_total, c(7L, 2L, 1L))
  expect_equal(out$n_valid, c(4L, 0L, 1L))
  expect_equal(out$median_R2, c(0.1, NA, 0.3))
  expect_equal(out$q25_R2, c(-0.05, NA, 0.3))
  expect_equal(out$q75_R2, c(0.25, NA, 0.3))
  expect_false("Method" %in% names(summarize_r2(x)))
})

test_that("summarize_r2 handles empty input and validates labels and columns", {
  empty <- data.frame(Technical = character(), partial_R2 = numeric())
  out <- summarize_r2(empty, "Raw")
  expect_equal(nrow(out), 0L)
  expect_type(out$Method, "character")
  expect_type(out$n_valid, "integer")
  expect_type(out$median_R2, "double")
  expect_error(summarize_r2(data.frame(Technical = "scanner")), "numeric `partial_R2`")
  expect_error(summarize_r2(data.frame(Technical = "scanner", partial_R2 = "0.1")), "numeric `partial_R2`")
  expect_error(summarize_r2(data.frame(Technical = NA_character_, partial_R2 = 0)), "nonmissing")
  expect_error(summarize_r2(empty, c("Raw", "CALM")), "single nonempty")
  expect_error(summarize_r2(empty, NA_character_), "single nonempty")
})

test_that("relative heatmaps align features and summarize paired percentage-point changes", {
  raw <- data.frame(Feature = c("a", "b", "c"), Technical = "scanner",
                    partial_R2 = c(0.1, 0.8, 0.9))
  after <- raw
  after$partial_R2 <- c(0.2, 0.2, 0.95)
  after <- after[c(3, 1, 2), ]
  plot <- plot_technical_r2(raw, after, labels = c("Raw", "CALM"))
  expect_equal(plot$data$median_delta_r2, 5)
  expect_equal(plot$data$change_label, "+5.00")
  expect_equal(plot$data$n_valid, 3L)
  expect_equal(as.character(plot$data$Method), "CALM")
  expect_equal(plot$scales$scales[[1]]$limits, c(-5, 5))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$data[[1]]), 1L)
  expect_equal(built$data[[2]]$label, "+5.00")
  expect_equal(as.character(plot$data$Technical), "scanner")
})

test_that("aggregate heatmaps support multiple methods as columns", {
  raw <- data.frame(Feature = c("a", "b", "a", "b"),
                    Technical = rep(c("scanner", "field"), each = 2),
                    partial_R2 = c(0.2, 0.4, 0.6, 0.8))
  calm <- raw
  calm$partial_R2 <- raw$partial_R2 - 0.1
  combat <- raw
  combat$partial_R2 <- raw$partial_R2 + 0.05
  plot <- plot_technical_r2(raw, list(CALM = calm, ComBat = combat))
  expect_equal(levels(plot$data$Method), c("CALM", "ComBat"))
  expect_equal(plot$data$median_delta_r2[plot$data$Method == "CALM"], c(-10, -10))
  expect_equal(plot$data$median_delta_r2[plot$data$Method == "ComBat"], c(5, 5))
  expect_equal(plot$scales$scales[[1]]$limits, c(-10, 10))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(length(unique(built$data[[1]]$x)), 2L)
  expect_equal(length(unique(built$data[[1]]$y)), 2L)
  expect_equal(nrow(built$layout$layout), 1L)
})

test_that("relative heatmaps retain unmatched and nonfinite pairs in counts", {
  raw <- data.frame(Feature = c("a", "b", "c", "a"),
                    Technical = c("scanner", "scanner", "scanner", "field"),
                    partial_R2 = c(0.4, NA, 0.3, 0.1))
  after <- data.frame(Feature = c("a", "b", "d", "a"),
                      Technical = c("scanner", "scanner", "scanner", "field"),
                      partial_R2 = c(0.2, 0.1, 0.3, Inf))
  plot <- plot_technical_r2(raw, after, limit = 5)
  scanner <- plot$data[plot$data$Technical == "scanner", ]
  field <- plot$data[plot$data$Technical == "field", ]
  expect_equal(scanner$median_delta_r2, -20)
  expect_equal(scanner$n_total, 4L)
  expect_equal(scanner$n_valid, 1L)
  expect_equal(scanner$n_failed, 3L)
  expect_true(is.na(field$median_delta_r2))
  expect_equal(field$change_label, "NA")
  expect_equal(plot$scales$scales[[1]]$limits, c(-5, 5))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[2]]$label[1], "-20.00")
  expect_equal(built$data[[1]]$fill[2], "grey85")
  expect_false(built$data[[1]]$fill[1] == "grey85")
  absent <- after
  absent$partial_R2[] <- NA_real_
  unavailable <- plot_technical_r2(raw, absent)
  expect_equal(unavailable$scales$scales[[1]]$limits, c(-1, 1))
  expect_true(all(unavailable$data$n_valid == 0L))
})

test_that("relative heatmaps reject ambiguous pairs and invalid parameters", {
  raw <- data.frame(Feature = "a", Technical = "scanner", partial_R2 = 0.2)
  expect_error(plot_technical_r2(rbind(raw, raw), raw), "unique row")
  expect_error(plot_technical_r2(raw, list(raw, raw)), "distinct method names")
  expect_error(plot_technical_r2(raw, raw, type = "summary"), "arg")
  expect_error(plot_technical_r2(raw, raw, confounding = "scenario"), "unused argument")
  expect_error(plot_technical_r2(raw, raw, limit = 0), "positive finite")
  missing_key <- raw
  missing_key$Feature <- NA_character_
  expect_error(plot_technical_r2(raw, missing_key), "nonmissing")
  empty <- raw[FALSE, ]
  plot <- plot_technical_r2(empty, empty)
  expect_equal(nrow(plot$data), 0L)
  expect_equal(plot$scales$scales[[1]]$limits, c(-1, 1))
})

test_that("feature boxplots include raw and absolute partial R2 distributions", {
  raw <- data.frame(Feature = c("b", "a", "c", "b", "a", "c"),
                    Technical = rep(c("scanner", "field"), each = 3),
                    partial_R2 = c(0.2, 0.5, 0.1, 0.3, 0.4, 0.2))
  calm <- raw
  calm$partial_R2 <- raw$partial_R2 + c(-0.3, 0.2, NA, 0, -0.3, 0.05)
  calm <- calm[c(4, 2, 6, 1, 5, 3), ]
  combat <- raw
  combat$partial_R2 <- raw$partial_R2 - 0.05
  plot <- plot_technical_r2(raw, list(CALM = calm, ComBat = combat), type = "feature")
  expect_equal(levels(plot$data$Technical), c("scanner", "field"))
  expect_equal(levels(plot$data$Method), c("Before", "CALM", "ComBat"))
  expect_equal(nrow(plot$data), 18L)
  expect_equal(plot$data$partial_R2[plot$data$Method == "Before"], raw$partial_R2)
  expect_s3_class(plot$layers[[1]]$geom, "GeomBoxplot")
  expect_equal(length(plot$layers), 1L)
  expect_false(any(c("delta_r2", "raw_R2", "harmonized_R2") %in% names(plot$data)))
  row <- with(plot$data, Feature == "b" & Technical == "scanner" & Method == "CALM")
  expect_equal(plot$data$partial_R2[row], -0.1)
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$data[[1]]), 6L)
  expect_equal(nrow(built$layout$layout), 2L)
  expect_equal(built$data[[1]]$middle, c(0.2, 0.3, 0.15, 0.3, 0.25, 0.25))
  expect_equal(plot$labels$y, "Partial R2")
})

test_that("feature boxplots use finite absolute estimates without requiring pairs", {
  raw <- data.frame(Feature = c("a", "b", "c"), Technical = "scanner",
                    partial_R2 = c(NA, 0.4, Inf))
  after <- data.frame(Feature = c("a", "d", "e"), Technical = "scanner",
                      partial_R2 = c(0.1, -0.2, NaN))
  plot <- plot_technical_r2(raw, after, type = "feature", labels = c("Raw", "CALM"), limit = 0.01)
  expect_equal(levels(plot$data$Method), c("Raw", "CALM"))
  expect_equal(nrow(plot$data), 6L)
  expect_true(any(is.infinite(plot$data$partial_R2)))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1]]$middle, c(0.4, -0.05))
  expect_true(built$layout$panel_params[[1]]$y.range[1] < -0.2)
  after$partial_R2[] <- NA_real_
  unavailable <- plot_technical_r2(raw, after, type = "feature")
  unavailable_built <- ggplot2::ggplot_build(unavailable)
  expect_equal(nrow(unavailable_built$data[[1]]), 1L)
  expect_equal(unavailable_built$layout$panel_params[[1]]$x$get_labels(), c("Before", "After"))
})

test_that("joint mode tests each factor against biology and batch alone", {
  set.seed(81)
  n <- 60L
  covar <- data.frame(age = rnorm(n), diagnosis = factor(rep(c("A", "B"), 30)))
  technical <- data.frame(scanner = rnorm(n), field = rnorm(n))
  bat <- rep(1:3, each = 20)
  y <- covar$age + technical$scanner + technical$field + rnorm(n)
  calls <- 0L
  counted_lm <- function(...) {
    calls <<- calls + 1L
    stats::lm(...)
  }
  out <- technical_partial_r2(cbind(a = y, b = y + 1), technical, covar, bat,
                              formula = y ~ age * diagnosis, model = counted_lm, test = "joint")
  frame <- data.frame(y, covar, technical, batch = factor(bat))
  reduced <- lm(y ~ age * diagnosis + batch, frame)
  scanner_full <- lm(y ~ age * diagnosis + batch + scanner, frame)
  field_full <- lm(y ~ age * diagnosis + batch + field, frame)
  expected <- c(1 - sum(residuals(scanner_full)^2) / sum(residuals(reduced)^2),
                1 - sum(residuals(field_full)^2) / sum(residuals(reduced)^2))
  expect_equal(nrow(out), 4L)
  expect_equal(out$Technical, rep(names(technical), 2))
  expect_equal(out$partial_R2, rep(expected, 2))
  expect_equal(out$rank_added, rep(1L, 4))
  expect_equal(calls, 8L)
  expect_equal(attr(out, "test"), "joint")
  expect_equal(attr(out, "technical_factors"), names(technical))
  formulas <- attr(out, "formula")
  expect_equal(names(formulas), names(technical))
  expect_false("field" %in% all.vars(formulas$scanner))
  expect_false("scanner" %in% all.vars(formulas$field))
  default <- technical_partial_r2(cbind(y), technical, covar, bat, formula = y ~ age)
  separate <- technical_partial_r2(cbind(y), technical, covar, bat, formula = y ~ age, test = "separate")
  expect_equal(default, separate)
  expect_error(technical_partial_r2(cbind(y), technical, bat = bat, test = "marginal"), "arg")
})

test_that("joint mode includes associations shared with other technical factors", {
  set.seed(82)
  scanner <- rnorm(40)
  technical <- data.frame(scanner, scanner_duplicate = scanner)
  data <- cbind(a = 2 * scanner + rnorm(40))
  separate <- technical_partial_r2(data, technical, adjust_batch = FALSE)
  joint <- technical_partial_r2(data, technical, adjust_batch = FALSE, test = "joint")
  expect_true(all(separate$status == "non_estimable"))
  expect_false(any(joint$rank_deficient_full))
  expect_true(all(joint$rank_added == 1L))
  expect_true(all(joint$estimable))
  expect_true(all(joint$partial_R2 > 0.5))
  expect_equal(joint$partial_R2[1], joint$partial_R2[2])
  bat <- factor(rep(1:2, each = 20))
  confounded <- technical_partial_r2(data, data.frame(scanner = bat), bat = bat, test = "joint")
  expect_equal(confounded$status, "non_estimable")
  expect_equal(confounded$rank_added, 0L)
  expect_true(is.na(confounded$partial_R2))
})

test_that("joint models select complete cases per tested factor without batch adjustment", {
  set.seed(83)
  n <- 40L
  covar <- data.frame(age = rnorm(n), unused = NA_real_)
  technical <- data.frame(scanner = rnorm(n), field = rnorm(n))
  technical$field[1] <- NA_real_
  y <- rnorm(n)
  y[2] <- NA_real_
  out <- technical_partial_r2(cbind(y), technical, covar, bat = rep(NA, n),
                              formula = y ~ age, adjust_batch = FALSE, test = "joint")
  for (name in names(technical)) {
    frame <- data.frame(y, age = covar$age, tested = technical[[name]])
    frame <- frame[complete.cases(frame), ]
    full <- lm(y ~ age + tested, frame)
    reduced <- lm(y ~ age, frame)
    expect_equal(out$partial_R2[out$Technical == name],
                 1 - sum(residuals(full)^2) / sum(residuals(reduced)^2))
  }
  expect_equal(out$n_obs, c(39L, 38L))
  one <- technical[, "scanner", drop = FALSE]
  joint_one <- technical_partial_r2(cbind(y), one, adjust_batch = FALSE, test = "joint")
  separate_one <- technical_partial_r2(cbind(y), one, adjust_batch = FALSE)
  expect_equal(joint_one$partial_R2, separate_one$partial_R2)
  expect_equal(joint_one$rank_added, separate_one$rank_added)
  all_missing <- technical
  all_missing$field[] <- NA_real_
  mixed <- technical_partial_r2(cbind(y), all_missing, adjust_batch = FALSE, test = "joint")
  expect_true(is.finite(mixed$partial_R2[1]))
  expect_equal(mixed$status[2], "fit_error")
  expect_equal(mixed$n_obs[2], 0L)
  broken <- technical_partial_r2(cbind(y), technical, adjust_batch = FALSE, test = "joint",
                                 model = function(...) stop("Joint fitting error"))
  expect_true(all(broken$status == "fit_error"))
  expect_match(broken$error[1], "Joint fitting error")
})

test_that("joint GAM comparisons and visualizations retain individual technical names", {
  skip_if_not_installed("mgcv")
  set.seed(84)
  n <- 60L
  covar <- data.frame(age = runif(n))
  technical <- data.frame(scanner = rnorm(n), field = rnorm(n))
  y <- sin(6 * covar$age) + technical$scanner + technical$field + rnorm(n, sd = 0.2)
  out <- technical_partial_r2(cbind(a = y, b = y + rnorm(n, sd = 0.1)), technical, covar,
                              adjust_batch = FALSE, test = "joint", model = mgcv::gam,
                              formula = y ~ s(age, k = 5), method = "REML")
  frame <- data.frame(y, covar, technical)
  full <- mgcv::gam(y ~ s(age, k = 5) + scanner, data = frame, method = "REML")
  reduced <- mgcv::gam(y ~ s(age, k = 5), data = frame, method = "REML")
  expect_equal(out$partial_R2[1], 1 - sum(residuals(full)^2) / sum(residuals(reduced)^2))
  summary <- summarize_r2(out, "Raw")
  expect_equal(summary$Technical, names(technical))
  expect_equal(summary$n_valid, c(2L, 2L))
  aggregate <- plot_technical_r2(out, list(CALM = out))
  expect_equal(as.character(aggregate$data$Technical), names(technical))
  expect_equal(aggregate$data$median_delta_r2, c(0, 0))
  expect_equal(nrow(ggplot2::ggplot_build(aggregate)$data[[1]]), 2L)
  boxplot <- plot_technical_r2(out, list(CALM = out), type = "feature")
  expect_equal(levels(boxplot$data$Method), c("Before", "CALM"))
  expect_equal(nrow(ggplot2::ggplot_build(boxplot)$data[[1]]), 4L)
})
