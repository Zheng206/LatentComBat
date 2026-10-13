test_that("fixed effects match direct OLS coefficients and term-wise BH", {
  set.seed(71)
  n <- 60L
  covar <- data.frame(age = rnorm(n), sex = factor(rep(c("F", "M"), n / 2)))
  data <- cbind(a = 2 * covar$age + rnorm(n), b = rnorm(n), c = rnorm(n))
  out <- extract_fixed_effects(data, covar, formula = y ~ age + sex, method = "CALM")
  expect_equal(nrow(out), 6L)
  expect_setequal(out$Term, c("age", "sexM"))
  expect_false("(Intercept)" %in% out$Term)
  expect_true(all(out$Method == "CALM"))
  expect_true(all(out$N_obs == n))
  for (feature in colnames(data)) {
    tab <- summary(lm(y ~ age + sex, data.frame(y = data[, feature], covar)))$coefficients
    for (term in c("age", "sexM")) {
      row <- out[out$Feature == feature & out$Term == term, ]
      expect_equal(row$Estimate, unname(tab[term, 1]))
      expect_equal(row$SE, unname(tab[term, 2]))
      expect_equal(row$Statistic, unname(tab[term, 3]))
      expect_equal(row$P_value, unname(tab[term, 4]))
    }
  }
  for (term in unique(out$Term)) {
    rows <- out$Term == term
    expect_equal(out$Q_value[rows], p.adjust(out$P_value[rows], "BH"))
  }
  filtered <- extract_fixed_effects(data, covar, terms = "sexM")
  expect_equal(filtered$Q_value, out$Q_value[out$Term == "sexM"])
  expect_equal(unique(filtered$Term), "sexM")
})

test_that("missing data, failed fits, and aliases remain visible", {
  set.seed(72)
  n <- 40L
  covar <- data.frame(age = rnorm(n), unused = NA_real_)
  data <- cbind(a = rnorm(n), b = rnorm(n), failed = NA_real_)
  data[1:2, "b"] <- NA_real_
  out <- extract_fixed_effects(data, covar, formula = y ~ age)
  expect_equal(out$N_obs, c(40L, 38L, 0L))
  expect_equal(out$Status, c("ok", "ok", "fit_error"))
  expect_match(out$Error[3], "No complete observations")
  expect_true(is.na(out$Estimate[3]))
  expect_equal(out$Q_value[1:2], p.adjust(out$P_value[1:2], "BH", n = 3))
  expect_equal(nrow(attr(out, "feature_diagnostics")), 3L)
  covar$duplicate <- covar$age
  aliased <- extract_fixed_effects(data[, 1:2], covar, formula = y ~ age + duplicate)
  expect_true(all(aliased$Status[aliased$Term == "duplicate"] == "non_estimable"))
  expect_true(all(is.na(aliased$Estimate[aliased$Term == "duplicate"])))
  expect_true(all(is.finite(aliased$Estimate[aliased$Term == "age"])))
  broken <- function(...) stop("Specific fitting failure")
  failed <- extract_fixed_effects(data, covar, formula = y ~ age, terms = "age", model = broken)
  expect_true(all(failed$Status == "fit_error"))
  expect_match(failed$Error[1], "Specific fitting failure")
  expect_error(extract_fixed_effects(data[, 1:2], covar, formula = y ~ age, model = broken), "Specific fitting failure")
  subset <- extract_fixed_effects(data[, 1:2], covar, formula = y ~ age, terms = "age", subset = 1:10)
  expect_true(all(subset$Status == "fit_error"))
  expect_match(subset$Error[1], "analysis sample")
})

test_that("factor contrasts stay consistent across feature-specific samples", {
  set.seed(73)
  covar <- data.frame(sex = factor(rep(c("F", "M"), each = 15)))
  data <- cbind(a = rnorm(30), b = rnorm(30))
  data[1:15, "b"] <- NA_real_
  out <- extract_fixed_effects(data, covar)
  expect_equal(out$Term, c("sexM", "sexM"))
  expect_equal(out$Status, c("ok", "fit_error"))
  expect_match(out$Error[2], "contrasts would change")
})

test_that("GAM extraction returns parametric terms and excludes smooth coefficients", {
  skip_if_not_installed("mgcv")
  set.seed(74)
  covar <- data.frame(age = runif(60), sex = factor(rep(c("F", "M"), 30)))
  data <- cbind(a = sin(6 * covar$age) + (covar$sex == "M") + rnorm(60, sd = 0.2))
  out <- extract_fixed_effects(data, covar, model = mgcv::gam,
                              formula = y ~ s(age, k = 5) + sex, method = "CALM")
  fit <- mgcv::gam(y ~ s(age, k = 5) + sex, data = data.frame(y = data[, 1], covar))
  expect_equal(out$Term, "sexM")
  expect_equal(out$Estimate, unname(summary(fit)$p.table["sexM", 1]))
  expect_equal(out$P_value, unname(summary(fit)$p.table["sexM", 4]))
  smooth_only <- extract_fixed_effects(data, covar, model = mgcv::gam, formula = y ~ s(age, k = 5))
  expect_equal(nrow(smooth_only), 0L)
  expect_equal(attr(smooth_only, "feature_diagnostics")$Status, "ok")
})

test_that("fixed-effect extraction validates formulas, names, and requested terms", {
  data <- cbind(a = 1:20)
  covar <- data.frame(age = rnorm(20))
  expect_error(extract_fixed_effects(data, covar[1:10, , drop = FALSE]), "same number")
  expect_error(extract_fixed_effects(data, setNames(covar, "y")), "reserved")
  expect_error(extract_fixed_effects(data, covar, formula = log(y) ~ age), "response")
  expect_error(extract_fixed_effects(data, covar, formula = y ~ missing), "columns of `covar`")
  expect_error(extract_fixed_effects(data, covar, terms = "missing"), "Unknown coefficient")
  expect_error(extract_fixed_effects(data, covar, terms = "(Intercept)"), "excluding the intercept")
  expect_error(extract_fixed_effects(data, covar, method = NA_character_), "nonempty")
  expect_error(extract_fixed_effects(data, covar, model = function(...) structure(list(), class = "other")), "Only lm and gam")
  covar <- setNames(covar, "age years")
  out <- extract_fixed_effects(data, covar, formula = y ~ `age years`)
  expect_equal(out$Term, "`age years`")
})

fixed_plot_fixture <- function() {
  raw <- data.frame(Feature = c("a", "b", "c", "a", "b", "c"),
                    Term = rep(c("age", "sexM"), each = 3), Estimate = 1:6,
                    Q_value = c(0.01, 0.05, NA, NA, NA, NA), Method = "Raw")
  calm <- raw
  calm$Estimate <- raw$Estimate + c(0.2, -0.1, 0, -0.2, 0.1, 0)
  calm$Q_value <- c(0.01, 0.02, 0.2, NA, NA, NA)
  calm$Method <- "CALM"
  list(raw = raw, calm = calm)
}

test_that("fixed-effect estimate and change plots use paired keys", {
  fixture <- fixed_plot_fixture()
  fixture$calm <- fixture$calm[c(6, 4, 2, 1, 3, 5), ]
  estimate <- plot_fixed_effect_comparison(fixture$raw, fixture$calm)
  expect_equal(levels(estimate$data$Method), "CALM")
  expect_equal(nrow(estimate$data), 6L)
  row <- estimate$data$Feature == "a" & estimate$data$Term == "age"
  expect_equal(estimate$data$Raw_Estimate[row], 1)
  expect_equal(estimate$data$Estimate[row], 1.2)
  expect_equal(estimate$data$Delta[row], 0.2)
  built <- ggplot2::ggplot_build(estimate)
  expect_equal(built$data[[1]]$slope, c(1, 1))
  change <- plot_fixed_effect_comparison(fixture$raw, list(CALM = fixture$calm), "change")
  expect_s3_class(change$layers[[2]]$geom, "GeomBoxplot")
  expect_equal(nrow(ggplot2::ggplot_build(change)$data[[2]]), 2L)
  extra <- fixture$calm[1, ]
  extra$Feature <- "unmatched"
  fixture$calm <- rbind(fixture$calm, extra)
  paired <- plot_fixed_effect_comparison(fixture$raw, fixture$calm, "change")
  expect_equal(nrow(paired$data), 7L)
  expect_equal(nrow(ggplot2::ggplot_build(paired)$data[[2]]), 2L)
  expect_true(is.na(paired$data$Raw_Estimate[paired$data$Feature == "unmatched"]))
})

test_that("significance plots show valid denominators and unavailable groups", {
  fixture <- fixed_plot_fixture()
  plot <- plot_fixed_effect_comparison(fixture$raw, list(CALM = fixture$calm), "significance")
  age <- plot$data[plot$data$Term == "age", ]
  expect_equal(age$significant, c(0.5, 2 / 3))
  expect_equal(age$n_valid, c(2L, 3L))
  expect_equal(age$n_total, c(3L, 3L))
  expect_equal(age$n_significant, c(1L, 2L))
  expect_equal(age$count_label, c("1/2", "2/3"))
  sex <- plot$data[plot$data$Term == "sexM", ]
  expect_true(all(is.na(sex$significant)))
  expect_equal(sex$n_failed, c(3L, 3L))
  expect_equal(sex$count_label, c("0/0", "0/0"))
  expect_equal(levels(plot$data$Method), c("Raw", "CALM"))
  built <- ggplot2::ggplot_build(plot)
  expect_true(all(c("1/2", "2/3", "0/0") %in% built$data[[2]]$label))
  fixture$calm <- fixture$calm[-1, ]
  missing <- plot_fixed_effect_comparison(fixture$raw, fixture$calm, "significance")
  expect_equal(missing$data$n_total[missing$data$Term == "age"], c(3L, 3L))
  expect_equal(missing$data$n_valid[missing$data$Term == "age"], c(2L, 2L))
})

test_that("fixed-effect comparison validates names, keys, and significance thresholds", {
  fixture <- fixed_plot_fixture()
  expect_error(plot_fixed_effect_comparison(fixture$raw, list(fixture$calm)), "named list")
  expect_error(plot_fixed_effect_comparison(fixture$raw, list(Raw = fixture$calm)), "other than `Raw`")
  expect_error(plot_fixed_effect_comparison(rbind(fixture$raw, fixture$raw), fixture$calm), "Duplicate")
  expect_error(plot_fixed_effect_comparison(fixture$raw, fixture$calm, alpha = 1), "between zero and one")
  bad <- fixture$calm
  bad$Q_value[1] <- -1
  expect_error(plot_fixed_effect_comparison(fixture$raw, bad), "between zero and one")
})

test_that("significance annotations agree with zero, half, and all significant bars", {
  raw <- data.frame(Feature = rep(paste0("f", 1:60), 2),
                    Term = rep(c("age_z", "sex1"), each = 60),
                    Estimate = 1, Q_value = c(rep(0.01, 60), rep(0.5, 60)))
  calm <- raw
  calm$Q_value[61:90] <- 0.01
  plot <- plot_fixed_effect_comparison(raw, list(ComBat = raw, CALM = calm), "significance")
  age <- plot$data[plot$data$Term == "age_z", ]
  sex <- plot$data[plot$data$Term == "sex1", ]
  expect_equal(age$count_label, rep("60/60", 3))
  expect_equal(age$significant, rep(1, 3))
  expect_equal(sex$count_label, c("0/60", "0/60", "30/60"))
  expect_equal(sex$significant, c(0, 0, 0.5))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(built$data[[1]]$y, plot$data$significant[order(plot$data$Term, plot$data$Method)])
  expect_equal(built$data[[2]]$label, plot$data$count_label)
  expect_equal(built$layout$panel_params[[1]]$y.range, c(0, 1.08))
  expect_equal(built$layout$panel_params[[1]]$y$get_breaks(), seq(0, 1, 0.25))
})
