make_fake_cmdstan_fit <- function(model = c("latent", "univariate"), bad = FALSE) {
  model <- match.arg(model)

  vars <- if (model == "latent") {
    c("mu_g[1]", "mu_g[2]", "tau_b[1]", "sigma0[1]", "sigma_logdelta_b[1]")
  } else {
    c("gamma_i[1]", "sigma_psi[1]", "mu_ig[1]", "sigma_delta[1]")
  }

  summary_fun <- function(variables = NULL) {
    keep <- if (is.null(variables)) {
      rep(TRUE, length(vars))
    } else {
      vapply(vars, function(v) any(startsWith(v, paste0(variables, "[")) | v %in% variables), logical(1))
    }
    if (!any(keep)) return(data.frame())

    data.frame(
      variable = vars[keep],
      mean = 0,
      median = 0,
      sd = 1,
      mad = 1,
      q5 = -1,
      q95 = 1,
      rhat = if (bad) rep(1.08, sum(keep)) else rep(1.001, sum(keep)),
      ess_bulk = if (bad) rep(40, sum(keep)) else rep(500, sum(keep)),
      ess_tail = if (bad) rep(45, sum(keep)) else rep(450, sum(keep)),
      stringsAsFactors = FALSE
    )
  }

  sampler_fun <- function(format = "draws_array") {
    d <- array(
      0,
      dim = c(20, 2, 3),
      dimnames = list(NULL, NULL, c("divergent__", "treedepth__", "energy__"))
    )
    d[, , "energy__"] <- cbind(rep(c(0, 1), 10), rep(c(1, 0), 10))
    if (bad) {
      d[1:2, 1, "divergent__"] <- 1
      d[1:3, 1, "treedepth__"] <- 12
    }
    d
  }

  out <- list(summary = summary_fun, sampler_diagnostics = sampler_fun)
  class(out) <- c("CmdStanMCMC", "list")
  out
}

test_that("Bayesian diagnostic helpers classify numerical edge cases", {
  expect_true(is.na(LatentComBat:::.safe_max(c(NA, Inf))))
  expect_true(is.na(LatentComBat:::.safe_min(c(NA, Inf))))
  expect_true(is.na(LatentComBat:::.safe_prop(c(NA, NA))))
  expect_equal(LatentComBat:::.safe_prop(c(TRUE, FALSE, NA)), 0.5)

  expect_equal(LatentComBat:::.worst_status("PASS", "WARNING"), "WARNING")
  expect_equal(LatentComBat:::.worst_status("WARNING", "FAIL"), "FAIL")
  expect_equal(LatentComBat:::.parameter_status(1.001, 200, 200), "PASS")
  expect_equal(LatentComBat:::.parameter_status(1.02, 200, 200), "WARNING")
  expect_equal(LatentComBat:::.parameter_status(1.06, 200, 200), "FAIL")
})

test_that("Bayesian model type and core variables are detected", {
  latent <- make_fake_cmdstan_fit("latent")
  uni <- make_fake_cmdstan_fit("univariate")

  expect_equal(LatentComBat:::.detect_bayes_model(latent), "latent")
  expect_equal(LatentComBat:::.detect_bayes_model(uni), "univariate")
  expect_true("mu_g" %in% LatentComBat:::.default_core_variables(latent))
  expect_true("gamma_i" %in% LatentComBat:::.default_core_variables(uni))
})

test_that("Bayesian diagnostics summarize a healthy fake latent fit", {
  fit <- make_fake_cmdstan_fit("latent", bad = FALSE)
  out <- bayes_diagnostics(fit)

  expect_s3_class(out, "latentcombat_bayes_diagnostics")
  expect_equal(out$model_type, "latent")
  expect_equal(out$parameter_status, "PASS")
  expect_equal(out$sampling_status, "PASS")
  expect_equal(out$status, "PASS")
  expect_equal(out$sampler$n_chains, 2L)
  expect_equal(out$sampler$n_draws_total, 40L)

  sm <- summary(out)
  expect_s3_class(sm, "summary_latentcombat_bayes_diagnostics")
  expect_true(nrow(sm$group_summary) >= 1L)
  expect_output(print(out), "Overall status")
  expect_output(print(sm), "Diagnostic summary|diagnostic summary")
})

test_that("Bayesian diagnostics flag poor fake sampling", {
  fit <- make_fake_cmdstan_fit("univariate", bad = TRUE)
  out <- bayes_diagnostics(fit)

  expect_equal(out$model_type, "univariate")
  expect_equal(out$parameter_status, "FAIL")
  expect_true(out$sampling_status %in% c("WARNING", "FAIL"))
  expect_equal(out$status, "FAIL")
  expect_gt(out$convergence$n_rhat_bad, 0L)
})

test_that("Bayesian diagnostic fit extraction supports wrapper lists", {
  fit <- make_fake_cmdstan_fit("latent")
  expect_identical(LatentComBat:::.get_bayes_stan_fit(fit), fit)
  expect_identical(LatentComBat:::.get_bayes_stan_fit(list(fit = fit)), fit)
  expect_identical(LatentComBat:::.get_bayes_stan_fit(list(stan_result = list(fit = fit))), fit)
  expect_error(LatentComBat:::.get_bayes_stan_fit(list()), "Could not locate")
})
