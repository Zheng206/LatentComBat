# Internal primary-batch ComBat/CovBat stage.
.com_harm_combat <- function(
    data,
    bat,
    covar,
    batch_result,
    model,
    formula,
    eb = TRUE,
    stan = FALSE,
    robust.LS = FALSE,
    ref.batch = NULL,
    cov = FALSE,
    var_thresh = 0.95,
    min_rblock = 1,
    max_rblock = Inf,
    adapt_delta = 0.98,
    iter_sampling = 1000,
    iter_warmup = 1000,
    max_treedepth = 12,
    ...
) {
  feature <- colnames(data)

  fitted_model <- model_fitting(
    data = data,
    batch = batch_result$batch_matrix,
    covar = covar,
    model = model,
    formula = formula,
    batch_include = TRUE,
    ...
  )

  data_stand_result <- standardize_data(fitted_model, batch_result, robust.LS = robust.LS)
  type <- combat_type(type = "univariate")

  if (!stan) {
    fit <- eb_algorithm(type, data_stand_result, batch_result, eb = eb, robust.LS = robust.LS)
    class(fit) <- type
    gamma_star <- fit$gamma_star
    delta_star <- fit$delta_star
    engine <- "eb"
  } else {
    stan_data <- stan_data_prep(type, data_stand_result, bat)
    fit <- stan_algorithm(type, stan_data, batch_names = levels(bat), adapt_delta = adapt_delta, iter_sampling = iter_sampling,
                          iter_warmup = iter_warmup, max_treedepth = max_treedepth)
    gamma_star <- fit$gamma_star
    delta_star <- fit$delta_star
    engine <- "stan"
  }

  data_nb <- data_stand_result$data_stand
  batch_levels <- rownames(gamma_star)

  for (b in batch_levels) {
    idx <- batch_result$batch_index[[b]]
    data_nb[idx, ] <- sweep(data_nb[idx, , drop = FALSE], 2, gamma_star[b, ], "-")
    data_nb[idx, ] <- sweep(data_nb[idx, , drop = FALSE], 2, sqrt(delta_star[b, ]), "/")
  }

  data_combat <- data_nb * data_stand_result$sd_mat + data_stand_result$stand_mean

  if (!is.null(batch_result$ref)) {data_combat[batch_result$ref, ] <- data[batch_result$ref, ]}

  if (isTRUE(cov)) {
    if (!exists("covbat", mode = "function", inherits = TRUE)) {
      stop("`cov = TRUE` requires a `covbat()` implementation, but none is ", "available in the package namespace.", call. = FALSE)
    }
    R <- data_combat - data_stand_result$stand_mean
    R_star <- covbat(
      R,
      site = bat,
      center = TRUE,
      scale_scores = TRUE,
      var_thresh = var_thresh,
      min_rblock = min_rblock,
      max_rblock = max_rblock,
      ref.batch = ref.batch
    )
    data_combat <- R_star + data_stand_result$stand_mean
  }

  list(
    harm_data = data_combat,
    fit = fit,
    engine = engine,
    resid = data_nb,
    data_stand_result = data_stand_result
  )
}


.check_latentcombat_stan <- function() {

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      paste0(
        "The Bayesian method requires `cmdstanr`.\n",
        "Run `setup_latentcombat(install = TRUE)` to configure the Stan backend."
      ),
      call. = FALSE
    )
  }

  version <- tryCatch(
    cmdstanr::cmdstan_version(error_on_NA = FALSE),
    error = function(e) NULL
  )

  if (is.null(version) || all(is.na(version))) {
    stop(
      paste0(
        "CmdStan is not available.\n",
        "Run `setup_latentcombat(install = TRUE)` to configure the Stan backend."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}
