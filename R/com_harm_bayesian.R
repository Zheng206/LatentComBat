#' Joint Bayesian Latent ComBat
#'
#' Fits batch effects and latent structure jointly using a Bayesian latent
#' ComBat model. Biological covariates and optional additional protected
#' covariates are included in a protected projection space, while the primary
#' batch variable is excluded from that space.
#'
#' @inheritParams com_harm
#'
#' @param posterior_draws Number of posterior draws of the harmonized data to
#'   retain for uncertainty-aware downstream inference. Set to `0` to retain
#'   only posterior means. Defaults to `0`.
#' @param posterior_seed Integer random seed used when subsampling retained
#'   posterior draws. Defaults to `1`.
#' @param adapt_delta Target acceptance probability for the NUTS sampler.
#'   Larger values use smaller leapfrog step sizes and can reduce divergent
#'   transitions at the cost of additional computation. Defaults to `0.98`.
#' @param iter_sampling Number of post-warmup MCMC iterations retained per
#'   chain. Defaults to `1000`.
#' @param iter_warmup Number of warmup iterations per chain used for NUTS
#'   adaptation. Warmup draws are not retained for posterior inference.
#'   Defaults to `1000`.
#' @param max_treedepth Maximum tree depth allowed for each NUTS trajectory.
#'   Larger values permit longer trajectories when the posterior geometry is
#'   difficult to explore, at the cost of additional computation. Defaults to
#'   `12`.
#' @param keep_inference_data Logical indicating whether the original inputs
#'   and model settings should be retained for optional downstream inference.
#'   Defaults to `TRUE`.
#'
#' @return An object of class `"latentcombat_bayesian"` containing:
#' \describe{
#'   \item{harm_data}{
#'     Harmonized data matrix.
#'   }
#'   \item{stan_result}{
#'     Results returned by the joint Bayesian Stan model.
#'   }
#'   \item{latent}{
#'     Information on latent-dimension selection, including `K_target`,
#'     `K_bayes`, and the latent-factor detection results.
#'   }
#'   \item{protection}{
#'     Information describing the protected covariate subspace, including the
#'     raw and reduced design matrices, orthonormal basis, rank, and whether
#'     batch was included in the protected space.
#'   }
#'   \item{data_stand_result}{
#'     Standardization results used to construct and rescale the Bayesian
#'     harmonization model.
#'   }
#'   \item{posterior_harmonization}{
#'     Optional original-scale posterior harmonized-data draws retained when
#'     `posterior_draws > 0`.
#'   }
#'   \item{inference_data}{
#'     Original inputs and settings retained for optional downstream inference
#'     when `keep_inference_data = TRUE`.
#'   }
#'   \item{method}{
#'     Character string `"bayesian"`.
#'   }
#' }
#'
#' @details
#' The data are first standardized after fitting the observed mean structure.
#' A protected covariate space is then constructed from the modeled biological
#' terms and any variables supplied through `protected_covar`. The primary
#' batch variable is intentionally excluded from this protected space.
#'
#' The target latent dimension is estimated from the standardized residual
#' matrix, and the Bayesian latent dimension is expanded by `K_buffer` up to
#' `K_max`. Batch effects and latent structure are then estimated jointly in
#' Stan, subject to the protected projection constraint.
#'
#' For longitudinal or repeated-measures data, `subj` is passed to the
#' corresponding longitudinal Stan model.
#'
#' @export
com_harm.bayesian <- function(
    bat,
    data,
    covar,
    model = lm,
    formula = NULL,
    ref.batch = NULL,
    robust.LS = FALSE,
    K_max = 10,
    K_buffer = 2,
    subj = NULL,
    protected_covar = NULL,
    posterior_draws = 0L,
    posterior_seed = 1L,
    adapt_delta = 0.98,
    iter_sampling = 1000,
    iter_warmup = 1000,
    max_treedepth = 12,
    keep_inference_data = TRUE,
    ...
) {
  .check_latentcombat_stan()
  batch_result <- batch_matrix(bat, ref.batch = ref.batch)

  fitted_model <- model_fitting(
    data = data,
    batch = batch_result$batch_matrix,
    covar = covar,
    model = model,
    formula = formula,
    batch_include = TRUE,
    ...
  )

  data_stand_result <- standardize_data(
    fitted_model,
    batch_result,
    robust.LS = robust.LS,
    bayes_sva = TRUE
  )

  type <- combat_type(type = "sva_combat")
  Y <- as.matrix(data)
  n <- nrow(Y)
  p <- ncol(Y)
  R <- as.matrix(data_stand_result$data_stand)
  sig_hat <- rep(1, p)

  protect <- build_Z_preserve(
    model = model,
    formula = formula,
    covar = covar,
    bat = bat,
    preserve_batch = FALSE,
    protected_covar = protected_covar,
    qr_tol = 1e-8,
    return_Q = TRUE
  )

  Z_bio_raw <- protect$Z_raw
  Z_bio <- protect$Z
  Q_bio <- protect$Q

  Z_K <- matrix(1, nrow = n, ncol = 1L)
  det <- detect_K_one(R, Z = Z_K, B = 50, alpha = 0.05)
  K_target <- det$K
  K_bayes <- min(K_max, max(1L, K_target + as.integer(K_buffer)))

  stan_data <- stan_data_prep(
    type = type,
    R = R,
    bat = bat,
    K_bayes = K_bayes,
    K_target = K_target,
    subj = subj,
    sigma_hat = sig_hat
  )
  stan_data$Qz <- as.integer(ncol(Q_bio))
  stan_data$Q_protect <- Q_bio

  stan_result <- stan_algorithm(
    type,
    stan_data,
    posterior_draws = posterior_draws,
    posterior_seed = posterior_seed,
    adapt_delta = adapt_delta,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    max_treedepth = max_treedepth
  )
  R_clean_stan <- stan_result$R_adj
  data_combat <- R_clean_stan * data_stand_result$sd_mat + data_stand_result$stand_mean

  posterior_harmonization <- NULL
  if (!is.null(stan_result$R_adj_draws)) {
    S <- dim(stan_result$R_adj_draws)[1L]
    harm_draws <- array(
      NA_real_, dim = c(S, n, p),
      dimnames = list(
        draw = seq_len(S),
        observation = rownames(Y),
        feature = colnames(Y)
      )
    )

    for (s in seq_len(S)) {
      harm_draws[s, , ] <- stan_result$R_adj_draws[s, , ] * data_stand_result$sd_mat + data_stand_result$stand_mean
    }

    posterior_harmonization <- list(
      harm_data_draws = harm_draws,
      n_draws = S,
      draw_ids = stan_result$draw_ids
    )
  }

  inference_data <- NULL
  if (isTRUE(keep_inference_data)) {
    inference_data <- list(
      data = as.matrix(data),
      bat = bat,
      covar = covar,
      model = model,
      formula = formula,
      ref.batch = ref.batch,
      robust.LS = robust.LS,
      K_max = K_max,
      K_buffer = K_buffer,
      subj = subj,
      protected_covar = protected_covar,
      adapt_delta = adapt_delta,
      iter_sampling = iter_sampling,
      iter_warmup = iter_warmup,
      max_treedepth = max_treedepth,
      extra_args = list(...)
    )
  }

  out <- list(
    harm_data = data_combat,
    stan_result = stan_result,
    latent = list(
      K_target = K_target,
      K_bayes = K_bayes,
      K_detection = det
    ),
    protection = list(
      Z_raw = Z_bio_raw,
      Z = Z_bio,
      Q = Q_bio,
      rank = protect$rank,
      preserve_batch = FALSE
    ),
    data_stand_result = data_stand_result,
    posterior_harmonization = posterior_harmonization,
    inference_data = inference_data,
    method = "bayesian"
  )
  class(out) <- c("latentcombat_bayesian", "latentcombat_fit", "list")
  message("Harmonization Complete.")
  out
}
