#' Load and Compile a Stan Model
#'
#' Locates a Stan model included with the package and compiles it using
#' [cmdstanr::cmdstan_model()]. Models are compiled when first requested.
#'
#' @param name Character string specifying the Stan model name, without the
#'   `.stan` extension.
#'
#' @return A compiled `CmdStanModel` object.
#'
#' @export
get_stan_model <- function(name = "my_model") {
  model_path <- system.file("stan", paste0(name, ".stan"), package = "LatentComBat")
  if (!file.exists(model_path)) {
    stop("Model not found: ", name)
  }
  cmdstanr::cmdstan_model(model_path)
}


#' Run a Stan Harmonization Algorithm
#'
#' S3 generic for fitting Stan-based harmonization models.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return The fitted Stan harmonization result.
#'
#' @export
stan_algorithm <- function(type, ...) {
  UseMethod("stan_algorithm")
}

#' Prepare Data for Stan Harmonization
#'
#' S3 generic for constructing data inputs required by Stan harmonization
#' models.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return A list containing data formatted for the selected Stan model.
#'
#' @export
stan_data_prep <- function(type, ...) {
  UseMethod("stan_data_prep")
}

#' @rdname stan_data_prep
#'
#' @param data_stand_result Standardized univariate data object.
#' @param bat Factor or vector specifying batch membership, with one
#'   value per observation.
#'
#' @return A list containing the standardized observations, feature and batch
#'   indices, and prior hyperparameters required by the univariate Stan model.
#'
#' @method stan_data_prep univariate
#' @export
stan_data_prep.univariate <- function(
    type,
    data_stand_result,
    bat,
    ...
) {

  data_mat <- as.matrix(data_stand_result$data_stand)

  n <- nrow(data_mat)
  p <- ncol(data_mat)

  if (length(bat) != n) {
    stop("`bat` must have one entry per observation.", call. = FALSE)
  }

  bat <- droplevels(as.factor(bat))

  observation <- rep(seq_len(n), times = p)
  feature <- rep(seq_len(p), each = n)

  batch_id <- as.integer(bat[observation])

  list(
    N = n * p,
    G = p,
    I = nlevels(bat),
    y = as.vector(data_mat),
    g = feature,
    i = batch_id,
    nu_psi = 3,
    nu_tau = 4,
    nu_delta = 3,
    A_psi = 2.5,
    A_tau = 2.5,
    A_delta = 2.5
  )
}

#' Reshape Univariate Batch-Effect Estimates
#'
#' Converts long-format posterior summaries of feature-wise batch effects into
#' a batch-by-feature matrix.
#'
#' @param gamma_star Data frame containing posterior batch-effect summaries.
#'
#' @return A numeric matrix with batches in rows and features in columns.
#'
#' @keywords internal
gamma_star_reshape <- function(gamma_star){
  result <- gamma_star %>% pivot_wider(names_from = .data[["feature"]], values_from = .data[["mean"]]) %>% dplyr::select(-1) %>% as.matrix()
  return(result)
}

#' @rdname stan_algorithm
#'
#' @param stan_data List of data and hyperparameters prepared for the
#'   univariate Stan model.
#' @param chains Number of MCMC chains. Defaults to `3`.
#' @param parallel_chains Number of chains run in parallel. Defaults to `3`.
#' @param batch_names Character vector containing the batch labels used to
#'   assign row names to the returned batch-effect matrices.
#' @param adapt_delta Target acceptance probability for the NUTS sampler.
#'   Larger values use smaller leapfrog step sizes and can reduce divergent
#'   transitions at the cost of additional computation. Defaults to `0.98`.
#' @param iter_sampling Number of post-warmup MCMC iterations retained per
#'   chain. Defaults to `1000`.
#' @param iter_warmup Number of warmup iterations per chain used for sampler
#'   adaptation. Warmup draws are not retained for posterior inference.
#'   Defaults to `1000`.
#' @param max_treedepth Maximum tree depth allowed for each NUTS trajectory.
#'   Larger values permit longer trajectories when the posterior geometry is
#'   difficult to explore, at the cost of additional computation. Defaults to
#'   `12`.
#' @param ... Additional arguments reserved for method compatibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma_star}{
#'     Matrix of posterior mean additive batch-effect estimates, with batches
#'     in rows and features in columns.
#'   }
#'   \item{delta_star}{
#'     Matrix of posterior mean multiplicative batch-effect variance estimates,
#'     with batches in rows and features in columns.
#'   }
#' }
#'
#'
#' @method stan_algorithm univariate
#' @export
stan_algorithm.univariate <- function(type, stan_data, chains = 3, parallel_chains = 3, batch_names, adapt_delta = 0.98, iter_sampling = 1000,
                                      iter_warmup = 1000, max_treedepth = 12, ...){
  message("Starting univariate Stan algorithm...")
  start_time <- Sys.time()
  message("[1/4] Loading Stan model...")
  stan_model <- get_stan_model("univariate_model")
  message("[2/4] Running MCMC sampling (this may take a while)...")
  fit <- stan_model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = parallel_chains,
    adapt_delta = adapt_delta,
    iter_sampling = iter_sampling,
    iter_warmup = iter_warmup,
    show_messages = TRUE
  )
  message("[3/4] Processing gamma_star results...")
  gamma_star <- fit$summary(variable = "mu_ig", "mean") %>%
    tidyr::extract(
      .data[["variable"]],
      into = c("feature", "batch"),
      regex = "mu_ig\\[([0-9]+),([0-9]+)\\]",
      convert = TRUE  # Converts to integers
    )
  message("[4/4] Processing delta_star results...")
  delta_star <- fit$summary(variable = "sigma_delta", "mean") %>%
    dplyr::mutate(mean = .data[["mean"]]^2) %>%
    tidyr::extract(
      .data[["variable"]],
      into = c("feature", "batch"),
      regex = "sigma_delta\\[([0-9]+),([0-9]+)\\]",
      convert = TRUE  # Converts to integers
    )
  message("Reshaping results...")
  gamma_star <- gamma_star_reshape(gamma_star)
  rownames(gamma_star) <- batch_names


  delta_star <- gamma_star_reshape(delta_star)
  rownames(delta_star) <- batch_names

  end_time <- Sys.time()
  message(sprintf("\nCompleted in %s", format(end_time - start_time)))
  return(list(gamma_star = gamma_star, delta_star = delta_star, fit = fit))
}

#' @rdname stan_data_prep
#'
#' @param R Numeric residual matrix with observations in rows and features in
#'   columns.
#' @param bat Factor specifying batch membership.
#' @param K_bayes Number of latent factors included in the Bayesian model.
#' @param K_target Target latent dimension used by the latent-factor prior.
#' @param subj Optional factor identifying subjects for longitudinal or
#'   repeated-measures data.
#' @param sigma_hat Feature-specific residual scale estimates used to anchor
#'   the baseline variance parameters.
#' @param sigma_logsigma0 Prior scale for baseline log-standard deviations.
#' @param sigma_orth Prior scale controlling latent-factor orthogonality.
#' @param sigma_H0 Prior scale for subject-level latent factors.
#' @param sigma_u_prior Prior scale for observation-level latent deviations.
#' @param a_alpha,b_alpha Shape parameters controlling longitudinal latent
#'   variance priors.
#' @param sigma_logdelta_prior Prior scale for batch-specific log-scale effects.
#' @param sigma_logs_prior Prior scale for cross-sectional latent
#'   heteroscedasticity.
#' @param sigma_u_batch_subj Prior scale for batch-by-subject latent variation
#'   in longitudinal models.
#' @param sigma_orth_HplusU Prior scale controlling orthogonality of combined
#'   subject- and observation-level latent effects.
#'
#' @return A list containing data, indices, and prior hyperparameters required
#'   by either the cross-sectional or longitudinal SVA-ComBat Stan model.
#'
#' @method stan_data_prep sva_combat
#' @export
stan_data_prep.sva_combat <- function(
    type,
    R,
    bat,
    K_bayes,
    K_target,
    subj = NULL,
    sigma_hat = NULL,
    sigma_logsigma0 = 0.25,
    sigma_orth = 0.2,
    sigma_H0  = 1,
    sigma_u_prior = 0.15,
    a_alpha   = 2,
    b_alpha   = 2,
    sigma_logdelta_prior = 0.25,
    sigma_logs_prior = 0.25,
    sigma_u_batch_subj = 0.10,
    sigma_orth_HplusU  = 0.20,
    ...
) {

  R <- as.matrix(R)
  N <- nrow(R)
  P <- ncol(R)

  bat <- droplevels(bat)
  B   <- nlevels(bat)

  K_bayes  <- as.integer(K_bayes)
  K_target <- as.integer(K_target)

  # sigma_hat is required for sigma0 models
  if (is.null(sigma_hat)) stop("sigma_hat is required as an anchor when using sigma0.")
  sigma_hat <- as.numeric(sigma_hat)
  if (length(sigma_hat) != P) stop("sigma_hat must have length P.")
  if (any(!is.finite(sigma_hat)) || any(sigma_hat <= 0)) {
    stop("sigma_hat must be finite and > 0.")
  }

  # ==================================================
  # LONGITUDINAL / REPEATED MEASURES CASE (sigma0 Stan)
  # ==================================================
  if (!is.null(subj)) {
    subj <- droplevels(subj)
    S    <- nlevels(subj)

    batch_int <- as.integer(bat)
    subj_int  <- as.integer(subj)

    # unique subjects per batch
    subj_by_batch <- lapply(seq_len(B), function(b) unique(subj_int[batch_int == b]))
    S_uniq_b <- vapply(subj_by_batch, length, integer(1))
    max_S_uniq_b <- max(S_uniq_b)
    if (max_S_uniq_b < 1) max_S_uniq_b <- 1L

    subj_uniq_b <- matrix(1L, nrow = B, ncol = max_S_uniq_b)
    for (b in seq_len(B)) {
      m <- S_uniq_b[b]
      if (m > 0) subj_uniq_b[b, seq_len(m)] <- subj_by_batch[[b]]
    }

    # per-(batch, subject) visit indices
    idx_list_bs <- vector("list", B)
    N_bs_mat    <- matrix(0L, nrow = B, ncol = max_S_uniq_b)

    max_N_bs <- 0L
    for (b in seq_len(B)) {
      mS <- S_uniq_b[b]
      idx_list_bs[[b]] <- vector("list", mS)
      if (mS > 0) {
        for (j in seq_len(mS)) {
          s_id <- subj_uniq_b[b, j]
          idx  <- which(batch_int == b & subj_int == s_id)
          idx_list_bs[[b]][[j]] <- idx
          N_bs_mat[b, j] <- length(idx)
          if (length(idx) > max_N_bs) max_N_bs <- length(idx)
        }
      }
    }
    if (max_N_bs < 1) max_N_bs <- 1L

    idx_bs_arr <- array(1L, dim = c(B, max_S_uniq_b, max_N_bs))
    for (b in seq_len(B)) {
      mS <- S_uniq_b[b]
      if (mS > 0) {
        for (j in seq_len(mS)) {
          idx <- idx_list_bs[[b]][[j]]
          m   <- length(idx)
          if (m > 0) idx_bs_arr[b, j, seq_len(m)] <- as.integer(idx)
        }
      }
    }

    return(list(
      N = N, P = P, B = B, K = K_bayes, S = S,
      R = R, batch = batch_int, subj = subj_int,

      sigma_hat       = sigma_hat,
      sigma_logsigma0 = sigma_logsigma0,

      sigma_mu = 10,
      a_gamma  = 2,
      b_gamma  = 2,

      a_alpha  = a_alpha,
      b_alpha  = b_alpha,

      a_kappa_base = 2,
      b_kappa_base = 2,
      K_target     = K_target,

      sigma_H0      = sigma_H0,
      sigma_u_prior = sigma_u_prior,
      sigma_orth    = sigma_orth,

      sigma_logdelta_prior = sigma_logdelta_prior,

      S_uniq_b     = as.integer(S_uniq_b),
      max_S_uniq_b = as.integer(max_S_uniq_b),
      subj_uniq_b  = array(as.integer(subj_uniq_b), dim = c(B, max_S_uniq_b)),

      max_N_bs = as.integer(max_N_bs),
      N_bs     = array(as.integer(N_bs_mat), dim = c(B, max_S_uniq_b)),
      idx_bs   = idx_bs_arr,

      sigma_u_batch_subj = sigma_u_batch_subj,
      sigma_orth_HplusU  = sigma_orth_HplusU
    ))
  }

  # ==================================================
  # CROSS-SECTIONAL CASE (sigma0 Stan)
  # ==================================================
  return(list(
    N = N, P = P, B = B, K = K_bayes,
    R = R, batch = as.integer(bat),

    sigma_hat       = sigma_hat,
    sigma_logsigma0 = sigma_logsigma0,

    sigma_mu = 1,
    a_gamma  = 2,
    b_gamma  = 2,

    a_kappa_base = 2,
    b_kappa_base = 2,
    K_target     = K_target,

    sigma_orth = sigma_orth,
    sigma_logdelta_prior = sigma_logdelta_prior,

    # keep ONLY if your CS sigma0 Stan still has sigma_logs / log_s_raw
    sigma_logs_prior = sigma_logs_prior
  ))
}

#' @rdname stan_algorithm
#'
#' @param stan_data List returned by [stan_data_prep()] for the joint
#'   SVA-ComBat model.
#' @param chains Number of MCMC chains. Defaults to `3`.
#' @param parallel_chains Number of chains run in parallel. Defaults to `3`.
#' @param posterior_draws Number of posterior draws of `R_adj` to retain for
#'   uncertainty-aware downstream inference. Set to `0` to retain only the
#'   posterior mean. Defaults to `0`.
#' @param posterior_seed Integer random seed used when subsampling retained
#'   posterior draws. Defaults to `1`.
#' @param adapt_delta Target acceptance probability for the NUTS sampler.
#'   Larger values use smaller leapfrog step sizes and can reduce divergent
#'   transitions at the cost of additional computation. Defaults to `0.98`.
#' @param iter_sampling Number of post-warmup MCMC iterations retained per
#'   chain. Defaults to `1000`.
#' @param iter_warmup Number of warmup iterations per chain used for sampler
#'   adaptation. Warmup draws are not retained for posterior inference.
#'   Defaults to `1000`.
#' @param max_treedepth Maximum tree depth allowed for each NUTS trajectory.
#'   Larger values permit longer trajectories when the posterior geometry is
#'   difficult to explore, at the cost of additional computation. Defaults to
#'   `12`.
#' @param ... Additional arguments reserved for method compatibility.
#'
#' @return A list containing:
#' \describe{
#'   \item{R_adj}{
#'     Posterior mean of the batch-adjusted residual matrix, with dimensions
#'     `N x P`.
#'   }
#'   \item{R_adj_draws}{
#'     Optional array of retained posterior draws of the adjusted residual
#'     matrix, with dimensions `S x N x P`, where `S` is the number of retained
#'     draws. `NULL` when `posterior_draws = 0`.
#'   }
#'   \item{draw_ids}{
#'     Integer indices identifying the posterior draws retained in
#'     `R_adj_draws`. Empty when `posterior_draws = 0`.
#'   }
#'   \item{H_hat}{
#'     Posterior mean estimate of the latent-factor scores. For
#'     cross-sectional models this is the posterior mean of `H`; for
#'     longitudinal models it is reconstructed from the subject-level and
#'     observation-level latent components.
#'   }
#'   \item{latent_hat}{
#'     Posterior mean of the feature-level latent contribution, with dimensions
#'     `N x P`.
#'   }
#'   \item{fit}{
#'     The fitted [cmdstanr::CmdStanMCMC] object, retained for convergence
#'     diagnostics and additional posterior summaries.
#'   }
#' }
#'
#' @details
#' MCMC sampling is performed using Stan's NUTS sampler. The `adapt_delta`
#' argument controls the target acceptance probability used during step-size
#' adaptation, while `max_treedepth` controls the maximum length of the NUTS
#' trajectory. Increasing either setting can improve sampling for difficult
#' posterior geometries, but may increase computation time.
#'
#' The `iter_warmup` argument controls the number of warmup iterations used for
#' sampler adaptation, and `iter_sampling` controls the number of retained
#' post-warmup draws per chain.
#'
#' Posterior draws of `R_adj` are only retained when `posterior_draws > 0`,
#' which can substantially reduce memory usage when uncertainty propagation is
#' not required.
#'
#' @method stan_algorithm sva_combat
#' @export
stan_algorithm.sva_combat <- function(type,
                                      stan_data,
                                      chains = 3,
                                      parallel_chains = 3,
                                      posterior_draws = 0L,
                                      posterior_seed = 1L,
                                      adapt_delta = 0.98,
                                      iter_sampling = 1000,
                                      iter_warmup = 1000,
                                      max_treedepth = 12,
                                      ...) {

  message("Starting SVA + ComBat Stan algorithm...")
  start_time <- Sys.time()

  is_longitudinal <- all(c("S", "subj") %in% names(stan_data))
  stan_model_name <- if (is_longitudinal) {"combat_sva_long"} else {"combat_sva_cross"}

  message("[1/4] Loading Stan model: ", stan_model_name, " ...")
  stan_model <- get_stan_model(stan_model_name)

  message("[2/4] Running MCMC sampling...")
  fit <- stan_model$sample(
    data            = stan_data,
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  message("[3/4] Extracting posterior summaries...")

  N <- stan_data$N
  P <- stan_data$P
  K <- stan_data$K

  posterior_draws <- as.integer(posterior_draws)
  if (length(posterior_draws) != 1L || is.na(posterior_draws) || posterior_draws < 0L) {
    stop("`posterior_draws` must be a single non-negative integer.", call. = FALSE)}

  posterior_seed <- as.integer(posterior_seed)
  if ( length(posterior_seed) != 1L || is.na(posterior_seed) ) {
    stop("`posterior_seed` must be a single integer.", call. = FALSE)
    }

  draws_R_adj <- fit$draws("R_adj", format = "draws_matrix")
  nm <- colnames(draws_R_adj)
  match_idx <- regexec("^R_adj\\[([0-9]+),([0-9]+)\\]$", nm)
  parsed <- regmatches(nm, match_idx)
  keep <- lengths(parsed) == 3L

  if (sum(keep) != N * P) {
    stop("Unable to reconstruct posterior `R_adj` draws from Stan variable names.", call. = FALSE )
    }

  row_idx <- as.integer(vapply(parsed[keep], `[[`, character(1), 2L))
  col_idx <- as.integer(vapply(parsed[keep], `[[`, character(1), 3L))
  cols <- which(keep)
  linear_idx <- row_idx + (col_idx - 1L) * N
  ord <- order(linear_idx)

  if (!identical(as.integer(linear_idx[ord]), seq_len(N * P))) {
    stop("Unexpected indexing of posterior `R_adj` elements.", call. = FALSE)
  }

  cols <- cols[ord]
  R_adj_flat_mean <- colMeans(draws_R_adj[, cols, drop = FALSE])
  R_adj_mean <- matrix(R_adj_flat_mean, nrow = N, ncol = P)
  rm(R_adj_flat_mean)

  R_adj_draws <- NULL
  draw_ids <- integer(0)

  if (posterior_draws > 0L) {
    n_available <- nrow(draws_R_adj)
    if (posterior_draws > n_available) {
      warning("Requested ", posterior_draws, " posterior draws, but only ", n_available, " are available. Retaining all available draws.", call. = FALSE)
      posterior_draws <- n_available
      }
    set.seed(posterior_seed)
    draw_ids <- if (posterior_draws == n_available) {seq_len(n_available)} else {sort(sample.int(n_available, posterior_draws, replace = FALSE))}
    selected_draws <- draws_R_adj[draw_ids, cols, drop = FALSE]
    R_adj_draws <- array(selected_draws, dim = c(length(draw_ids), N, P))
    rm(selected_draws)
  }
  rm(draws_R_adj, parsed, match_idx, keep, nm, cols, row_idx, col_idx, ord)

  draws_latent <- fit$draws("latent_hat", format = "draws_matrix")
  latent_flat_mean <- colMeans(draws_latent)
  stopifnot(length(latent_flat_mean) == N * P)
  latent_hat <- matrix(latent_flat_mean, nrow = N, ncol = P, byrow = FALSE)
  rm(draws_latent, latent_flat_mean)

  if (K == 0L) {
    H_hat <- matrix(numeric(0), nrow = N, ncol = 0L)
  } else if (is_longitudinal) {
    # H0: S x K
    S <- stan_data$S
    draws_H0 <- fit$draws("H0", format = "draws_matrix")
    H0_flat <- colMeans(draws_H0)
    stopifnot(length(H0_flat) == S * K)
    H0_hat <- matrix(H0_flat, nrow = S, ncol = K, byrow = FALSE)
    rm(draws_H0, H0_flat)

    # U_raw: N x K
    draws_U <- fit$draws("U_raw", format = "draws_matrix")
    U_flat <- colMeans(draws_U)
    stopifnot(length(U_flat) == N * K)
    U_raw_hat <- matrix(U_flat, nrow = N, ncol = K, byrow = FALSE)
    rm(draws_U, U_flat)

    # sigma_u: K
    draws_su <- fit$draws("sigma_u", format = "draws_matrix")
    sigma_u_hat <- as.numeric(colMeans(draws_su))
    stopifnot(length(sigma_u_hat) == K)
    rm(draws_su)

    # Construct per-observation Hn
    U_hat <- sweep(U_raw_hat, 2, sigma_u_hat, `*`)
    subj_int <- stan_data$subj
    H_hat <- H0_hat[subj_int, , drop = FALSE] + U_hat
    rm(H0_hat, U_raw_hat, U_hat, sigma_u_hat)
  } else {
    draws_H <- fit$draws("H", format = "draws_matrix")
    H_flat <- colMeans(draws_H)
    stopifnot(length(H_flat) == N * K)
    H_hat <- matrix(H_flat, nrow = N, ncol = K, byrow = FALSE)
    rm(draws_H, H_flat)
  }

  message("[4/4] Finished in ",
          round(difftime(Sys.time(), start_time, units = "mins"), 2),
          " minutes.")

  return(list(
    R_adj = R_adj_mean,
    R_adj_draws = R_adj_draws,
    draw_ids = draw_ids,
    H_hat = H_hat,
    latent_hat = latent_hat,
    fit = fit
  ))
}
