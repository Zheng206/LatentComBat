# Internal latent variance stabilization used by the sequential pipeline.
.stabilize_latent_variance <- function(
    data,
    bat,
    covar,
    batch_result,
    model,
    formula,
    ...
) {
  EM_max <- 100L
  tol_ell <- 5e-3
  damp_ell <- 0.5
  damp_lat <- 0.1
  ridge <- 0.1
  eps_var <- NULL
  K_var <- 3L
  cap_lat <- 2.0
  best_obj <- Inf
  no_improve <- 0L
  patience_obj <- 8L
  min_improve <- 1e-4
  tol_latent <- 1e-3
  patience_lat <- 10L
  latent_stable <- 0L
  k_max_auto <- 8L
  B_perm_auto <- 30L
  quant_auto <- 0.95
  k_cap <- 5L

  Y_em <- as.matrix(data)
  N <- nrow(Y_em)
  P <- ncol(Y_em)

  message(sprintf("Starting variance stabilization (N=%d, P=%d)...", N, P))


  Xb <- stats::model.matrix(~ bat - 1)
  B <- ncol(Xb)
  b_int <- as.integer(bat)
  w_b <- as.numeric(table(bat)) / N

  message("Pre-pass for K_var selection ", "(unweighted mean fit)...")

  fitted_model0 <- model_fitting(
    data = Y_em,
    batch = batch_result$batch_matrix,
    covar = covar,
    model = model,
    formula = formula,
    batch_include = TRUE,
    weights = NULL,
    ...
  )

  mu_hat0 <- get_mu_hat(fitted_model0)
  e_hat0 <- Y_em - mu_hat0

  med0 <- stats::median(as.numeric(e_hat0^2), na.rm = TRUE)

  if (!is.finite(med0) || med0 <= 0) {med0 <- 1}

  eps_var <- 1e-8 * med0
  Z0 <- log(e_hat0^2 + eps_var)

  coef_bg0 <- qr.solve(Xb, Z0)

  if (nrow(coef_bg0) != B) {coef_bg0 <- t(coef_bg0)}

  c_g0 <- as.numeric(t(w_b) %*% coef_bg0)
  theta_bg0 <- sweep(coef_bg0, 2, c_g0, "-")
  ell_global0 <- matrix(rep(c_g0, each = N), N, P, byrow = FALSE)
  ell_batch0 <- theta_bg0[b_int,,drop = FALSE]
  Z_resid0 <- Z0 - (ell_global0 + ell_batch0)

  sel <- select_k_var_parallel(
    Z_resid0,
    k_max = k_max_auto,
    B_perm = B_perm_auto,
    quant = quant_auto,
    seed = 1L
  )

  K_var <- sel$k
  K_var <- max(0L, min(K_var, k_cap))

  message(sprintf(paste0("Auto-selected K_var = %d ", "(parallel analysis, q=%.2f, B=%d)."), K_var, quant_auto, B_perm_auto))

  ell_global <- ell_global0
  ell_batch <- ell_batch0

  if (K_var == 0L) {
    message("K_var = 0: no latent variance factors detected; ", "skipping latent variance stabilization.")

    ell_latent <- matrix(0, N, P)
    W <- matrix(0, N, 0L)
    Psi <- matrix(0, 0L, P)
  } else {
    sv0 <- svd(Z_resid0, nu = K_var, nv = K_var)
    W <- sv0$u[, seq_len(K_var), drop = FALSE] %*% diag(sv0$d[seq_len(K_var)], K_var, K_var)
    Psi <- t(sv0$v[,seq_len(K_var),drop = FALSE])

    for (bb in seq_len(B)) {
      idx <- which(b_int == bb)
      if (length(idx) > 0) {
        W[idx,] <- sweep(W[idx,,drop = FALSE], 2, colMeans(W[idx,,drop = FALSE]), "-")
      }
    }

    ell_latent <- W %*% Psi
    ell_latent <- pmin(pmax(ell_latent,-cap_lat),cap_lat)
  }

  ell_prev <- ell_global + ell_batch + ell_latent
  latent_prev <- norm(ell_latent, "F")

  for (it in seq_len(EM_max)) {
    ell_total <- ell_global + ell_batch + ell_latent
    W_weights <- exp(-ell_total)

    fitted_model_em <- model_fitting(
      data = Y_em,
      batch = batch_result$batch_matrix,
      covar = covar,
      model = model,
      formula = formula,
      batch_include = TRUE,
      weights = W_weights,
      ...
    )

    mu_hat <- get_mu_hat(fitted_model_em)
    e_hat <- Y_em - mu_hat

    Z <- log(e_hat^2 + eps_var)
    coef_bg <- qr.solve(Xb, Z)

    if (nrow(coef_bg) != B) {
      coef_bg <- t(coef_bg)
    }

    c_g <- as.numeric(t(w_b) %*% coef_bg)
    theta_bg <- sweep(coef_bg, 2, c_g, "-")
    ell_global_new <- matrix(rep(c_g, each = N), N, P, byrow = FALSE)

    ell_batch_new <- theta_bg[b_int,,drop = FALSE]

    if (K_var == 0L) {
      ell_latent_new <- matrix(0, N, P)
      W_new <- W
      Psi_new <- Psi
    } else {
      Z_resid <- Z - (ell_global_new + ell_batch_new)
      PsiPsiT <- Psi %*% t(Psi) + ridge * diag(K_var)
      W_new <- Z_resid %*% t(Psi) %*% solve(PsiPsiT)

      for (bb in seq_len(B)) {
        idx <- which(b_int == bb)
        if (length(idx) > 0) {
          W_new[idx,] <- sweep(W_new[idx,,drop = FALSE], 2, colMeans(W_new[idx,,drop = FALSE]), "-")
        }
      }

      WtW <- t(W_new) %*% W_new + ridge * diag(K_var)
      Psi_new <- solve(WtW) %*% t(W_new) %*% Z_resid
      ell_latent_new <- W_new %*% Psi_new
    }

    ell_latent <- (1 - damp_lat) * ell_latent + damp_lat * ell_latent_new

    ell_latent <- pmin(pmax(ell_latent,-cap_lat),cap_lat)
    W <-(1 - damp_lat) * W + damp_lat * W_new
    Psi <- (1 - damp_lat) * Psi + damp_lat * Psi_new
    ell_global <- (1 - damp_ell) * ell_global + damp_ell * ell_global_new
    ell_batch <- (1 - damp_ell) * ell_batch + damp_ell * ell_batch_new
    ell_total <- ell_global + ell_batch + ell_latent
    obj <- 0.5 * sum(ell_total + e_hat^2 * exp(-ell_total))


    if (obj < best_obj - min_improve) {
      best_obj <- obj
      no_improve <- 0L
    } else {
      no_improve <- no_improve + 1L
    }

    d_ell <- norm(ell_total - ell_prev, "F") / max(1e-8, norm(ell_prev, "F"))
    ell_prev <- ell_total

    latent_now <- norm(ell_latent, "F")
    rel_latent <- abs(latent_now - latent_prev) / max(1e-8, latent_prev)
    latent_prev <- latent_now

    if (rel_latent < tol_latent) {
      latent_stable <- latent_stable + 1L
    } else {
      latent_stable <- 0L
    }

    message(sprintf(paste0("  EM iter %03d | ", "NLL=%.4f | ", "d_ell=%.2e | ", "frob_WPsi=%.2f | ", "d_latent=%.2e | ", "stab=%d | ", "K_var=%d"),
        it, obj, d_ell, latent_now, rel_latent, latent_stable, K_var))

    if (latent_stable >= patience_lat && d_ell < tol_ell) {
      message(sprintf("  OK: stopping (latent stabilized) at iter %d", it))
      break
    }

    if (no_improve >= patience_obj && d_ell < tol_ell) {
      message(sprintf("  OK: stopping (plateau) at iter %d", it))
      break
    }
  }


  sigma_stabilize <- exp(ell_latent / 2)

  if (!exists("mu_hat")) {
    mu_hat <- mu_hat0
  }

  data <- mu_hat + (Y_em - mu_hat) / sigma_stabilize

  if (!exists("theta_bg")) {
    theta_bg <- theta_bg0
  }

  log_delta_bg <- 0.5 * theta_bg
  message("Variance stabilization complete ", "(latent variance only).")
  message("Proceeding to Stage 2 (ComBat)...")

  list(
    data = data,
    K_var = K_var,
    W = W,
    Psi = Psi,
    ell_global = ell_global,
    ell_batch = ell_batch,
    ell_latent = ell_latent,
    log_delta_bg = log_delta_bg,
    selection = sel
  )
}
