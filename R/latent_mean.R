# Internal sequential latent-mean correction.
.com_harm_sva_mean <- function(
    data,
    bat,
    covar,
    fitted_model,
    model,
    formula,
    protected_covar = NULL,
    B = 50L,
    alpha = 0.05
) {
  Y_raw <- as.matrix(data)
  p <- ncol(Y_raw)

  resid_fits <- sapply(seq_len(p), function(g) {stats::resid(fitted_model$fits[[g]])})
  resid_fits <- as.matrix(resid_fits)

  protect <- build_Z_preserve(
    model = model,
    formula = formula,
    covar = covar,
    bat = bat,
    preserve_batch = TRUE,
    protected_covar = protected_covar,
    qr_tol = 1e-8,
    return_Q = TRUE
  )

  Z_pres <- protect$Z_raw
  sva_result <- sva_per_measurement_np(
    data_nb = resid_fits,
    Z = Z_pres,
    B = B,
    alpha = alpha,
    preserve_Z = TRUE
  )

  H_svs <- sva_result$H
  H_orth <- NULL
  beta_H <- NULL
  latent_hat <- NULL
  data_clean <- Y_raw

  if (!is.null(H_svs) && ncol(H_svs) > 0L) {
    message("Stage 1A (SVA): found ", ncol(H_svs), " hidden factors. Removing latent variation while preserving known terms + batch...")
    H_orth <- residualize_one(Y = H_svs, Z = Z_pres)
    qrH <- qr(H_orth)
    beta_H <- qr.coef(qrH, Y_raw)
    beta_H[is.na(beta_H)] <- 0
    latent_hat <- H_orth %*% beta_H
    data_clean <- Y_raw - latent_hat
  } else {
    message("Stage 1A (SVA): found 0 hidden factors; skipping latent mean removal.")
  }

  list(
    data = data_clean,
    latent = list(
      H = H_svs,
      H_orth = H_orth,
      beta = beta_H,
      fitted = latent_hat,
      sva_result = sva_result
    ),
    protection = list(
      Z_raw = protect$Z_raw,
      Z = protect$Z,
      Q = protect$Q,
      rank = protect$rank,
      preserve_batch = TRUE
    )
  )
}
