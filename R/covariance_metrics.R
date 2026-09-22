#' Frobenius Norm Error
#'
#' Computes the Frobenius norm of the difference between a reference covariance
#' matrix and an estimated covariance matrix.
#'
#' @param Sigma Reference covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return A numeric value giving the Frobenius norm error.
#'
#' @export
frobenius_norm <- function(Sigma, Sigma_hat) {
  sqrt(sum((Sigma - Sigma_hat)^2))
}

#' Mean Squared Covariance Error
#'
#' Computes the mean squared element-wise difference between a reference
#' covariance matrix and an estimated covariance matrix.
#'
#' @param Sigma Reference covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return A numeric value giving the mean squared covariance error.
#'
#' @export
mse_cov <- function(Sigma, Sigma_hat) {
  mean((Sigma - Sigma_hat)^2)
}

#' Spectral Norm Error
#'
#' Computes the spectral norm of the difference between a reference covariance
#' matrix and an estimated covariance matrix.
#'
#' @param Sigma Reference covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return A numeric value giving the spectral norm error.
#'
#' @export
spectral_norm <- function(Sigma, Sigma_hat) {
  vals <- eigen(Sigma - Sigma_hat, symmetric = TRUE, only.values = TRUE)$values
  max(abs(vals))
}

#' Eigenvalue Error
#'
#' Computes the mean absolute difference between the eigenvalues of a reference
#' covariance matrix and an estimated covariance matrix.
#'
#' @param Sigma Reference covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return A numeric value giving the mean absolute eigenvalue error.
#'
#' @export
eigen_error <- function(Sigma, Sigma_hat) {
  true_eigen <- eigen(Sigma, symmetric = TRUE)$values
  est_eigen <- eigen(Sigma_hat, symmetric = TRUE)$values
  mean(abs(true_eigen - est_eigen))
}

#' Evaluate Covariance Recovery
#'
#' Evaluates agreement between a reference covariance matrix and an estimated
#' covariance matrix using several complementary error measures.
#'
#' @param Sigma Reference covariance matrix.
#' @param Sigma_hat Estimated covariance matrix.
#'
#' @return A list containing Frobenius norm error, mean squared error,
#'   spectral norm error, and mean absolute eigenvalue error.
#'
#' @export
evaluate_cov_recovery <- function(Sigma, Sigma_hat) {
  list(
    Frobenius = frobenius_norm(Sigma, Sigma_hat),
    MSE = mse_cov(Sigma, Sigma_hat),
    Spectral = spectral_norm(Sigma, Sigma_hat),
    EigenError = eigen_error(Sigma, Sigma_hat)
  )
}
