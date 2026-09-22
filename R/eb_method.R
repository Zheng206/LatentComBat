#' Create a ComBat Type Object
#'
#' Creates an object used for S3 dispatch within the ComBat estimation
#' framework.
#'
#' @param type Character string specifying the ComBat implementation type,
#'   such as `"univariate"`.
#'
#' @return An object whose class is set to `type`.
#'
#' @export
combat_type <- function(type){
  result <- type
  class(result) <- type
  return(result)
}

#' Run Empirical-Bayes Estimation
#'
#' S3 generic for estimating batch-specific location and scale effects using
#' empirical-Bayes or unshrunk estimates.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return Batch-specific location and scale estimates.
#'
#' @export
eb_algorithm <- function(type, ...) {
  UseMethod("eb_algorithm")
}

#' Estimate Additive Batch Effects
#'
#' S3 generic for estimating batch-specific additive effects from standardized
#' data.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return Batch-specific additive-effect estimates.
#'
#' @export
gamm_hat_gen <- function(type, ...){
  UseMethod("gamm_hat_gen")
}

#' Estimate Multiplicative Batch Effects
#'
#' S3 generic for estimating batch-specific scale or variance effects from
#' standardized data.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return Batch-specific scale-effect estimates.
#'
#' @export
delta_hat_gen <- function(type, ...){
  UseMethod("delta_hat_gen")
}

#' Estimate Empirical-Bayes Hyperparameters
#'
#' S3 generic for estimating method-of-moments hyperparameters used in
#' empirical-Bayes shrinkage.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return Estimated empirical-Bayes hyperparameters.
#'
#' @export
mom_calculation <- function(type, ...){
  UseMethod("mom_calculation")
}

#' Perform One Empirical-Bayes Update
#'
#' S3 generic for performing one iterative empirical-Bayes update of
#' batch-specific location and scale effects.
#'
#' @param type Object used for S3 method dispatch.
#' @param ... Additional arguments passed to the corresponding method.
#'
#' @return Updated batch-effect estimates and a convergence measure.
#'
#' @export
eb_one_iteration <- function(type, ...){
  UseMethod("eb_one_iteration")
}


#' @rdname gamm_hat_gen
#'
#' @param data_stand_result Standardized data object.
#' @param batch_result Batch-information object containing batch labels and
#'   indices.
#'
#' @return A matrix of batch-by-feature additive-effect estimates.
#'
#' @method gamm_hat_gen univariate
#' @export
gamm_hat_gen.univariate <- function(type, data_stand_result, batch_result, ...){
  data_stand <- data_stand_result$data_stand
  batch_vector <- batch_result$batch_vector
  gamma_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, mean)))
  rownames(gamma_hat) <- levels(batch_vector)
  return(gamma_hat)
}

#' @rdname delta_hat_gen
#'
#' @param data_stand_result Standardized data object.
#' @param batch_result Batch-information object containing batch labels and
#'   indices.
#' @param robust.LS Logical indicating whether a robust variance estimator is
#'   used instead of the sample variance.
#'
#' @return A matrix of batch-by-feature scale-effect estimates.
#'
#' @method delta_hat_gen univariate
#' @export
delta_hat_gen.univariate <- function(type, data_stand_result, batch_result, robust.LS = FALSE, ...){
  data_stand <- data_stand_result$data_stand
  batch_vector <- batch_result$batch_vector
  if(robust.LS){
    delta_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, biweight_midvar)))
  }else{
    delta_hat <- Reduce(rbind, by(data_stand, batch_vector, function(x) apply(x, 2, var)))
  }
  rownames(delta_hat) <- levels(batch_vector)
  return(delta_hat)
}

#' @rdname mom_calculation
#'
#' @param gamma_hat Vector of additive batch-effect estimates.
#' @param delta_hat Vector of multiplicative batch-effect estimates.
#'
#' @return A list containing method-of-moments estimates for the location and
#'   scale hyperparameters.
#'
#' @method mom_calculation univariate
#' @export
mom_calculation.univariate <- function(type, gamma_hat, delta_hat, ...){
  g_bar <- mean(gamma_hat)
  g_var <- var(gamma_hat)
  d_bar <- mean(delta_hat)
  d_var <- var(delta_hat)
  d_a <- (2 * d_var + d_bar^2)/d_var
  d_b <- (d_bar * d_var + d_bar^3)/d_var
  return(list("g_bar" = g_bar, "g_var" = g_var, "d_bar" = d_bar, "d_var" = d_var, "d_a" = d_a, "d_b" = d_b))
}


#' @rdname eb_one_iteration
#'
#' @param bdat Standardized data matrix for one batch.
#' @param g_orig Initial additive batch-effect estimates.
#' @param g_old Additive-effect estimates from the previous iteration.
#' @param d_old Scale-effect estimates from the previous iteration.
#' @param mom List of empirical-Bayes hyperparameters returned by
#'   [mom_calculation()].
#'
#' @return A list containing updated additive effects (`g_new`), updated scale
#'   effects (`d_new`), and the maximum relative change (`change`).
#'
#' @method eb_one_iteration univariate
#' @export
eb_one_iteration.univariate <- function(type, bdat, g_orig, g_old, d_old, mom, ...){
  n_b <- nrow(bdat)
  g_new <- (n_b*mom$g_var*g_orig + d_old*mom$g_bar)/(n_b*mom$g_var + d_old)
  sum2   <- colSums(sweep(bdat, 2, g_new)^2)
  d_new <- (sum2/2 + mom$d_b)/(n_b/2 + mom$d_a - 1)
  change <- max(abs(g_new - g_old)/g_old, abs(d_new - d_old)/d_old)
  return(list("g_new" = g_new, "d_new" = d_new, "change" = change))
}


#' @rdname eb_algorithm
#'
#' @param data_stand_result Standardized data object.
#' @param batch_result Batch-information object containing batch labels,
#'   indices, and batch sizes.
#' @param eb Logical indicating whether empirical-Bayes shrinkage is applied.
#'   If `FALSE`, the initial batch-effect estimates are returned.
#' @param robust.LS Logical indicating whether robust variance estimates are
#'   used when estimating multiplicative batch effects.
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma_star}{Final additive batch-effect estimates.}
#'   \item{delta_star}{Final multiplicative batch-effect estimates.}
#'   \item{gamma_hat}{Initial additive batch-effect estimates.}
#'   \item{delta_hat}{Initial multiplicative batch-effect estimates.}
#'   \item{mom}{Estimated empirical-Bayes hyperparameters, or `NULL` when
#'   `eb = FALSE`.}
#' }
#'
#' @method eb_algorithm univariate
#' @export
eb_algorithm.univariate <- function(type, data_stand_result, batch_result, eb = TRUE, robust.LS = FALSE, ...){
  data_stand <- data_stand_result$data_stand
  gamma_hat <- gamm_hat_gen(type, data_stand_result, batch_result)
  delta_hat <- delta_hat_gen(type, data_stand_result, batch_result, robust.LS = robust.LS)
  if(eb){
    gamma_star <- NULL
    delta_star <- NULL
    batch_level <- levels(batch_result$batch_vector)
    eb_result <- lapply(batch_level, function(b){
      n_b <- batch_result$n_batches[b]
      mom <- mom_calculation.univariate(type, gamma_hat[b, ], delta_hat[b, ])
      # adjust within batch
      bdat <- data_stand[batch_result$batch_index[[b]],]
      g_orig <- gamma_hat[b,]
      g_old  <- gamma_hat[b,]
      d_old  <- delta_hat[b,]
      change_old <- 1
      change <- 1
      count  <- 0
      while(change > 10e-5){
        eb_one_result <- eb_one_iteration(type, bdat, g_orig, g_old, d_old, mom)
        change <- eb_one_result$change
        if (count > 30) {
          if (change > change_old) {
            warning("Empirical Bayes step failed to converge after 30 iterations,
    	            using estimate before change between iterations increases.")
            break
          }
        }
        g_old <- eb_one_result$g_new
        d_old <- eb_one_result$d_new
        change_old <- eb_one_result$change
        count <- count+1
      }
      gamma_star <- data.frame(t(eb_one_result$g_new))
      delta_star <- data.frame(t(eb_one_result$d_new))
      return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "mom" = mom))
    })
    gamma_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$gamma_star) %>% bind_rows() %>% as.matrix()
    delta_star <- lapply(1:length(batch_level), function(i) eb_result[[i]]$delta_star) %>% bind_rows() %>% as.matrix()
    mom <- lapply(1:length(batch_level), function(i) eb_result[[i]]$mom)
    rownames(gamma_star) <- rownames(delta_star) <- batch_level
  }else{
    gamma_star <- gamma_hat
    delta_star <- delta_hat
    mom <- NULL
  }
  return(list("gamma_star" = gamma_star, "delta_star" = delta_star, "gamma_hat" = gamma_hat, "delta_hat" = delta_hat, "mom" = mom))
}

