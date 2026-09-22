#' Run Univariate Batch-Effect Diagnostics
#'
#' Fits feature-wise diagnostic models and summarizes evidence of additive and
#' variance-related batch effects using a collection of univariate tests.
#'
#' @param bat Factor specifying batch or site membership.
#' @param data Numeric matrix or data frame with observations in rows and
#'   features in columns.
#' @param covar Optional data frame of observed covariates included in the
#'   diagnostic models.
#' @param model Model-fitting function, such as [stats::lm()], [mgcv::gam()],
#'   or [lme4::lmer()].
#' @param formula Optional model formula.
#' @param ref.batch Optional reference batch.
#' @param test_param List of additional test-specific settings. For mixed-effects
#'   models, `random` may be used to specify the random-effect grouping factor.
#' @param ... Additional arguments passed to [diag_model_gen()] and the
#'   underlying model-fitting function.
#'
#' @return A one-row data frame containing the percentage of features with
#'   significant batch effects for the applicable diagnostic tests. Depending
#'   on the fitted model class, the table may include:
#' \describe{
#'   \item{ANOVA}{
#'     Percentage of features showing a significant batch effect from the
#'     model-based comparison.
#'   }
#'   \item{Kruskal}{
#'     Percentage of features showing significant differences in additive
#'     residual distributions across batches.
#'   }
#'   \item{KenwardRoger}{
#'     Percentage of features showing a significant batch effect using the
#'     Kenward--Roger test for mixed-effects models.
#'   }
#'   \item{Levene}{
#'     Percentage of features showing significant residual-variance differences
#'     across batches using Levene's test.
#'   }
#'   \item{Bartlett}{
#'     Percentage of features showing significant residual-variance differences
#'     across batches using Bartlett's test.
#'   }
#'   \item{FlignerKilleen}{
#'     Percentage of features showing significant residual-variance differences
#'     across batches using the Fligner--Killeen test.
#'   }
#' }
#'
#' @details
#' Additive batch effects are assessed using the model-based batch test and,
#' for non-mixed models, the Kruskal--Wallis test applied to additive residuals.
#' Variance-related batch effects are assessed using Levene, Bartlett, and
#' Fligner--Killeen tests applied to residuals from the fitted mean model.
#'
#' For mixed-effects models, the random-effect grouping factor is inferred from
#' the fitted model when possible. It can be supplied explicitly through
#' `test_param = list(random = ...)`. The Kenward--Roger test is used for
#' mixed-effects models, whereas Kruskal--Wallis and Bartlett tests are omitted.
#'
#' @export

uni_test <- function(
    bat,
    data,
    covar,
    model,
    formula = NULL,
    ref.batch = NULL,
    test_param = list(random = NULL),
    ...
) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }
  bat <- as.factor(bat)
  if (length(bat) != nrow(data)) {
    stop("`bat` must have the same number of observations as `data`.", call. = FALSE)
  }

  diag_model <- diag_model_gen(
    bat = bat,
    data = data,
    covar = covar,
    model = model,
    formula = formula,
    ref.batch = ref.batch,
    ...
  )

  first_fit <- diag_model$fits[[1L]]
  is_lmer <- inherits(first_fit, "lmerMod")
  is_gam  <- inherits(first_fit, "gam")
  is_lm   <- inherits(first_fit, "lm") && !is_gam

  random <- test_param$random
  if (is_lmer) {
    if (is.null(random)) {
      random_terms <- lme4::findbars(stats::formula(first_fit))
      random <- unique(unlist(lapply(random_terms,function(x) all.vars(x[[3L]]))))

      if (length(random) == 0L) {
        stop("Unable to determine the random-effect grouping factor. ", "Please specify it through `test_param = list(random = ...)`.")
      }
    }

  } else {
    random <- NULL
  }

  if (is_lmer) {
    diag_summary <- diag_model_summary(diag_model, random = random)
  } else {
    diag_summary <- diag_model_summary(diag_model)
  }

  anova_result <- anova_test(diag_model)$perc.sig

  if (is_lmer) {
    kenward_result <- kenward_test(diag_model)$perc.sig
    kruskal_result <- NA_real_
  } else {
    kenward_result <- NA_real_
    kruskal_result <- kruskal_test(diag_summary$resid_add,bat)$perc.sig
  }

  levene_result <- lv_test(diag_summary$resid_mul, bat)$perc.sig
  fligner_result <- fk_test(diag_summary$resid_mul, bat)$perc.sig

  if (is_lmer) {
    bartlett_result <- NA_real_
  } else {
    bartlett_result <- bl_test(diag_summary$resid_mul, bat)$perc.sig
  }

  result_table <- data.frame(
    ANOVA = anova_result,
    Kruskal = kruskal_result,
    KenwardRoger = kenward_result,
    Levene = levene_result,
    Bartlett = bartlett_result,
    FlignerKilleen = fligner_result,
    row.names = NULL,
    check.names = FALSE
  )

  result_table <- result_table[, colSums(!is.na(result_table)) > 0, drop = FALSE]
  return(result_table)
}
