test_that("diagnostic model generation and lm summary return expected residuals", {
  fx <- make_lc_fixture(n_per_batch = 15L, p = 5L)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )
  sm <- diag_model_summary(dm)

  expect_s3_class(dm, "lm_model")
  expect_equal(dim(sm$resid_add), dim(fx$data))
  expect_equal(dim(sm$resid_mul), dim(fx$data))
})

test_that("batch-only diagnostic summaries retain the observed data", {
  fx <- make_lc_fixture(n_per_batch = 8L, p = 3L)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = NULL,
    model = stats::lm
  )
  sm <- diag_model_summary(dm)

  expect_equal(sm$resid_add, fx$data)
  expect_equal(dim(sm$resid_mul), dim(fx$data))
})

test_that("PCA preparation forwards batch controls and plot ellipse option", {
  fx <- make_lc_fixture(n_per_batch = 10L, p = 4L)
  adjusted <- pca_prep(
    fx$bat, fx$data, fx$covar, stats::lm, fx$formula,
    ref.batch = "B"
  )
  unadjusted <- pca_prep(
    fx$bat, fx$data, fx$covar, stats::lm, fx$formula,
    bat_adjust = FALSE
  )

  expect_false(isTRUE(all.equal(adjusted$F_t, unadjusted$F_t)))

  plot_with_ellipse <- pca_plot(adjusted, ellipse = TRUE)
  plot_without_ellipse <- pca_plot(adjusted, ellipse = FALSE)
  is_ellipse <- function(layer) inherits(layer$stat, "StatEllipse")
  expect_equal(sum(vapply(plot_with_ellipse$layers, is_ellipse, logical(1))), 1L)
  expect_equal(sum(vapply(plot_without_ellipse$layers, is_ellipse, logical(1))), 0L)
})

test_that("feature-wise nonparametric variance diagnostics return valid tables", {
  set.seed(8)
  bat <- factor(rep(c("A", "B"), each = 20))
  R <- cbind(
    c(rnorm(20, sd = 1), rnorm(20, sd = 2)),
    c(rnorm(20), rnorm(20)),
    c(rnorm(20, sd = 0.5), rnorm(20, sd = 1.5))
  )
  colnames(R) <- paste0("f", 1:3)

  for (fun in list(kruskal_test, lv_test, bl_test, fk_test)) {
    out <- fun(R, bat)
    expect_equal(nrow(out$test_table), 3L)
    expect_match(out$perc.sig, "%$")
  }
})

test_that("anova_test detects the batch term in diagnostic models", {
  fx <- make_lc_fixture(n_per_batch = 20L, p = 5L, seed = 88)
  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )
  out <- anova_test(dm)

  expect_equal(nrow(out$test_table), ncol(fx$data))
  expect_true(all(as.numeric(out$test_table$p_value) >= 0))
  expect_match(out$perc.sig, "%$")
})

test_that("uni_test validates input dimensions", {
  fx <- make_lc_fixture(n_per_batch = 5L, p = 3L)
  expect_error(
    uni_test(fx$bat[-1], fx$data, fx$covar, stats::lm, fx$formula),
    "same number of observations"
  )
  expect_error(
    uni_test(fx$bat, as.list(fx$data), fx$covar, stats::lm, fx$formula),
    "data.frame or matrix"
  )
})

test_that("uni_test returns applicable lm diagnostics", {
  fx <- make_lc_fixture(n_per_batch = 18L, p = 4L)
  out <- uni_test(
    bat = fx$bat, data = fx$data, covar = fx$covar,
    model = stats::lm, formula = fx$formula
  )

  expect_s3_class(out, "data.frame")
  expect_true(all(c("ANOVA", "Kruskal", "Levene", "Bartlett", "FlignerKilleen") %in% names(out)))
  expect_false("KenwardRoger" %in% names(out))
  expect_equal(nrow(out), 1L)
})

test_that("GAM diagnostic model and summary return valid residuals", {
  skip_if_not_installed("mgcv")

  fx <- make_lc_fixture(
    n_per_batch = 15L,
    p = 3L,
    seed = 201
  )

  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat,
    data = fx$data,
    covar = fx$covar,
    model = mgcv::gam,
    formula = y ~ x + sex
  )

  sm <- diag_model_summary(dm)

  expect_s3_class(dm, "gam_model")
  expect_s3_class(sm, "diag_summary")

  expect_equal(
    dim(sm$resid_add),
    dim(fx$data)
  )

  expect_equal(
    dim(sm$resid_mul),
    dim(fx$data)
  )

  expect_true(all(is.finite(sm$resid_add)))
  expect_true(all(is.finite(sm$resid_mul)))
})


test_that("diag_model_gen supports models without batch adjustment", {
  fx <- make_lc_fixture(
    n_per_batch = 10L,
    p = 4L
  )

  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat,
    data = fx$data,
    covar = fx$covar,
    model = stats::lm,
    formula = y ~ x + sex,
    bat_adjust = FALSE
  )

  expect_s3_class(dm, "lm_model")
  expect_length(dm$fits, ncol(fx$data))
  expect_identical(
    names(dm$fits),
    colnames(fx$data)
  )
})


test_that("diag_model_gen supports no covariates without batch adjustment", {
  fx <- make_lc_fixture(
    n_per_batch = 8L,
    p = 3L
  )

  dm <- LatentComBat:::diag_model_gen(
    bat = fx$bat,
    data = fx$data,
    covar = NULL,
    model = stats::lm,
    formula = y ~ 1,
    bat_adjust = FALSE
  )

  expect_s3_class(dm, "lm_model")
  expect_length(dm$fits, 3L)
})

test_that("t-SNE preparation and plotting return valid objects", {
  skip_if_not_installed("Rtsne")

  fx <- make_lc_fixture(
    n_per_batch = 50L,
    p = 6L,
    seed = 202
  )

  ts <- tsne_prep(
    bat = fx$bat,
    data = fx$data,
    covar = fx$covar,
    model = stats::lm,
    formula = fx$formula
  )

  expect_s3_class(ts, "data.frame")

  expect_named(
    ts,
    c("dim1", "dim2", "batch")
  )

  expect_equal(
    nrow(ts),
    nrow(fx$data)
  )

  expect_true(all(is.finite(ts$dim1)))
  expect_true(all(is.finite(ts$dim2)))

  expect_equal(
    ts$batch,
    fx$bat
  )

  p1 <- tsne_plot(
    ts,
    ellipse = TRUE
  )

  p2 <- tsne_plot(
    ts,
    ellipse = FALSE
  )

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")

  is_ellipse <- function(layer) {
    inherits(layer$stat, "StatEllipse")
  }

  expect_equal(
    sum(vapply(p1$layers, is_ellipse, logical(1))),
    1L
  )

  expect_equal(
    sum(vapply(p2$layers, is_ellipse, logical(1))),
    0L
  )
})

test_that("t-SNE preparation rejects unsupported data types", {
  fx <- make_lc_fixture()

  expect_error(
    tsne_prep(
      bat = fx$bat,
      data = as.list(fx$data),
      covar = fx$covar,
      model = stats::lm,
      formula = fx$formula
    ),
    "data.frame, matrix"
  )
})

test_that("PCA preparation rejects unsupported data types", {
  fx <- make_lc_fixture()

  expect_error(
    pca_prep(
      bat = fx$bat,
      data = as.list(fx$data),
      covar = fx$covar,
      model = stats::lm,
      formula = fx$formula
    ),
    "data.frame, matrix"
  )
})



test_that("assoc_one returns expected association metadata", {
  set.seed(203)

  x <- rnorm(30)
  y <- 2 * x + rnorm(30, sd = 0.2)

  out <- assoc_one(
    y = y,
    x = x,
    var_name = "age",
    panel_name = "Covariates"
  )

  expect_equal(nrow(out), 1L)
  expect_identical(out$Variable, "age")
  expect_identical(out$Panel, "Covariates")
  expect_true(is.finite(out$R2))
  expect_true(out$R2 >= 0 && out$R2 <= 1)
})


test_that("sv_assoc_table combines batch and biological associations", {
  set.seed(204)

  n <- 30L

  H <- cbind(
    SV1 = rnorm(n),
    SV2 = rnorm(n)
  )

  batch_df <- data.frame(
    site = factor(rep(c("A", "B"), each = 15))
  )

  cov_df <- data.frame(
    age = rnorm(n),
    sex = factor(rep(c("F", "M"), 15))
  )

  out <- sv_assoc_table(
    H = H,
    batch_df = batch_df,
    cov_df = cov_df
  )

  expect_equal(
    nrow(out),
    2L * (1L + 2L)
  )

  expect_setequal(
    unique(out$Panel),
    c("Batch variables", "Covariates")
  )

  expect_setequal(
    unique(out$SV),
    c("SV1", "SV2")
  )

  expect_true(
    all(out$R2 >= 0 & out$R2 <= 1)
  )

  expect_s3_class(
    plot_sv_assoc_heatmap(out),
    "ggplot"
  )
})


test_that("sv_assoc_table creates SV names when absent", {
  set.seed(205)

  H <- matrix(
    rnorm(40),
    20,
    2
  )

  out <- sv_assoc_table(
    H,
    cov_df = data.frame(x = rnorm(20))
  )

  expect_setequal(
    unique(out$SV),
    c("SV1", "SV2")
  )
})


test_that("cca_diagnostic rejects constant surrogate variables", {
  H <- matrix(
    1,
    nrow = 20,
    ncol = 2
  )

  X <- data.frame(
    x = rnorm(20)
  )

  expect_error(
    cca_diagnostic(H, X),
    "No non-constant columns remain in `H`"
  )
})

test_that("cca_diagnostic rejects constant explanatory block", {
  H <- cbind(
    rnorm(20),
    rnorm(20)
  )

  X <- data.frame(
    constant = rep(1, 20)
  )

  expect_error(
    cca_diagnostic(H, X),
    "No non-constant columns remain in `X_block`"
  )
})

test_that("cca_diagnostic removes linearly dependent columns", {
  set.seed(206)

  x <- rnorm(30)
  h <- rnorm(30)

  H <- cbind(
    h1 = h,
    h2 = 2 * h,
    h3 = rnorm(30)
  )

  X <- data.frame(
    x1 = x,
    x2 = 3 * x,
    x3 = rnorm(30)
  )

  out <- cca_diagnostic(H, X)

  expect_lt(out$rank_H, ncol(H))
  expect_lt(out$rank_X, 3L)

  expect_true(
    all(out$cor >= 0 & out$cor <= 1)
  )
})


test_that("make_model_matrix_block returns zero-column matrix for constant metadata", {
  df <- data.frame(
    a = rep(1, 10),
    b = rep("x", 10)
  )

  X <- make_model_matrix_block(df)

  expect_equal(
    dim(X),
    c(10L, 0L)
  )
})


test_that("pc_assoc_table handles character metadata and default panel", {
  fx <- make_lc_fixture(
    n_per_batch = 10L,
    p = 5L
  )

  metadata <- data.frame(
    site = as.character(fx$bat),
    age = fx$covar$x,
    stringsAsFactors = FALSE
  )

  out <- pc_assoc_table(
    fx$data,
    metadata,
    n_pc = 2
  )

  expect_equal(nrow(out), 4L)
  expect_identical(
    unique(out$Panel),
    "Metadata"
  )

  expect_true(
    all(is.finite(out$R2))
  )
})
