test_that("residualize_one removes design-space components", {
  set.seed(1)

  x <- seq(-1, 1, length.out = 30)
  Z <- cbind(1, x)

  Y <- cbind(
    2 + 3 * x,
    -1 + 0.5 * x
  ) + matrix(rnorm(60, sd = 0.05), 30, 2)

  R <- LatentComBat:::residualize_one(Y, Z)

  expect_equal(dim(R), dim(Y))
  expect_lt(max(abs(crossprod(Z, R))), 1e-10)
})


test_that("residualize_one rejects missing values", {
  Y <- matrix(rnorm(20), 10, 2)
  Z <- cbind(1, rnorm(10))

  Y[1, 1] <- NA

  expect_error(
    LatentComBat:::residualize_one(Y, Z),
    "NA values"
  )
})


test_that("within-subject permutation preserves each subject-feature multiset", {
  set.seed(5)

  R <- matrix(seq_len(24), 8, 3)
  subj <- rep(c("s1", "s2", "s3", "s4"), each = 2)

  Rp <- LatentComBat:::permute_within_subject_featurewise(
    R,
    subj
  )

  expect_equal(dim(Rp), dim(R))

  for (s in unique(subj)) {
    idx <- which(subj == s)

    for (g in seq_len(ncol(R))) {
      expect_setequal(
        Rp[idx, g],
        R[idx, g]
      )
    }
  }
})


test_that("within-subject permutation leaves singleton subjects unchanged", {
  R <- matrix(seq_len(15), 5, 3)

  subj <- c(
    "s1",
    "s1",
    "s2",
    "s3",
    "s3"
  )

  out <- LatentComBat:::permute_within_subject_featurewise(
    R,
    subj
  )

  expect_equal(
    out[subj == "s2", , drop = FALSE],
    R[subj == "s2", , drop = FALSE]
  )

  expect_equal(dim(out), dim(R))
})


test_that("detect_K_one returns zero for a zero residual matrix", {
  R <- matrix(0, 12, 4)
  Z <- matrix(1, 12, 1)

  out <- LatentComBat:::detect_K_one(
    R,
    Z,
    B = 3L,
    alpha = 0.05
  )

  expect_identical(out$K, 0L)
  expect_true(all(is.na(out$pvals)))
})


test_that("detect_K_one returns zero when no residual degrees of freedom remain", {
  set.seed(10)

  R <- matrix(rnorm(20), 5, 4)

  # Full-rank design: rank(Z) = n
  Z <- diag(5)

  out <- LatentComBat:::detect_K_one(
    R = R,
    Z = Z,
    B = 3L,
    alpha = 0.10
  )

  expect_identical(out$K, 0L)
  expect_true(is.na(out$pvals))
  expect_true(is.list(out$svd))
})


test_that("detect_K_one returns a valid sequential p-value sequence", {
  set.seed(6)

  h <- rnorm(24)

  load <- c(
    2,
    1.5,
    1,
    0.5,
    0.2
  )

  R <- h %o% load +
    matrix(
      rnorm(24 * 5, sd = 0.25),
      24,
      5
    )

  Z <- matrix(1, 24, 1)

  out <- LatentComBat:::detect_K_one(
    R,
    Z,
    B = 8L,
    alpha = 0.2
  )

  expect_true(out$K >= 0L)
  expect_true(out$K <= 5L)

  expect_true(
    all(diff(out$pvals) >= -1e-12)
  )

  expect_true(
    all(out$pvals >= 0 & out$pvals <= 1)
  )
})


test_that("construct_svs_one returns requested number of factors", {
  set.seed(7)

  R <- matrix(
    rnorm(30 * 12),
    30,
    12
  )

  s <- svd(R)

  H <- LatentComBat:::construct_svs_one(
    R,
    K = 2L,
    svd_R = s
  )

  expect_equal(
    dim(H),
    c(30L, 2L)
  )

  expect_identical(
    colnames(H),
    c("SV1", "SV2")
  )

  expect_true(
    all(is.finite(H))
  )
})


test_that("construct_svs_one handles constant feature correlations", {
  set.seed(11)

  R <- cbind(
    constant = rep(1, 30),
    matrix(rnorm(30 * 5), 30, 5)
  )

  s <- svd(R)

  H <- suppressWarnings(
    LatentComBat:::construct_svs_one(
      R = R,
      K = 1L,
      svd_R = s
    )
  )

  expect_equal(
    dim(H),
    c(30L, 1L)
  )

  expect_true(
    all(is.finite(H))
  )
})


test_that("remove_sv_one removes variation explained by surrogate variables", {
  set.seed(20)

  n <- 30L

  H <- cbind(
    h1 = rnorm(n),
    h2 = rnorm(n)
  )

  A <- matrix(
    c(
      2, -1,
      0.5, 1.5
    ),
    nrow = 2L
  )

  Y <- H %*% A

  out <- LatentComBat:::remove_sv_one(
    Y = Y,
    H = H
  )

  expect_equal(
    dim(out),
    dim(Y)
  )

  expect_lt(
    max(abs(out)),
    1e-10
  )
})


test_that("remove_sv_one handles rank-deficient surrogate variables", {
  set.seed(21)

  n <- 25L

  h <- rnorm(n)

  H <- cbind(
    h1 = h,
    h2 = h
  )

  Y <- cbind(
    3 * h,
    -2 * h
  )

  out <- LatentComBat:::remove_sv_one(
    Y = Y,
    H = H
  )

  expect_equal(
    dim(out),
    dim(Y)
  )

  expect_true(
    all(is.finite(out))
  )
})


test_that("remove_sv_one_preserve_Z removes H while preserving Z", {
  set.seed(22)

  n <- 40L

  z <- seq(
    -1,
    1,
    length.out = n
  )

  Z <- cbind(
    "(Intercept)" = 1,
    z = z
  )

  h <- rnorm(n)

  H <- matrix(
    h,
    ncol = 1L
  )

  beta <- matrix(
    c(2, 3),
    ncol = 1L
  )

  alpha <- matrix(
    4,
    nrow = 1L
  )

  Y <- Z %*% beta +
    H %*% alpha

  out <- LatentComBat:::remove_sv_one_preserve_Z(
    Y = Y,
    Z = Z,
    H = H
  )

  expected <- Z %*% beta

  expect_equal(
    out,
    expected,
    tolerance = 1e-10
  )
})


test_that("remove_sv_one_preserve_Z returns conformable finite output", {
  set.seed(23)

  n <- 30L

  Z <- cbind(
    1,
    seq(-1, 1, length.out = n)
  )

  H <- cbind(
    rnorm(n),
    rnorm(n)
  )

  Y <- matrix(
    rnorm(n * 4),
    n,
    4
  )

  out <- LatentComBat:::remove_sv_one_preserve_Z(
    Y = Y,
    Z = Z,
    H = H
  )

  expect_equal(
    dim(out),
    dim(Y)
  )

  expect_true(
    all(is.finite(out))
  )
})


test_that("sva_per_measurement_np validates inputs", {
  set.seed(24)

  Y <- matrix(
    rnorm(40),
    10,
    4
  )

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      Z = matrix(1, 9, 1)
    ),
    "same number of rows"
  )

  Y_na <- Y
  Y_na[1, 1] <- NA

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y_na
    ),
    "Missing values"
  )

  Z_na <- matrix(
    1,
    10,
    1
  )

  Z_na[1, 1] <- NA

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      Z = Z_na
    ),
    "Missing values"
  )

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      B = 0
    ),
    "positive integer"
  )

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      alpha = 0
    ),
    "between 0 and 1"
  )

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      alpha = 1
    ),
    "between 0 and 1"
  )

  expect_error(
    LatentComBat:::sva_per_measurement_np(
      Y,
      alpha = "0.1"
    ),
    "between 0 and 1"
  )
})


test_that("sva_per_measurement_np returns original data when K is zero", {
  set.seed(25)

  Y <- matrix(
    rnorm(30),
    10,
    3
  )

  fake_detect <- function(R, Z, B, alpha) {
    list(
      K = 0L,
      pvals = c(0.8, 0.9),
      svd = svd(R)
    )
  }

  testthat::local_mocked_bindings(
    detect_K_one = fake_detect,
    .package = "LatentComBat"
  )

  out <- suppressMessages(
    LatentComBat:::sva_per_measurement_np(
      data_nb = Y,
      Z = NULL,
      B = 3L
    )
  )

  expect_equal(
    out$data_nb_sva,
    Y
  )

  expect_null(
    out$H
  )

  expect_identical(
    out$K,
    0L
  )

  expect_equal(
    out$pvals,
    c(0.8, 0.9)
  )
})


test_that("sva_per_measurement_np uses unprotected removal when preserve_Z is FALSE", {
  set.seed(26)

  n <- 20L

  Y <- matrix(
    rnorm(n * 4),
    n,
    4
  )

  Z <- cbind(
    1,
    seq(
      -1,
      1,
      length.out = n
    )
  )

  H_fake <- matrix(
    seq(
      -1,
      1,
      length.out = n
    ),
    ncol = 1L
  )

  fake_detect <- function(R, Z, B, alpha) {
    list(
      K = 1L,
      pvals = 0.01,
      svd = svd(R)
    )
  }

  fake_construct <- function(R, K, svd_R) {
    H_fake
  }

  testthat::local_mocked_bindings(
    detect_K_one = fake_detect,
    construct_svs_one = fake_construct,
    .package = "LatentComBat"
  )

  out <- suppressMessages(
    LatentComBat:::sva_per_measurement_np(
      data_nb = Y,
      Z = Z,
      B = 3L,
      preserve_Z = FALSE
    )
  )

  expected <- LatentComBat:::remove_sv_one(
    Y,
    H_fake
  )

  expect_identical(
    out$K,
    1L
  )

  expect_equal(
    out$H,
    H_fake
  )

  expect_equal(
    out$data_nb_sva,
    expected
  )

  expect_equal(
    dim(out$data_nb_sva),
    dim(Y)
  )
})


test_that("sva_per_measurement_np uses protected removal when preserve_Z is TRUE", {
  set.seed(27)

  n <- 20L

  z <- seq(
    -1,
    1,
    length.out = n
  )

  Z <- cbind(
    1,
    z
  )

  H <- matrix(
    rnorm(n),
    ncol = 1L
  )

  beta <- matrix(
    c(1, 2),
    ncol = 1L
  )

  alpha <- matrix(
    3,
    nrow = 1L
  )

  Y <- Z %*% beta +
    H %*% alpha

  fake_detect <- function(R, Z, B, alpha) {
    list(
      K = 1L,
      pvals = 0.01,
      svd = svd(R)
    )
  }

  fake_construct <- function(R, K, svd_R) {
    H
  }

  testthat::local_mocked_bindings(
    detect_K_one = fake_detect,
    construct_svs_one = fake_construct,
    .package = "LatentComBat"
  )

  out <- suppressMessages(
    LatentComBat:::sva_per_measurement_np(
      data_nb = Y,
      Z = Z,
      B = 3L,
      preserve_Z = TRUE
    )
  )

  expected <- Z %*% beta

  expect_identical(
    out$K,
    1L
  )

  expect_equal(
    out$H,
    H
  )

  expect_equal(
    out$data_nb_sva,
    expected,
    tolerance = 1e-8
  )
})


test_that("sva_per_measurement_np emits detected-K message", {
  set.seed(28)

  Y <- matrix(
    rnorm(30),
    10,
    3
  )

  fake_detect <- function(R, Z, B, alpha) {
    list(
      K = 0L,
      pvals = NA_real_,
      svd = svd(R)
    )
  }

  testthat::local_mocked_bindings(
    detect_K_one = fake_detect,
    .package = "LatentComBat"
  )

  expect_message(
    LatentComBat:::sva_per_measurement_np(
      data_nb = Y,
      B = 3L
    ),
    "Detected latent factors: K = 0"
  )
})


test_that("select_k_var_parallel returns valid component selection", {
  set.seed(29)

  n <- 40L

  h <- rnorm(n)

  Z <- cbind(
    3 * h + rnorm(n, sd = 0.2),
    2 * h + rnorm(n, sd = 0.2),
    rnorm(n),
    rnorm(n)
  )

  out <- LatentComBat:::select_k_var_parallel(
    Z_resid = Z,
    k_max = 3L,
    B_perm = 10L,
    quant = 0.90,
    seed = 1L
  )

  expect_true(
    out$k >= 0L
  )

  expect_true(
    out$k <= 3L
  )

  expect_length(
    out$d_obs,
    3L
  )

  expect_length(
    out$thr,
    3L
  )

  expect_true(
    all(is.finite(out$d_obs))
  )

  expect_true(
    all(is.finite(out$thr))
  )
})


test_that("select_k_var_parallel supports batch residualization", {
  set.seed(30)

  n <- 30L

  batch <- factor(
    rep(
      c("A", "B"),
      each = 15L
    )
  )

  X_batch <- stats::model.matrix(
    ~ batch
  )

  Z <- matrix(
    rnorm(n * 4),
    n,
    4
  )

  out <- LatentComBat:::select_k_var_parallel(
    Z_resid = Z,
    X_batch = X_batch,
    k_max = 2L,
    B_perm = 5L,
    quant = 0.90,
    seed = 2L
  )

  expect_true(
    out$k %in% 0:2
  )

  expect_length(
    out$d_obs,
    2L
  )

  expect_length(
    out$thr,
    2L
  )

  expect_true(
    all(is.finite(out$thr))
  )
})


test_that("select_k_var_parallel respects k_max", {
  set.seed(31)

  Z <- matrix(
    rnorm(40),
    10,
    4
  )

  out <- LatentComBat:::select_k_var_parallel(
    Z_resid = Z,
    k_max = 2L,
    B_perm = 5L,
    quant = 0.90,
    seed = 3L
  )

  expect_length(
    out$d_obs,
    2L
  )

  expect_length(
    out$thr,
    2L
  )

  expect_true(
    out$k >= 0L &&
      out$k <= 2L
  )
})


test_that("select_k_var_parallel handles zero candidate dimension", {
  Z <- matrix(
    numeric(0),
    nrow = 10L,
    ncol = 0L
  )

  out <- LatentComBat:::select_k_var_parallel(
    Z_resid = Z,
    k_max = 0L,
    B_perm = 3L,
    seed = 1L
  )

  expect_identical(
    out$k,
    0L
  )

  expect_identical(
    out$d_obs,
    numeric(0)
  )

  expect_identical(
    out$thr,
    numeric(0)
  )
})
