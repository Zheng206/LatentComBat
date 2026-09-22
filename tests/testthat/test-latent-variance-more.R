test_that("latent variance stabilization handles K_var = 0", {
  fx <- make_lc_fixture(n_per_batch = 6L, p = 3L)
  br <- batch_matrix(fx$bat)

  fake_select <- function(Z_resid, k_max, B_perm, quant, seed) {
    list(k = 0L, observed = numeric(), threshold = numeric())
  }

  testthat::local_mocked_bindings(
    select_k_var_parallel = fake_select,
    .package = "LatentComBat"
  )

  out <- suppressMessages(LatentComBat:::.stabilize_latent_variance(
    data = fx$data,
    bat = fx$bat,
    covar = fx$covar,
    batch_result = br,
    model = stats::lm,
    formula = fx$formula
  ))

  expect_equal(out$K_var, 0L)
  expect_equal(dim(out$data), dim(fx$data))
  expect_equal(dim(out$ell_latent), dim(fx$data))
  expect_true(all(out$ell_latent == 0))
  expect_equal(dim(out$W), c(nrow(fx$data), 0L))
  expect_equal(dim(out$Psi), c(0L, ncol(fx$data)))
  expect_true(all(is.finite(out$data)))
})
