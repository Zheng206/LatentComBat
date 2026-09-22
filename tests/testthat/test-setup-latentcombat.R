test_that("setup_latentcombat returns a structured dependency report without Stan", {
  out <- setup_latentcombat(stan = FALSE, quiet = TRUE)

  expect_s3_class(out, "latentcombat_setup")
  expect_named(
    out,
    c("r_packages_ok", "missing_packages", "cmdstanr_available", "toolchain_ok",
      "cmdstan_ok", "cmdstan_version", "stan_ready")
  )
  expect_true(is.logical(out$r_packages_ok) && length(out$r_packages_ok) == 1L)
  expect_true(is.character(out$missing_packages))
  expect_true(is.na(out$stan_ready))
})
