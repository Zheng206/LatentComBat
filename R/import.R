#' Package imports
#'
#' Centralized import declarations for functions used without explicit
#' `package::function()` qualification throughout LatentComBat.
#'
#' @importFrom broom tidy
#' @importFrom car leveneTest
#' @importFrom dplyr %>% bind_rows filter mutate pull
#' @importFrom ggplot2 aes coord_equal element_blank element_line element_rect
#'   element_text facet_grid geom_line geom_point geom_text geom_tile ggplot labs
#'   margin scale_color_brewer scale_fill_gradientn scale_y_continuous
#'   stat_ellipse theme theme_bw theme_minimal facet_wrap vars scale_fill_gradient2
#' @importFrom grid unit
#' @importFrom pbkrtest KRmodcomp
#' @importFrom purrr map map_dfr
#' @importFrom rlang .data
#' @importFrom stats bartlett.test coef fligner.test kruskal.test p.adjust prcomp cor lm model.matrix pt var median qnorm
#'   predict resid update
#' @importFrom tidyr pivot_wider
#'
#' @keywords internal
"_PACKAGE"




#' Check and Set Up LatentComBat Dependencies
#'
#' Checks whether the R packages and system components required by
#' LatentComBat are available. When requested, missing R packages and CmdStan
#' are installed automatically where possible.
#'
#' @param install Logical indicating whether missing dependencies should be
#'   installed automatically. Defaults to `FALSE`.
#' @param stan Logical indicating whether Stan-specific dependencies should be
#'   checked. Defaults to `TRUE`.
#' @param cores Number of CPU cores used when installing CmdStan.
#' @param quiet Logical indicating whether setup messages should be suppressed.
#'
#' @return Invisibly returns a list describing the availability of required
#'   dependencies.
#'
#' @details
#' The function checks required R packages first. If `stan = TRUE`, it also
#' checks for `cmdstanr`, a working C++ toolchain, and a CmdStan installation.
#'
#' System-level compiler requirements may require manual installation on some
#' operating systems and cannot always be installed automatically from R.
#'
#' @export
setup_latentcombat <- function(
    install = FALSE,
    stan = TRUE,
    cores = 2L,
    quiet = FALSE
) {

  required_packages <- c(
    "broom",
    "car",
    "dplyr",
    "ggplot2",
    "grid",
    "lme4",
    "mgcv",
    "pbkrtest",
    "purrr",
    "rlang",
    "Rtsne",
    "scales",
    "tidyr"
  )

  if (isTRUE(stan)) {required_packages <- c(required_packages, "cmdstanr")}

  installed <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  missing_packages <- required_packages[!installed]

  if (length(missing_packages) > 0L) {
    if (!quiet) {
      message("Missing R packages: ", paste(missing_packages, collapse = ", "))
    }

    if (isTRUE(install)) {
      utils::install.packages(missing_packages, dependencies = TRUE)
      installed <- vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
      missing_packages <- required_packages[!installed]
    }
  }

  cmdstanr_available <- requireNamespace("cmdstanr", quietly = TRUE)

  toolchain_ok <- NA
  cmdstan_ok <- NA
  cmdstan_version <- NA_character_

  if (isTRUE(stan) && cmdstanr_available) {

    toolchain_ok <- tryCatch(
      {
        cmdstanr::check_cmdstan_toolchain(
          fix = isTRUE(install) && .Platform$OS.type == "windows",
          quiet = quiet
        )

        TRUE
      },
      error = function(e) {

        if (!quiet) {
          message("CmdStan C++ toolchain check failed:\n", conditionMessage(e))
        }

        FALSE
      }
    )

    version <- tryCatch(
      cmdstanr::cmdstan_version(
        error_on_NA = FALSE
      ),
      error = function(e) NULL
    )

    cmdstan_ok <- !is.null(version) && length(version) > 0L && !all(is.na(version))

    if (!isTRUE(cmdstan_ok) && isTRUE(install) && isTRUE(toolchain_ok)) {

      if (!quiet) {
        message("CmdStan not found. Installing CmdStan...")
      }

      tryCatch(
        {
          cmdstanr::install_cmdstan(cores = cores, quiet = quiet)
        },
        error = function(e) {
          if (!quiet) {
            message("CmdStan installation failed:\n", conditionMessage(e))
          }
        }
      )

      version <- tryCatch(
        cmdstanr::cmdstan_version(error_on_NA = FALSE),
        error = function(e) NULL
      )

      cmdstan_ok <- !is.null(version) && length(version) > 0L && !all(is.na(version))
    }

    if (isTRUE(cmdstan_ok)) {cmdstan_version <- paste(version, collapse = ".")}
  }

  r_packages_ok <- length(missing_packages) == 0L

  stan_ready <- if (isTRUE(stan)) {isTRUE(cmdstanr_available) && isTRUE(toolchain_ok) && isTRUE(cmdstan_ok)} else {NA}

  out <- list(
    r_packages_ok = r_packages_ok,
    missing_packages = missing_packages,
    cmdstanr_available = cmdstanr_available,
    toolchain_ok = toolchain_ok,
    cmdstan_ok = cmdstan_ok,
    cmdstan_version = cmdstan_version,
    stan_ready = stan_ready
  )

  class(out) <- "latentcombat_setup"

  if (!quiet) {

    message("")
    message("LatentComBat dependency check")
    message("-----------------------------")

    message("R packages:        ", if (r_packages_ok) "OK" else "missing dependencies")

    if (isTRUE(stan)) {
      message("cmdstanr:          ", if (cmdstanr_available) "OK" else "not installed")
      message("C++ toolchain:     ", if (isTRUE(toolchain_ok)) "OK" else "not available")
      message("CmdStan:           ", if (isTRUE(cmdstan_ok)) {paste0("OK (", cmdstan_version, ")")} else {"not installed"})
      message("Stan backend:      ", if (isTRUE(stan_ready)) "READY" else "NOT READY")
    }
  }

  invisible(out)
}


utils::globalVariables(c("density"))
