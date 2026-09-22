# LatentComBat

[![R-CMD-check](https://github.com/Zheng206/LatentComBat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zheng206/LatentComBat/actions/workflows/R-CMD-check.yaml) [![codecov](https://codecov.io/gh/Zheng206/LatentComBat/branch/main/graph/badge.svg)](https://codecov.io/gh/Zheng206/LatentComBat)

`LatentComBat` is an R package for harmonizing high-dimensional data across observed batches while accounting for latent technical variation.

The package supports two complementary workflows:

- **Sequential LatentComBat**: optionally estimates and removes latent mean structure with surrogate variables, optionally stabilizes latent variance structure, and then applies ComBat with optional covariance harmonization.
- **Joint Bayesian LatentComBat**: jointly estimates latent structure and batch effects in Stan while protecting modeled biological variation.

LatentComBat also includes pre/post harmonization diagnostics, covariance-structure diagnostics, Bayesian MCMC diagnostics, and optional uncertainty-aware downstream inference.

## When should I use each workflow?

| Goal | Suggested workflow |
|---|---|
| Standard observed-batch harmonization | Sequential, `sva = FALSE`, `var_stable = FALSE` |
| Residual latent mean shifts are suspected | Sequential, `sva = TRUE` |
| Latent heteroskedasticity is suspected | Sequential, `var_stable = TRUE` |
| Both latent mean and variance structure may remain | Sequential, `sva = TRUE`, `var_stable = TRUE` |
| Cross-feature covariance differs across batches | Add `cov = TRUE` after checking covariance diagnostics |
| Strong overlap between latent technical and biological variation is a concern | Joint Bayesian workflow |
| Harmonization uncertainty should be propagated into downstream inference | `harm_inference()` |

The choice should be guided by the scientific setting and diagnostic evidence rather than by turning on every adjustment by default.

## Installation and setup

Install LatentComBat from the package source, then load it with:

```r
library(LatentComBat)
```

For the sequential empirical-Bayes workflow, check the standard R dependencies with:

```r
setup_latentcombat(stan = FALSE)
```

The Bayesian workflow requires `cmdstanr`, a working C++ toolchain, and CmdStan:

```r
setup_latentcombat(stan = TRUE)
```

To allow LatentComBat to install missing R dependencies and CmdStan where possible:

```r
setup_latentcombat(install = TRUE, stan = TRUE)
```

## Quick start

The examples below use a small simulated dataset with three batches, biological covariates, and feature-wise batch effects.

```r
set.seed(1)

n <- 120
p <- 20

site <- factor(rep(c("A", "B", "C"), each = n / 3))
age <- scale(rnorm(n, mean = 50, sd = 12))[, 1]
sex <- factor(rbinom(n, 1, 0.5), labels = c("F", "M"))

covar <- data.frame(age = age, sex = sex)

Y <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(Y) <- paste0("feature_", seq_len(p))

# Add biological signal.
Y[, 1:5] <- Y[, 1:5] + 0.6 * age
Y[, 6:10] <- Y[, 6:10] + 0.4 * (sex == "M")

# Add observed batch mean and scale effects.
Y[site == "B", ] <- 1.15 * Y[site == "B", ] + 0.7
Y[site == "C", ] <- 0.85 * Y[site == "C", ] - 0.5
```

### 1. Sequential harmonization

```r
fit_seq <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential"
)

fit_seq
summary(fit_seq)

Y_harm <- fit_seq$harm_data
```

The harmonized matrix is stored in `fit_seq$harm_data`.

### 2. Add latent mean adjustment

Set `sva = TRUE` when diagnostics or study design suggest that important technical mean structure is not explained by the observed primary batch variable.

```r
fit_sva <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential",
  sva = TRUE
)
```

Latent-factor information is stored in:

```r
fit_sva$latent
```

Additional covariates can be protected from latent removal:

```r
bio1 <- rnorm(n)

fit_protected <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential",
  sva = TRUE,
  protected_covar = data.frame(bio1 = bio1)
)
```

### 3. Stabilize latent variance structure

Latent variance stabilization is controlled independently of latent mean adjustment:

```r
fit_var <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential",
  var_stable = TRUE
)
```

Both stages can be used together:

```r
fit_latent <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential",
  sva = TRUE,
  var_stable = TRUE
)
```

When both are enabled, latent mean adjustment is performed first, followed by latent variance stabilization and then primary-batch ComBat.

## Covariance diagnostics and CovBat

Batch effects may remain in cross-feature covariance even after location/scale correction. LatentComBat provides tools to examine this directly.

```r
plot_batch_covariance(
  data = Y,
  bat = site,
  type = "correlation"
)

batch_covariance_distance(
  data = Y,
  bat = site,
  type = "correlation",
  metric = "frobenius"
)
```

After harmonization, compare the same structure before and after adjustment:

```r
plot_batch_distance(
  raw_data = Y,
  harm_data = fit_seq$harm_data,
  bat = site,
  type = "correlation",
  metric = "frobenius"
)

compare_batch_structure(
  raw_data = Y,
  harm_data = fit_seq$harm_data,
  bat = site
)
```

If cross-feature covariance heterogeneity is an important target of harmonization, enable covariance adjustment:

```r
fit_cov <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "sequential",
  cov = TRUE,
  var_thresh = 0.95
)
```

`var_thresh` controls the cumulative variance explained by the principal-component subspace used for covariance harmonization.

## Joint Bayesian LatentComBat

The Bayesian workflow jointly estimates latent technical structure and observed batch effects rather than estimating them in separate stages.

```r
fit_bayes <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "bayesian",
  K_max = 10,
  K_buffer = 2,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.98
)
```

The fitted object reports the target and fitted latent dimensions:

```r
fit_bayes$latent$K_target
fit_bayes$latent$K_bayes
```

For Bayesian analyses, always inspect MCMC diagnostics before interpreting harmonized results:

```r
diag <- bayes_diagnostics(fit_bayes)

diag
summary(diag)

plot(diag, type = "rhat")
plot(diag, type = "ess")
plot(diag, type = "trace")
```

## Uncertainty-aware downstream inference

Harmonization is usually followed by a downstream scientific model. LatentComBat can optionally propagate harmonization uncertainty into feature-wise inference.

### Sequential workflow: full-pipeline bootstrap

```r
inf_seq <- harm_inference(
  fit_seq,
  formula = y ~ age + sex,
  B = 200,
  seed = 1
)

summary(inf_seq)
plot(inf_seq, type = "forest")
```

Each bootstrap replicate resamples observations within batch, reruns the sequential harmonization pipeline, and refits the downstream model.

The built-in sequential bootstrap currently assumes independent observational units. Longitudinal or repeated-measures analyses require cluster-aware resampling.

### Bayesian workflow: posterior propagation

Retain harmonized posterior draws during fitting:

```r
fit_bayes_draws <- com_harm(
  bat = site,
  data = Y,
  covar = covar,
  formula = y ~ age + sex,
  method = "bayesian",
  posterior_draws = 100
)

inf_bayes <- harm_inference(
  fit_bayes_draws,
  formula = y ~ age + sex
)

summary(inf_bayes)
plot(inf_bayes, type = "forest")
```

For each retained harmonized-data draw, the downstream model is refit. LatentComBat combines conditional downstream variance with between-draw variation to quantify propagated uncertainty. This is posterior propagation of harmonization uncertainty; it is not a fully Bayesian downstream model.

## Diagnostics

LatentComBat includes several diagnostic layers.

### Univariate batch diagnostics

```r
uni_test(
  bat = site,
  data = Y,
  covar = covar,
  model = stats::lm,
  formula = y ~ age + sex
)
```

These tests summarize evidence of batch-associated location and scale differences across features.

### PCA and multivariate structure

```r
pca_obj <- pca_prep(
  bat = site,
  data = Y,
  covar = covar,
  model = stats::lm,
  formula = y ~ age + sex,
  bat_adjust = FALSE
)

pca_plot(pca_obj)
```

Additional utilities include CCA diagnostics, PC/batch association summaries, covariance heatmaps, covariance-distance metrics, Box's M, and MANOVA-based batch-structure tests.

## Main fitted objects

`com_harm()` returns a fitted object with a common interface across workflows.

Common components include:

- `harm_data`: harmonized observations-by-features matrix.
- `method`: harmonization strategy used.
- `latent`: latent-factor information when applicable.
- `protection`: information on protected biological structure when applicable.
- `data_stand_result`: standardization information used during harmonization.
- `inference_data`: stored inputs used by optional uncertainty-aware inference.

Sequential fits may additionally contain `variance`, `eb_result`, or `stan_result`. Bayesian fits contain `stan_result` and may contain retained posterior harmonized-data draws.

## Recommended workflow

A practical analysis generally follows this sequence:

1. Specify the observed batch variable and biological covariates that should be preserved.
2. Examine pre-harmonization location, scale, multivariate, and covariance diagnostics.
3. Decide whether latent mean adjustment, latent variance stabilization, or covariance harmonization is scientifically justified.
4. Choose sequential or joint Bayesian harmonization.
5. For Bayesian models, verify MCMC convergence and sampler diagnostics.
6. Re-run the same diagnostics after harmonization to assess batch removal and biological preservation.
7. If downstream inferential uncertainty matters, use `harm_inference()` or a custom uncertainty-propagation procedure.

Harmonization should be treated as an iterative model-checking process rather than as a single automatic preprocessing step.

## Tutorials

The package includes focused tutorials for the main workflows:

- **Getting started with sequential LatentComBat**
- **Latent mean and variance adjustment**
- **Covariance diagnostics and CovBat**
- **Joint Bayesian LatentComBat and MCMC diagnostics**
- **Uncertainty-aware downstream inference**

See the package vignettes for complete examples.

## Package architecture

The main implementation is organized around the following modules:

```text
R/
├── com_harm.R               # public dispatcher
├── com_harm_sequential.R    # sequential workflow
├── com_harm_bayesian.R      # joint Bayesian workflow
├── latent_mean.R            # sequential SVA latent-mean adjustment
├── latent_variance.R        # latent variance stabilization
├── latent_sva.R             # SVA estimation / latent-dimension utilities
├── bio_protect.R            # protected biological subspace
├── combat_core.R            # final primary-batch ComBat/CovBat stage
├── covariance_diagnostics.R # covariance and multivariate diagnostics
├── bayesian_diagnostics.R   # MCMC diagnostics
├── inference.R              # uncertainty-aware downstream inference
└── inference-plots.R        # downstream inference visualization
```

## Development status

LatentComBat is under active development. Interfaces and defaults may change as the methodology, diagnostics, and inference tools are refined.
