data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=2> B;
  int<lower=0> K;

  matrix[N, P] R;
  array[N] int<lower=1, upper=B> batch;
  int<lower=0> Qz;
  matrix[N, Qz] Q_protect;
  vector<lower=0>[P] sigma_hat;
  real<lower=0> sigma_logsigma0;

  real<lower=0> sigma_mu;
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;
  real<lower=0> a_kappa_base;
  real<lower=0> b_kappa_base;
  int<lower=0> K_target;

  // Retained for backward compatibility with the current R interface.
  // The revised model uses exact within-batch centering of latent scores,
  // so the former soft-centering penalty controlled by sigma_orth is no
  // longer needed.
  real<lower=0> sigma_orth;

  real<lower=0> sigma_logdelta_prior;
  real<lower=0> sigma_logs_prior;
}

transformed data {
  array[B] int N_b;
  int max_N_b = 0;

  for (b in 1:B) N_b[b] = 0;
  for (n in 1:N) N_b[batch[n]] += 1;
  for (b in 1:B) if (N_b[b] > max_N_b) max_N_b = N_b[b];

  array[B, max_N_b] int idx_b;
  {
    array[B] int counter;
    for (b in 1:B) counter[b] = 1;
    for (n in 1:N) {
      int bb = batch[n];
      idx_b[bb, counter[bb]] = n;
      counter[bb] += 1;
    }
  }

  vector[B] w_b;
  row_vector[B] w_row;
  for (b in 1:B) w_b[b] = N_b[b] / (1.0 * N);
  w_row = to_row_vector(w_b);
}

parameters {
  // Mean / observed-batch model
  vector[P] mu_g;
  matrix[B, P] gamma_raw;
  vector<lower=0>[P] tau_g;

  // Latent mean structure
  matrix[N, K] H_raw;
  matrix[P, K] Lambda;
  vector<lower=0>[K] kappa;

  // Non-centered batch residual-scale hierarchy
  matrix[B, P] z_log_delta;
  real<lower=0> sigma_logdelta;

  // Non-centered latent batch-scale hierarchy
  matrix[B, K] z_log_s;
  real<lower=0> sigma_logs;

  // Feature residual scales
  vector<lower=0>[P] sigma0;
}

transformed parameters {
  matrix[B, P] gamma;

  matrix[B, P] log_delta_raw;
  matrix[B, P] log_delta;
  matrix[B, P] delta;

  matrix[B, K] log_s_raw;
  matrix[B, K] log_s;
  matrix[B, K] s;

  matrix[N, K] H_centered;
  matrix[N, K] H;

  // Weighted sum-to-zero constraints for observed batch effects.
  {
    row_vector[P] gamma_means = w_row * gamma_raw;
    gamma = gamma_raw - rep_matrix(gamma_means, B);
  }

  // Non-centered residual batch scaling followed by the same weighted
  // sum-to-zero identifiability constraint as before.
  log_delta_raw = sigma_logdelta * z_log_delta;
  {
    row_vector[P] logdelta_means = w_row * log_delta_raw;
    log_delta = log_delta_raw - rep_matrix(logdelta_means, B);
    delta = exp(log_delta);
  }

  if (K > 0) {
    // Non-centered latent batch scaling.
    log_s_raw = sigma_logs * z_log_s;
    {
      row_vector[K] logs_means = w_row * log_s_raw;
      log_s = log_s_raw - rep_matrix(logs_means, B);
      s = exp(log_s);
    }

    // Exact within-batch centering of the latent scores.
    // This removes the location tradeoff between H*Lambda', mu_g, and gamma.
    H_centered = H_raw;
    for (b in 1:B) {
      int n_count = N_b[b];
      for (k in 1:K) {
        real hbar = mean(col(H_raw[idx_b[b, 1:n_count]], k));
        for (j in 1:n_count)
          H_centered[idx_b[b, j], k] -= hbar;
      }
    }

    // Apply batch-specific latent scaling after centering.  Because each
    // centered column has batch mean zero, scaling preserves that constraint.
    for (b in 1:B) {
      int n_count = N_b[b];
      H[idx_b[b, 1:n_count]] =
        H_centered[idx_b[b, 1:n_count]] .* rep_matrix(s[b], n_count);
    }
  }
}

model {
  // Mean / observed batch effects
  mu_g ~ normal(0, sigma_mu);
  tau_g ~ gamma(a_gamma, b_gamma);
  for (b in 1:B)
    gamma_raw[b] ~ normal(0, tau_g);

  // Latent subject scores
  if (K > 0)
    to_vector(H_raw) ~ std_normal();

  // Latent loading shrinkage.
  // IMPORTANT: buffer factors now truly receive stronger shrinkage.
  // In the previous code both gamma shape and rate were multiplied by w,
  // which left E[kappa] unchanged.  Here only the shape is multiplied.
  if (K > 0) {
    for (k in 1:K) {
      real w = (k <= K_target) ? 1.0 : 5.0;
      kappa[k] ~ gamma(a_kappa_base * w, b_kappa_base);
      Lambda[, k] ~ normal(0, inv_sqrt(kappa[k]));
    }
  }

  // Non-centered hierarchical scale priors.
  sigma_logs ~ normal(0, sigma_logs_prior);
  if (K > 0)
    to_vector(z_log_s) ~ std_normal();

  sigma_logdelta ~ normal(0, sigma_logdelta_prior);
  to_vector(z_log_delta) ~ std_normal();

  sigma0 ~ lognormal(log(sigma_hat), sigma_logsigma0);

  // Likelihood
  {
    row_vector[P] base_sd = sigma0';

    for (b in 1:B) {
      int n_count = N_b[b];

      matrix[n_count, P] Rb = R[idx_b[b, 1:n_count]];

      matrix[n_count, P] Mub = rep_matrix(mu_g', n_count);
      Mub += rep_matrix(gamma[b], n_count);
      if (K > 0)
        Mub += H[idx_b[b, 1:n_count]] * Lambda';

      row_vector[P] sigma_row = delta[b] .* base_sd;
      matrix[n_count, P] Sigb = rep_matrix(sigma_row, n_count);

      to_vector(Rb) ~ normal(to_vector(Mub), to_vector(Sigb));
    }
  }
}

generated quantities {
  matrix[N, P] R_adj;
  matrix[N, P] latent_hat;
  matrix[N, P] latent_protect;
  matrix[N, P] latent_remove;
  row_vector[P] delta_bar;

  latent_hat     = rep_matrix(0, N, P);
  latent_protect = rep_matrix(0, N, P);
  latent_remove  = rep_matrix(0, N, P);

  if (K > 0) {
    latent_hat = H * Lambda';

    if (Qz > 0) {
      latent_protect = Q_protect * (Q_protect' * latent_hat);
      latent_remove  = latent_hat - latent_protect;
    } else {
      latent_remove = latent_hat;
    }
  }

  for (g in 1:P) {
    real tmp = 0;
    for (b in 1:B)
      tmp += w_b[b] * log(delta[b, g]);
    delta_bar[g] = exp(tmp);
  }

  for (b in 1:B) {
    int n_count = N_b[b];

    matrix[n_count, P] Rb = R[idx_b[b, 1:n_count]];
    matrix[n_count, P] HL = latent_hat[idx_b[b, 1:n_count]];
    matrix[n_count, P] LP = latent_protect[idx_b[b, 1:n_count]];

    for (i in 1:n_count) {
      row_vector[P] resid = Rb[i] - mu_g' - gamma[b] - HL[i];

      R_adj[idx_b[b, i]] =
        (resid ./ delta[b]) .* delta_bar + mu_g' + LP[i];
    }
  }
}
