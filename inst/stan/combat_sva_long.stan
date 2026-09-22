// combat_sva_resid_long_gaussian_H0plusUn_centered_sigma0_fastCP_weightedCenter_Ucentered_UbatchSUBJ_HplusUorth_COMBATGAMMA.stan
// Longitudinal Joint ComBat + SVA (H0 + U) with sigma0 anchoring + weighted-centering
//
// KEY CHANGE (per your recommendation):
//   - Make gamma prior ComBat-like (feature-wise shrinkage), consistent with the cross-sectional joint model.
//   - REMOVE batch-level hyper-mean mu_gamma[b] and batch-level precision/variance tau_gamma[b].
//   - NEW: tau_g[g] shared across batches for each feature g:
//         gamma_raw[b,g] ~ Normal(0, tau_g[g])
//         tau_g[g] ~ Gamma(a_gamma, b_gamma)
//
// RATIONALE:
//   - H0/U * Lambda already provides low-rank cross-feature sharing.
//   - A batch-level mean term mu_gamma[b] competes with HLambda and can absorb low-rank structure,
//     harming identifiability + leaving residual batch signal.
//
// Everything else is unchanged:
//   - sigma0[g] learned, anchored to sigma_hat[g] (log-normal)
//   - log_delta_raw[b,g] weighted-centered -> delta[b,g] has weighted geometric mean 1 per feature
//   - U within-subject centered
//   - optional orthogonality penalties (H0, U batch-subj, H0+U batch-subj)
//   - fast likelihood sliced by batch

data {
  int<lower=1> N;
  int<lower=1> P;
  int<lower=2> B;
  int<lower=0> K;
  int<lower=1> S;

  matrix[N, P] R;
  array[N] int<lower=1, upper=B> batch;
  array[N] int<lower=1, upper=S> subj;

  // ============================ CHANGED ============================
  // Orthonormal basis for the protected covariate space in observation
  // space. This may contain the original biological design, primary batch,
  // and optional extra protected covariates.
  int<lower=0> Qz;
  matrix[N, Qz] Q_protect;
  // ========================== END CHANGED ==========================

  vector<lower=0>[P] sigma_hat;

  // Prior: log(sigma0[g]) ~ Normal(log(sigma_hat[g]), sigma_logsigma0)
  real<lower=0> sigma_logsigma0;

  // NOTE: sigma_mu kept for backward compatibility (unused here)
  real<lower=0> sigma_mu;

  //   tau_g[g] ~ Gamma(a_gamma, b_gamma)
  //   gamma_raw[b,g] ~ Normal(0, tau_g[g])
  real<lower=0> a_gamma;
  real<lower=0> b_gamma;

  // subject intercept shrinkage
  real<lower=0> a_alpha;
  real<lower=0> b_alpha;

  // ARD on loadings
  real<lower=0> a_kappa_base;
  real<lower=0> b_kappa_base;
  int<lower=0> K_target;

  real<lower=0> sigma_H0;
  real<lower=0> sigma_u_prior;

  // Existing soft orthogonality strength (H0-only; optional)
  real<lower=0> sigma_orth;

  // log-delta shrinkage prior
  real<lower=0> sigma_logdelta_prior;

  // ---- Unique subjects per batch ----
  array[B] int<lower=0> S_uniq_b;
  int<lower=1> max_S_uniq_b;
  array[B, max_S_uniq_b] int<lower=1, upper=S> subj_uniq_b;

  // ---- Visit indices for each (batch, unique subject) ----
  int<lower=1> max_N_bs;
  array[B, max_S_uniq_b] int<lower=0> N_bs;
  array[B, max_S_uniq_b, max_N_bs] int<lower=1, upper=N> idx_bs;

  // (1) penalty strength for subject-weighted batch mean of U
  real<lower=0> sigma_u_batch_subj;

  // (2) penalty strength for subject-weighted batch mean of (H0 + U)
  real<lower=0> sigma_orth_HplusU;
}

transformed data {
  // ---- Precompute indices per batch for fast likelihood slicing ----
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

  // Subject indices aligned to idx_b
  array[B, max_N_b] int subj_b;
  for (b in 1:B)
    for (i in 1:N_b[b])
      subj_b[b, i] = subj[idx_b[b, i]];

  // ---- N-weighted centering weights (by observation counts) ----
  vector[B] w_b;
  row_vector[B] w_row;
  for (b in 1:B) w_b[b] = N_b[b] / (1.0 * N);
  w_row = to_row_vector(w_b);

  // ---- Precompute indices per subject for within-subject centering of U ----
  array[S] int N_s;
  int max_N_s = 0;

  for (s in 1:S) N_s[s] = 0;
  for (n in 1:N) N_s[subj[n]] += 1;
  for (s in 1:S) if (N_s[s] > max_N_s) max_N_s = N_s[s];

  array[S, max_N_s] int idx_s;
  {
    array[S] int counter_s;
    for (s in 1:S) counter_s[s] = 1;
    for (n in 1:N) {
      int ss = subj[n];
      idx_s[ss, counter_s[ss]] = n;
      counter_s[ss] += 1;
    }
  }
}

parameters {
  // Random intercepts
  matrix[S, P] alpha;
  // Multiplicative batch effects
  vector<lower=0>[P] tau_alpha;
  matrix[B, P] gamma_raw;
  vector<lower=0>[P] tau_g;

  // Subject-level latent
  matrix[S, K] H0;

  // Visit-level latent deviations
  matrix[N, K] U_raw_free;
  vector<lower=0>[K] sigma_u;

  // Loadings
  matrix[P, K] Lambda;
  vector<lower=0>[K] kappa;

  // Multiplicative batch effects
  matrix[B, P] log_delta_raw;
  real<lower=0> sigma_logdelta;

  // Baseline per-feature residual SD (learned)
  vector<lower=0>[P] sigma0;
}

transformed parameters {
  matrix[B, P] gamma;
  matrix[B, P] log_delta;
  matrix[B, P] delta;

  matrix[N, K] U_raw;   // within-subject centered

  // Weighted centering across batches per feature
  {
    row_vector[P] gamma_means    = w_row * gamma_raw;
    row_vector[P] logdelta_means = w_row * log_delta_raw;

    gamma     = gamma_raw     - rep_matrix(gamma_means, B);
    log_delta = log_delta_raw - rep_matrix(logdelta_means, B);
    delta     = exp(log_delta);
  }

  // Center U within subject
  if (K > 0) {
    U_raw = U_raw_free;
    for (s in 1:S) {
      int m = N_s[s];
      if (m > 0) {
        for (k in 1:K) {
          real mu_sk = mean(col(U_raw[idx_s[s, 1:m]], k));
          for (i in 1:m) {
            int n_idx = idx_s[s, i];
            U_raw[n_idx, k] -= mu_sk;
          }
        }
      }
    }
  } else {
    U_raw = rep_matrix(0.0, N, K);
  }
}

model {
  // alpha
  tau_alpha ~ gamma(a_alpha, b_alpha);
  for (s in 1:S)
    for (g in 1:P)
      alpha[s, g] ~ normal(0, inv_sqrt(tau_alpha[g]));

  // UPDATED gamma prior (ComBat-like)
  tau_g ~ gamma(a_gamma, b_gamma);
  for (b in 1:B)
    gamma_raw[b] ~ normal(0, tau_g);

  // H0
  if (K > 0) to_vector(H0) ~ normal(0, sigma_H0);

  // U
  if (K > 0) {
    to_vector(U_raw_free) ~ normal(0, 1);
    sigma_u ~ normal(0, sigma_u_prior);
  }

  // Lambda (ARD)
  if (K > 0) {
    for (k in 1:K) {
      real w = (k <= K_target) ? 1.0 : 5.0;
      kappa[k] ~ gamma(a_kappa_base * w, b_kappa_base * w);
    }
    for (g in 1:P)
      for (k in 1:K)
        Lambda[g, k] ~ normal(0, inv_sqrt(kappa[k]));
  }

  // log-delta
  sigma_logdelta ~ normal(0, sigma_logdelta_prior);
  to_vector(log_delta_raw) ~ normal(0, sigma_logdelta);

  // sigma0 anchored to sigma_hat (log-normal)
  target += normal_lpdf(log(sigma0) | log(sigma_hat), sigma_logsigma0);

  // Existing soft orthogonality on H0
  if (K > 0 && sigma_orth > 0) {
    for (k in 1:K) {
      for (b in 1:B) {
        int m = S_uniq_b[b];
        if (m > 0) {
          real h_mean = mean(col(H0[subj_uniq_b[b, 1:m]], k));
          target += normal_lpdf(h_mean | 0, sigma_orth);
        }
      }
    }
  }


  // SUBJECT-WEIGHTED batch-mean-zero penalty for U

  if (K > 0 && sigma_u_batch_subj > 0) {
    for (k in 1:K) {
      for (b in 1:B) {
        int mS = S_uniq_b[b];
        if (mS > 0) {
          real acc = 0;
          int  cnt = 0;
          for (j in 1:mS) {
            int m = N_bs[b, j];
            if (m > 0) {
              real ubar = mean(col(U_raw[idx_bs[b, j, 1:m]], k));
              acc += ubar;
              cnt += 1;
            }
          }
          if (cnt > 0) {
            real u_mean_subj = acc / cnt;
            target += normal_lpdf(u_mean_subj | 0, sigma_u_batch_subj);
          }
        }
      }
    }
  }

  // (2) SUBJECT-WEIGHTED batch-orthogonality on (H0 + U)

  if (K > 0 && sigma_orth_HplusU > 0) {
    for (k in 1:K) {
      for (b in 1:B) {
        int mS = S_uniq_b[b];
        if (mS > 0) {
          real acc = 0;
          int  cnt = 0;
          for (j in 1:mS) {
            int s_id = subj_uniq_b[b, j];
            int m = N_bs[b, j];
            if (m > 0) {
              real ubar = mean(col(U_raw[idx_bs[b, j, 1:m]], k));
              real hplus = H0[s_id, k] + sigma_u[k] * ubar;
              acc += hplus;
              cnt += 1;
            }
          }
          if (cnt > 0) {
            real hplus_mean = acc / cnt;
            target += normal_lpdf(hplus_mean | 0, sigma_orth_HplusU);
          }
        }
      }
    }
  }

  // Vectorized Likelihood by batch (fast)
  {
    row_vector[P] base_sd = sigma0';   // baseline, learned

    for (b in 1:B) {
      int n_c = N_b[b];
      array[n_c] int idx   = idx_b[b, 1:n_c];
      array[n_c] int s_idx = subj_b[b, 1:n_c];

      matrix[n_c, P] Rb  = R[idx];
      matrix[n_c, P] Mub = alpha[s_idx] + rep_matrix(gamma[b], n_c);

      if (K > 0) {
        matrix[n_c, K] U_sc = U_raw[idx] * diag_matrix(sigma_u);
        Mub += (H0[s_idx] + U_sc) * Lambda';
      }

      // sd: sigma0[g] * delta[b,g]
      matrix[n_c, P] Sigb = rep_matrix(delta[b] .* base_sd, n_c);

      to_vector(Rb) ~ normal(to_vector(Mub), to_vector(Sigb));
    }
  }
}

generated quantities {
  matrix[N, P] resid_raw;
  matrix[N, P] resid_adj;
  matrix[N, P] adj_plus_alpha;
  matrix[N, P] adj_final;
  matrix[N, P] R_adj;
  matrix[N, P] latent_hat;

  // ============================ CHANGED ============================
  // Decompose the fitted longitudinal latent surface into:
  //
  //   latent_protect = P_Z L = Q_protect Q_protect' L
  //   latent_remove  = (I - P_Z) L
  //
  // where L is the observation-level latent surface generated by
  // H0[subj] + U_raw * diag_matrix(sigma_u).
  //
  // Projection is performed on the full N x P latent surface rather than
  // on H0 or U separately. This is rotation-invariant and lets the
  // protected design contain both subject-level and visit-level covariates.
  matrix[N, P] latent_protect;
  matrix[N, P] latent_remove;
  // ========================== END CHANGED ==========================

  row_vector[P] delta_bar;

  // weighted geometric mean of delta per feature
  // (should be ~1 due to weighted centering)
  for (g in 1:P) {
    real tmp = 0;
    for (b in 1:B)
      tmp += (N_b[b] / (1.0 * N)) * log(delta[b, g]);
    delta_bar[g] = exp(tmp);
  }

  latent_hat     = rep_matrix(0.0, N, P);
  latent_protect = rep_matrix(0.0, N, P);  // CHANGED
  latent_remove  = rep_matrix(0.0, N, P);  // CHANGED

  // Construct observation-level latent surface:
  //
  // H_obs[n,k] =
  //   H0[subj[n],k] + sigma_u[k] * U_raw[n,k]
  //
  // L = H_obs * Lambda'
  if (K > 0) {
    matrix[N, K] H_obs;

    H_obs =
      H0[subj] +
      U_raw * diag_matrix(sigma_u);

    latent_hat =
      H_obs * Lambda';

    // ============================ CHANGED ==========================
    // Project the entire fitted latent surface onto the protected
    // observation-space design.
    //
    // Because Q_protect has orthonormal columns,
    // Q_protect * Q_protect' is the orthogonal projector onto span(Z).
    if (Qz > 0) {
      latent_protect =
        Q_protect * (Q_protect' * latent_hat);

      latent_remove =
        latent_hat - latent_protect;
    } else {
      // No protected design => recover original aggressive latent removal.
      latent_remove = latent_hat;
    }
    // ========================== END CHANGED ========================
  }

  // ================================================================
  // CHANGED reconstruction:
  //
  // 1. Remove the FULL fitted latent surface before batch scaling,
  //    just as in the original aggressive model.
  //
  // 2. Remove additive batch gamma.
  //
  // 3. Scale only the residual component by delta.
  //
  // 4. Restore alpha[subj] UN-SCALED.
  //
  // 5. Restore latent_protect AFTER scaling, so protected latent
  //    variation is not altered by multiplicative batch correction.
  //
  // Therefore:
  //
  //   R_adj_protected
  //     = R_adj_aggressive + latent_protect
  //
  // while preserving the longitudinal random intercept alpha.
  // ================================================================
  for (n in 1:N) {
    for (g in 1:P) {

      resid_raw[n, g] =
        R[n, g]
        - alpha[subj[n], g]
        - gamma[batch[n], g]
        - latent_hat[n, g];

      resid_adj[n, g] =
        resid_raw[n, g]
        / delta[batch[n], g]
        * delta_bar[g];

      adj_plus_alpha[n, g] =
        resid_adj[n, g]
        + alpha[subj[n], g];

      adj_final[n, g] =
        adj_plus_alpha[n, g]
        + latent_protect[n, g];

      R_adj[n, g] = adj_final[n, g];
    }
  }
}
