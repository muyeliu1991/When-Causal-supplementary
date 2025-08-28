library(MASS)
library(MatchIt)
library(marginaleffects)
library(grf)
library(PSweight)
library(dplyr)

# ----- 1. Define Functions -----

# This function computes the ATE via a stabilized IPW approach.
compute_ATE <- function(A, Y, e){
  w1 <- sum(A / e)
  w2 <- sum((1 - A) / (1 - e))
  # Y_star is the IPW transformation
  Y_star <- A * Y / (w1 * e) - (1 - A) * Y / ((1 - e) * w2)
  ate <- sum(Y_star)
  return(ate)
}

# This function runs one simulation:
#  - Generates data
#  - Splits into RCT data (n=100) and Obs data (n=2000)
#  - Performs bootstrap to get dist. of diff = tau_rct - tau_obs
#  - Returns the vector of 'diff' from bootstrap
run_one_simulation <- function(N = 2100, n_rct = 100, B = 1000, delta_a, beta_u) {
  # 1) Generate full data
  X <- rnorm(N, 2, 0.5)
  A_obs <- rbinom(N, 1, 0.2)
  eps <- rnorm(N, 0, 0.1)
  sigma_u <- 0.5
  nu  <- rnorm(N, 0, sigma_u)
  beta_0 <- 1
  beta_a <- 2
  beta_x <- 0.5
  delta_0 <- 1
  delta_x <- 1
  U <- delta_0 + delta_a * A_obs + delta_x*X + nu  # latent/confounder
  
  # 2) Select RCT subset
  sample_indices <- sample(seq_len(N), n_rct)
  
  # RCT data (initially copy from full)
  rct_data <- data.frame(
    X   = X[sample_indices],
    A   = A_obs[sample_indices],
    eps = eps[sample_indices],
    nu  = nu[sample_indices],
    U   = U[sample_indices]
  )
  # Overwrite A with new Bernoulli(0.5)
  rct_data$A <- rbinom(n_rct, 1, 0.5)
  
  # 3) Construct outcome for RCT & Observational
  
  # RCT outcome
  rct_data$Y <- beta_0 +
    beta_a * rct_data$A +
    beta_x * rct_data$X +
    beta_u * rct_data$U +
    rct_data$eps
  
  # Observational data
  obs_indices <- setdiff(seq_len(N), sample_indices)
  obs_data <- data.frame(
    X   = X[obs_indices],
    A   = A_obs[obs_indices],  # keep the original A ~ Bernoulli(0.2)
    eps = eps[obs_indices],
    nu  = nu[obs_indices],
    U   = U[obs_indices]
  )
  obs_data$Y <- beta_0 +
    beta_a * obs_data$A +
    beta_x * obs_data$X +
    beta_u * obs_data$U +
    obs_data$eps
  
  # 4) Do bootstrap
  diff_vec <- numeric(B)
  ps.formula <- A ~ X
  
  for (b in 1:B) {
    # Resample each dataset
    rct_boot <- rct_data[sample(n_rct, replace = TRUE), ]
    obs_boot <- obs_data[sample(nrow(obs_data), replace = TRUE), ]
    
    # Estimate propensity scores
    msstat_rct <- SumStat(ps.formula, trtgrp = NULL, data = rct_boot, weight = "IPW")
    p_rct <- msstat_rct$propensity[, 2]
    
    msstat_obs <- SumStat(ps.formula, trtgrp = NULL, data = obs_boot, weight = "IPW")
    p_obs <- msstat_obs$propensity[, 2]
    
    # Compute IPW-based ATE for each
    tau_rct <- compute_ATE(rct_boot$A, rct_boot$Y, p_rct)
    tau_obs <- compute_ATE(obs_boot$A, obs_boot$Y, p_obs)
    
    # Store difference
    diff_vec[b] <- tau_rct - tau_obs
  }
  
  return(diff_vec)
}

# ----- 2. Run Multiple Simulations & Check Coverage -----

# Number of total simulations
M <- 1000           # e.g., run 50 times (you can make this larger)
alpha <- 0.1     # significance level
B <- 1000         # bootstrap reps within each simulation

cover_count <- 0
all_ci_lower <- numeric(M)
all_ci_upper <- numeric(M)

for (m in 1:M) {
  cat("Current m:", m, "\n")
  # 1) Run one simulation => get bootstrap distribution of diff
  diff_vec <- run_one_simulation(N = 2100, n_rct = 100, B = B, beta_u=2, delta_a = 0)
  
  # 2) Build an empirical (1 - alpha) = 95% CI via percentiles
  ci_lower <- quantile(diff_vec, probs = alpha / 2)
  ci_upper <- quantile(diff_vec, probs = 1 - alpha / 2)
  
  all_ci_lower[m] <- ci_lower
  all_ci_upper[m] <- ci_upper
  
  # 3) Check if 0 is in this interval
  if (ci_lower <= 0 && 0 <= ci_upper) {
    cover_count <- cover_count + 1
  }
}

# Final coverage proportion
coverage_proportion <- cover_count / M
cat("Out of", M, "simulations, 0 was in the (1-alpha) CI", cover_count, "times.\n")
cat("Type I error:", 1 - coverage_proportion, "\n", "at alpha:", alpha , "\n")