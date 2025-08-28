library(MASS)
library(MatchIt)
library(marginaleffects)
library(grf)
library(PSweight)
library(dplyr)

set.seed(123)

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

# Modified function: now includes beta_u as an argument.
run_one_simulation <- function(N = 1100, n_rct = 100, B = 400, beta_u, delta_a, sigma_nu) {
  # 1) Generate full data
  X <- rnorm(N, 2, 0.5)
  A_obs <- rbinom(N, 1, 0.2)
  eps <- rnorm(N, 0, 0.1)
  nu  <- rnorm(N, 0, sigma_nu)
  beta_0 <- 1
  beta_a <- 2
  beta_x <- 0.5*beta_u
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

# ----- 2. Run Simulations for Different beta_u Values and Plot -----

# Simulation parameters
M <- 200          # Number of simulations per beta_u
alpha <- 0.05       # Significance level
B <- 500           # Bootstrap reps within each simulation

# Vary beta_u from 1 to 10 by 0.5
beta_u_seq <- seq(1,10, by = 0.5)
sigma_nu_seq <- c(2,3,5)
reject_rate_vec <- numeric(length(beta_u_seq))
reject_rate_matrix_u <- matrix(nrow = 3, ncol = length(beta_u_seq))

for (s in 1:3){
  sigma_nu <- sigma_nu_seq[s]
  cat("Current sigma_nu:", sigma_nu, "\n")
for (i in seq_along(beta_u_seq)) {
  beta_u <- beta_u_seq[i]
  cat("Current beta_u:", beta_u, "\n")  # Report the current value of beta_u
  
  cover_count <- 0
  all_ci_lower <- numeric(M)
  all_ci_upper <- numeric(M)
  
  for (m in 1:M) {
    # Run one simulation with the current beta_u
    cat("Current m:", m, "\n")
    diff_vec <- run_one_simulation(N = 1100, n_rct = 100, B = B, beta_u = beta_u, delta_a = 1, sigma_nu = sigma_nu)
    
    # Build an empirical (1 - alpha) CI via percentiles
    ci_lower <- quantile(diff_vec, probs = alpha / 2)
    ci_upper <- quantile(diff_vec, probs = 1 - alpha / 2)
    all_ci_lower[m] <- ci_lower
    all_ci_upper[m] <- ci_upper
    
    # Check if 0 is in this interval
    if (ci_lower <= 0 && 0 <= ci_upper) {
      cover_count <- cover_count + 1
    }
  }
  
  coverage_proportion <- cover_count / M
  reject_rate <- 1 - coverage_proportion
  reject_rate_vec[i] <- reject_rate
  
  cat("Reject proportion for beta_u", beta_u, ":", reject_rate, "\n")
}
  reject_rate_matrix_u[s,] <- reject_rate_vec
}

###################################
M <- 200          # Number of simulations per beta_u
alpha <- 0.05       # Significance level (for a 90% CI)
B <- 500           # Bootstrap reps within each simulation

# Vary delta_a from 1 to 10 by 0.5
delta_a_seq <- seq(1,10, by = 0.5)
sigma_nu_seq <- c(2,3,5)
reject_rate_vec <- numeric(length(delta_a_seq))
reject_rate_matrix_a <- matrix(nrow = 3, ncol = length(delta_a_seq))

for (s in 1:3){
  sigma_nu <- sigma_nu_seq[s]
  cat("Current sigma_nu:", sigma_nu, "\n")
  for (i in seq_along(delta_a_seq)) {
    delta_a <- delta_a_seq[i]
    cat("Current delta_a:", delta_a, "\n")  # Report the current value of delta_a
    
    cover_count <- 0
    all_ci_lower <- numeric(M)
    all_ci_upper <- numeric(M)
    
    for (m in 1:M) {
      # Run one simulation with the current beta_u
      cat("Current m:", m, "\n")
      diff_vec <- run_one_simulation(N = 1100, n_rct = 100, B = B, beta_u = 1, delta_a = delta_a, sigma_nu = sigma_nu)
      
      # Build an empirical (1 - alpha) CI via percentiles
      ci_lower <- quantile(diff_vec, probs = alpha / 2)
      ci_upper <- quantile(diff_vec, probs = 1 - alpha / 2)
      all_ci_lower[m] <- ci_lower
      all_ci_upper[m] <- ci_upper
      
      # Check if 0 is in this interval
      if (ci_lower <= 0 && 0 <= ci_upper) {
        cover_count <- cover_count + 1
      }
    }
    
    coverage_proportion <- cover_count / M
    reject_rate <- 1 - coverage_proportion
    reject_rate_vec[i] <- reject_rate
    
    cat("Reject proportion for delta_a", delta_a, ":", reject_rate, "\n")
  }
  reject_rate_matrix_a[s,] <- reject_rate_vec
}

result <- reject_rate_matrix_u
result2 <- reject_rate_matrix_a

#Figure 5
plot(beta_u_seq, result[3,], type = "l", lty = 3, col = "red",lwd = 2,
     ylim = c(0, 1),           
     xlab = expression(beta[U]),
     ylab = "Reject Proportion",
     main = "", xaxt = "n")

lines(beta_u_seq, result[2,], type = "l", lty = 2, col = "blue", lwd = 2)
lines(beta_u_seq, result[1,], type = "l", lty = 1, col ="darkgreen",  lwd = 1)
axis(side = 1, at = c(1,3,5,7,9), labels = c(1,3,5,7,9))

legend("topright",
       legend = c(
         expression(sigma[U] == 5),
         expression(sigma[U] == 3),
         expression(sigma[U] == 2)
       ),
       col      = c("red", "blue", "darkgreen"),
       lty      = c(3, 2, 1),
       bg       = "white",
       box.lty  = 1,
       cex      = 1.1,                          
       y.intersp= 1,                          
       x.intersp= 1.5,                          
       seg.len  = 2,                           
       text.width = strwidth(expression(sigma[nu] == 50))*1.2
)

#Figure 6
plot(delta_a_seq, result2[3,], type = "l", lty = 3, col = "red",lwd = 2,
     ylim = c(0, 1),          
     xlab = expression(delta[A]),
     ylab = "Reject Proportion",
     main = "", xaxt = "n")

lines(delta_a_seq, result2[2,], type = "l", lty = 2, col = "blue", lwd = 2)
lines(delta_a_seq, result2[1,], type = "l", lty = 1, col ="darkgreen",  lwd = 1)
axis(side = 1, at = c(1,3,5,7,9), labels = c(1,3,5,7,9))

legend("bottomright",
       legend = c(
         expression(sigma[U] == 5),
         expression(sigma[U] == 3),
         expression(sigma[U] == 2)
       ),
       col      = c("red", "blue", "darkgreen"),
       lty      = c(3, 2, 1),
       bg       = "white",
       box.lty  = 1,
       cex      = 1.1,                          
       y.intersp= 1,                      
       x.intersp= 1.5,                      
       seg.len  = 2,                           
       text.width = strwidth(expression(sigma[nu] == 50))*1.2
)