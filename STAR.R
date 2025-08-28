library(haven)
library(AER)
library(dplyr)
library(tidyr)
library(MASS)
library(MatchIt)
library(marginaleffects)
library(grf)
library(PSweight)
library(glmnet)


# STAR Data
data("STAR")
keep <- c("gender","ethnicity","birth","star1","read1","math1","lunch1","school1","experience1")
df2 <- STAR[, keep, drop = FALSE]

df2[] <- lapply(df2, function(x) if (is.character(x)) na_if(x, "") else x)
df2 <- df2[complete.cases(df2), ]
df2 <- subset(df2, star1 %in% c("small", "regular"))
if (is.factor(df2$star1)) df2$star1 <- droplevels(df2$star1)

# Draw RCT sample
set.seed(123) 

idx <- sample.int(nrow(df2), size = 200, replace = FALSE)
rct <- df2[idx, , drop = FALSE]   
obs <- df2[-idx, , drop = FALSE] # keep the rest as OBS

# Define (X,A,Y).
rct <- rct %>% mutate(Y = math1)
obs <- obs %>% mutate(Y = math1)
# 1) School loaction included
X_vars_unc <- c("gender","ethnicity","lunch1","school1")
# 2) School loaction NOT included
X_vars_con <- c("gender","ethnicity","lunch1")

rct$A <- ifelse(rct$star1 == "small", 1L, 0L)
obs$A <- ifelse(obs$star1 == "small", 1L, 0L)

# Induce confounding in OBS
top_half <- function(df) {
  if (nrow(df) == 0) return(df[0, , drop = FALSE])
  thr <- median(df$Y, na.rm = TRUE)
  dplyr::filter(df, Y >= thr)
}
bottom_half <- function(df) {
  if (nrow(df) == 0) return(df[0, , drop = FALSE])
  thr <- median(df$Y, na.rm = TRUE)
  dplyr::filter(df, Y <= thr)
}

# urban/suburban
urban_sub <- obs %>% filter(school1 %in% c("urban","suburban"))
urban_keep <- bind_rows(
  urban_sub %>% filter(A == 0),
  urban_sub %>% filter(A == 1) %>% top_half()
)

# rural/inner-city
rural_inner <- obs %>% filter(school1 %in% c("rural","inner-city"))
rural_keep <- bind_rows(
  rural_inner %>% filter(A == 1),
  rural_inner %>% filter(A == 0) %>% bottom_half()
)

obs <- bind_rows(urban_keep, rural_keep) %>% droplevels()


# ------------Define functions -----------
# This function computes the ATE via a stabilized IPW approach.
compute_ATE <- function(A, Y, e){
  e <- pmin(pmax(e, 1e-6), 1 - 1e-6)
  w1 <- sum(A / e)
  w2 <- sum((1 - A) / (1 - e))
  # Y_star is the IPW transformation
  Y_star <- A * Y / (w1 * e) - (1 - A) * Y / ((1 - e) * w2)
  ate <- sum(Y_star)
  return(ate)
}

# Perform bootstrap to construct ci
bootstrap <- function(rct_data, obs_data, X_vars, B){
  diff_vec <- numeric(B)
  omega_obs <- numeric(B)
  omega_rct <- numeric(B)
  ps.formula <- reformulate(X_vars, response = "A")
  
  for (b in 1:B) {
    # Resample each dataset
    rct_boot <- rct_data[sample(nrow(rct_data), replace = TRUE), ]
    obs_boot <- obs_data[sample(nrow(obs_data), replace = TRUE), ]
    
    # Estimate propensity scores
    msstat_rct <- SumStat(ps.formula, trtgrp = NULL, data = rct_boot, weight = "IPW")
    p_rct <- msstat_rct$propensity[, 2]
    
    msstat_obs <- SumStat(ps.formula, trtgrp = NULL, data = obs_boot, weight = "IPW")
    p_obs <- msstat_obs$propensity[, 2]
    
    # Compute IPW-based ATE for each
    omega_rct[b] <- compute_ATE(rct_boot$A, rct_boot$Y, p_rct)
    omega_obs[b] <- compute_ATE(obs_boot$A, obs_boot$Y, p_obs)
    
    # Store difference
    diff_vec[b] <- omega_obs[b] - omega_rct[b]
  }
  return(list(diff_vec,omega_obs,omega_rct))
}

X_vars_con <- c("gender","ethnicity","birth", "lunch1")
diff_con <- bootstrap(rct, obs, B = 1000, X_vars = X_vars_con)

X_vars_unc <- c("gender","ethnicity","birth","lunch1","school1")
diff_unc <- bootstrap(rct, obs, B = 1000, X_vars = X_vars_unc)

# Draw histogram for CI
hist(diff_con[[1]], breaks = 25, col = rgb(1, 0, 0, 0.5), 
     xlim = c(-30, 50), ylim = c(0, 200),
     xlab = "Value", ylab = "Frequency", 
     main = "")

hist(diff_unc[[1]], breaks = 25, col = rgb(0, 0, 1, 0.5), add = TRUE)

legend("topright", legend = c("No Unmeasured", "Unmeasured"),
       fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)))