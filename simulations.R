# ============================================================
# Goal: make plot(mod1, "classification") show *more covariance structure*
# (i.e., class-specific ellipses / shapes, not just mean shifts).
#
# Fix: simulate cats/dogs from *different 8D covariance matrices* (not iid items).
# This increases visible variance/covariance differences in the mclust classification plot.
# ============================================================

library(tidyverse)
library(mclust)

set.seed(123)

# ---- 1) 8 features (same names as before) ----
feat <- c("size","obedience","agility","sociability",
          "independence","trainability","grooming_need","vocality")
p <- length(feat)

# ---- 2) MVN generator for p-dim without extra packages ----
rmvnorm_p <- function(n, mu, Sigma) {
  Z <- matrix(rnorm(n * length(mu)), nrow = n, ncol = length(mu))
  X <- Z %*% chol(Sigma)
  sweep(X, 2, mu, "+")
}

# ---- 3) Helpers to build *valid* correlation/covariance matrices ----
# Create a random SPD matrix with controllable "covariance richness"
random_spd <- function(p, richness = 1.2) {
  # richness > 1 increases spread of eigenvalues (more varied shapes)
  A <- matrix(rnorm(p * p), p, p)
  S <- crossprod(A)  # SPD
  # reshape eigen-spectrum
  e <- eigen(S)
  vals <- e$values
  vals2 <- (vals / mean(vals))^richness
  S2 <- e$vectors %*% diag(vals2) %*% t(e$vectors)
  # stabilize
  S2 + diag(1e-6, p)
}

cov_from_spd <- function(S_spd, sds) {
  # convert SPD "shape" to correlation, then to covariance with desired SDs
  Dinv <- diag(1 / sqrt(diag(S_spd)), p)
  R <- Dinv %*% S_spd %*% Dinv
  D <- diag(sds, p)
  D %*% R %*% D
}

# ---- 4) Means: keep earlier intent (some big separations, some overlap) ----
mu_cat <- c(
  size =  0.00,
  obedience = -1.10,
  agility =  0.90,
  sociability = -0.50,
  independence =  1.00,
  trainability = -0.60,
  grooming_need =  0.80,
  vocality =  0.30
)

mu_dog <- c(
  size =  0.35,
  obedience =  1.10,
  agility =  0.20,
  sociability =  0.80,
  independence = -0.20,
  trainability =  0.70,
  grooming_need =  0.10,
  vocality =  0.00
)

# ---- 5) SDs: allow some heteroscedasticity too (helps mclust see structure) ----
sd_cat <- c(1.10, 0.75, 0.85, 0.95, 0.80, 0.95, 0.90, 1.00)
sd_dog <- c(0.90, 0.70, 1.05, 0.85, 0.95, 0.80, 1.00, 0.85)

# ---- 6) Build class-specific covariance matrices with more varied correlation structure ----
# Richness knobs: increase for more dramatic covariance differences
S_cat <- random_spd(p, richness = 1.35)
S_dog <- random_spd(p, richness = 1.55)

Sigma_cat <- cov_from_spd(S_cat, sd_cat)
Sigma_dog <- cov_from_spd(S_dog, sd_dog)

# Optional: force a couple intuitive correlations to make the story coherent
# (and still SPD because we're editing gently)
# Example: dogs show size~obedience +, cats show size~obedience ~0 or slight -
Sigma_cat[1,2] <- Sigma_cat[2,1] <- 0.00 * sd_cat[1] * sd_cat[2]
Sigma_dog[1,2] <- Sigma_dog[2,1] <- 0.25 * sd_dog[1] * sd_dog[2]

# Re-stabilize after manual edits (ensure SPD)
make_spd <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  e$values[e$values < 1e-6] <- 1e-6
  e$vectors %*% diag(e$values) %*% t(e$vectors)
}
Sigma_cat <- make_spd(Sigma_cat)
Sigma_dog <- make_spd(Sigma_dog)

# ---- 7) Simulate data ----
n_cat <- 700
n_dog <- 700

X_cat <- rmvnorm_p(n_cat, mu = mu_cat[feat], Sigma = Sigma_cat)
X_dog <- rmvnorm_p(n_dog, mu = mu_dog[feat], Sigma = Sigma_dog)

d8 <- bind_rows(
  as_tibble(X_cat) %>% setNames(feat) %>% mutate(group_true = "cat"),
  as_tibble(X_dog) %>% setNames(feat) %>% mutate(group_true = "dog")
) %>% mutate(group_true = factor(group_true, levels = c("cat","dog")))

X <- as.matrix(d8[, feat])


# ---- 8) Fit 2-Gaussian mclust and plot classification ----
mod1 <- Mclust(X, G = 2)

# This is the plot you asked about: should now show more covariance variation
plot(mod1, what = "classification", colors = c("gray20", "gray80"))

# ---- 9) (Optional) If you want the covariance contrast even stronger ----
# - increase richness (e.g., 1.8 / 2.2)
# - increase SD differences between sd_cat and sd_dog
# - or manually inject a few higher correlations in one class and not the other
#
# Example quick tweaks:
# S_cat <- random_spd(p, richness = 1.2)
# S_dog <- random_spd(p, richness = 2.0)
# sd_dog <- sd_dog * c(1,1,1.2,1,1.1,1,1.2,1)
