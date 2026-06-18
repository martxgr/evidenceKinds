# libraries
library(tidyverse)
library(lme4)
library(psych)
library(corrplot)
library(mclust)
library(patchwork)
library(GGally)
library(ade4)
library(lavaan)
#for the notation
options(scipen=999)

#### Clean Data ####

DVs <- c("exi_1", "avl_1", "jus_1", "inv_1", "avd_1", "mot_1", "vul_1", "sen_1",
         "really_1", "consider_1", "predict_1", "context_1", 
         "surprise_1", "bet_1", "voluntarism_1", "double_1")

d <- read.csv("lb.p1.csv") %>% 
  slice(-c(1:2)) %>% 
  mutate(domain = case_when(superstition != "" ~ "superstition",
                            religious != "" ~ "religious",
                            superlative != "" ~ "superlative",
                            political != "" ~ "political")) %>% 
  filter(age != "") %>% 
  mutate(across(all_of(DVs), as.numeric)) %>% 
  rename_all(~ gsub("_1$", "", .)) %>% 
  dplyr::select(prolific, text, 
                sen, exi, avl, jus, inv, avd, mot, vul, 
                really, consider, predict, context,
                surprise, bet, voluntarism, double, 
                age, gen, race, domain)

#reverse score items so that high scores are more evidency/really believe
d$exi <- -1*(d$exi - 100)
d$avd <- -1*(d$avd - 100)
d$vul <- -1*(d$vul - 100)
d$double <- -1*(d$double - 100)
d$context <- -1*(d$context - 100)

# prepare for clustering
d.c <- d %>% 
  dplyr::select(really, consider, predict, context,
                surprise, bet, mot, jus, double, domain) %>% 
  slice_sample(prop = 0.50)

vars <- c("really", "consider", "predict", "context",
          "surprise", "bet", "mot", "jus", "double")

d.complete <- d.c[complete.cases(d.c[, c(vars, "domain")]), ]
residuals_df <- data.frame(matrix(ncol = length(vars), nrow = nrow(d.complete)))
colnames(residuals_df) <- vars
for(var in vars) {
  mod <- lmer(paste(var, "~ 1 + (1|domain)"), data = d.complete)
  residuals_df[[var]] <- residuals(mod)
}

#### GMM ####
# configure
d.cluster <- residuals_df
# d.cluster <- d.complete %>% 
#   select(-domain)

# how many clusters?
m <- mclustBIC(d.cluster)
summary(m)
plot(m)
ICL <- mclustICL(d.cluster)
summary(ICL)

# analysis
mod1 <- Mclust(d.cluster, x = m, G=2, modelNames = "VVE")
summary(mod1, parameters = TRUE)

# plot
plot(mod1, what = "classification")

##### 2-cluster #####

# Add cluster assignments to data
# d_with_clusters_gauss <- d.cluster
# d_with_clusters_gauss$cluster <- factor(m$classification,
#                                         labels=c("Group 1", "Group 2"))
# 
# # Beautiful comprehensive pairs plot
# ggpairs(d_with_clusters_gauss, 
#         columns = c("really", "consider", "predict", "context", 
#                     "surprise", "bet", "mot", "jus", "double"),
#         aes(color = cluster, alpha = 0.6),
#         upper = list(continuous = wrap("cor", size = 3)),  # Correlations in upper
#         lower = list(continuous = wrap("points", alpha = 0.4, size = 0.8)),  # Scatter in lower
#         diag = list(continuous = wrap("densityDiag", alpha = 0.7))) +  # Densities on diagonal
#   scale_color_manual(values = c("#E69F00", "#56B4E9")) +
#   scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
#   theme_minimal() +
#   theme(strip.text = element_text(size = 8),
#         axis.text = element_text(size = 6)) +
#   ggtitle("Gaussian Mixture Model (k=2): Complete Pairwise Analysis")

##### similarity of covariance matrices #####

# Extract covariance matrices for each cluster
cov1 <- mod1$parameters$variance$sigma[,,1]
cov2 <- mod1$parameters$variance$sigma[,,2]

# Full sample covariance (pooled across groups)
cov_full <- cov(d.cluster)

# Convert to correlations (removes scale/variance effects)
cor1 <- cov2cor(cov1)
cor2 <- cov2cor(cov2)
cor_full <- cov2cor(cov_full)

# Frobenius distances on correlations
frob_cor_between <- norm(cor1 - cor2, type = "F")
frob_cor_g1_full <- norm(cor1 - cor_full, type = "F")
frob_cor_g2_full <- norm(cor2 - cor_full, type = "F")

# Normalize by number of unique correlations
n_vars <- ncol(d.cluster)
n_unique_cors <- n_vars * (n_vars - 1) / 2
norm_factor <- sqrt(n_unique_cors)

structural_similarity <- tibble(
  Comparison = c("Group 1 vs Group 2", 
                 "Group 1 vs Full",
                 "Group 2 vs Full"),
  Frobenius_Distance = c(frob_cor_between, frob_cor_g1_full, frob_cor_g2_full),
  Normalized_Distance = c(frob_cor_between, frob_cor_g1_full, frob_cor_g2_full) / norm_factor,
  Similarity_Percent = 100 * (1 - c(frob_cor_between, frob_cor_g1_full, frob_cor_g2_full) / norm_factor)
)

print(structural_similarity)

#### Permutation test ####
library(clusterGeneration)

set.seed(42)
n_sims <- 10000
n_vars <- ncol(d.cluster)
n_unique_cors <- n_vars * (n_vars - 1) / 2
norm_factor <- sqrt(n_unique_cors)

# Observed quantities from GMM
cluster_assignments <- mod1$classification
n1 <- sum(cluster_assignments == 1)
n2 <- sum(cluster_assignments == 2)
N <- n1 + n2
w1 <- n1 / N
w2 <- n2 / N

# Observed within-class variances (as SD vectors for scaling)
sd1 <- sqrt(diag(cov1))
sd2 <- sqrt(diag(cov2))

# Observed mean difference
mu_diff <- mod1$parameters$mean[, 1] - mod1$parameters$mean[, 2]

# Observed normalized Frobenius distances
cor_full <- cov2cor(cov(d.cluster))
obs_dist_g1 <- norm(cov2cor(cov1) - cor_full, type = "F") / norm_factor
obs_dist_g2 <- norm(cov2cor(cov2) - cor_full, type = "F") / norm_factor

# Analytic pooled correlation from within-class covariances and mean separation
compute_pooled_cor <- function(cov1_sim, cov2_sim) {
  # Pooled covariance = weighted within-class + between-class mean term
  cov_pool <- w1 * cov1_sim + w2 * cov2_sim + w1 * w2 * (mu_diff %o% mu_diff)
  return(cov2cor(cov_pool))
}

# Run simulations
sim_dists <- matrix(NA, nrow = n_sims, ncol = 2)

for (i in 1:n_sims) {
  # Random valid correlation matrices
  cor_sim1 <- rcorrmatrix(n_vars)
  cor_sim2 <- rcorrmatrix(n_vars)
  
  # Scale to covariance using observed variances
  D1 <- diag(sd1)
  D2 <- diag(sd2)
  cov_sim1 <- D1 %*% cor_sim1 %*% D1
  cov_sim2 <- D2 %*% cor_sim2 %*% D2
  
  # Analytic pooled correlation
  cor_pool_sim <- compute_pooled_cor(cov_sim1, cov_sim2)
  
  # Normalized Frobenius distances
  sim_dists[i, 1] <- norm(cor_sim1 - cor_pool_sim, type = "F") / norm_factor
  sim_dists[i, 2] <- norm(cor_sim2 - cor_pool_sim, type = "F") / norm_factor
}

sim_dists_g1 <- sim_dists[, 1]
sim_dists_g2 <- sim_dists[, 2]

# Separate p-values
p_g1 <- mean(sim_dists_g1 <= obs_dist_g1)
p_g2 <- mean(sim_dists_g2 <= obs_dist_g2)

cat("Observed normalized Frobenius distances:\n")
cat(sprintf("  Group 1 vs Full: %.4f  (p = %.4f)\n", obs_dist_g1, p_g1))
cat(sprintf("  Group 2 vs Full: %.4f  (p = %.4f)\n\n", obs_dist_g2, p_g2))
cat(sprintf("Null distribution summaries (n = %d simulations):\n", n_sims))
cat(sprintf("  Group 1: median = %.4f, 5th pctile = %.4f\n", median(sim_dists_g1), quantile(sim_dists_g1, 0.05)))
cat(sprintf("  Group 2: median = %.4f, 5th pctile = %.4f\n", median(sim_dists_g2), quantile(sim_dists_g2, 0.05)))

# Plot both
par(mfrow = c(1, 2))

hist(sim_dists_g1, breaks = 50, col = "gray80", border = "white",
     main = "Group 1 vs Full",
     xlab = "Normalized Frobenius Distance",
     xlim = c(min(sim_dists_g1, obs_dist_g1) - 0.02, max(sim_dists_g1, obs_dist_g1) + 0.02))
abline(v = obs_dist_g1, col = "red", lwd = 2, lty = 2)
text(obs_dist_g1, par("usr")[4] * 0.9, sprintf("Obs = %.3f\np = %.4f", obs_dist_g1, p_g1),
     col = "red", pos = 4, cex = 0.85)

hist(sim_dists_g2, breaks = 50, col = "gray80", border = "white",
     main = "Group 2 vs Full",
     xlab = "Normalized Frobenius Distance",
     xlim = c(min(sim_dists_g2, obs_dist_g2) - 0.02, max(sim_dists_g2, obs_dist_g2) + 0.02))
abline(v = obs_dist_g2, col = "red", lwd = 2, lty = 2)
text(obs_dist_g2, par("usr")[4] * 0.9, sprintf("Obs = %.3f\np = %.4f", obs_dist_g2, p_g2),
     col = "red", pos = 4, cex = 0.85)

par(mfrow = c(1, 1))

##### Permutation test: Is there one thing or two things? #####

library(parallel)
library(MASS)
library(mclust)

set.seed(42)
n_sims <- 1000
n_vars <- ncol(d.cluster)
n_unique_cors <- n_vars * (n_vars - 1) / 2
norm_factor <- sqrt(n_unique_cors)
N <- nrow(d.cluster)

# Observed pooled mean and covariance (the "one thing")
mu_pool <- colMeans(d.cluster)
cov_pool <- cov(d.cluster)
cor_pool <- cor(d.cluster)

# Observed distances
cluster_assignments <- mod1$classification
cor_obs1 <- cor(d.cluster[cluster_assignments == 1, ])
cor_obs2 <- cor(d.cluster[cluster_assignments == 2, ])
obs_dist_g1 <- norm(cor_obs1 - cor_pool, type = "F") / norm_factor
obs_dist_g2 <- norm(cor_obs2 - cor_pool, type = "F") / norm_factor

# Set up cluster
n_cores <- 10
cl <- makeCluster(n_cores)
clusterExport(cl, c("N", "mu_pool", "cov_pool", "n_vars", "norm_factor"))
clusterEvalQ(cl, {
  library(MASS)
  library(mclust)
})

# Set different seeds on each worker
clusterSetRNGStream(cl, 42)

# Run simulations in parallel
sim_results <- parLapply(cl, 1:n_sims, function(i) {
  # Generate data from one population
  X_sim <- mvrnorm(n = N, mu = mu_pool, Sigma = cov_pool)
  
  # Fit 2-class GMM
  gmm_sim <- Mclust(X_sim, G = 2, modelNames = "VVE", verbose = FALSE)
  
  # Get within-class correlations
  assignments <- gmm_sim$classification
  cor_sim1 <- cor(X_sim[assignments == 1, ])
  cor_sim2 <- cor(X_sim[assignments == 2, ])
  cor_full_sim <- cor(X_sim)
  
  # Distances
  d1 <- norm(cor_sim1 - cor_full_sim, type = "F") / norm_factor
  d2 <- norm(cor_sim2 - cor_full_sim, type = "F") / norm_factor
  
  # Assign by group size (larger = g1)
  n1 <- sum(assignments == 1)
  n2 <- sum(assignments == 2)
  if (n1 >= n2) c(d1, d2) else c(d2, d1)
})

stopCluster(cl)

# Collect results
sim_dists <- do.call(rbind, sim_results)
sim_dists_g1 <- sim_dists[, 1]
sim_dists_g2 <- sim_dists[, 2]

# p-values: how often does a single-population GMM split produce distances as large as observed?
p_g1 <- mean(sim_dists_g1 >= obs_dist_g1)
p_g2 <- mean(sim_dists_g2 >= obs_dist_g2)

cat("Observed normalized Frobenius distances:\n")
cat(sprintf("  Group 1 (n=165, larger) vs Full: %.4f  (p = %.4f)\n", obs_dist_g1, p_g1))
cat(sprintf("  Group 2 (n=45, smaller) vs Full: %.4f  (p = %.4f)\n\n", obs_dist_g2, p_g2))
cat(sprintf("Null distribution summaries (n = %d simulations):\n", n_sims))
cat(sprintf("  Larger group:  median = %.4f, 95th pctile = %.4f\n", median(sim_dists_g1), quantile(sim_dists_g1, 0.95)))
cat(sprintf("  Smaller group: median = %.4f, 95th pctile = %.4f\n", median(sim_dists_g2), quantile(sim_dists_g2, 0.95)))

# Plot
par(mfrow = c(1, 2))

hist(sim_dists_g1, breaks = 50, col = "gray80", border = "white",
     main = "Larger Group vs Full",
     xlab = "Normalized Frobenius Distance",
     xlim = c(0, max(sim_dists_g1, obs_dist_g1) + 0.02))
abline(v = obs_dist_g1, col = "red", lwd = 2, lty = 2)
text(obs_dist_g1, par("usr")[4] * 0.9, sprintf("Obs = %.3f\np = %.4f", obs_dist_g1, p_g1),
     col = "red", pos = 4, cex = 0.85)

hist(sim_dists_g2, breaks = 50, col = "gray80", border = "white",
     main = "Smaller Group vs Full",
     xlab = "Normalized Frobenius Distance",
     xlim = c(0, max(sim_dists_g2, obs_dist_g2) + 0.02))
abline(v = obs_dist_g2, col = "red", lwd = 2, lty = 2)
text(obs_dist_g2, par("usr")[4] * 0.9, sprintf("Obs = %.3f\np = %.4f", obs_dist_g2, p_g2),
     col = "red", pos = 4, cex = 0.85)

par(mfrow = c(1, 1))

##### CFA #####

# Define single factor model
cfa_model <- '
  # Single latent factor
  general =~ really + consider + predict + context + surprise + bet + mot + jus + double
'

# Fit the model
fit_cfa <- cfa(cfa_model, data = d.cluster)

# Summary with fit indices and standardized loadings
summary(fit_cfa, fit.measures = TRUE, standardized = TRUE)

# Extract key fit indices
fit_indices <- fitMeasures(fit_cfa, c("chisq", "df", "pvalue", 
                                      "cfi", "tli", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                                      "srmr", "aic", "bic"))
print(fit_indices)

fa.parallel(d.cluster, fa = "fa")
fa(d.cluster, fm = "pa", rotate = "oblimin", nfactors = 2)


##### .5 SD of the mean of the factor scores (or something like that) and plot them on the individual item distributions

#### Factor Score Range Overlay Visualization ####

# Get factor scores
factor_scores <- predict(fit_cfa)
d_with_factors <- d.cluster %>%
  mutate(factor_score = factor_scores[,1])

# Identify individuals within ±0.5 SD of mean factor score
mean_fs <- mean(d_with_factors$factor_score)
sd_fs <- sd(d_with_factors$factor_score)

lower_bound <- mean_fs - 0.5 * sd_fs
upper_bound <- mean_fs + 0.5 * sd_fs

# Flag individuals in this range
d_with_factors <- d_with_factors %>%
  mutate(in_middle_range = factor_score >= lower_bound & factor_score <= upper_bound)

# Subset for middle range
d_middle <- d_with_factors %>% filter(in_middle_range)

# Calculate proportion for scaling
prop_middle <- nrow(d_middle) / nrow(d_with_factors)

# Create plots for each item
item_names <- colnames(d.cluster)
plots_list <- list()

for(item in item_names) {
  p <- ggplot() +
    # Full distribution (background)
    geom_density(data = d.cluster, aes(x = .data[[item]]), 
                 fill = "gray70", alpha = 0.4, color = "black", linewidth = 0.8) +
    # Middle range distribution (scaled overlay)
    geom_density(data = d_middle, aes(x = .data[[item]], y = after_stat(density) * prop_middle), 
                 fill = "#56B4E9", alpha = 0.7, color = "#0072B2", linewidth = 0.8) +
    theme_minimal() +
    labs(title = item, x = NULL, y = "Density")
  
  plots_list[[item]] <- p
}

# Combine
wrap_plots(plots_list, ncol = 3) +
  plot_annotation(
    title = "Item Distributions",
    subtitle = sprintf("",
                       sum(d_with_factors$in_middle_range),
                       nrow(d_with_factors),
                       100 * mean(d_with_factors$in_middle_range))
  )

# Get factor scores
factor_scores <- predict(fit_cfa)
d_with_factors <- d.cluster %>%
  mutate(factor_score = factor_scores[,1])

# Identify individuals within ±0.5 SD of mean factor score
mean_fs <- mean(d_with_factors$factor_score)
sd_fs <- sd(d_with_factors$factor_score)
lower_bound <- mean_fs - 0.5 * sd_fs
upper_bound <- mean_fs + 0.5 * sd_fs

# Flag individuals in this range
d_with_factors <- d_with_factors %>%
  mutate(in_middle_range = factor_score >= lower_bound & factor_score <= upper_bound)

# Subset for middle range
d_middle <- d_with_factors %>% filter(in_middle_range)

# Create plots for each item
item_names <- colnames(d.cluster)
plots_list <- list()

for(item in item_names) {
  p <- ggplot() +
    # Full distribution (background)
    geom_histogram(data = d.cluster, aes(x = .data[[item]]), 
                   fill = "gray70", color = "black", alpha = 0.4, 
                   bins = 8, position = "identity") +
    # Middle range distribution (overlay)
    geom_histogram(data = d_middle, aes(x = .data[[item]]), 
                   fill = "#56B4E9", color = "#0072B2", alpha = 0.7,
                   bins = 8, position = "identity") +
    theme_minimal() +
    labs(title = item, x = NULL, y = "Count")
  
  plots_list[[item]] <- p
}

# Combine
wrap_plots(plots_list, ncol = 3) +
  plot_annotation(
    title = "Item Distributions",
    subtitle = sprintf("",
                       sum(d_with_factors$in_middle_range),
                       nrow(d_with_factors),
                       100 * mean(d_with_factors$in_middle_range))
  )

#### Correlation plots ####

colnames(residuals_df) <- vars
for(var in vars) {
  mod <- lmer(paste(var, "~ 1 + (1|domain)"), data = d.complete)
  residuals_df[[var]] <- residuals(mod)
}
c <- corr.test(residuals_df, method = "pearson", alpha = 0.05, adjust = "holm")
corrplot(c$r, 
         p.mat = c$p, 
         method = "color", 
         type = "upper", 
         insig = "pch", 
         pch = "/",
         pch.col = "grey20",
         pch.cex = 2,
         addCoef.col = "black", 
         diag = TRUE,
         title = "All domains",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         cl.pos = "n",
         mar = c(0,0,2,0)
)

# broken by domains
superstitions <- d.complete %>% 
  filter(domain == "superstition") %>% 
  select(-domain)
m <- cor(d.e)
c <- corr.test(
  d.e,
  method = "pearson",
  alpha = 0.05,
  adjust = "holm"
)
corrplot(m, 
         p.mat = c$p, 
         method = "color", 
         type = "upper", 
         insig = "pch", 
         pch = "/",
         pch.col = "grey20",
         pch.cex = 2,
         addCoef.col = "black", 
         diag = T,
         title = "Superstitions",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         cl.pos = "n",
         mar=c(0,0,2,0)
)

religion <- d.complete %>% 
  filter(domain == "religious") %>% 
  select(-domain)
m <- cor(religion)
c <- corr.test(
  religion,
  method = "pearson",
  alpha = 0.05,
  adjust = "holm"
)
corrplot(m, 
         p.mat = c$p, 
         method = "color", 
         type = "upper", 
         insig = "pch", 
         pch = "/",
         pch.col = "grey20",
         pch.cex = 2,
         addCoef.col = "black", 
         diag = T,
         title = "Religious",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         cl.pos = "n",
         mar=c(0,0,2,0)
)

superlative <- d.complete %>% 
  filter(domain == "superlative") %>% 
  select(-domain)
m <- cor(superlative)
c <- corr.test(
  superlative,
  method = "pearson",
  alpha = 0.05,
  adjust = "holm"
)
corrplot(m, 
         p.mat = c$p, 
         method = "color", 
         type = "upper", 
         insig = "pch", 
         pch = "/",
         pch.col = "grey20",
         pch.cex = 2,
         addCoef.col = "black", 
         diag = T,
         title = "Superlative",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         cl.pos = "n",
         mar=c(0,0,2,0)
)

political <- d.complete %>% 
  filter(domain == "political") %>% 
  select(-domain)
m <- cor(political)
c <- corr.test(
  political,
  method = "pearson",
  alpha = 0.05,
  adjust = "holm"
)
corrplot(m, 
         p.mat = c$p, 
         method = "color", 
         type = "upper", 
         insig = "pch", 
         pch = "/",
         pch.col = "grey20",
         pch.cex = 2,
         addCoef.col = "black", 
         diag = T,
         title = "Political",
         tl.col = "black",
         tl.cex = 0.8,
         tl.srt = 45,
         cl.pos = "n",
         mar=c(0,0,2,0)
)

#### Distributions ####
# Distributions of items (for supplement)
p1 <- ggplot(d.c, aes(Evidence))+
  geom_histogram(aes(y = ..density..), bins = 10, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  geom_vline(xintercept = -0.3, linetype = "dashed", color = "red") +
  theme_minimal()

p2 <- ggplot(d.c, aes(really))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p3 <- ggplot(d.c, aes(consider))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p4 <- ggplot(d.c, aes(predict))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p5 <- ggplot(d.c, aes(context))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p6 <- ggplot(d.c, aes(surprise))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p7 <- ggplot(d.c, aes(bet))+
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p8 <- ggplot(d.c, aes(voluntarism)) +
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

p9 <- ggplot(d.c, aes(double)) +
  geom_histogram(aes(y = ..density..), bins = 5, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  theme_minimal()

# Combine plots
(p1 | p2 | p3 | p4) /
  (p5 | p6 | p7 | p8 | p9)

# Two factor solution
p8 <- ggplot(d.c, aes(Existence)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  geom_vline(xintercept = -0.35, linetype = "dashed", color = "red") +
  theme_minimal()

p9 <- ggplot(d.c, aes(Transigence)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", alpha = 0.7, boundary = 0) +
  geom_density(color = "salmon", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal()


p8+p9

# Evidence intercept
d.o <- d.e

#pretend participant who answered 50 for everything
new <- data.frame(matrix(50, nrow = 1, ncol = ncol(d.o)))
colnames(new) <- colnames(d.o)

# Add the new row to the dataset
d_INTERCEPT <- rbind(d.o, new)

factor.scores_INTERCEPT <- fa(d_INTERCEPT, nfactors = 2, fm = "pa", rotate = "oblimin", scores = "tenBerge")
INTERCEPT <- data.frame(factor.scores_INTERCEPT[["scores"]])
colnames(INTERCEPT)[colnames(INTERCEPT) == "PA1"] <- "Evidence"
colnames(INTERCEPT)[colnames(INTERCEPT) == "PA2"] <- "Transigence"

symbolic_intercept <- INTERCEPT[nrow(INTERCEPT), 1]
objectivity_intercept <- INTERCEPT[nrow(INTERCEPT), 2]