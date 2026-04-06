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
old <- c("jus_1", "mot_1", "consider_1", "predict_1", "context_1", "surprise_1", "bet_1", "double_1")

new <- c("dishonest1_1", "dishonest2_1", "serious1_1", "serious2_1", "sense_1", "single_1", "really_1", "deeper_1")

d <- read.csv("rb.p2.csv") %>% 
  slice(-c(1:2)) %>% 
  mutate(domain = case_when(geography != "" ~ "geography",
                            religious != "" ~ "religious",
                            number != "" ~ "number",
                            political != "" ~ "political",
                            maybe != "" ~ "maybe")) %>% 
  filter(age != "") %>% 
  mutate(across(all_of(old), as.numeric),
         across(all_of(new), as.numeric)) %>% 
  rename_all(~ gsub("_1$", "", .)) %>% 
  dplyr::select(prolific, jus, mot, consider, predict, context, surprise, bet, double,
                dishonest1, dishonest2, serious1, serious2, sense, single, really, deeper, domain) %>% 
  filter(prolific != "a")

# reverse score items so that high scores are more evidency/really believe
d$double <- -1*(d$double - 100)
d$context <- -1*(d$context - 100)
d$sense <- -1*(d$sense - 100)
d$dishonest1 <- -1*(d$dishonest1 - 100)
d$dishonest2 <- -1*(d$dishonest2 - 100)
d$serious2 <- -1*(d$serious2 - 100)

################################################################################

d <- d %>%
  mutate(
    dishonest = coalesce(dishonest1, dishonest2),
    serious   = coalesce(serious1, serious2)
  )

d %>% 
  pivot_longer(cols = all_of(c("dishonest","serious",
                               "sense","single","really","deeper")),
               names_to = "item",
               values_to = "value") %>% 
  ggplot(aes(x = value, color = item)) +
  geom_density() +
  facet_wrap(~ domain) +
  theme_minimal()

################################################################################



believe <- d %>% 
  filter(!(domain %in% c("geography", "number", "maybe"))) %>%
  dplyr::select(serious, really, deeper, dishonest, sense, single)

m <- cor(believe)
c <- corr.test(believe, method = "pearson", alpha = 0.05, adjust = "holm")
corrplot(m, p.mat = c$p, method = "number", type = "upper", insig = "pch", pch = "/",
         pch.col = "grey20", pch.cex = 2, addCoef.col = "black", diag = T, title = "",
         tl.col = "black", tl.cex = 0.8, tl.srt = 45, cl.pos = "n", mar=c(0,0,2,0))
fa.parallel(believe, fm = "ml", fa = "fa")
fa(believe, nfactors = 1, rotate = "oblimin", fm = "ml")
fa(believe, nfactors = 2, rotate = "oblimin", fm = "ml")


really_believe <- d %>% 
  filter(!(domain %in% c("geography", "number", "maybe"))) %>% 
  dplyr::select(serious, really, deeper)

cor(really_believe, use = "pairwise.complete.obs")
psych::alpha(really_believe)

just_believe <- d %>% 
  filter(!(domain %in% c("geography", "number", "maybe"))) %>% 
  dplyr::select(dishonest, sense, single)
cor(just_believe, use = "pairwise.complete.obs")
psych::alpha(just_believe)

test <- d %>% 
  mutate(
    really_believe = rowMeans(select(., serious, really, deeper), na.rm = TRUE),
    believe        = rowMeans(select(., dishonest, sense, single), na.rm = TRUE)
  ) %>% 
  mutate(mundane = case_when(domain %in% c("geography", "maybe", "number") ~ "mundane",
                             domain %in% c("religious", "political") ~ "not mundane"))

cor(test$believe, test$really_believe)



ggplot(test, aes(believe, really_believe, color = mundane)) +
  geom_point(size = 4, alpha = 0.2) + 
  theme_classic()

ggplot(test, aes(believe, really_believe)) +
  geom_bin2d(aes(fill = after_stat(count / sum(count))), bins = 4) +
  facet_wrap(~ mundane) +
  scale_fill_continuous(labels = scales::percent) +
  theme_classic()


#### Old items ####

# data with just religious/political
d_rich <- d %>%
  filter(!(domain %in% c("geography", "number", "maybe"))) %>%
  select(really, consider, predict, context, surprise, bet, mot, jus, double)

# Define single factor model
cfa_model <- '
  # Single latent factor
  general =~ really + consider + predict + context + surprise + bet + mot + jus + double
'

# Fit the model
fit_cfa <- cfa(cfa_model, data = d_rich)

# Summary with fit indices and standardized loadings
summary(fit_cfa, fit.measures = TRUE, standardized = TRUE)

fa.parallel(d_rich, fa = "fa")
fa(d_rich, fm = "pa", rotate = "oblimin", nfactors = 1)

fit_multi <- cfa(cfa_model,
                 data = d,
                 group = "domain")

summary(fit_multi)


################################################################################

library(tidyverse)
library(lavaan)
library(patchwork)

# ------------------------------------------------------------------
# 1. Reduced-domain CFA data
# ------------------------------------------------------------------

d_reduced <- d
  filter(!(domain %in% c("geography", "number", "maybe")))

d_rich <- d_reduced %>%
  select(really, consider, predict, context, surprise, bet, mot, jus, double)

cfa_model <- '
  general =~ really + consider + predict + context + surprise + bet + mot + jus + double
'

fit_cfa <- cfa(cfa_model, data = d_rich)

# factor scores
fs <- lavPredict(fit_cfa)
fs <- lavPredict(fit_cfa, newdata = d)

d_plot <- d %>%
  mutate(general_factor = as.numeric(fs[, 1]),
         really_believe_z      = as.numeric(scale(general_factor))) %>%
  mutate(
    factor_group = case_when(
      really_believe_z <= -0.5 ~ "Low (-0.5 SD or lower)",
      really_believe_z >=  0.5 ~ "High (+0.5 SD or higher)",
      TRUE            ~ "Intermediate"
    )
  ) 
  filter(domain %in% c("geography", "number"))

d_plot$factor_group <- factor(
  d_plot$factor_group,
  levels = c("Low (-0.5 SD or lower)", "Intermediate", "High (+0.5 SD or higher)")
)

# ------------------------------------------------------------------
# 2. Long data for new items (sans really)
# ------------------------------------------------------------------

plot_df <- d_plot %>%
  select(prolific, general_factor, really_believe_z, factor_group,
         dishonest, serious, sense, single, deeper) %>%
  pivot_longer(
    cols = c(dishonest, serious, sense, single, deeper),
    names_to = "item",
    values_to = "value"
  )

plot_df$item <- factor(
  plot_df$item,
  levels = c("dishonest", "sense", "single", "serious", "deeper")
)

# ------------------------------------------------------------------
# 3. Factor-score distribution with shaded regions
# ------------------------------------------------------------------

p_factor <- ggplot(d_plot, aes(x = really_believe_z)) +
  annotate("rect", xmin = -Inf, xmax = -0.5, ymin = -Inf, ymax = Inf,
           fill = "grey90", alpha = 0.9) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf,
           fill = "grey30", alpha = 0.35) +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "grey2", alpha = 0.75) +
  geom_density(linewidth = 1, color = "black", fill = NA) +
  geom_vline(xintercept = c(-0.5, 0.5),
             linetype = c("dashed", "dashed"),
             linewidth = c(0.7, 0.7),
             color = "black") +
  theme_classic() +
  labs(
    title = "",
    y = "Density",
    x = "Really believe factor"
  )

# ------------------------------------------------------------------
# 4. Density of new items, highlighted by factor-score group
# ------------------------------------------------------------------

p_items <- ggplot(plot_df, aes(x = value, fill = factor_group, color = factor_group)) +
  geom_density(alpha = 0.35, linewidth = 0.8, adjust = 1.1) +
  geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.8, color = "black") +
  facet_wrap(~ item, scales = "free_y") +
  scale_fill_manual(values = c(
    "Low (-0.5 SD or lower)"   = "grey80",
    "Intermediate"           = "grey50",
    "High (+0.5 SD or higher)" = "black"
  )) +
  scale_color_manual(values = c(
    "Low (-0.5 SD or lower)"   = "grey65",
    "Intermediate"           = "grey40",
    "High (+0.5 SD or higher)" = "black"
  )) +
  theme_classic() +
  labs(
    x = "Item response",
    y = "Density"
  )

# ------------------------------------------------------------------
# 5. Combined figure
# ------------------------------------------------------------------

p_factor / p_items + plot_layout(heights = c(1, 2))

plot_df2 <- plot_df %>%
  mutate(
    really_believe_z_j = jitter(really_believe_z, amount = 0.02),
    value_j     = jitter(value, amount = 0.02)
  )

ggplot(plot_df2, aes(x = really_believe_z_j, y = value_j)) +
  stat_density_2d_filled(contour_var = "ndensity") +
  scale_fill_grey(start = 0.9, end = 0.2) +
  geom_hline(yintercept = 50, linewidth = 0.8, color = "black", linetype = "dashed") +
  facet_grid(factor_group ~ item, scales = "free_y") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Really believe factor",
    y = "Item response"
  )

ggplot(plot_df2, aes(x = really_believe_z_j, y = value_j)) +
  geom_point(alpha = 0.12, size = 1.5, color = "black") +
  geom_hline(yintercept = 50, linewidth = 0.8, color = "black", linetype = "dashed") +
  facet_grid(factor_group ~ item, scales = "free_y") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Really believe factor",
    y = "Item response"
  )

################################################################################

# prepare for clustering
d.c <- d %>% 
  dplyr::select(jus, mot, consider, predict, context, surprise, bet, double, really, domain)

# ,dishonest, serious, sense, single, really, deeper, domain

vars <- c("really", "consider", "predict", "context",
          "surprise", "bet", "mot", "jus", "double", "really")

d.complete <- d.c[complete.cases(d.c[, c(vars, "domain")]), ]
residuals_df <- data.frame(matrix(ncol = length(vars), nrow = nrow(d.complete)))
colnames(residuals_df) <- vars
for(var in vars) {
  mod <- lmer(paste(var, "~ 1 + (1|domain)"), data = d.complete)
  residuals_df[[var]] <- residuals(mod)
}

#### GMM ####
# configure
# d.cluster <- residuals_df
d.cluster <- d.c %>% 
  filter((domain %in% c("religious", "political"))) %>% 
  dplyr::select(-domain)

# how many clusters?
m <- mclustBIC(d.cluster)
summary(m)
plot(m)

# analysis
mod1 <- Mclust(d.cluster, x = m, G=2, modelNames = "VEV")
summary(mod1, parameters = TRUE)

# plot
plot(mod1, what = "classification")

# Add cluster assignments to data
d_with_clusters_gauss <- d.cluster
d_with_clusters_gauss$cluster <- factor(mod1$classification,
                                        labels=c("Group 1", "Group 2"))

# Extract covariance matrices for each cluster
cov1 <- mod1$parameters$variance$sigma[,,1]
cov2 <- mod1$parameters$variance$sigma[,,2]

# Full sample covariance (pooled across groups)
cov_full <- cov(d.cluster)

# Convert to correlations (removes scale/variance effects)
cor1 <- cov2cor(cov1)
cor2 <- cov2cor(cov2)
cor_full <- cov2cor(cov_full)

library(corrplot)

# Calculate p-values for the full dataset (you already have pmat1 and pmat2)
pmat_full <- cor.mtest(d.cluster)$p

# p-values for each cluster based on assigned cases
pmat1 <- cor.mtest(d_with_clusters_gauss %>% 
                     filter(cluster == "Group 1") %>% 
                     select(-cluster))$p

pmat2 <- cor.mtest(d_with_clusters_gauss %>% 
                     filter(cluster == "Group 2") %>% 
                     select(-cluster))$p

# Set up a 1x3 plot layout
par(mfrow = c(1, 3))

# Plot Group 1
corrplot(cor1, method = "number", type = "upper", 
         tl.col = "black", tl.srt = 45,
         p.mat = pmat1, sig.level = 0.05, insig = "pch", pch = 47, pch.col = "black",
         mar = c(0,0,0,0))

# Plot Group 2
corrplot(cor2, method = "number", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         p.mat = pmat2, sig.level = 0.05, insig = "pch", pch = 47, pch.col = "black",
         mar = c(0,0,0,0))

# Plot Overall Sample
corrplot(cor_full, method = "number", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         p.mat = pmat_full, sig.level = 0.05, insig = "pch", pch = 47, pch.col = "black",
         mar = c(0,0,0,0))

# Reset plot layout
par(mfrow = c(1, 1))

# install.packages(c("ggplot2", "tidyr")) # Run if you don't have these installed
library(ggplot2)
library(tidyr)

# Use the dataframe you already created that contains the cluster labels
# Reshape the data to "long" format so ggplot can map over every variable
d_long <- d_with_clusters_gauss |> 
  pivot_longer(cols = -cluster, names_to = "Variable", values_to = "Value")

library(ggplot2)

# Assuming d_long is still in your environment from the previous step
ggplot(d_long, aes(x = Value, fill = cluster)) +
  # Use position="identity" to overlap, or "dodge" to put bars side-by-side
  geom_histogram(position = "identity", alpha = 0.5, bins = 5, color = "white") +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(title = "Histograms by Variable and Cluster",
       x = "Value",
       y = "Count",
       fill = "Cluster") +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

# install.packages("GGally") # Run this once if you don't have it
library(GGally)

library(GGally)
library(ggplot2)

# 1. Define a custom function for the scatterplots
custom_scatter <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    # Draw the points colored by group
    geom_point(alpha = 0.2) +
    
    # Draw the group-specific trendlines (inherits the color mapping)
    geom_smooth(method = "lm", se = FALSE, ...) +
    
    # Draw the full cohort trendline 
    # aes(group = 1) forces it to ignore the cluster colors and treat all data as one group
    geom_smooth(aes(group = 1), method = "lm", color = "black", linetype = "dashed", se = FALSE)
}

# 2. Feed the custom function into ggpairs
n_vars <- ncol(d.cluster)

ggpairs(d_with_clusters_gauss, 
        columns = 1:n_vars, 
        mapping = aes(color = cluster),
        # Apply our custom function to the continuous lower panels
        lower = list(continuous = custom_scatter),
        title = "") +
  theme_minimal()

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
fa(d.cluster, fm = "pa", rotate = "oblimin", nfactors = 1)


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

