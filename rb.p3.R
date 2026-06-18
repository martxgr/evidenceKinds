library(tidyverse)
library(mclust)
library(lme4)
library(psych)
library(corrplot)
library(GGally)

# ==========================================
# 1. CONFIGURATION & TOGGLES
# ==========================================
# Toggle 1: Filter 'really_believe' to +/- 1 SD from the mean?
use_intermediate_filter <- FALSE

# Toggle 2: Cluster on domain residuals instead of raw scores?
use_domain_residuals <- FALSE

# ==========================================
# 2. DATA CLEANING & PREP
# ==========================================
d <- read.csv("rb.p3.csv") %>% 
  filter(attention == 2 | attention == "") %>%
  mutate(domain = case_when(!(religious == "") ~ "religious",
                            !(superlative == "") ~ "superlative",
                            !(superstition == "") ~ "superstition",
                            !(political == "") ~ "political"),
         really_believe    = as.numeric(really_1), 
         confident_think   = as.numeric(y_t_1), 
         unconfident_think = as.numeric(n_t_1),
         uncertain_think   = as.numeric(u_t_1), 
         confident_say     = as.numeric(y_s_1), 
         unconfident_say   = as.numeric(n_s_1),
         uncertain_say     = as.numeric(u_s_1), 
         confident_do      = as.numeric(y_d_1), 
         unconfident_do    = as.numeric(n_d_1),
         uncertain_do      = as.numeric(u_d_1)) %>%
  filter(!is.na(domain))

vars <- c("confident_think", "unconfident_think", "uncertain_think",
          "confident_say", "unconfident_say", "uncertain_say",
          "confident_do", "unconfident_do", "uncertain_do")

# Create ONE unified dataset without NAs, adding Row_ID for distinct counts later
d_clean <- d %>% 
  dplyr::select(all_of(vars), really_believe, domain) %>%
  drop_na() %>%
  mutate(Row_ID = row_number())

# Apply the intermediate filter if toggled TRUE
if(use_intermediate_filter) {
  mean_val <- mean(d_clean$really_believe)
  sd_val <- sd(d_clean$really_believe)
  
  d_clean <- d_clean %>%
    filter(really_believe > (50 - sd_val) & 
             really_believe < (50 + sd_val))
}

hist(d_clean$really_believe)

# ==========================================
# 3. GMM CLUSTERING
# ==========================================
# Base dataset for clustering
d.cluster <- d_clean %>% 
  dplyr::select(all_of(vars))

# Apply the domain residualizer if toggled TRUE
if(use_domain_residuals) {
  for(var in vars) {
    mod <- lm(paste(var, "~ domain"), data = d_clean)
    d.cluster[[var]] <- residuals(mod)
  }
}

# how many clusters?
m <- mclustBIC(d.cluster)
summary(m)
plot(m)

# analysis
mod1 <- Mclust(d.cluster, x = m, G = 5, modelNames = "VII")
summary(mod1, parameters = TRUE)
plot(mod1, what = "classification")

# Append class assignments back to both our working datasets
d_clean$Class_Num <- mod1$classification
d.cluster$Class_Num <- mod1$classification

# ==========================================
# 4. PROFILE VISUALIZATION
# ==========================================
plot_data <- as.data.frame(mod1$parameters$mean) %>%
  rownames_to_column(var = "Item") %>%
  pivot_longer(cols = -Item, 
               names_to = "Class_Num", 
               values_to = "Mean") %>%
  mutate(Class = paste("Class", str_extract(Class_Num, "\\d+")),
         Class_Idx = as.integer(str_extract(Class_Num, "\\d+")))

class_counts <- as.numeric(table(mod1$classification))

# EII models store variance as a scalar, VII uses a vector, others use arrays. 
# Using length.out safely handles scalars vs vectors without over-repeating!
if (!is.null(mod1$parameters$variance$sigmasq)) {
  variances <- rep(mod1$parameters$variance$sigmasq, length.out = mod1$G)
} else {
  variances <- mod1$parameters$variance$sigma[1, 1, ]
}

se_df <- tibble(Class_Idx = 1:mod1$G,
                n = class_counts, 
                SE = sqrt(variances / class_counts))

plot_data <- plot_data %>%
  left_join(se_df, by = "Class_Idx") %>%
  mutate(Lower = Mean - SE,
         Upper = Mean + SE,
         Item = fct_inorder(Item),
         Class = paste0(Class, " (n = ", n, ")"))

# Dynamic X-axis label based on our toggle
x_axis_label <- ifelse(use_domain_residuals, "Deviation from Domain Mean (Residual)", "Slider Score")

ggplot(plot_data, aes(x = Mean, y = fct_rev(Item))) +
  geom_col(fill = "gray70", color = "black", alpha = 0.8) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), width = 0.25) +
  facet_wrap(~ Class, nrow = 1) + 
  theme_classic() +
  labs(title = "Average Belief Profiles by Class", x = x_axis_label, y = NULL) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90", color = "black"),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.y = element_text(size = 10, color = "black"))

# ==========================================
# 5. OUTCOME ANALYSIS: REALLY BELIEVE
# ==========================================
# 1. Establish statistical boundaries for intermediacy
rb_sd <- sd(d_clean$really_believe, na.rm = TRUE)
lower_bound <- 50 - rb_sd
upper_bound <- 50 + rb_sd

# 2. Calculate class-level summary statistics
rb_summary <- d_clean %>%
  group_by(Class_Num) %>%
  summarize(
    n = n(),
    mean_rb = mean(really_believe, na.rm = TRUE),
    sd_rb = sd(really_believe, na.rm = TRUE),
    se_rb = sd_rb / sqrt(n),
    .groups = "drop" 
  ) %>%
  mutate(
    Class_Label = factor(paste0("Class ", Class_Num, "\n(n = ", n, ")")),
    Lower = mean_rb - se_rb,
    Upper = mean_rb + se_rb
  )

# 3. Create master plotting dataframe combining raw data, metrics, and labels
d_plot <- d_clean %>%
  mutate(
    Class_Factor = as.factor(paste("Class", Class_Num)),
    intermediate = case_when(
      really_believe >= lower_bound & really_believe <= upper_bound ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  left_join(rb_summary, by = "Class_Num")

# 4. Pre-calculate text labels for Figure 5A
class_proportions <- d_plot %>%
  group_by(Class_Label) %>%
  summarize(
    prop_inter = mean(intermediate),
    label = paste0("Intermediate: ", round(prop_inter * 100, 1), "%"),
    .groups = "drop"
  )
class1_domain_summary <- d_plot %>%
  filter(Class_Num == 1) %>%
  group_by(domain) %>%
  summarize(Count = n(), .groups = "drop") %>%
  mutate(Percentage = round((Count / sum(Count)) * 100, 1)) %>%
  arrange(desc(Count))

print("--- Domain Breakdown for Class 1 ---")
print(class1_domain_summary)


# ------------------------------------------
# STEP 3: VISUALIZATION SUITE
# ------------------------------------------

# --- FIGURE 5A: Faceted Distribution with Proportions ---
ggplot(d_plot, aes(x = really_believe, fill = intermediate)) +
  geom_histogram(color = "black", bins = 15, alpha = 0.8) +
  geom_text(data = class_proportions, aes(x = 50, y = Inf, label = label), 
            vjust = 1.8, fontface = "bold", color = "black", size = 3.5, inherit.aes = FALSE) +
  facet_wrap(~ Class_Label) +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "gray20")) +
  theme_classic() +
  labs(
    title = "Distribution of 'Really Believe' Scores by Class",
    x = "Really Believe Score",
    y = "Count",
    fill = "Intermediate Believer"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 10)
  )

# --- FIGURE 5B: Continuous Spectrum & Bounds ---
ggplot(d_plot, aes(x = really_believe, fill = Class_Factor)) +
  geom_histogram(position = "stack", color = "black", bins = 20, alpha = 0.8) +
  geom_vline(xintercept = c(lower_bound, upper_bound), 
             color = "black", linetype = "dashed", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "Class Breakdown Across the 'Really Believe' Spectrum",
    subtitle = "Dashed lines indicate the bounds of Intermediacy (+/- 1 SD from 50)",
    x = "Really Believe Score",
    y = "Count",
    fill = "Assigned Class"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12))

# --- FIGURE 5C: Categorical Class Proportions ---
ggplot(d_plot, aes(x = intermediate, fill = Class_Factor)) +
  geom_bar(position = "fill", color = "black", alpha = 0.8, width = 0.5) +
  scale_y_continuous(labels = scales::percent) + 
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
  labs(
    title = "Class Proportions among Intermediate vs. Extreme Believers",
    x = "Is the Participant an Intermediate Believer?",
    y = "Proportion",
    fill = "Assigned Class"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold")
  )

# --- FIGURE 5D: Class 1 Specific Domain Distribution ---
ggplot(filter(d_plot, Class_Num == 1), aes(x = fct_infreq(domain), fill = domain)) +
  geom_bar(color = "black", alpha = 0.8, width = 0.6) +
  scale_fill_brewer(palette = "Pastel1") +
  theme_classic() +
  labs(
    title = "Distribution of Content Domains within Class 1",
    subtitle = "Showing which domains make up the Class 1 profile",
    x = "Domain",
    y = "Number of Observations"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 10),
    legend.position = "none"
  )

# --- FIGURE 5E: Domain Distributions Faceted by Class ---
ggplot(d_plot, aes(x = really_believe, fill = domain)) +
  geom_histogram(position = "stack", color = "black", bins = 20, alpha = 0.8) +
  geom_vline(xintercept = c(lower_bound, upper_bound), 
             color = "black", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ Class_Factor) +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    title = "Domain Distributions Across 'Really Believe' Spectrum by Latent Class",
    x = "Really Believe Score",
    y = "Count",
    fill = "Domain"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11)
  )

# --- FIGURE 5F: Item-by-Item Raw Distributions (Grand Finale with Density Shadows) ---
# 1. Define your explicit item order
desired_order <- c(
  "confident_think", "unconfident_think", "uncertain_think",
  "confident_say",   "unconfident_say",   "uncertain_say",
  "confident_do",    "unconfident_do",    "uncertain_do"
)

# 2. Extract the exact model profile parameters
model_means <- as.data.frame(mod1$parameters$mean) %>%
  rownames_to_column(var = "Item") %>%
  pivot_longer(cols = -Item, names_to = "Class_Str", values_to = "Mean") %>%
  mutate(Class_Num = as.integer(str_extract(Class_Str, "\\d+")))

class_counts <- as.numeric(table(mod1$classification))
if (!is.null(mod1$parameters$variance$sigmasq)) {
  variances <- rep(mod1$parameters$variance$sigmasq, length.out = mod1$G)
} else {
  variances <- mod1$parameters$variance$sigma[1, 1, ]
}

se_df <- tibble(
  Class_Num = 1:mod1$G,
  n = class_counts, 
  SE = sqrt(variances / class_counts)
)

model_profile <- model_means %>%
  left_join(se_df, by = "Class_Num") %>%
  mutate(
    Lower = Mean - SE,
    Upper = Mean + SE
  ) %>%
  select(Item, Class_Num, Mean, Lower, Upper)

# 3. Pivot the raw item variables long
item_plot_data <- d_plot %>%
  select(all_of(vars), Class_Label, Class_Num, domain, intermediate) %>%
  pivot_longer(cols = all_of(vars), names_to = "Item", values_to = "Score") %>%
  mutate(Item = factor(Item, levels = desired_order)) %>%
  left_join(model_profile, by = c("Item", "Class_Num"))

# 4. Handle the dynamic X-axis label
x_axis_label <- ifelse(exists("use_domain_residuals") && use_domain_residuals, 
                       "Deviation from Domain Mean (Residual)", 
                       "Slider Score")

# 5. Build the item-by-item distribution plot
ggplot(item_plot_data, aes(x = Score, y = fct_rev(Item))) +
  geom_violin(orientation = "y", fill = "gray85", color = NA, alpha = 0.35, scale = "width", width = 0.6) +
  geom_jitter(aes(color = intermediate, shape = domain), 
              width = 0,       
              height = 0.18,  
              alpha = 0.4,     
              size = 1.8) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), height = 0.25, color = "black", linewidth = 0.8, width= 0.1) +
  geom_point(aes(x = Mean), color = "black", size = 4, shape = 18) + 
  facet_wrap(~ Class_Label, nrow = 1) + 
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  labs(
    title = "",
    x = x_axis_label,
    y = NULL,
    color = "Intermediacy",
    shape = "Domain"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ==========================================
# 6. CORRELATION MATRICES BY CLASS
# ==========================================
num_cols <- ceiling(mod1$G / 1)
par(mfrow = c(1, num_cols)) 

for (i in 1:mod1$G) {
  class_data <- d.cluster %>%
    filter(Class_Num == i) %>%
    dplyr::select(all_of(vars))
  
  cormat <- cor(class_data, use = "pairwise.complete.obs")
  
  corrplot(cormat, 
           method = "number", 
           type = "upper", 
           tl.col = "black", 
           tl.cex = 0.8,
           title = paste("Class", i, "Correlations"),
           mar = c(0, 0, 2, 0))
}

par(mfrow = c(1, 1))

# ==========================================
# 7. OVERLAID REGRESSION MATRICES
# ==========================================
pair_data <- d.cluster %>%
  dplyr::select(all_of(vars), Class_Num) %>%
  mutate(Class = as.factor(paste("Class", Class_Num))) %>%
  dplyr::select(-Class_Num)

plot_overlaid_lm <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha = 0.15, size = 0.5) +
    geom_smooth(method = "lm", se = FALSE, size = 1, ...) +
    theme_classic()
}

ggpairs(
  pair_data, 
  columns = 1:length(vars),
  mapping = aes(color = Class, fill = Class),
  lower = list(continuous = plot_overlaid_lm),
  diag = list(continuous = wrap("densityDiag", alpha = 0.5)), 
  upper = list(continuous = wrap("cor", size = 3)) 
) +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    axis.text = element_blank(), 
    axis.ticks = element_blank()
  )

# ==========================================
# 8. MODERATION ANALYSIS (FULL SAMPLE)
# ==========================================
mod_data <- d.cluster %>%
  mutate(Class_Factor = as.factor(Class_Num))

var_pairs <- combn(vars, 2, simplify = FALSE)
moderation_results <- list()

for(i in seq_along(var_pairs)) {
  var_y <- var_pairs[[i]][1]
  var_x <- var_pairs[[i]][2]
  
  form <- as.formula(paste(var_y, "~", var_x, "* Class_Factor"))
  fit <- lm(form, data = mod_data)
  
  aov_res <- anova(fit)
  interaction_term <- paste0(var_x, ":Class_Factor")
  
  p_val <- aov_res[interaction_term, "Pr(>F)"]
  f_val <- aov_res[interaction_term, "F value"]
  
  moderation_results[[i]] <- tibble(
    Outcome = var_y,
    Predictor = var_x,
    F_value = round(f_val, 3),
    p_value = p_val
  )
}

mod_df <- bind_rows(moderation_results) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"), 
    Significant = ifelse(p_adj < 0.05, "Yes", "No")
  ) %>%
  arrange(p_adj) 

print(head(mod_df, 10))


# -------------------------------------------------------------------------
# =========================================================================
#                    PART II: INTERMEDIATE BELIEFS DEEP DIVE
# =========================================================================
# -------------------------------------------------------------------------

target_classes <- c(1, 3, 4)

# ==========================================
# 9. ITEM DISTRIBUTIONS (INTERMEDIATE ONLY)
# ==========================================
inter_raw <- d_plot %>% filter(intermediate == TRUE, Class_Num %in% target_classes)

inter_plot_data <- inter_raw %>%
  select(all_of(vars), Class_Num, domain, Row_ID) %>%
  pivot_longer(cols = all_of(vars), names_to = "Item", values_to = "Score") %>%
  mutate(Item = factor(Item, levels = desired_order)) %>%
  group_by(Class_Num) %>%
  mutate(Class_Label_Inter = paste0("Class ", Class_Num, "\n(n = ", n_distinct(Row_ID), ")")) %>%
  ungroup()

inter_summary <- inter_plot_data %>%
  group_by(Class_Label_Inter, Item) %>%
  summarize(Mean = mean(Score, na.rm = TRUE), SE = sd(Score, na.rm = TRUE)/sqrt(n()), .groups = "drop")

ggplot(inter_plot_data, aes(x = Score, y = fct_rev(Item))) +
  geom_violin(fill = "gray85", color = NA, alpha = 0.35, scale = "width", width = 0.6) +
  geom_jitter(aes(color = domain, shape = domain), width = 0, height = 0.18, alpha = 0.6, size = 2.2) +
  geom_errorbar(data = inter_summary, aes(xmin = Mean-SE, xmax = Mean+SE, y = fct_rev(Item)), height = 0.25, color = "black", linewidth = 0.8, width = 0, inherit.aes = FALSE) +
  geom_point(data = inter_summary, aes(x = Mean, y = fct_rev(Item)), color = "black", size = 3.5, shape = 18, inherit.aes = FALSE) +
  facet_wrap(~ Class_Label_Inter, nrow = 1) + 
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
  labs(title = "Item-by-Item Intermediate Belief Distributions", x = "Score", y = NULL, color = "Domain", shape = "Domain") +
  theme(plot.title = element_text(face="bold"), strip.background = element_rect(fill="grey90"), strip.text = element_text(face="bold"))

# ==========================================
# 10. CORRELATION MATRICES (INTERMEDIATE)
# ==========================================
cor_long <- bind_rows(lapply(target_classes, function(c) {
  as.data.frame(as.table(cor(inter_raw %>% filter(Class_Num == c) %>% select(all_of(vars)), use = "pairwise"))) %>%
    mutate(Class_Label = paste("Class", c))
})) %>%
  mutate(Var1 = factor(Var1, levels = desired_order), Var2 = factor(Var2, levels = rev(desired_order)))

ggplot(cor_long, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Freq, 2)), color = "black", size = 3) +
  facet_wrap(~ Class_Label, nrow = 1) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() + labs(title = "", x = NULL, y = NULL, fill = "r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face="bold"), panel.grid = element_blank())

# ==========================================
# 11. PARALLEL ANALYSIS
# ==========================================
par(mfrow = c(1, length(target_classes)))
cat("\n--- PARALLEL ANALYSIS RECOMMENDATIONS ---\n")

for (c_num in target_classes) {
  pa <- fa.parallel(inter_raw %>% filter(Class_Num == c_num) %>% select(all_of(desired_order)), fm = "minres", fa = "fa", main = paste("Class", c_num), plot = TRUE)
  cat("Class", c_num, "Recommended Factors:", pa$nfact, "\n")
}
par(mfrow = c(1, 1))

# Inspect extracted models based on parallel analysis hits
fa(inter_raw %>% filter(Class_Num == 1) %>% select(all_of(desired_order)), nfactors = 1, fm = "minres")
fa(inter_raw %>% filter(Class_Num == 3) %>% select(all_of(desired_order)), nfactors = 1, fm = "minres")
fa(inter_raw %>% filter(Class_Num == 4) %>% select(all_of(desired_order)), nfactors = 2, fm = "minres", rotate = "oblimin")

# ==========================================
# 12. OMNIBUS MODERATION: GRID VISUALIZATION
# ==========================================
pairs_df <- tribble(
  ~Comparison,         ~Action, ~X_Var,              ~Y_Var,
  "Conf vs Unconf",    "Think", "confident_think",   "unconfident_think",
  "Conf vs Unconf",    "Say",   "confident_say",     "unconfident_say",
  "Conf vs Unconf",    "Do",    "confident_do",      "unconfident_do",
  "Conf vs Uncert",    "Think", "confident_think",   "uncertain_think",
  "Conf vs Uncert",    "Say",   "confident_say",     "uncertain_say",
  "Conf vs Uncert",    "Do",    "confident_do",      "uncertain_do",
  "Unconf vs Uncert",  "Think", "unconfident_think", "uncertain_think",
  "Unconf vs Uncert",  "Say",   "unconfident_say",   "uncertain_say",
  "Unconf vs Uncert",  "Do",    "unconfident_do",    "uncertain_do"
)

mod_plot_df <- bind_rows(lapply(1:nrow(pairs_df), function(i) {
  r <- pairs_df[i, ]
  p_int <- anova(lm(reformulate(c(r$X_Var, "Class_Factor"), r$Y_Var), data=inter_raw),
                 lm(reformulate(c(r$X_Var, "Class_Factor", paste0(r$X_Var, ":Class_Factor")), r$Y_Var), data=inter_raw))$`Pr(>F)`[2]
  
  inter_raw %>%
    select(Class_Factor, X_Val = all_of(r$X_Var), Y_Val = all_of(r$Y_Var)) %>%
    mutate(Comparison = r$Comparison, Action = r$Action, p_label = sprintf("Interaction p = %.3f", p_int))
})) %>%
  mutate(Action = factor(Action, levels = c("Think", "Say", "Do")), 
         Comparison = factor(Comparison, levels = c("Conf vs Unconf", "Conf vs Uncert", "Unconf vs Uncert")))

ggplot(mod_plot_df, aes(x = X_Val, y = Y_Val)) +
  geom_jitter(aes(color = Class_Factor), alpha = 0.25, width = 1.5, height = 1.5, size = 1) +
  geom_smooth(aes(color = Class_Factor), method = "lm", se = FALSE, linewidth = 1.2) +
  geom_text(data = mod_plot_df %>% distinct(Comparison, Action, p_label), aes(x = -Inf, y = Inf, label = p_label), hjust = -0.1, vjust = 1.5, fontface = "italic", size = 3.5) +
  facet_grid(Comparison ~ Action, scales = "free") +
  scale_color_manual(values = c("Class 1" = "green", "Class 3" = "grey", "Class 4" = "black")) +
  theme_bw() + labs(title = "", x = "", y = "", color = "Class") +
  theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"), legend.position = "bottom")

# ==========================================
# 13. CORRELATION GRIDS (TARGETED PAIRS)
# ==========================================
cor_grid_df <- bind_rows(lapply(target_classes, function(c_num) {
  df_sub <- inter_raw %>% filter(Class_Num == c_num)
  bind_rows(lapply(1:nrow(pairs_df), function(i) {
    tibble(Class_Label = paste0("Class ", c_num, "\n(n = ", nrow(df_sub), ")"),
           Comparison = pairs_df$Comparison[i], Action = pairs_df$Action[i],
           Correlation = cor(df_sub[[pairs_df$X_Var[i]]], df_sub[[pairs_df$Y_Var[i]]], use = "pairwise.complete.obs"))
  }))
})) %>%
  mutate(Action = factor(Action, levels = c("Think", "Say", "Do")),
         Comparison = factor(Comparison, levels = rev(c("Conf vs Unconf", "Conf vs Uncert", "Unconf vs Uncert"))),
         Class_Label = fct_inorder(Class_Label))

ggplot(cor_grid_df, aes(x = Action, y = Comparison, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 4.5, fontface = "bold") +
  facet_wrap(~ Class_Label, nrow = 1) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() + labs(title = "", x = "", y = "", fill = "r") +
  theme(strip.background = element_rect(fill = "grey90"), strip.text = element_text(face = "bold"), 
        panel.grid = element_blank(), legend.position = "none")