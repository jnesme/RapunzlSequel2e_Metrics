#!/usr/bin/env Rscript
# PacBio Sequel IIe Run Analysis
# Analysis of loading metrics (P0, P1, P2) vs yield with insert size effects

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(broom)

# Set theme for all plots
theme_set(theme_bw(base_size = 12))

# Read data
cat("Reading data...\n")
data <- read_csv("ALL_260205.csv", show_col_types = FALSE)

cat(sprintf("Loaded %d runs\n", nrow(data)))
cat(sprintf("  CCS/HiFi: %d\n", sum(data$run_type == "CCS/HiFi", na.rm = TRUE)))
cat(sprintf("  CLR: %d\n", sum(data$run_type == "CLR", na.rm = TRUE)))
cat("\n")

# Filter to complete cases for main analysis
data_complete <- data %>%
  filter(!is.na(yield_gb) & !is.na(p0_percent) & !is.na(p1_percent) & 
         !is.na(p2_percent) & !is.na(insert_size))

cat(sprintf("Complete cases for analysis: %d runs\n\n", nrow(data_complete)))

# Create output directory
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ============================================================================
# 1. CORRELATION ANALYSIS
# ============================================================================
cat("="*70, "\n")
cat("CORRELATION ANALYSIS: Loading Metrics vs Yield\n")
cat("="*70, "\n\n")

# Calculate correlations
cor_p0 <- cor(data_complete$p0_percent, data_complete$yield_gb, use = "complete.obs")
cor_p1 <- cor(data_complete$p1_percent, data_complete$yield_gb, use = "complete.obs")
cor_p2 <- cor(data_complete$p2_percent, data_complete$yield_gb, use = "complete.obs")

cat(sprintf("Pearson correlations with yield:\n"))
cat(sprintf("  P0 (Empty):   r = %.3f\n", cor_p0))
cat(sprintf("  P1 (Single):  r = %.3f\n", cor_p1))
cat(sprintf("  P2 (Multi):   r = %.3f\n", cor_p2))
cat("\n")

# Test significance
cor_test_p0 <- cor.test(data_complete$p0_percent, data_complete$yield_gb)
cor_test_p1 <- cor.test(data_complete$p1_percent, data_complete$yield_gb)
cor_test_p2 <- cor.test(data_complete$p2_percent, data_complete$yield_gb)

cat("Correlation tests (p-values):\n")
cat(sprintf("  P0: p = %.2e\n", cor_test_p0$p.value))
cat(sprintf("  P1: p = %.2e\n", cor_test_p1$p.value))
cat(sprintf("  P2: p = %.2e\n", cor_test_p2$p.value))
cat("\n")

# ============================================================================
# 2. VISUALIZATION: P0, P1, P2 vs Yield
# ============================================================================
cat("Creating scatter plots...\n")

# P1 vs Yield
p1_plot <- ggplot(data_complete, aes(x = p1_percent, y = yield_gb)) +
  geom_point(aes(color = run_type), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P1 (Single Loading) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p1, cor_test_p1$p.value),
    x = "P1 - Single Loading (%)",
    y = "Yield (Gb)",
    color = "Run Type"
  ) +
  scale_color_manual(values = c("CCS/HiFi" = "#E41A1C", "CLR" = "#377EB8")) +
  theme(legend.position = "bottom")

ggsave("figures/01_p1_vs_yield.png", p1_plot, width = 8, height = 6, dpi = 300)

# P0 vs Yield
p0_plot <- ggplot(data_complete, aes(x = p0_percent, y = yield_gb)) +
  geom_point(aes(color = run_type), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P0 (Empty ZMWs) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p0, cor_test_p0$p.value),
    x = "P0 - Empty ZMWs (%)",
    y = "Yield (Gb)",
    color = "Run Type"
  ) +
  scale_color_manual(values = c("CCS/HiFi" = "#E41A1C", "CLR" = "#377EB8")) +
  theme(legend.position = "bottom")

ggsave("figures/02_p0_vs_yield.png", p0_plot, width = 8, height = 6, dpi = 300)

# P2 vs Yield
p2_plot <- ggplot(data_complete, aes(x = p2_percent, y = yield_gb)) +
  geom_point(aes(color = run_type), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P2 (Multi Loading) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p2, cor_test_p2$p.value),
    x = "P2 - Multi Loading (%)",
    y = "Yield (Gb)",
    color = "Run Type"
  ) +
  scale_color_manual(values = c("CCS/HiFi" = "#E41A1C", "CLR" = "#377EB8")) +
  theme(legend.position = "bottom")

ggsave("figures/03_p2_vs_yield.png", p2_plot, width = 8, height = 6, dpi = 300)

# Combined panel
combined_plot <- ggarrange(p1_plot, p0_plot, p2_plot, 
                           ncol = 2, nrow = 2,
                           common.legend = TRUE, 
                           legend = "bottom")

ggsave("figures/04_combined_loading_vs_yield.png", combined_plot, 
       width = 12, height = 10, dpi = 300)

cat("  ✓ Saved scatter plots to figures/\n\n")

# ============================================================================
# 3. LINEAR MODELS: Effect of Loading Metrics on Yield
# ============================================================================
cat("="*70, "\n")
cat("LINEAR REGRESSION MODELS\n")
cat("="*70, "\n\n")

# Model 1: P1 only
model_p1 <- lm(yield_gb ~ p1_percent, data = data_complete)

cat("Model 1: Yield ~ P1\n")
cat(sprintf("  R² = %.4f\n", summary(model_p1)$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary(model_p1)$adj.r.squared))
cat(sprintf("  P1 coefficient: %.3f (p = %.2e)\n", 
            coef(model_p1)[2], 
            summary(model_p1)$coefficients[2, 4]))
cat("\n")

# Model 2: All loading metrics
model_all_loading <- lm(yield_gb ~ p0_percent + p1_percent + p2_percent, 
                        data = data_complete)

cat("Model 2: Yield ~ P0 + P1 + P2\n")
cat(sprintf("  R² = %.4f\n", summary(model_all_loading)$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary(model_all_loading)$adj.r.squared))
print(summary(model_all_loading)$coefficients)
cat("\n")

# Model 3: Including insert size
model_with_insert <- lm(yield_gb ~ p0_percent + p1_percent + p2_percent + insert_size, 
                        data = data_complete)

cat("Model 3: Yield ~ P0 + P1 + P2 + Insert_Size\n")
cat(sprintf("  R² = %.4f\n", summary(model_with_insert)$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary(model_with_insert)$adj.r.squared))
print(summary(model_with_insert)$coefficients)
cat("\n")

# Model 4: With interaction between P1 and insert size
model_interaction <- lm(yield_gb ~ p1_percent * insert_size + p0_percent + p2_percent, 
                        data = data_complete)

cat("Model 4: Yield ~ P1 * Insert_Size + P0 + P2\n")
cat(sprintf("  R² = %.4f\n", summary(model_interaction)$r.squared))
cat(sprintf("  Adjusted R² = %.4f\n", summary(model_interaction)$adj.r.squared))
print(summary(model_interaction)$coefficients)
cat("\n")

# Compare models with ANOVA
cat("Model Comparison (ANOVA):\n")
anova_result <- anova(model_p1, model_all_loading, model_with_insert, model_interaction)
print(anova_result)
cat("\n")

# ============================================================================
# 4. VISUALIZATION: Insert Size Effect
# ============================================================================
cat("Creating insert size effect plots...\n")

# Bin insert sizes for visualization
data_complete <- data_complete %>%
  mutate(insert_bin = cut(insert_size, 
                          breaks = c(0, 3000, 6000, 10000, Inf),
                          labels = c("<3kb", "3-6kb", "6-10kb", ">10kb")))

# P1 vs Yield colored by insert size
p1_insert_plot <- ggplot(data_complete, aes(x = p1_percent, y = yield_gb)) +
  geom_point(aes(color = insert_bin), alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
  geom_smooth(aes(color = insert_bin), method = "lm", se = FALSE, linewidth = 0.8) +
  labs(
    title = "P1 vs Yield: Effect of Insert Size",
    x = "P1 - Single Loading (%)",
    y = "Yield (Gb)",
    color = "Insert Size"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")

ggsave("figures/05_p1_vs_yield_by_insert_size.png", p1_insert_plot, 
       width = 9, height = 6, dpi = 300)

# Box plot: Yield by insert size bins
yield_by_insert <- ggplot(data_complete, aes(x = insert_bin, y = yield_gb)) +
  geom_boxplot(aes(fill = insert_bin), alpha = 0.7, outlier.shape = 1) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1.5) +
  stat_compare_means(method = "anova", label.y = max(data_complete$yield_gb, na.rm = TRUE) * 1.05) +
  labs(
    title = "Yield Distribution by Insert Size",
    x = "Insert Size",
    y = "Yield (Gb)"
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")

ggsave("figures/06_yield_by_insert_size.png", yield_by_insert, 
       width = 8, height = 6, dpi = 300)

cat("  ✓ Saved insert size plots to figures/\n\n")

# ============================================================================
# 5. RESIDUAL DIAGNOSTICS
# ============================================================================
cat("Creating diagnostic plots...\n")

# Diagnostic plots for best model
png("figures/07_model_diagnostics.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))
plot(model_with_insert, main = "Model 3: With Insert Size")
dev.off()

cat("  ✓ Saved diagnostic plots\n\n")

# ============================================================================
# 6. SUMMARY TABLE
# ============================================================================
cat("Creating summary table...\n")

# Model comparison table
model_comparison <- tibble(
  Model = c("P1 only", "All loading (P0+P1+P2)", 
            "With insert size", "P1*Insert interaction"),
  R_squared = c(summary(model_p1)$r.squared,
                summary(model_all_loading)$r.squared,
                summary(model_with_insert)$r.squared,
                summary(model_interaction)$r.squared),
  Adj_R_squared = c(summary(model_p1)$adj.r.squared,
                    summary(model_all_loading)$adj.r.squared,
                    summary(model_with_insert)$adj.r.squared,
                    summary(model_interaction)$adj.r.squared),
  AIC = c(AIC(model_p1), AIC(model_all_loading), 
          AIC(model_with_insert), AIC(model_interaction))
)

print(model_comparison)
write_csv(model_comparison, "results/model_comparison.csv")
cat("\n")

# Summary statistics by run type
summary_by_type <- data_complete %>%
  group_by(run_type) %>%
  summarise(
    n_runs = n(),
    mean_yield = mean(yield_gb, na.rm = TRUE),
    sd_yield = sd(yield_gb, na.rm = TRUE),
    mean_p0 = mean(p0_percent, na.rm = TRUE),
    mean_p1 = mean(p1_percent, na.rm = TRUE),
    mean_p2 = mean(p2_percent, na.rm = TRUE),
    mean_insert = mean(insert_size, na.rm = TRUE),
    .groups = "drop"
  )

cat("Summary by Run Type:\n")
print(summary_by_type)
write_csv(summary_by_type, "results/summary_by_run_type.csv")
cat("\n")

# ============================================================================
# 7. KEY FINDINGS
# ============================================================================
cat("="*70, "\n")
cat("KEY FINDINGS\n")
cat("="*70, "\n\n")

best_model <- if (summary(model_with_insert)$adj.r.squared > 
                   summary(model_interaction)$adj.r.squared) {
  model_with_insert
} else {
  model_interaction
}

cat(sprintf("1. P1 (single loading) shows %s correlation with yield (r=%.3f)\n",
            ifelse(cor_p1 > 0, "positive", "negative"), cor_p1))
cat(sprintf("2. P0 (empty ZMWs) shows %s correlation with yield (r=%.3f)\n",
            ifelse(cor_p0 > 0, "positive", "negative"), cor_p0))
cat(sprintf("3. Insert size %s a significant predictor of yield\n",
            ifelse(summary(model_with_insert)$coefficients["insert_size", "Pr(>|t|)"] < 0.05,
                   "IS", "is NOT")))
cat(sprintf("4. Best model explains %.1f%% of variance in yield\n",
            summary(best_model)$adj.r.squared * 100))
cat(sprintf("5. Average P1 is %.1f%% (target: 30-40%%)\n",
            mean(data_complete$p1_percent, na.rm = TRUE)))
cat(sprintf("6. Average yield is %.1f Gb per run\n",
            mean(data_complete$yield_gb, na.rm = TRUE)))
cat("\n")

cat("Analysis complete! Results saved to:\n")
cat("  - figures/     : All plots\n")
cat("  - results/     : CSV summaries\n")
cat("\n")
