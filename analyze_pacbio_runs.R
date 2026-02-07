# PacBio Sequel IIe Run Analysis
# Analysis of loading metrics (P0, P1, P2) vs yield at SMRT Cell level
# Each row = one independent SMRT Cell

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(broom)

# Set theme for all plots
theme_set(theme_bw(base_size = 12))

# Read data
data <- read_csv("ALL_260207.csv", show_col_types = FALSE)

# Cleanup: remove CLR runs, diagnostic runs, and non-HiFi BAM references
data <- data %>%
  filter(run_type == "CCS/HiFi") %>%
  filter(!str_detect(run_name, regex("Diag", ignore_case = TRUE))) %>%
  filter(is_hifi_bam == TRUE)  # Only include runs with .hifi_reads.bam (not .reads.bam)

# Convert total_length to Gb
data <- data %>%
  mutate(yield_gb = total_length / 1e9)

# Filter to complete cases for main analysis
data_complete <- data %>%
  filter(!is.na(yield_gb) & !is.na(p0_percent) & !is.na(p1_percent) &
         !is.na(p2_percent) & !is.na(insert_size))

# Create output directory
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ============================================================================
# 1. CORRELATION ANALYSIS
# ============================================================================

# Calculate correlations
cor_p0 <- cor(data_complete$p0_percent, data_complete$yield_gb, use = "complete.obs")
cor_p1 <- cor(data_complete$p1_percent, data_complete$yield_gb, use = "complete.obs")
cor_p2 <- cor(data_complete$p2_percent, data_complete$yield_gb, use = "complete.obs")

# Test significance
cor_test_p0 <- cor.test(data_complete$p0_percent, data_complete$yield_gb)
cor_test_p1 <- cor.test(data_complete$p1_percent, data_complete$yield_gb)
cor_test_p2 <- cor.test(data_complete$p2_percent, data_complete$yield_gb)

# ============================================================================
# 2. DISTRIBUTIONS: Histograms of Loading Metrics
# ============================================================================

# P1 distribution
p1_hist <- gghistogram(data_complete, x = "p1_percent",
                       bins = 50,
                       fill = "#E41A1C",
                       add = "mean") +
  geom_vline(xintercept = c(30, 40), linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  annotate("rect", xmin = 30, xmax = 40, ymin = 0, ymax = Inf,
           alpha = 0.1, fill = "green") +
  labs(
    title = "P1 (Single Loading) Distribution",
    subtitle = "Green shaded area = optimal range (30-40%)",
    x = "P1 - Single Loading (%)",
    y = "Count"
  )

ggsave("figures/00a_p1_distribution.png", p1_hist, width = 8, height = 6, dpi = 300)

# P0 distribution
p0_hist <- gghistogram(data_complete, x = "p0_percent",
                       bins = 50,
                       fill = "#E41A1C",
                       add = "mean") +
  labs(
    title = "P0 (Empty ZMWs) Distribution",
    x = "P0 - Empty ZMWs (%)",
    y = "Count"
  )

ggsave("figures/00b_p0_distribution.png", p0_hist, width = 8, height = 6, dpi = 300)

# P2 distribution
p2_hist <- gghistogram(data_complete, x = "p2_percent",
                       bins = 50,
                       fill = "#E41A1C",
                       add = "mean") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red", linewidth = 0.8) +
  annotate("text", x = 10, y = Inf, label = "Target: <10%",
           vjust = 1.5, hjust = -0.1, color = "red", size = 3.5) +
  labs(
    title = "P2 (Multi Loading) Distribution",
    subtitle = "Lower is better (target <10%)",
    x = "P2 - Multi Loading (%)",
    y = "Count"
  )

ggsave("figures/00c_p2_distribution.png", p2_hist, width = 8, height = 6, dpi = 300)

# Yield distribution
yield_hist <- gghistogram(data_complete, x = "yield_gb",
                          bins = 50,
                          fill = "#E41A1C",
                          add = "mean") +
  labs(
    title = "Yield Distribution",
    x = "Yield (Gb)",
    y = "Count"
  )

ggsave("figures/00d_yield_distribution.png", yield_hist, width = 8, height = 6, dpi = 300)

# Combined loading metrics panel
loading_combined <- ggarrange(p1_hist, p0_hist, p2_hist,
                              ncol = 3, nrow = 1)

ggsave("figures/00e_loading_distributions_combined.png", loading_combined,
       width = 15, height = 6, dpi = 300)

# P1 summary statistics
p1_stats <- data_complete %>%
  summarise(
    mean_p1 = mean(p1_percent, na.rm = TRUE),
    median_p1 = median(p1_percent, na.rm = TRUE),
    sd_p1 = sd(p1_percent, na.rm = TRUE),
    min_p1 = min(p1_percent, na.rm = TRUE),
    max_p1 = max(p1_percent, na.rm = TRUE),
    n_optimal = sum(p1_percent >= 30 & p1_percent <= 40, na.rm = TRUE),
    pct_optimal = 100 * n_optimal / n()
  )

# ============================================================================
# 3. VISUALIZATION: P0, P1, P2 vs Yield
# ============================================================================

# P1 vs Yield
p1_plot <- ggplot(data_complete, aes(x = p1_percent, y = yield_gb)) +
  geom_point(alpha = 0.6, size = 2, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P1 (Single Loading) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p1, cor_test_p1$p.value),
    x = "P1 - Single Loading (%)",
    y = "Yield (Gb)"
  )

ggsave("figures/01_p1_vs_yield.png", p1_plot, width = 8, height = 6, dpi = 300)

# P0 vs Yield
p0_plot <- ggplot(data_complete, aes(x = p0_percent, y = yield_gb)) +
  geom_point(alpha = 0.6, size = 2, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P0 (Empty ZMWs) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p0, cor_test_p0$p.value),
    x = "P0 - Empty ZMWs (%)",
    y = "Yield (Gb)"
  )

ggsave("figures/02_p0_vs_yield.png", p0_plot, width = 8, height = 6, dpi = 300)

# P2 vs Yield
p2_plot <- ggplot(data_complete, aes(x = p2_percent, y = yield_gb)) +
  geom_point(alpha = 0.6, size = 2, color = "#E41A1C") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linewidth = 1) +
  labs(
    title = "P2 (Multi Loading) vs Yield",
    subtitle = sprintf("r = %.3f, p = %.2e", cor_p2, cor_test_p2$p.value),
    x = "P2 - Multi Loading (%)",
    y = "Yield (Gb)"
  )

ggsave("figures/03_p2_vs_yield.png", p2_plot, width = 8, height = 6, dpi = 300)

# Combined panel
combined_plot <- ggarrange(p1_plot, p0_plot, p2_plot,
                           ncol = 2, nrow = 2)

ggsave("figures/04_combined_loading_vs_yield.png", combined_plot,
       width = 12, height = 10, dpi = 300)

# ============================================================================
# 4. LINEAR MODELS: Effect of Loading Metrics on Yield
# ============================================================================

# Model 1: P1 only
model_p1 <- lm(yield_gb ~ p1_percent, data = data_complete)

# Model 2: All loading metrics
model_all_loading <- lm(yield_gb ~ p0_percent + p1_percent + p2_percent,
                        data = data_complete)

# Model 3: Including insert size
model_with_insert <- lm(yield_gb ~ p0_percent + p1_percent + p2_percent + insert_size,
                        data = data_complete)

# Model 4: With interaction between P1 and insert size
model_interaction <- lm(yield_gb ~ p1_percent * insert_size + p0_percent + p2_percent,
                        data = data_complete)

# Compare models with ANOVA
anova_result <- anova(model_p1, model_all_loading, model_with_insert, model_interaction)

# ============================================================================
# 5. VISUALIZATION: Insert Size Effect
# ============================================================================

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

# ============================================================================
# 6. RESIDUAL DIAGNOSTICS
# ============================================================================

# Diagnostic plots for best model
png("figures/07_model_diagnostics.png", width = 10, height = 8, units = "in", res = 300)
par(mfrow = c(2, 2))
plot(model_with_insert, main = "Model 3: With Insert Size")
dev.off()

# ============================================================================
# 7. SUMMARY TABLE
# ============================================================================

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

write_csv(model_comparison, "results/model_comparison.csv")

# Summary statistics
summary_stats <- data_complete %>%
  summarise(
    n_cells = n(),
    mean_yield = mean(yield_gb, na.rm = TRUE),
    sd_yield = sd(yield_gb, na.rm = TRUE),
    mean_p0 = mean(p0_percent, na.rm = TRUE),
    mean_p1 = mean(p1_percent, na.rm = TRUE),
    mean_p2 = mean(p2_percent, na.rm = TRUE),
    mean_insert = mean(insert_size, na.rm = TRUE)
  )

write_csv(summary_stats, "results/summary_stats.csv")

# ============================================================================
# 8. KEY FINDINGS (stored in variables for interactive inspection)
# ============================================================================

best_model <- if (summary(model_with_insert)$adj.r.squared >
                   summary(model_interaction)$adj.r.squared) {
  model_with_insert
} else {
  model_interaction
}

