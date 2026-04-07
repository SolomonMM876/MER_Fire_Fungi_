# Load libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(car)
library(broom)

# Load data
load('HMSC_MER/results/mcmc_output.RData')  # loads Beta1 (PA model)
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')

# Prepare main beta and probability tables
Beta_estimates_PA <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Severity)

prob_x_greater_0 <- Beta1$pos %>%
  dplyr::select(OTU = Species, prob_pos = Severity)

prob_x_less_0 <- Beta1$neg %>%
  dplyr::select(OTU = Species, prob_neg = Severity)

# Taxonomy
soil_OTU_phy <- myco_tax_soil %>%
  distinct(OTU, phylum)

# Thresholds to test
thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9)

# Containers for results
combined_plot_data <- list()
anova_results_list <- list()

# Loop across thresholds
for (thresh in thresholds) {
  
  # 1. Filter significant positive and negative effects
  sig_pos <- Beta_estimates_PA %>%
    left_join(prob_x_greater_0, by = "OTU") %>%
    filter(prob_pos >= thresh)
  
  sig_neg <- Beta_estimates_PA %>%
    left_join(prob_x_less_0, by = "OTU") %>%
    filter(prob_neg >= thresh)
  
  # 2. Combine and join taxonomy
  sig_all <- bind_rows(sig_pos, sig_neg) %>%
    inner_join(soil_OTU_phy, by = "OTU") %>%
    mutate(threshold = thresh)
  
  # 3. Store for plotting
  combined_plot_data[[as.character(thresh)]] <- sig_all
  
  # 4. Run ANOVA if more than 1 phylum present
  if (n_distinct(sig_all$phylum) > 1 && nrow(sig_all) > 3) {
    lm_model <- lm(Beta ~ phylum, data = sig_all)
    smry <- summary(lm_model)
    anov <- Anova(lm_model, type = "II")
    
    anova_tbl <- tibble(
      Threshold = thresh,
      Term = rownames(anov)[1],
      Slope = coef(smry)[2, 1],
      SE = coef(smry)[2, 2],
      F = anov$`F value`[1],
      DF_num = anov$Df[1],
      DF_denom = anov$Df[2],
      P = coef(smry)[2, 4],
      R2 = smry$r.squared
    )
    
    anova_results_list[[as.character(thresh)]] <- anova_tbl
  }
}

# Combine all data for plotting
combined_plot_df <- bind_rows(combined_plot_data) %>%
  mutate(threshold = factor(threshold, levels = thresholds))

# Combine ANOVA results
anova_results_df <- bind_rows(anova_results_list)

# Show ANOVA results
print(anova_results_df)

# Plot
ggplot(combined_plot_df, aes(x = phylum, y = Beta, group=interaction(phylum,threshold))) +
  geom_boxplot(
    aes(alpha = as.numeric(as.character(threshold)), fill = phylum),
    outlier.shape = NA,
    width = 0.6,
    alpha=.3
  ) +
  geom_jitter(
    aes(alpha = as.numeric(as.character(threshold)), color = phylum),
    size = 1.5,
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.6)
  ) +
  scale_alpha_continuous(range = c(0.1, 1),  guide="none") +
  scale_fill_brewer(palette = "Set2", guide="none") +
  scale_color_brewer(palette = "Dark2", guide="none") +
  labs(
    title = "Effect of Fire on Fungal OTUs by Phylum and Certainty Threshold",
    subtitle = "Higher certainty thresholds (0.9) shown with higher opacity",
    x = "Fungal Phylum",
    y = "Fire Effect Size (Beta coefficient)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )




# Plot
ggplot(combined_plot_df, aes(x = phylum, y = Beta, group=interaction(phylum,threshold))) +
  geom_boxplot(
    aes(alpha = as.numeric(as.character(threshold)), fill = phylum),
    outlier.shape = NA,
    width = 0.6
  ) +
  geom_jitter(
    aes(alpha = as.numeric(as.character(threshold)), color = phylum),
    size = 1.5,
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.6)
  ) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  scale_fill_brewer(palette = "Set2", guide="none") +
  scale_color_brewer(palette = "Dark2", guide="none") +
  labs(
    title = "Effect of Fire on Fungal OTUs by Phylum and Certainty Threshold",
    subtitle = "Higher certainty thresholds (0.9) shown with higher opacity",
    x = "Fungal Phylum",
    y = "Fire Effect Size (Beta coefficient)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )




# Summarize min/max for bars
range_summary <- combined_plot_df %>%
  group_by(threshold, phylum) %>%
  summarise(
    ymin = min(Beta, na.rm = TRUE),
    ymax = max(Beta, na.rm = TRUE),
    .groups = "drop"
  )
# Create plot: range bars per phylum per threshold
ggplot() + 
  geom_linerange(
    data = range_summary,
    aes(
      x = phylum, 
      ymin = ymin, 
      ymax = ymax, 
      color = phylum, 
      alpha = as.numeric(as.character(threshold))
    ), 
    size = 4, 
    position = position_dodge(width = 0)
  ) +
  geom_jitter(
    data = combined_plot_df,
    aes(
      x = phylum, 
      y = Beta, 
      alpha = as.numeric(as.character(threshold)), 
      color = phylum
    ),
    size = 1.5,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7)
  ) +
  scale_color_brewer(palette = "Dark2", guide="none") +
  scale_alpha_continuous(range = c(0.1, 1)) +
  labs(
    title = "Range of Fire Effect Size (Beta) by Fungal Phylum",
    subtitle = "Alpha by increasing certainty thresholds",
    x = "Phylum",
    y = "Fire Effect Size (Beta)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )
