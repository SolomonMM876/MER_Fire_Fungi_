# ----------------------------
# Load libraries
# ----------------------------
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(performance)
library(ggpubr)
library(emmeans)
library(broom.mixed)
library(patchwork)
library(cowplot)

# ----------------------------
# Load data
# ----------------------------
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
load('HMSC_MER/results/Beta1.RData')  # Beta1 model
load('HMSC_MER/Output/Processed_Data/Vegitation_alpha_metrics_raw.RDS')

# Beta estimates from HMSC model
Beta_estimates_PA <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Severity)

# ----------------------------
# Prepare data
# ----------------------------
# Merge soil OTUs with Beta, taxonomy, and alpha metrics
Site_OTU <- myco_dat_soil %>%
  distinct(Site, Plot, OTU, .keep_all = TRUE) %>%
  left_join(Beta_estimates_PA, by = "OTU") %>%
  left_join(myco_tax_soil) %>%
  left_join(Veg_richness_plot)

# Select alpha metrics
alpha <- c("Shannon", "Simpson", "Chao1", "Observed", "Pielou")

# ------------------------------------------------------------
# A. Interaction models: Beta ~ alpha * guild2 + (1|Site/Plot)
# ------------------------------------------------------------
interaction_outputs <- purrr::map(alpha, function(alpha) {
  
  # Formula: Beta ~ alpha_metric * guild2 + (1|Site/Plot)
  formula <- as.formula(paste0("Beta ~ ", alpha, " * guild2 + (1|Site/Plot)"))
  
  model <- lmer(formula, data = Site_OTU)
  anov <- Anova(model, test = "F")
  tidy_model <- tidy(model, effects = "fixed")
  
  # Extract R²
  R2_vals <- performance::r2(model)
  
  # Create ANOVA summary
  anova_tbl <- anov %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(
      alpha_metric = alpha,
      R2_marginal = R2_vals$R2_marginal,
      R2_conditional = R2_vals$R2_conditional
    )
  
  # Get p-value for interaction term
  interaction_p <- anova_tbl %>%
    filter(str_detect(Factor, paste0(alpha, ":guild2"))) %>%
    pull(`Pr(>F)`)
  
  # Build plot with interaction
  plot_inter <- ggplot(Site_OTU, aes(x = .data[[alpha]], y = Beta, color = guild2)) +
    geom_point(size = 4, alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 2) +
    labs(
      title = paste0("Interaction: ", alpha, " × guild2"),
      y = expression(beta[fire]),
      x = paste0("Alpha Diversity (", alpha, ")")
    ) +
    annotate("text",
             x = min(Site_OTU[[alpha]], na.rm = TRUE),
             y = max(Site_OTU$Beta, na.rm = TRUE),
             label = paste0("P (interaction) = ", signif(interaction_p, 3)),
             hjust = 0, size = 6) +
    theme_classic(base_size = 18)
  
  return(list(
    anova = anova_tbl,
    plot = plot_inter
  ))
})

# Save ANOVA + plots for interaction models
interaction_anova <- bind_rows(purrr::map(interaction_outputs, "anova"))
interaction_plots <- purrr::map(interaction_outputs, "plot")

# Combine first four plots nicely
interaction_plot_grid <- wrap_plots(interaction_plots, ncol = 2)


#ggsave("plots/Interaction_Alpha_Guild_Beta.png", interaction_plot_grid, width = 30, height = 30, dpi = 300)

# ------------------------------------------------------------
# B. Independent models: Beta ~ alpha + (1|Site/Plot)
# ------------------------------------------------------------
main_outputs <- purrr::map(alpha, function(alpha) {
  
  formula <- as.formula(paste0("Beta ~ ", alpha, " + (1|Site/Plot)"))
  model <- lmer(formula, data = Site_OTU)
  anov <- Anova(model, test = "F")
  tidy_model <- tidy(model, effects = "fixed")
  
  R2_vals <- performance::r2(model)
  
  anova_tbl <- anov %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(
      alpha_metric = alpha,
      Estimate = tidy_model$estimate[tidy_model$term == alpha],
      Std_Error = tidy_model$std.error[tidy_model$term == alpha],
      R2_marginal = R2_vals$R2_marginal,
      R2_conditional = R2_vals$R2_conditional
    )
  
  # Build plot
  plot_main <- ggplot(Site_OTU, aes(x = .data[[alpha]], y = Beta)) +
    geom_point(size = 4, alpha = 0.6, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 2, color = "purple4") +
    labs(
      title = paste0("Main Effect: ", alpha),
      y = expression(beta[fire]),
      x = paste0("Alpha Diversity (", alpha, ")")
    ) +
    theme_classic(base_size = 18)
  
  return(list(
    anova = anova_tbl,
    plot = plot_main
  ))
})

# Save ANOVA + plots for main effects
main_anova <- bind_rows(purrr::map(main_outputs, "anova"))
main_plots <- purrr::map(main_outputs, "plot")
main_plot_grid <- wrap_plots(main_plots, ncol = 2)
#ggsave("plots/Main_Alpha_Beta.png", main_plot_grid, width = 30, height = 30, dpi = 300)

# Save outputs
#save(interaction_anova, file = "HMSC_MER/Output/Processed_Data/Interaction_Alpha_Beta_ANOVA.RDS")
#save(main_anova, file = "HMSC_MER/Output/Processed_Data/Main_Alpha_Beta_ANOVA.RDS")
