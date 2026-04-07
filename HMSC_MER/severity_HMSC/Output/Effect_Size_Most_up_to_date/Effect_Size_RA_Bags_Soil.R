library(tidyverse)
library(car)
library(performance)
library(broom)

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#load most up to date PA model
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

# ============================================================
# 1. Prepare Data
# ============================================================

Beta_estimates_PA <- Beta1$mean

Beta_estimates_PA <- Beta_estimates_PA %>%
  dplyr::select(OTU = Species, Beta = Severity)

# Step 1: Total counts across all samples
hyph_count <- myco_dat_hyph %>% summarise(total = sum(count, na.rm = TRUE)) %>% pull()
soil_count <- myco_dat_soil %>% summarise(total = sum(count, na.rm = TRUE)) %>% pull()

# Step 2: Calculate RA across entire dataset (not per sample)
hyph_RA <- myco_dat_hyph %>%
  group_by(OTU) %>%
  summarise(RA_hyph = sum(count, na.rm = TRUE) / hyph_count, .groups = "drop")

soil_RA <- myco_dat_soil %>%
  group_by(OTU) %>%
  summarise(RA_soil = sum(count, na.rm = TRUE) / soil_count, .groups = "drop")

RA_joined <- full_join(hyph_RA, soil_RA, by = "OTU") %>%
  mutate(RA_hyph = replace_na(RA_hyph, 0))

# Step 3b: Calculate the min non-zero ratio
min_ratio <- RA_joined %>%
  filter(RA_hyph > 0 & RA_soil > 0) %>%
  mutate(ratio = RA_hyph / RA_soil) %>%
  summarise(min_ratio = min(ratio, na.rm = TRUE)) %>%
  mutate(pseudocount = round(min_ratio / 2, 4)) %>%
  pull(pseudocount)

# Step 4: Join, calculate log-ratio, and merge with Beta_estimates
RA_ratio_PA <- RA_joined %>%
  mutate(
    ratio = RA_hyph / RA_soil,
    log_ratio = log10(ratio + min_ratio)
  ) %>%
  inner_join(Beta_estimates_PA, by = "OTU")

RA_ratio_PA <- RA_ratio_PA %>%
  mutate(Abs_Beta = abs(Beta))

# ============================================================
# 2. Helper Function
# ============================================================

analyze_beta_type <- function(data, y_var, beta_type_label, plot_color) {
  # Fit model
  model <- lm(as.formula(paste(y_var, "~ log_ratio")), data = data)
  model_check <- check_model(model)
  smry <- summary(model)
  anov <- Anova(model, test.statistic = "F")
  
  # ANOVA table
  anova_tbl <- tibble(
    Type = beta_type_label,
    sum_sq = anov$`Sum Sq`[1],
    sum_sq_resid = anov$`Sum Sq`[2],
    DF_num = anov$Df[1],
    DF_denom = anov$Df[2],
    F = anov$`F value`[1],
    P = anov$`Pr(>F)`[1],
    R2 = smry$r.squared,
    adj_R2 = smry$adj.r.squared
  )
  
  anova_tbl <- left_join(
    tidy(model) %>%
      filter(term == "log_ratio") %>%
      mutate(Type = beta_type_label) %>%
      dplyr::select(Type, Estimate = estimate, Std_Error = std.error),
    anova_tbl,
    by = "Type"
  )
  
  # Label text (render-safe for ggplot >= 3.5)
  if (anova_tbl$P > 0.10) {
    label_text_p <- "ns"
    parse_flag <- FALSE
  } else {
    label_text_p <- paste0(
      "italic(p):", round(anova_tbl$P, 2)    )
    parse_flag <- TRUE
  }
  # Label text (render-safe for ggplot >= 3.5)
  if (anova_tbl$P > 0.10) {
    label_text_R <- ""
    parse_flag <- FALSE
  } else {
    label_text_R <- paste0(
      "R^2: ", round(anova_tbl$R2, 2)
    )
    parse_flag <- TRUE
  }
  
  # Plot
  plot <- ggplot(data, aes(x = log_ratio, y = !!sym(y_var))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1, alpha = 0.6) +
    geom_point(size = 4, alpha = 0.6, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = plot_color, linewidth = 2) +
    theme_classic() +
    labs(
      y = eval(parse(text = beta_type_label)),
      x = expression(log[10]~"(RA Hyphae / RA Soil)")
    ) +
    annotate(
      "text",
      x = min(data$log_ratio, na.rm = TRUE) + 0.4,
      y = max(data[[y_var]], na.rm = TRUE) - 0.03,
      label = label_text_p,
      parse = parse_flag,
      hjust = 0, size = 8, color = "black"
    ) +  
    annotate(
      "text",
      x = min(data$log_ratio, na.rm = TRUE) + 0.4,
      y = max(data[[y_var]], na.rm = TRUE) - 0.055,
      label = label_text_R,
      parse = parse_flag,
      hjust = 0, size = 8, color = "black"
    ) +
    theme(
      axis.text.x = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.line = element_line(linewidth = 2, colour = "black")
    )
  
  return(list(anova = anova_tbl, plot = plot, model_check = model_check))
}

# ============================================================
# 3. Run Analyses
# ============================================================

relevant_otus <- myco_tax_soil %>%
  filter(guild2 == "EctoMycorrhizal") %>%
  pull(OTU)

# Optionally filter relevant OTUs
# RA_ratio_PA <- RA_ratio_PA %>% filter(OTU %in% relevant_otus)

res_abs <- analyze_beta_type(RA_ratio_PA, "Abs_Beta", "expression(paste('Absolute likelihood of occurrence post-fire (', beta[fire], ')'))", "green")
res_raw <- analyze_beta_type(RA_ratio_PA, "Beta", "expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')'))", "purple4")
res_pos <- analyze_beta_type(RA_ratio_PA %>% filter(Beta > 0), "Beta", "expression(paste('Positive likelihood of occurrence post-fires (', beta[fire], '>0)'))", "blue")
res_neg <- analyze_beta_type(RA_ratio_PA %>% filter(Beta < 0), "Beta", "expression(paste('Negative likelihood of occurrence post-fire (', beta[fire], '<0)'))", "red")

# Combine ANOVA tables
anova_table_all <- bind_rows(
  res_abs$anova,
  res_raw$anova,
  res_pos$anova,
  res_neg$anova
)

# Output
print(anova_table_all)

save(anova_table_all, file = "HMSC_MER/Output/Processed_Data/Effect_Size_RA_Beta.RDS")

# Optional: view plots
print(res_abs$plot)
print(res_raw$plot)
print(res_pos$plot)
print(res_neg$plot)
library(patchwork)

# Tag and style each plot
plot_a <- res_raw$plot +
  labs(tag = "a)") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

plot_b <- res_abs$plot +
  labs(tag = "b)") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

plot_c <- res_pos$plot +
  labs(tag = "c)") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

plot_d <- res_neg$plot +
  labs(tag = "d)") +
  theme(plot.tag = element_text(size = 16, face = "bold"))

# Arrange plots
combined_plot <- (plot_a + plot_b) /
  (plot_c + plot_d)

final_plot<-combined_plot

# View plot
final_plot


# Arrange the plots using patchwork
res_abs$model_check
res_raw$model_check
res_pos$model_check
res_neg$model_check


ggsave("plots/RA_OTU_Effect_Size_Plot.png", final_plot, width = 14, height = 14, dpi = 300)
