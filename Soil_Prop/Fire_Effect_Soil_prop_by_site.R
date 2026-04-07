library(tidyverse)
library(purrr)
library(broom)

# Read in data
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)

MER_nute <- readRDS('Processed_data/nutrients/MER_All_Soil.rds')

Nute_Site <- left_join(All_Sites, MER_nute) %>% 
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")))

vars <- c('PO4', 'Soil_pH', 'NH4', 'NO3')

# -----------------------------
# 1. OVERALL MODELS
# -----------------------------
model_results_soil <- map_dfr(vars, function(var) {
  formula <- as.formula(paste(var, '~ Fire_Treatment'))
  model <- lm(formula, data = Nute_Site)
  summary_model <- summary(model)
  tidy_model <- broom::tidy(model)
  
  tibble(
    variable = var,
    estimate = tidy_model$estimate[2],
    std_error = tidy_model$std.error[2],
    p_value = tidy_model$p.value[2],
    adj_r2 = summary_model$adj.r.squared
  )
}) %>%
  mutate(
    sig_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

plot_data_overall <- Nute_Site %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  left_join(model_results_soil %>% select(variable, sig_label), by = "variable") %>%
  group_by(variable) %>%
  mutate(y_pos = max(value, na.rm = TRUE) * 1.05)

p_overall <- ggplot(plot_data_overall, aes(x = Fire_Treatment, y = value, fill = Fire_Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_fill_manual(values = c(B = "#9F2121", U = "#5497B6"))+
  geom_text(aes(x = 1.5, y = y_pos, label = sig_label), inherit.aes = FALSE) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(x = "Fire Treatment", y = "Value (mg/kg or pH)", title = "Overall soil property differences") +
  theme_bw()

p_overall




# -----------------------------
# 2. PER-SITE MODELS
# -----------------------------
nested_data <- Nute_Site %>%
  filter(!is.na(Fire_Treatment)) %>%
  group_by(Site) %>%
  nest()

model_results_by_site <- nested_data %>%
  mutate(models = map(data, function(df) {
    map_dfr(vars, function(var) {
      if (all(is.na(df[[var]]))) return(tibble())
      
      model_formula <- as.formula(paste(var, "~ Fire_Treatment"))
      model <- tryCatch(lm(model_formula, data = df), error = function(e) NULL)
      if (is.null(model)) return(tibble())
      
      summary_model <- summary(model)
      tidy_model <- tidy(model)
      
      tibble(
        variable = var,
        estimate = tidy_model$estimate[2],
        p_value = tidy_model$p.value[2]
      )
    })
  })) %>%
  select(Site, models) %>%
  unnest(models) %>%
  mutate(
    sig_label = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

plot_data_by_site <- Nute_Site %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "value") %>%
  left_join(model_results_by_site %>% select(Site, variable, sig_label), by = c("Site", "variable")) %>%
  group_by(Site, variable) %>%
  mutate(y_pos = max(value, na.rm = TRUE) * 1.05)

p_by_site <- ggplot(plot_data_by_site, aes(x = Fire_Treatment, y = value, fill = Fire_Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  geom_text(aes(x = 1.5, y = y_pos, label = sig_label), inherit.aes = FALSE) +
  facet_grid(variable ~ Site, scales = "free_y") +
  labs(x = "Fire Treatment", y = "Value", title = "Per-site nutrient differences") +
  theme_bw()

# -----------------------------
# Print plots
# -----------------------------
p_overall
p_by_site
