library(tidyverse)


AFA_invassives<-readRDS('processed_data/metadata/veg/AFA_native_invasive_freq_MER.Rds')



#read in all site locations
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% select(Site,Plot,Fire_Treatment)


veg_Site<-inner_join(All_Sites,AFA_invassives) %>% 
  mutate(log_rel_freq=log10(rel_freq)
  )



library(ggplot2)

ggplot(veg_Site, aes(x = Plot, y = rel_freq, fill = unified_status)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Plot",
    y = "Relative Frequency",
    fill = "Status",
    title = "Proportion of Species Categories per Plot"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    panel.spacing = unit(1, "lines")
  ) +
  guides(fill = guide_legend(title = "Species Status")) +
  facet_grid(. ~ Site , scales = "free_x", space = "free_x")

library(broom)

vars <- c("doubtfully_naturalised", "native", "native_colonising", "naturalised")

# Recode to avoid issues with column names
veg_wide <- veg_Site %>%
  mutate(unified_status = str_replace_all(unified_status, " ", "_")) %>%
  select(Site, Plot, Fire_Treatment, unified_status, rel_freq) %>%
  pivot_wider(
    names_from = unified_status,
    values_from = rel_freq,
    values_fill = 0  # assume 0 if status is absent in a plot
  )  %>% 
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")))



model_results_veg <- map_dfr(vars, function(var) {
  formula <- as.formula(paste(var, '~ Fire_Treatment + (1|Site)'))
  model <- lmer(formula, data = veg_wide)
  anov <- Anova(model, test = "F")
  tidy_model <- tidy(model, effects = "fixed")%>% 
    filter(term != "(Intercept)")
  
  # Compute R²
  R2_marg <- performance::r2(model)$R2_marginal
  R2_cond <- performance::r2(model)$R2_conditional
  # Return ANOVA and plot
  anova_tbl <- anov %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(
      Estimate = tidy_model$estimate,
      Std_Error = tidy_model$std.error,
      R2_marginal = R2_marg,
      R2_conditional= R2_cond
    ) %>%
    select(Factor,Estimate,Std_Error, Df,Df_res=Df.res, F, P = `Pr(>F)`,R2_marginal,R2_conditional)
  
}) %>%
  arrange(desc(R2_marginal))

model_results_veg



#nest by site, but needs work to work


nested_data <- veg_wide %>%
  filter(!is.na(Fire_Treatment)) %>%
  group_by(Site) %>%
  nest()

model_results_by_site <- nested_data %>%
  mutate(models = map(data, function(df) {
    map_dfr(vars, function(var) {
      if (all(is.na(df[[var]])) || length(unique(df[[var]])) <= 1) return(tibble())
      if (length(unique(df$Fire_Treatment)) <= 1) return(tibble())
      
      model_formula <- as.formula(paste(var, "~ Fire_Treatment"))
      model <- tryCatch(lm(model_formula, data = df), error = function(e) NULL)
      
      if (is.null(model)) return(tibble())
      
      summary_model <- summary(model)
      tidy_model <- tidy(model)
      
      tibble(
        variable = var,
        df1 = summary_model$fstatistic[2],
        df2 = summary_model$fstatistic[3],
        estimate = tidy_model$estimate[2],
        std_error = tidy_model$std.error[2],
        t_value = tidy_model$statistic[2],
        p_value = tidy_model$p.value[2],
        f_statistic = summary_model$fstatistic[1],
        r2 = summary_model$r.squared
      )
    })
  })) %>%
  select(Site, models) %>%
  unnest(models) %>%
  arrange(Site, variable, desc(p_value))

signif_veg_fire <- model_results_by_site %>%
  filter(p_value < 0.11) %>%
  arrange(Site)

veg_long <- veg_Site %>%
  mutate(unified_status = str_replace_all(unified_status, " ", "_")) %>%
  filter(Site %in% signif_veg_fire$Site, unified_status %in% signif_veg_fire$variable)

ggplot(veg_long, aes(x = Fire_Treatment, y = rel_freq, fill = unified_status)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Site, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Significant Sites: Species Status Relative Frequency by Fire Treatment",
    x = "Fire Treatment",
    y = "Relative Frequency",
    fill = "Species Status"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



