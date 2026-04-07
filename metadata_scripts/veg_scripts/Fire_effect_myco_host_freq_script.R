library(tidyverse)
library(purrr)
library(broom.mixed)
library(lme4)

#read in site level data
#myco_freq_final_sum_site<-readRDS('Processed_data/metadata/veg/myco_freq_Site_level_df.Rdata')

#plot level data
myco_freq_final_sum_plot<-readRDS('Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata')




#read in all site locations
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% select(Site,Plot,Fire_Treatment)


veg_Site<-right_join(All_Sites,myco_freq_final_sum_plot) %>% 
  mutate(log_Freq_non_host=log10(freq_non_host+1E-6)
  ) %>% 
  mutate(Fire_treatment = factor(Fire_Treatment, levels = c("U", "B")))


vars<-c( 'freq_total_myco','freq_AM','freq_ECM','log_Freq_non_host')

hist(log10(veg_Site$freq_ECM))

# Fit a univariate model for each explanatory variable
model_results_veg <- map_dfr(vars, function(var) {
  formula <- as.formula(paste(var, '~ Fire_Treatment + (1|Site)'))
  model <- lmer(formula, data = veg_Site)
  #summary_model <- summary(model)
  anov <- Anova(model, test = "F")
  tidy_model <- broom.mixed::tidy(model, effects = "fixed") %>% 
    filter(term == "Fire_TreatmentB") 
  
  # 5. Combine into summary table
  tibble(
    variable = var,
    Factor = rownames(anov)[1],
    Estimate = tidy_model$estimate,
    Std_Error = tidy_model$std.error,
    DF_num = anov$Df,
    DF_denom = anov$Df.res,
    F = anov$F,
    P = anov$`Pr(>F)`[1]
  )
  
}) %>%
  arrange(desc(P))

model_results_veg


levels(veg_Site$Fire_Treatment)
vars


# Nest data by Site
nested_data <- veg_Site %>%
  group_by(Site) %>%
  nest()
# Fit univariate models for each variable within each Site
model_results_by_site <- nested_data %>%
  mutate(models = map(data, function(df) {
    map_dfr(vars, function(var) {
      if (all(is.na(df[[var]])) || length(unique(df[[var]])) <= 1) return(tibble())
      if (length(unique(df$Fire_Treatment)) <= 1) return(tibble())
      
      model_formula <- as.formula(paste0(var, "~ Fire_Treatment + (1|Site)"))
      model <- tryCatch(lmer(model_formula, data = df), error = function(e) NULL)
      
      if (is.null(model)) return(tibble())
      
      anov <- Anova(model, test = "F")
      #summary_model <- summary(model)
      tidy_model <- broom.mixed::tidy(model, effects = "fixed") %>% 
        filter(term == "Fire_TreatmentB")
      
      tibble(
        variable = var,
        Factor = rownames(anov)[1],
        Estimate = tidy_model$estimate,
        Std_Error = tidy_model$std.error,
        DF_num = anov$Df,
        DF_denom = anov$Df.res,
        F = anov$F,
        P = anov$`Pr(>F)`[1]
      )
    })
  })) %>%
  select(Site, models) %>%
  unnest(models) %>%
  arrange(Site,variable, desc(p_value))

# View result
model_results_by_site

signif_veg_fire<-model_results_by_site %>% 
  filter(p_value<0.11) %>% arrange(Site)

signif_veg_fire


veg_long <- veg_Site %>%
  filter(Site %in% signif_veg_fire$Site) %>%
  pivot_longer(
    cols = c(freq_total_myco, freq_AM, freq_ECM),
    names_to = "Myco_Type",
    values_to = "Frequency"
  )

ggplot(veg_long, aes(x = Fire_Treatment, y = Frequency, fill = Myco_Type)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, color = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              alpha = 0.5, size = 1, color = "black") +
  facet_wrap(~ Site, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Significant Sites: Mycorrhizal Frequencies by Fire Treatment",
    x = "Fire Treatment",
    y = "Frequency of Mycorrhizal Hosts",
    fill = "Mycorrhizal Type"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

