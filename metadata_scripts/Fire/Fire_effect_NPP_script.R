library(tidyverse)
library(purrr)

#plot level NPP data
NPP_MER<-readRDS('Processed_data/metadata/NPP_MER_by_plot.Rdata') 

#read in all site locations
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% select(Site,Plot,Fire_Treatment)


NPP_Site<-left_join(All_Sites,NPP_MER)

vars<-c('mean_NDVI')

hist((NPP_Site$mean_NDVI))

# Fit a univariate model for each explanatory variable
model_results_NPP<- map_dfr(vars, function(var) {
  formula <- as.formula(paste(var, '~ Fire_Treatment'))
  model <- lm(formula, data = NPP_Site)
  summary_model <- summary(model)
  tidy_model <- broom::tidy(model)
  
  tibble(
    variable = var,
    df1 = summary_model$fstatistic[2],
    df2 = summary_model$fstatistic[3],
    estimate = tidy_model$estimate[2],
    std_error = tidy_model$std.error[2],
    p_value = tidy_model$p.value[2],
    adj_r2 = summary_model$adj.r.squared,
    f_statistic = summary_model$fstatistic[1],
    r2 = summary_model$r.squared,
  )
}) %>%
  arrange(desc(adj_r2))

model_results_NPP




# Nest data by Site
nested_data <- NPP_Site %>%
  filter(!is.na(Fire_Treatment)) %>%
  group_by(Site) %>%
  nest()
# Fit univariate models for each variable within each Site
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
  arrange(Site,variable, desc(p_value))

# View result
model_results_by_site

signif_NPP_fire<-model_results_by_site %>% 
  filter(p_value<0.11) %>% arrange(Site)

signif_NPP_fire

# Filter significant sites and compute mean and SE if needed
ndvi_summary <- NPP_Site %>%
  filter(Site %in% signif_NPP_fire$Site) 

# Plot
ggplot(ndvi_summary, aes(x = Fire_Treatment, y = mean_NDVI)) +
  geom_boxplot(aes(fill = Fire_Treatment), width = 0.6, color = "black", show.legend = FALSE) +
  geom_point()+
  facet_wrap(~ Site, scales = "free_y") +
  theme_bw() +
  labs(
    title = "NDVI  by Fire Treatment at Significant Sites",
    x = "Fire Treatment",
    y = "Mean NDVI"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

