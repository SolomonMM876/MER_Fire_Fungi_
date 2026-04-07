library(tidyverse)
library(purrr)
library(broom)

#import soil and hyphal permanovas
permanova_hyph<-readRDS('Processed_data/Seq_dat/permANOVA_data/hyph_perm.Rdata')
permanova_soil<-readRDS('Processed_data/Seq_dat/permANOVA_data/soil_perm.Rdata')

#read in NPP data
NPP_MER<-readRDS('Processed_data/metadata/NPP_MER_by_site.Rdata') 
#read in NPP data since last_burn
NPP_MER_last_burn<-readRDS('Processed_data/metadata/NPP_MER_since_fire.Rdata') 

#read in aridity data
aridity_MER<-readRDS('Processed_data/metadata/Aridity_MER.Rdata')

#read in Tmax data
tmax_MER<-readRDS('Processed_data/metadata/Tmax_MER.Rdata') 

#read in Tmin data
tmin_MER<-readRDS('Processed_data/metadata/Tmin_MER.Rdata')

#read in fire regime data
fire_regime_MER<-readRDS('Processed_data/metadata/fire/Fire_regime_MER.Rdata')

#read in myco host fre df
myco_freq_final_sum<-readRDS('Processed_data/metadata/veg/myco_freq_Site_level_df.Rdata')

#dif in host freq
myco_host_diff_final<-readRDS('Processed_data/metadata/veg/myco_freq_Site_differences.Rdata')

#unburnt host freq
host_freq_unburnt<-readRDS('Processed_data/metadata/veg/myco_freq_Site_unburnt_plots.Rdata')


soil_perm<-permanova_soil %>% 
  left_join(aridity_MER) %>% 
 # left_join(NPP_MER) %>% 
  left_join(NPP_MER_last_burn) %>% 
  left_join(tmax_MER) %>% 
  left_join(tmin_MER) %>% 
  left_join(fire_regime_MER)


# Define your explanatory variables
vars <- c("aridity","NPP_since_fire", "tmax", "tmin", "severity", "year_of_most_recent_fire","fire_interval")

# Fit a univariate model for each explanatory variable
model_results_soil <- map_dfr(vars, function(var) {
  formula <- as.formula(paste("R2 ~", var))
  model <- lm(formula, data = soil_perm)
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
    AIC = AIC(model)
  )
}) %>%
  arrange(desc(adj_r2))

model_results_soil
##################################################
#veg characteristics 
soil_perm_veg<-permanova_soil %>% 
  left_join(myco_host_diff_final) %>% 
  left_join(myco_freq_final_sum) %>% 
  left_join(host_freq_unburnt)


ggplot(soil_perm_veg, aes(x = freq_AM_U, y = R2)) +
  #geom_point() +
  geom_text(aes(label=Site))+
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "freq_AM hosts in unburnt",
       y = expression(R^2),
       title = "Relationship number of AM hosts and R²") +
  theme_minimal()

#plotly::ggplotly(t)
colnames(host_freq_unburnt)


vars_veg <- c("Site" ,"freq_AM_U" ,"freq_ECM_U", "freq_total_myco_U", "freq_non_host_U",  
  "diff_total_myco" ,     "diff_non_host"    ,    "diff_AM"       ,       "diff_ECM" ,
              "freq_total_myco",'total_non_myco_hosts', "freq_AM","freq_ECM")

# Fit a univariate model for each explanatory variable
model_results_soil_veg <- map_dfr(vars_veg, function(var) {
  formula <- as.formula(paste("R2 ~", var))
  model <- lm(formula, data = soil_perm_veg)
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
    AIC = AIC(model)
  )
}) %>%
  arrange(desc(adj_r2))

model_results_soil_veg


#saveRDS(model_results_soil_veg,)














library(ggcorrplot) # For plotting

# Select only numeric explanatory variables
explan_vars <- soil_perm %>%
  select(all_of(vars),-fire_interval) %>%
  distinct()

# Compute correlation matrix
cor_matrix <- cor(explan_vars, use = "pairwise.complete.obs")

# View plain matrix
print(round(cor_matrix, 2))
