library(tidyverse)
library(lme4)
library(car)
library(emmeans)  # optional for post-hocs
library(lmerTest) # if you want p-values in summary()
library(ggplot2)
library(performance)



#load Beta1 and 2

myco_tax<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')


explo_RA <- myco_dat %>%
  filter(!is.na(exploration_type)) %>%
  group_by(Site,Plot,Fire_Treatment, exploration_type) %>%
  summarise(explo_type_count = sum(count), .groups = "drop") %>%
  group_by(Plot) %>%
  mutate(explo_Plot_count = sum(explo_type_count),
         RA_explo = explo_type_count / explo_Plot_count,
         log_RA_explo = log10(RA_explo)) %>%
  ungroup() 



explo_RA %>% 
  ggplot(aes(x=Fire_Treatment, y=RA_explo))+
  geom_boxplot()+
  geom_jitter(width=.2, color='grey')+
  facet_grid(~exploration_type)


##############
#one way to test explo type

explo_model_fire<-lmer((log_RA_explo)~ exploration_type *(Fire_Treatment) + 
                         (1|Site/Plot),
                       data=explo_RA)


check_model(explo_model_fire)

summary(explo_model_fire)
Anova_resin<-round(Anova(explo_model_fire,test='F'), 2) 
Anova_resin
###########################
#another way to test explo type by Ra of each in samples



# Filter to one exploration strategy
contact_df <- explo_RA %>% 
  filter(exploration_type == "contact")


# Fit the model 
contact_model <- lmer( log_RA_explo ~ Fire_Treatment + (1 | Site), 
  data = contact_df
)

check_model(contact_model)

summary(contact_model)
Anova_resin<-round(Anova(contact_model,test='F'), 2) 
Anova_resin





exploration_types <- explo_RA %>%
  filter(!is.na(exploration_type)) %>%
  distinct(exploration_type) %>%
  pull(exploration_type)

# Loop through types and run models
models_list <- exploration_types %>%
  set_names() %>%
  map(function(type) {
    df <- explo_RA %>%
      filter(exploration_type == type, !is.na(RA_explo))
    
    model <- lmer(
      log_RA_explo ~ Fire_Treatment + (1 | Site),
      data = df
    )
    
    # Optional: show diagnostic plots as you go
    plot(sim_res, main = type)
    qqnorm(resid(model), main = paste("QQ Plot:", type))
    qqline(resid(model))
    
    # Return ANOVA results and summary
    list(
      model = model,
      anova = round(Anova(model, test = "F"), 2),
      summary = summary(model)
    )
  })

# View ANOVA for a specific exploration type:
models_list$`short-distance_coarse`$anova

# View model summary:
models_list$`short-distance_coarse`$summary
############Make a table

library(dplyr)
library(purrr)

anova_table <- map2_dfr(models_list, names(models_list), function(mod, label) {
  mod$anova %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    mutate(exploration_type = label, .before = 1)
})
