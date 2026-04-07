library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)

#load Beta1 and 2 from soil data HMSC
load('HMSC_MER/results/mcmc_output.RData')

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')


#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean
Beta_estimates_COP<-Beta2$mean


Beta_estimates_PA<-Beta_estimates_PA %>% 
  select(OTU=Species,Fire=Severity)

Beta_estimates_COP<-Beta_estimates_COP %>% 
  select(OTU=Species,Fire=Severity)

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

# Step 3: Join, calculate log-ratio, and merge with Beta_estimates
Effect_Index_PA <- hyph_RA %>%
  inner_join(soil_RA, by = "OTU") %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Fire = abs(Fire),
         Effect_Index_hyph= Fire*RA_hyph,
         Effect_Index_soil= Fire*RA_soil,
         Effect_Index_hyph_abs= abs_Fire*RA_hyph,
         Effect_Index_soil_abs= abs_Fire*RA_soil
      )

# Pivot to long format
Effect_Index_PA_long <- Effect_Index_PA %>%
  pivot_longer(cols = c(Effect_Index_hyph, Effect_Index_soil, Effect_Index_hyph_abs, Effect_Index_soil_abs),
               names_to = "Effect_Type",
               values_to = "Effect_Value")

# Plot
ggplot(Effect_Index_PA_long, aes(x = Effect_Type, y = Effect_Value)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA) + 
  geom_jitter(width = 0.2, alpha = 0.3, color = "darkblue") +
  theme_minimal() +
  labs(x = "Effect Index Type",
       y = "Effect Index Value(RA x Beta coefficent per OTU)",
       title = "Effect Index Values by Type")
