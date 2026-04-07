library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(ggpubr)
library(emmeans)
library(lmerTest)

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
  dplyr::select(OTU=Species,Beta=Severity)

Beta_estimates_COP<-Beta_estimates_COP %>% 
  dplyr::select(OTU=Species,Beta=Severity)


# Calculate Beta across entire dataset

soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site,Plot,OTU) 

#read in NPP data
ndvi_plot_post_burn<-readRDS('Processed_data/metadata/NPP_MER_by_plot.Rdata')


##############################
#######PA MODEL###############
##########################


# Join and merge with Beta_estimates
Beta_PA_Site_soil<- soil_Site_OTU %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta=abs(Beta)) %>% 
  left_join(ndvi_plot_post_burn)

colnames(ndvi_plot_post_burn)


# Plot for pH#######
ggplot(Beta_PA_Site_soil, aes(x = mean_NDVI , y = abs_Beta)) +
  geom_point( alpha = 0.9, size = 1 ) +
  geom_smooth(method='lm')+
  labs(y = "Beta values per OTU", x = "mean_NDVI", title = "Beta values by mean_NDVI ") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

hist(Beta_PA_Site_soil$mean_NDVI )
#PA model all
model <- lmer(abs_Beta ~ mean_NDVI + (1|Site/Plot), data = Beta_PA_Site_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
anova(model)

Anova(model, test = "F")

##############################################################
############interaction model#############
##################NOW testing the RA model#######################################
#########################################

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
RA_ratio_PA <- hyph_RA %>%
  inner_join(soil_RA, by = "OTU") %>%
  mutate(
    ratio = (RA_hyph) / (RA_soil),  # pseudocount to avoid log(0)
    log_ratio = log10(ratio)
  ) %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta=abs(Beta)) %>% 
  left_join(soil_Site_OTU)%>% 
  left_join(ndvi_plot_post_burn)


RA_ratio_PA %>% 
  ggplot(aes(x = log_ratio, y = Beta)) +
  geom_point() +
  geom_hline(yintercept =0)+
  theme_minimal() +
  labs(
    x = "log10(Hyphal rel. abundance / Soil rel. abundance)",
    y = "|Fire effect size|",
    title = "Relationship between OTU enrichment in hyphae vs soil and Fire effect size"
  )

#Absolute effect size
hist((RA_ratio_PA$abs_Beta))

# Additive model
mod_add <- lmer(abs_Beta ~ log_ratio + mean_NDVI  + (1|Site/Plot) , data = RA_ratio_PA)

# Interaction model
mod_int <- lmer(abs_Beta ~ log_ratio * mean_NDVI  + (1|Site/Plot) , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison

model <- lmer(abs_Beta ~ log_ratio * mean_NDVI  + (1|Site/Plot) , data = RA_ratio_PA)


model_summary <- summary(model)
check_model(model)
model_summary
Anova_model<-anova(model)
Anova_model


library(ggeffects)










#######################
######COP MODEL###############
##############################
# Join and merge with Beta_estimates
Beta_COP_Site_soil<- soil_Site_OTU %>%
  filter(OTU %in% myco_tax_hyph$OTU) %>% 
  left_join(Beta_estimates_COP, by = "OTU") %>% 
  mutate(    Type= "Soil")

Beta_COP_Site_hyph<- hyph_Site_OTU %>%
  inner_join(Beta_estimates_COP, by = "OTU") %>% 
  mutate(    Type= "Hyph")

Beta_Site_COP_both <- bind_rows(Beta_COP_Site_soil, Beta_COP_Site_hyph)


# Plot for both
ggplot(Beta_Site_COP_both, aes(x = Site, y = Beta, fill = Type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta values per OTU", x = "Site", title = "Beta values by Site") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


#COP model
model <- lm(Beta ~ Type * Site, data = Beta_Site_COP_both)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)

anova(model)



######### NOW just plot for soil OTUs
Beta_COP_Site_soil_all<- soil_Site_OTU %>%
  #filter(OTU %in% hyph_Site_OTU$OTU) %>% 
  left_join(Beta_estimates_COP, by = "OTU") %>% 
  mutate(    Type= "Soil")

# Plot for both
ggplot(Beta_COP_Site_soil_all, aes(x = Site, y = Beta, fill = Site)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta values per OTU", x = "Site", title = "Beta values by Site for COP model") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))


#COP model all
model <- lm(Beta ~ Site, data = Beta_COP_Site_soil_all)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)

anova(model)
