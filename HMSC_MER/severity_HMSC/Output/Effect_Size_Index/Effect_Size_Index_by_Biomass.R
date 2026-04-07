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

#bag biomass
All_biomass<-readRDS('processed_data/Bag_data/MER_biomass.RDS')

##############################
#######PA MODEL###############
##########################


# Join and merge with Beta_estimates
Beta_PA_Site_soil<- soil_Site_OTU %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta=abs(Beta)) %>% 
  left_join(All_biomass)

colnames(All_biomass)


# Plot for pH#######
ggplot(Beta_PA_Site_soil, aes(x = biomass_kg_ha_day, y = abs_Beta)) +
  geom_point( alpha = 0.9, size = 1 ) +
  geom_smooth(method='lm')+
  labs(y = "Beta values per OTU", x = "biomass_kg_ha_day", title = "Beta values by biomass_kg_ha_day") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

hist(Beta_PA_Site_soil$log_biomass_kg_ha_day)
#PA model all
model <- lmer(Beta ~ log_biomass_kg_ha_day+ (1|Site/Plot), data = Beta_PA_Site_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)

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


RA_joined <- full_join(hyph_RA, soil_RA, by = "OTU") %>%
  mutate(
    RA_hyph = replace_na(RA_hyph, 0),
    RA_soil = replace_na(RA_soil, 0)
  )


# Step 3b: Calculate the min non-zero ratio
min_ratio <- RA_joined %>%
  filter(RA_hyph > 0 & RA_soil > 0) %>%  # ensure denominator isn't 0
  mutate(ratio = RA_hyph / RA_soil) %>%
  summarise(min_ratio = min(ratio, na.rm = TRUE)) %>%
  mutate(pseudocount = min_ratio / 2) %>%
  pull(pseudocount)

# Step 3: Join, calculate log-ratio, and merge with Beta_estimates
RA_ratio_PA <- RA_joined %>% 
  mutate(
    ratio = (RA_hyph + min_ratio) / (RA_soil + min_ratio),
    log_ratio = log10(ratio)
  ) %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta = abs(Beta)) %>% 
  left_join(soil_Site_OTU)%>% 
  left_join(All_biomass)


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
hist((RA_ratio_PA$Beta))

# Additive model
mod_add <- lmer(Beta ~ log_ratio + log_biomass_kg_ha_day + (1|Site/Plot) , data = RA_ratio_PA)

# Interaction model
mod_int <- lmer(Beta ~ log_ratio * log_biomass_kg_ha_day + (1|Site/Plot) , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison

model <- lmer(Beta ~ log_ratio * log_biomass_kg_ha_day + (1|Site/Plot) , data = RA_ratio_PA)


model_summary <- summary(model)
check_model(model)
model_summary
Anova_model<-anova(model)
Anova_model


library(ggeffects)


preds <- ggpredict(model, terms = c("log_ratio", "log_biomass_kg_ha_day"))

preds %>% 
  mutate()
ggplot( aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), 
              alpha = 0.2, color = NA) +
  labs(
    title = "Interaction between log_ratio and log_biomass",
    x = "log_ratio",
    y = "Predicted abs_Beta",
    color = "log_biomass",
    fill = "log_biomass"
  ) +
  theme_minimal()



# Get quantiles of biomass to represent "low", "medium", and "high"
biomass_vals <- quantile(RA_ratio_PA$log_biomass_kg_ha_day, probs = c(0.1, 0.5, 0.9), na.rm = TRUE)

# Create a grid of values to predict
new_data <- expand.grid(
  log_ratio = seq(min(RA_ratio_PA$log_ratio, na.rm = TRUE), 
                  max(RA_ratio_PA$log_ratio, na.rm = TRUE), 
                  length.out = 100),
  log_biomass_kg_ha_day = biomass_vals
)

# Predict using the model (fixed effects only)
new_data$predicted <- predict(model, newdata = new_data, re.form = NA)



# Label biomass levels
biomass_labels <- c("Low biomass", "Medium biomass", "High biomass")
new_data$biomass_group <- factor(new_data$log_biomass_kg_ha_day, 
                                 labels = biomass_labels)

# Plot
ggplot(new_data, aes(x = log_ratio, y = predicted, color = biomass_group)) +
  geom_line(size = 1.2) +
  labs(
    title = "Interaction: Slope of log_ratio across biomass levels",
    x = "log_ratio",
    y = "Predicted abs_Beta",
    color = "Biomass level"
  ) +
  theme_minimal()














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
