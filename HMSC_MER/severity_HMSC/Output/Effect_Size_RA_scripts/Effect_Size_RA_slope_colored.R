library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(performance)

#load Beta1 and 2 from soil data HMSC
load('HMSC_MER/results/Beta.RData')

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')


#read in df extracted from Brundrett 2008 website
final_repo_strat<- readRDS('Processed_data/Seq_dat/datasets_external/repro_strat_df.Rds')

repo_strat_OTU<- myco_tax_hyph%>% 
  left_join(final_repo_strat) %>% 
  mutate(
    repo_strategy = case_when(
      phylum == "Glomeromycota" ~ "Glomeromycota",  # or "microscopic", "none", etc.
      TRUE ~ repo_strategy
    )
  ) 

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
RA_ratio_PA <- hyph_RA %>%
  inner_join(soil_RA, by = "OTU") %>%
  mutate(
    ratio = (RA_hyph) / (RA_soil),  # pseudocount to avoid log(0)
    log_ratio = log10(ratio)
  ) %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Fire = abs(Fire)) %>% 
  left_join(myco_tax_hyph) 

RA_ratio_PA %>% 
  ggplot(aes(x = log_ratio, y = Fire)) +
  geom_point() +
  geom_hline(yintercept =0)+
  theme_minimal() +
  labs(
    x = "log10(Hyphal rel. abundance / Soil rel. abundance)",
    y = "|Fire effect size|",
    title = "Relationship between OTU enrichment in hyphae vs soil and Fire effect size"
  )


####################################
##############phylum###############
#####################################


#Absolute effect size
hist((RA_ratio_PA$abs_Fire))

# Additive model
mod_add <- lm(abs_Fire ~ log_ratio + phylum , data = RA_ratio_PA)

# Interaction model
mod_int <- lm(abs_Fire ~ log_ratio * phylum , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison


model <- lm(abs_Fire ~ log_ratio * phylum , data = RA_ratio_PA)


model_summary <- summary(model)
check_model(model)
model_summary

# Get the coefficients, standard error, p-value, and R²
slope <- coef(model_summary)[2, 1]  # Slope
se_slope <- coef(model_summary)[2, 2]  # Standard error of the slope
p_value <- coef(model_summary)[2, 4]  # P-value for the slope
r2 <- model_summary$r.squared  # R-squared value


# Create the ggplot
RA_ratio_PA %>% 
  ggplot( aes(x = log_ratio, y = abs_Fire)) +
  geom_point(aes(color=phylum ),size=3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "lm(abs_Fire ~ log_ratio + phylum )", y = "|Fire effect size|", x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal() +
  annotate("text", x = max(RA_ratio_PA$log_ratio) * 0.7, y = max(RA_ratio_PA$abs_Fire) * 0.9,
           label = paste("Slope: ", round(slope, 3), "\n",
                         "SE: ", round(se_slope, 3), "\n",
                         "p-value: ", round(p_value,3), "\n",
                         "R²: ", round(r2, 3)),
           hjust = 0, size = 4, color = "black")














#############################
########Repo_strat############
########################
# Step 3: Join, calculate log-ratio, and merge with Beta_estimates
RA_ratio_PA <- hyph_RA %>%
  inner_join(soil_RA, by = "OTU") %>%
  mutate(
    ratio = (RA_hyph) / (RA_soil),  # pseudocount to avoid log(0)
    log_ratio = log10(ratio)
  ) %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Fire = abs(Fire)) %>% 
  left_join(repo_strat_OTU) %>% 
  filter(!is.na(repo_strategy))



#Absolute effect size
hist((RA_ratio_PA$abs_Fire))

# Additive model
mod_add <- lm(abs_Fire ~ log_ratio + repo_strategy , data = RA_ratio_PA)

# Interaction model
mod_int <- lm(abs_Fire ~ log_ratio * repo_strategy , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison


model_summary <- summary(mod_add)
check_model(mod_add)
model_summary

# Get the coefficients, standard error, p-value, and R²
slope <- coef(model_summary)[2, 1]  # Slope
se_slope <- coef(model_summary)[2, 2]  # Standard error of the slope
p_value <- coef(model_summary)[2, 4]  # P-value for the slope
r2 <- model_summary$r.squared  # R-squared value


# Create the ggplot
RA_ratio_PA %>% 
ggplot( aes(x = log_ratio, y = abs_Fire)) +
  geom_point(aes(color=repo_strategy ),size=3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "lm(abs_Fire ~ log_ratio + repo_strategy )", y = "|Fire effect size|", x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal() +
  annotate("text", x = max(RA_ratio_PA$log_ratio) * 0.7, y = max(RA_ratio_PA$abs_Fire) * 0.9,
           label = paste("Slope: ", round(slope, 3), "\n",
                         "SE: ", round(se_slope, 3), "\n",
                         "p-value: ", round(p_value,3), "\n",
                         "R²: ", round(r2, 3)),
           hjust = 0, size = 4, color = "black")


##################################
######explo type##############
######################################

# Additive model
mod_add <- lm(abs_Fire ~ log_ratio + exploration_type, data = RA_ratio_PA)

# Interaction model
mod_int <- lm(abs_Fire ~ log_ratio * exploration_type, data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison


model_summary <- summary(mod_add)
check_model(mod_add)
model_summary

# Get the coefficients, standard error, p-value, and R²
slope <- coef(model_summary)[2, 1]  # Slope
se_slope <- coef(model_summary)[2, 2]  # Standard error of the slope
p_value <- coef(model_summary)[2, 4]  # P-value for the slope
r2 <- model_summary$r.squared  # R-squared value


# Create the ggplot
RA_ratio_PA %>% 
  ggplot( aes(x = log_ratio, y = abs_Fire)) +
  geom_point(aes(color=exploration_type),size=3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(title = "lm(abs_Fire ~ log_ratio + exploration_type)", y = "|Fire effect size|", x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal() +
  annotate("text", x = max(RA_ratio_PA$log_ratio) * 0.7, y = max(RA_ratio_PA$abs_Fire) * 0.9,
           label = paste("Slope: ", round(slope, 3), "\n",
                         "SE: ", round(se_slope, 3), "\n",
                         "p-value: ", round(p_value,3), "\n",
                         "R²: ", round(r2, 3)),
           hjust = 0, size = 4, color = "black")

########################################
##########Repeat above, but for Abundance models####################
############################################################################

# Step 3: Join, calculate log-ratio, and merge with Beta_estimates
RA_ratio_COP <- hyph_RA %>%
  inner_join(soil_RA, by = "OTU") %>%
  mutate(
    ratio = (RA_hyph) / (RA_soil),  # pseudocount to avoid log(0)
    log_ratio = log10(ratio)
  ) %>%
  inner_join(Beta_estimates_COP, by = "OTU") %>% 
  mutate(abs_Fire = abs(Fire)) %>% 
  left_join(myco_tax_soil)


RA_ratio_COP %>% 
  ggplot(aes(x = log_ratio, y = Fire)) +
  geom_point() +
  geom_hline(yintercept =0)+
  theme_minimal() +
  labs(
    x = "log10(Hyphal rel. abundance / Soil rel. abundance)",
    y = "|Fire effect size|",
    title = "Relationship between OTU enrichment in hyphae vs soil and Fire effect size COP model"
  )

library(performance)
#Absolute effect size
hist(log10(RA_ratio_COP$abs_Fire))
model_abs <- lm((abs_Fire) ~ log_ratio+repo_strategy , data = RA_ratio_COP)
#model_abs <- lm((abs_Fire) ~ log_ratio, data = RA_ratio_COP)


summary(model_abs)
check_model(model_abs)


ggplot(RA_ratio_COP, aes(x = log_ratio, y = abs_Fire)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Absolute Fire Effect Size vs. log10(Hyphae/Soil RA)", y = "|Fire effect size|", x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal()

