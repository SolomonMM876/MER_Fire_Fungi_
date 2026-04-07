library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)

#load Beta1 and 2 from soil data HMSC
load('HMSC_MER/results/mcmc_output.RData')

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#load most up to date PA model
load('HMSC_MER/results/Beta1.RData')

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean
Beta_estimates_COP<-Beta2$mean


Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)

Beta_estimates_COP<-Beta_estimates_COP %>% 
  dplyr::select(OTU=Species,Beta=Severity)

# Select OTUs associated with family

soil_OTU_phy <- myco_tax_soil %>%
  distinct(OTU,phylum,order,family,genus) 

###############################################
###PA MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_PA_order_soil<- soil_OTU_phy %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta=abs(Beta))



# Plot
ggplot(Beta_PA_order_soil, aes(x = order, y = Beta, fill=order)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
) +
  labs(y = "Beta estimates", x = "order", title = "Beta values by order for PA model") +
  theme_minimal()


#PA model
model <- lm(Beta ~ 0+  order, data = Beta_PA_phy_soil)



# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model,test= 'F')

# Post-hoc pairwise comparisons for interaction
emm <- emmeans(model, ~ order)

# Tukey pairwise contrasts: Type within each order
contrasts_type_within_order <- emmeans(model, pairwise ~ order, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within_order$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.1) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ paste0('p=',round(p.value,3))
    ),
    xmin = str_trim(str_extract(contrast, "^[^-]+")),         # before "-"
    xmax = str_trim(str_extract(contrast, "(?<=-)\\s*.*$"))   # after "-"
  )


 order_max <- Beta_PA_phy_soil %>%
  group_by(order) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

sig_results <- sig_labels %>%
  left_join(order_max, by = c("xmin"="order") )



# Main boxplot
Beta_PA_phy_soil %>% 
  ggplot( aes(x = order, y =Beta, fill=order)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.4, # wider dodge increases sPAcing
  ) +
  geom_jitter( aes(colour=order),alpha = 0.9, size = 2, width = 0.05 # must match boxplot dodge
  ) +
  labs(y = "Beta coefficent", x = "order", title = "Beta values by order PA Model all taxa") +
  theme_minimal() +
  geom_bracket(
    data = sig_results,
    aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
    inherit.aes = FALSE,
    tip.length = 0.01,
    size = 1.1
  )

#plot
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
  left_join(myco_tax_hyph)

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
mod_add <- lm(abs_Beta ~ log_ratio + order , data = RA_ratio_PA)

# Interaction model
mod_int <- lm(abs_Beta ~ log_ratio * order , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison

model <- lm(Beta ~ log_ratio + order , data = RA_ratio_PA)


model_summary <- summary(model)
check_model(model)
model_summary
Anova(model)
# Get the coefficients, standard error, p-value, and R²
# slope <- coef(model_summary)[2, 1]  # Slope
# se_slope <- coef(model_summary)[2, 2]  # Standard error of the slope
# p_value <- coef(model_summary)[2, 4]  # P-value for the slope
# r2 <- model_summary$r.squared  # R-squared value


# Tukey pairwise contrasts: Type within each order
contrasts_type_within_order <- emmeans(model, pairwise ~ order*log_ratio, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within_order$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ "+"
    ),
    xmin = str_trim(str_extract(contrast, "(?<=\\().*?(?= log_ratio)")),  # extract between ( and log_ratio
    xmax = str_trim(str_extract(contrast, "(?<=- \\().*?(?= log_ratio)")) # same for second group after the dash
  )



order_max <- RA_ratio_PA %>%
  group_by(order) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

sig_results <- sig_labels %>%
  left_join(order_max, by = c("xmax"="order") )




ggplot(RA_ratio_PA, aes(x = log_ratio, y = Beta, color = order)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  labs(
    title = "Interaction model: Beta ~ log_ratio * order",
    y = "Beta fire effect",
    x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

#######################################################################
#####################testing for heteroskedasticity
########################
library(olsrr)

model <- lm(Beta ~ log_ratio * order, data = RA_ratio_PA)

# Breusch-Pagan Test
ols_test_breusch_pagan(model)

# White's Test
ols_test_f(model)

plot(model)






###############################################
###COP MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_COP_phy_soil<- soil_OTU_phy %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  inner_join(Beta_estimates_COP, by = "OTU") 



# Plot
ggplot(Beta_COP_phy_soil, aes(x = order, y = Beta, fill=order)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta estimates", x = "order", title = "Beta values by order for COP model") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


#COP model
model <- lm(Beta ~  order, data = Beta_COP_phy_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)

anova(model)

