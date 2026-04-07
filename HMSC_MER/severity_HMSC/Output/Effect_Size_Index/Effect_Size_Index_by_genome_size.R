library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(ggpubr)

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


#read in df extracted from Brundrett 2008 website
tax_gs<- readRDS('Processed_data/Seq_dat/datasets_external/tax_genome_size_MER.Rds')



# Select OTUs associated with mean_gs

soil_OTU_genome_size <- myco_tax_soil %>%
  left_join(tax_gs) %>%
  filter(!is.na(mean_gs)) %>% 
  distinct(OTU,mean_gs) 

###############################################
###PA MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_PA_genome_size_soil<- soil_OTU_genome_size %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  inner_join(Beta_estimates_PA, by = "OTU") 



# Plot
ggplot(Beta_PA_genome_size_soil, aes(x = mean_gs, y = Beta)) +
  geom_smooth()+
  geom_jitter( alpha = 0.9, size = 1  ) +
  labs(y = "Beta estimates", x = "mean_gs", title = "Beta values by mean_gs for hyph taxa PA model") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(palette = "Set2")

hist(Beta_PA_genome_size_soil$Beta)

#PA model
model <- lm(Beta ~  mean_gs, data = Beta_PA_genome_size_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)

anova(model)

anova_p <- anova(model)$`Pr(>F)`[1]  # First term in the model
# Get significance level
p_sig <- case_when(
  anova_p <= 0.001 ~ "***",
  anova_p <= 0.01 ~ "**",
  anova_p <= 0.05 ~ "*",
  anova_p <= 0.1 ~ "+",
  TRUE ~ "ns"
)

# Optional: also extract R-squared
r2 <- round(model_performance(model)$R2, 2)

# Create annotation label
model_label <- paste0("Fruiting Body Type ", p_sig, "\nR² = ", r2)

# Post-hoc pairwise comparisons for interaction
emm <- emmeans(model, ~ mean_gs)

# Tukey pairwise contrasts: Type within each mean_gs
contrasts_type_within_mean_gs <- emmeans(model, pairwise ~ mean_gs, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within_mean_gs$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ "+"
    ),
    xmin = str_trim(str_extract(contrast, "(^[^-]+)")),         # before "-"
    xmax = str_trim(ifelse(
      str_detect(contrast, "\\(.*\\)"),
      str_extract(contrast, "(?<=\\().*?(?=\\))"),   # extract inside parentheses
      str_extract(contrast, "(?<=-)\\s*.*$")         # fallback: extract after "-"
    )))



mean_gs_max <- Beta_PA_genome_size_soil %>%
  group_by(mean_gs) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

sig_results <- sig_labels %>%
  left_join(mean_gs_max, by = c("xmax"="mean_gs") )



# Main boxplot
Beta_PA_genome_size_soil %>% 
  ggplot( aes(x = mean_gs, y = Beta, fill=mean_gs)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.4,     position = position_dodge(width = 0.6)  # wider dodge increases sPAcing
  ) +
  geom_jitter( alpha = 0.9, size = 1,  width = 0.1 # must match boxplot dodge
  ) +
  labs(y = "Beta coefficent", x = "Fruiting body type", title = "Beta values by fruiting body type PA Model all taxa") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5))+
  geom_bracket(
    data = sig_results,
    aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
    inherit.aes = FALSE,
    tip.length = 0.01,
    size = 1.1
  )+
  annotate("text", x = 0, y = max(Beta_PA_genome_size_soil$Beta+0.05, na.rm = TRUE),
               label = model_label, hjust = 0, vjust = 1, size = 4.5)

###################################
##################NOW testing the RA model#######################################
#########with interaction
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
  left_join(soil_OTU_genome_size) %>% 
  filter(!is.na(mean_gs))

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
mod_add <- lm(abs_Beta ~ log_ratio + mean_gs , data = RA_ratio_PA)

# Interaction model
mod_int <- lm(abs_Beta ~ log_ratio * mean_gs , data = RA_ratio_PA)

anova(mod_add, mod_int)  # Likelihood Ratio Test
AIC(mod_add, mod_int)    # Model comparison

model <- lm(abs_Beta ~ log_ratio * mean_gs , data = RA_ratio_PA)


model_summary <- summary(model)
check_model(model)
model_summary

# Get the coefficients, standard error, p-value, and R²
slope <- coef(model_summary)[2, 1]  # Slope
se_slope <- coef(model_summary)[2, 2]  # Standard error of the slope
p_value <- coef(model_summary)[2, 4]  # P-value for the slope
r2 <- model_summary$r.squared  # R-squared value


# Tukey pairwise contrasts: Type within each mean_gs
contrasts_type_within_mean_gs <- emmeans(model, pairwise ~ mean_gs*log_ratio, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within_mean_gs$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.05) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ "+"
    ),
    xmin = str_trim(str_extract(contrast, "(^[^-]+)")),         # before "-"
    xmax = str_trim(ifelse(
      str_detect(contrast, "\\(.*\\)"),
      str_extract(contrast, "(?<=\\().*?(?=\\))"),   # extract inside parentheses
      str_extract(contrast, "(?<=-)\\s*.*$")         # fallback: extract after "-"
    )))



mean_gs_max <- Beta_PA_genome_size_soil %>%
  group_by(mean_gs) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

sig_results <- sig_labels %>%
  left_join(mean_gs_max, by = c("xmax"="mean_gs") )
 



ggplot(RA_ratio_PA, aes(x = log_ratio, y = abs_Beta, color = mean_gs)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  labs(
    title = "Interaction model: abs_Beta ~ log_ratio * mean_gs",
    y = "|Beta fire effect|",
    x = "log10(RA Hyphae / RA Soil)") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())
  









###############################################
###COP MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_COP_repo_strat_soil<- soil_OTU_repo_strat %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  inner_join(Beta_estimates_COP, by = "OTU") 



# Plot
ggplot(Beta_COP_repo_strat_soil, aes(x = mean_gs, y = Beta, fill=mean_gs)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta estimates", x = "mean_gs", title = "Beta values by mean_gs for COP model all taxa") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(palette = "Set2")


#COP model
model <- lm(Beta ~  mean_gs, data = Beta_COP_repo_strat_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)

anova(model)

