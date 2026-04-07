# ==============================================================================
# SCRIPT 03: Alpha Diversity Analysis 
# DESCRIPTION: Calculates alpha diversity metrics and fits LMMs testing Fire Treatment.
# INPUT: Reads clean CSVs from 'data/' folder.
# OUTPUT: ANOVA table saved to 'processed_data/'
# ==============================================================================

# 1. SETUP & LIBRARIES
library(tidyverse)
library(fossil)
library(vegan)
library(lme4)
library(car)
library(broom.mixed)

# Set Working Directory to Project Folder
setwd("xxx")

# Create output folder if not exists
dir.create("processed_data", showWarnings = FALSE)

# 2. LOAD DATA (FROM CLEAN CSVs)
wide_myco_soil <- read_csv("data/input_soil_community.csv", show_col_types = FALSE)
site_metadata  <- read_csv("data/input_site_metadata.csv", show_col_types = FALSE)

# 3. CALCULATE ALPHA DIVERSITY
# Calculates Shannon, Simpson, Chao1, Observed, and Pielou's Evenness
alpha_div_soil <- wide_myco_soil %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    # chao1 requires a numeric vector of abundances
    Chao1 = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall")) > 0), # Sum of species with count > 0
    Pielou = ifelse(Observed > 1, Shannon / log(Observed), NA)
  ) %>%
  ungroup() %>%
  dplyr::select(Plot, Shannon, Simpson, Chao1, Observed, Pielou)

# 4. JOIN WITH METADATA
alpha_diversity_soil <- site_metadata %>%
  dplyr::select(Site, Plot, Fire_Treatment) %>% 
  left_join(alpha_div_soil, by = "Plot") %>%
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B"))) %>%
  filter(!is.na(Shannon)) # Remove samples with no diversity data

# 5. FIT LMMs & EXTRACT STATISTICS
metrics <- c("Shannon", "Simpson", "Chao1", "Pielou")

anova_table <- map_dfr(metrics, function(metric) {
  # Formula: Metric ~ Fire_Treatment + (1|Site)
  formula <- as.formula(paste0(metric, " ~ Fire_Treatment + (1|Site)"))
  model <- lmer(formula, data = alpha_diversity_soil)
  
  # ANOVA (Type II Wald F tests with Kenward-Roger df is standard, 
  # but 'test="F"' in car::Anova approximates this behavior)
  anov <- Anova(model, test = "F")
  
  # Extract Fixed Effects (Fire Treatment)
  tidy_model <- tidy(model, effects = "fixed") %>%
    filter(term == "Fire_TreatmentB")
  
  # Compile Stats
  tibble(
    Metric = metric,
    Factor = rownames(anov)[1],
    Estimate = tidy_model$estimate,
    Std_Error = tidy_model$std.error,
    DF_num = anov$Df,
    DF_denom = anov$Df.res,
    F = anov$F,
    P = anov$`Pr(>F)`[1]
  )
})

# 6. FORMAT & SAVE
anova_table_formatted_soil <- anova_table %>%
  dplyr::select(Metric, Factor, Estimate, Std_Error, DF_num, DF_denom, F, P) %>% 
  mutate(Type = 'Soil')

# Print to console
print(anova_table_formatted_soil)

# Save to CSV
write_csv(anova_table_formatted_soil, "processed_data/Alpha_Diversity_ANOVA.csv")

print("Alpha Diversity Analysis Complete. Results saved to 'processed_data/Alpha_Diversity_ANOVA.csv'")