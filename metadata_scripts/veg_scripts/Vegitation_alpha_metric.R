library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(purrr)
library(APCalign)
library(fossil)
library(broom.mixed)


#read in Auscribe year 1 data already extract from veg_import_clean
auscribe_veg<-readRDS('Processed_data/metadata/veg/auscribe_veg_yr1.Rdata')

#read in all site locations
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% select(Site,Plot,Fire_Treatment)

#load for austraits
tax_resources <- load_taxonomic_resources()

# Extract distinct species names as a character vector
spp_names <- auscribe_veg %>%
  distinct(herbarium_determination) %>%
  filter(!is.na(herbarium_determination)) %>% 
  pull(herbarium_determination) %>%  # extract as a character vector
  iconv(from = "", to = "UTF-8", sub = "byte")  # Convert to valid UTF-8


# Align the species names using APC resources
aligned_spp <- align_taxa(spp_names, resources = tax_resources)

upd_spp_list<-update_taxonomy(aligned_spp, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
  select(original_name,suggested_name,aligned_name,genus:taxon_rank)

# Step 1: Join vegetation with mycorrhizal classifications, remove 'Unknown' types
auscribe_alpha<- auscribe_veg %>% 
  left_join(upd_spp_list %>% select(original_name,suggested_name,genus,family), by =c('herbarium_determination'='original_name'))%>% 
  select(-plot, Plot = plot_join) 


#-------------------------
# 1. Prepare species-by-plot matrix (wide format)
#-------------------------
veg_comm <- auscribe_alpha %>%
  filter(!is.na(suggested_name)) %>% 
  count(Plot, herbarium_determination) %>%
  pivot_wider(names_from = herbarium_determination, values_from = n, values_fill = 0)

#-------------------------
# 2. Calculate alpha diversity
#-------------------------
alpha_div_veg <- veg_comm %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(-Plot), index = "shannon"),
    Simpson = diversity(c_across(-Plot), index = "simpson"),
    Chao1   = chao1(c_across(-Plot)),
    Observed = sum(c_across(-Plot) > 0),
    Pielou  = ifelse(Observed > 1, Shannon / log(Observed), NA)
  ) %>%
  ungroup()

#-------------------------
# 3. Join with site metadata
#-------------------------
alpha_diversity_veg <- All_Sites %>%
  select(Site, Plot, Fire_Treatment) %>%
  left_join(alpha_div_veg, by = "Plot") 


Plots_not_sampled<-alpha_diversity_veg%>% 
  filter(is.na(Chao1)) %>% distinct(Site,Plot,Fire_Treatment) %>%
  group_by(Site,Fire_Treatment) %>% 
  summarise(n()) %>%
  rename(Count = `n()`) %>%
  pivot_wider(
    names_from = Fire_Treatment,
    values_from = Count,
    values_fill = 0
  ) %>%
  rename(Burnt = B, Unburnt = U)%>% 
  mutate(Type='Vegetation')

Plots_not_sampled 
save(Plots_not_sampled, file='Processed_data/metadata/veg/Plots_w_no_veg_data.rds')

alpha_diversity_veg<-alpha_diversity_veg%>% 
  filter(!is.na(Chao1)) 

#set factor order
alpha_diversity_veg <- alpha_diversity_veg %>%
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")))

Veg_richness_plot<-alpha_diversity_veg %>% 
  select(Site:Fire_Treatment,Shannon:Pielou)

save(Veg_richness_plot, file='HMSC_MER/Output/Processed_Data/Vegitation_alpha_metrics_raw.RDS')
#-------------------------
# 4. Fit LMMs and extract ANOVA
#-------------------------


# Define diversity metrics
metrics <- c("Shannon", "Simpson", "Chao1", "Pielou")

anova_table_all_veg <- map_dfr(metrics, function(metric) {
  
  # 1. Fit model
  formula <- as.formula(paste0(metric, " ~ Fire_Treatment + (1|Site)"))
  model <- lmer(formula, data = alpha_diversity_veg)
  
  # 2. ANOVA
  anov <- Anova(model, test = "F")
  
  # 3. R²
  #smry<-summary(model)
  
  # 4. Estimate and Std. Error for Fire_Treatment effect
  tidy_model <- tidy(model, effects = "fixed") %>%
    filter(term == "Fire_TreatmentB")  # or partial match if multiple levels
  
  # 5. Combine into summary table
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

#-------------------------
# 5. Final formatted ANOVA table
#-------------------------
# Print result
print(anova_table_all_veg)

save(anova_table_all_veg, file='HMSC_MER/Output/Processed_Data/Vegitation_alpha_metrics.RDS')
####################################################
library(ggpubr)   # for geom_bracket

#-------------------------
# 1. Reshape alpha_diversity_veg to long format
#-------------------------
alpha_long <- alpha_diversity_veg %>%
  pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>% 
  select(Site,Plot,Fire_Treatment, Metric,Value)

#-------------------------
# 2. Add p-values and labels from anova_table_formatted
#-------------------------
# Create a lookup table for brackets
# Prepare bracket data with y.position
bracket_df <- anova_table_all_veg %>%
  filter(Factor == "Fire_Treatment") %>%
  mutate(
    p_label = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    xmin = 1,
    xmax = 2
  ) %>%
  left_join(
    alpha_long %>%
      group_by(Metric) %>%
      summarise(y.position = max(Value, na.rm = TRUE) * .95),
    by = "Metric"
  ) %>%
  select(Metric, xmin, xmax, y.position, p_label)

# Now use geom_bracket with data = bracket_df
ggplot(alpha_long, aes(x = Fire_Treatment, y = Value)) +
  geom_boxplot(aes(fill = Fire_Treatment), outliers=FALSE,alpha = 0.6, linewidth=1) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 3) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_text(
    data = bracket_df,
    aes(x = .5, y = y.position, label = p_label),  # position label near top-left
    hjust = 0, vjust = 0,
    size = 9
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    x = "Fire Treatment",
    y = "Alpha Diversity Value",
    title = "Alpha Diversity Metrics by Fire Treatment for all vegitation"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        strip.text = element_text(size = 27, face = "bold"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 3),
        axis.line = element_line(linewidth = 2, colour = "black")
        )


unique(alpha_long$Site)
