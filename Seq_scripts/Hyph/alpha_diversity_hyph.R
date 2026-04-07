
# Load libraries
library(tidyverse)
library(fossil)
library(vegan)
library(lme4)
library(car)
library(broom.mixed)

#-------------------------
# 1. Load data
#-------------------------
wide_myco_hyph<-readRDS('Processed_data/Seq_dat/Hyph/wide_myco.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS")

#-------------------------
# 2. Calculate alpha diversity
#-------------------------
alpha_div_hyph <- wide_myco_hyph %>%
  rowwise() %>%
  mutate(
    Shannon = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    Chao1 = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall"))),
    Pielou = ifelse(Observed > 1, Shannon / log(Observed), NA)
  ) %>%
  ungroup()

#-------------------------
# 3. Join with site metadata
#-------------------------
alpha_diversity_hyph <- All_Sites %>%
  select(Site,Plot,Fire_Treatment) %>% 
  left_join(alpha_div_hyph %>% 
              select(Plot,Shannon,Simpson,Chao1,Observed,Pielou), by = "Plot")

#-------------------------
# 4. Fit LMMs and extract ANOVA
#-------------------------
# Define metrics
metrics <- c("Shannon", "Simpson", "Chao1", "Pielou")


#set factor order
alpha_diversity_hyph <- alpha_diversity_hyph %>%
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")))

# Fit models and collect ANOVA tables
anova_table <- map_dfr(metrics, function(metric) {
  formula <- as.formula(paste0(metric, " ~ Fire_Treatment + (1|Site)"))
  model <- lmer(formula, data = alpha_diversity_hyph)
  
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
# 5. Final ANOVA result table
#-------------------------
anova_table_formatted_bag <- anova_table %>%
  select(Metric, Factor,Estimate, Std_Error, DF_num,DF_denom, F, P) %>% 
  mutate(Type='Ingrowth bag')

print(anova_table_formatted_bag)

save(anova_table_formatted_bag, file='HMSC_MER/Output/Processed_Data/bag_myco_alpha_metrics.RDS')



####################################################
library(ggpubr)   # for geom_bracket

#-------------------------
# 1. Reshape alpha_diversity_veg to long format
#-------------------------
alpha_long <- alpha_diversity_hyph %>%
  pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>% 
  select(Site,Plot,Fire_Treatment, Metric,Value) %>% 
  filter(!is.na(Value))

#-------------------------
# 2. Add p-values and labels from anova_table_formatted
#-------------------------
# Create a lookup table for brackets
# Prepare bracket data with y.position
bracket_df <- anova_table_formatted_bag %>%
  mutate(
    p_label = case_when(
      P <= 0.0001 ~ "p <0.0001",
      P <= 0.001 ~ "p <0.001",
      P <= 0.01  ~ "p <0.01",
      P <= 0.05  ~ paste0("p = ", round(P, 3)),
      P <= 0.1 ~ "+",
      TRUE ~ "ns"
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
hyph<-ggplot(alpha_long, aes(x = Fire_Treatment, y = Value)) +
  geom_boxplot(aes(fill = Fire_Treatment),outliers=FALSE,alpha = 0.6, linewidth=2) +
  scale_fill_manual(values = c(B = "#9F2121", U = "#5497B6"))+
  geom_jitter(width = 0.15, alpha = 0.5, size = 6) +
  #stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_text(
    data = bracket_df,
    aes(x = .5, y = y.position, label = p_label),  # position label near top-left
    hjust = 0, vjust = 0, 
    size = 15
  ) +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(
    x = "Fire Treatment",
    y = "Alpha Diversity Value"
  ) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text( hjust = 0.5, size = 35, face = "bold"),
        axis.text.y = element_text(size = 35, face = "bold"),
        axis.title.x = element_text(size = 50, face = "bold"),
        axis.title.y = element_text(size = 45, face = "bold"),
        strip.text = element_text(size = 45, face = "bold"),
        axis.line = element_line(linewidth = 3, colour = "black"),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 3),
  )

hyph  


#ggsave("plots/alpha_diversity_hyph.png", hyph, width = 30, height = 30, dpi = 300)


library(patchwork)

alpha<-soil+labs( tag="a)")+
  theme(plot.tag = element_text(size = 65, face = "bold"))+
  hyph+labs(tag = "b)")+
  theme(plot.tag = element_text(size = 65, face = "bold"))

alpha

ggsave("plots/alpha_diversity_both.png", alpha, width = 30, height = 15, dpi = 300)
