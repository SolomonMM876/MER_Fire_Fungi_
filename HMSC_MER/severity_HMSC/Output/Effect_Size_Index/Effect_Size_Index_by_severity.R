# ----------------------------
# Load libraries
# ----------------------------
library(tidyverse)
library(ggplot2)
library(car)
library(readxl)

# ----------------------------
# Load data
# ----------------------------
myco_tax_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds")
myco_dat_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_RA_soil.rds")

TSF <- readRDS("Processed_data/metadata/fire/Fire_regime_MER.Rdata") %>%
  distinct(Site, severity, TSF_years)

MER_Site_description <- read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx",
                                   sheet = "SP suggestions (highlighted)")

Mer_veg <- MER_Site_description %>%
  left_join(TSF) %>%
  select(Site, Site_nickname, Last_fire_severity, severity, TSF_years,
         Veg_type, Fire_interval, Min_TFI = `Min TFI`,
         Simple_fire_response_dominants , Fire_response)

# ----------------------------
# Load HMSC beta estimates
# ----------------------------
load("HMSC_MER/results/Beta1.RData")
Beta_estimates_PA <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Severity)

# ----------------------------
# Combine OTU + Beta + Taxonomy + Veg metadata
# ----------------------------
soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site, Plot, OTU) %>%
  left_join(Mer_veg)

Beta_PA_Veg_soil <- soil_Site_OTU %>%
  left_join(Beta_estimates_PA, by = "OTU") %>%
  left_join(myco_tax_soil, by = "OTU") %>%
  filter(grepl("Mycorrhizal", guild2)) %>%  # keep all mycorrhizal OTUs
  mutate(Type = "Soil") %>%
  drop_na(Beta, severity)

# ----------------------------
# Fit model for all mycorrhizal OTUs
# ----------------------------
model_all <- lm(Beta ~ severity, data = Beta_PA_Veg_soil)

# Get ANOVA table
anova_all <- Anova(model_all,  type = 3, test.statistic = "F")
summary_model <- summary(model_all)

# Create tidy ANOVA results table
anova_table <- tibble(
  Term = rownames(anova_all)[1],
  Sum_Sq = anova_all$`Sum Sq`[1],
  DF_num = anova_all$Df[1],
  DF_denom = anova_all$Df[2],
  F_value = anova_all$`F value`[1],
  R2 = summary_model$r.squared,
  Adj_R2 = summary_model$adj.r.squared,
  P_value = anova_all$`Pr(>F)`[1]
)

print(anova_table)

# ----------------------------
# Create single plot
# ----------------------------http://127.0.0.1:46287/graphics/plot_zoom_png?width=1184&height=861
label_p <- paste0("p = ", signif(anova_table$P_value, 2))
label_r <- paste0("R²: ", signif(anova_table$R2, 2))

final_plot <- ggplot(Beta_PA_Veg_soil, aes(x = severity, y = Beta)) +
  geom_point(size = 3, alpha = 0.5, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "purple4", linewidth = 2.5) +
  geom_hline(yintercept = 0, color='grey')+
  annotate("text",
           x = min(Beta_PA_Veg_soil$severity, na.rm = TRUE),
           y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE) + 0.2,
           label = label_p, size = 5, fontface = "bold", hjust = 0) +
  annotate("text",
           x = min(Beta_PA_Veg_soil$severity, na.rm = TRUE),
           y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE) + 0.15,
           label = label_r, size = 5, fontface = "bold", hjust = 0) +
  labs(
    x = "Average fire severity of burnt plots",
    y = expression(paste("Effect size (", beta[fire], ")"))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# Display the plot
print(final_plot)

# ----------------------------
# Save outputs
# ----------------------------
#saveRDS(anova_table, "HMSC_MER/Output/Processed_Data/fire_severity_ANOVA_output_all.RDS")
ggsave("plots/Fire_severity_Effect_Size_All_Mycorrhizal.png", final_plot,
       width = 10, height = 8, dpi = 300)
















# Create a dataset with group labels for analysis
Beta_PA_Veg_soil_grouped <- Beta_PA_Veg_soil %>%
  mutate(Group = case_when(
    guild2 == "EctoMycorrhizal" ~ "EctoMycorrhizal",
    guild2 == "Arbuscular Mycorrhizal" ~ "Arbuscular Mycorrhizal",
    TRUE ~ "Other"
  )) %>%
  filter(Group %in% c("EctoMycorrhizal", "Arbuscular Mycorrhizal")) %>% 
  drop_na(Beta, severity)


# Fit linear model with interaction
interaction_model <- lm(Beta ~ severity * Group, data = Beta_PA_Veg_soil_grouped)

# Get ANOVA table for interaction significance
anova_interaction <- Anova(interaction_model, type = "II")  # Type III to test main effects properly
anova_interaction

summary(interaction_model)




