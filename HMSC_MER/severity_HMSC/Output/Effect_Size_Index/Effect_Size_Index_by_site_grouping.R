# ----------------------------
# Load libraries
# ----------------------------
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(ggpubr)
library(emmeans)
library(ggsignif)
library(ggpp)
library(ggforce)
library(readxl)

# ----------------------------
# Load data
# ----------------------------
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
myco_dat_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

TSF <- readRDS('Processed_data/metadata/fire/Fire_regime_MER.Rdata') %>%
  distinct(Site, severity, TSF_years)

MER_Site_description<-read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx",
           sheet = "SP suggestions (highlighted)")

colnames(MER_Site_description)
Mer_veg <- MER_Site_description %>%
  left_join(TSF) %>% 
  dplyr::select(Site, Site_nickname, Last_fire_severity, severity,TSF_years,
         Veg_type,Veg_type_alt,Veg_type_alt_2, Fire_interval, Min_TFI = `Min TFI`,
         Simple_fire_response_dominants, Fire_response) %>% 
  mutate(Min_TFI=as.numeric(Min_TFI))

load('HMSC_MER/results/Beta1.RData')
load('HMSC_MER/results/Beta1_host_freq.RData')


Beta_estimates_PA <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Severity)

soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site, Plot, OTU) %>%
  left_join(Mer_veg)

Beta_PA_Veg_soil <- soil_Site_OTU %>%
  left_join(Beta_estimates_PA, by = "OTU") %>%
  mutate(Type = "Soil")

# ----------------------------
# Variables to analyze
# ----------------------------
vars <- c('severity', 'Veg_type','Veg_type_alt', 'Fire_interval','TSF_years',
          'Min_TFI', 'Simple_fire_response_dominants', 'Fire_response')

# Output directory for plots
output_dir <- "plots/fire_regime"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------
# Function to analyze each variable
# ----------------------------
analyze_var <- function(var) {
  message("Running analysis for: ", var)
  
  # Fit model dynamically
  formula <- as.formula(paste("Beta ~", var))
  model <- lm(formula, data = Beta_PA_Veg_soil)
  
  # Summaries
  smry <- summary(model)
  anov <- Anova(model, test.statistic="F")
  
  # Build ANOVA table
  anova_tbl <- tibble(
    Variable = var,
    Term = rownames(anov)[1],
    sum_sq = anov$`Sum Sq`[1],
    sum_sq_resid = anov$`Sum Sq`[2],
    F_value = anov$`F value`[1],
    DF_num = anov$Df[1],
    DF_denom = anov$Df[2],
    R2 = smry$r.squared,
    adj_R2 = smry$adj.r.squared,
    P = anov$`Pr(>F)`[1]
  )
  
  # Check if variable is categorical or numeric
  is_factor <- is.factor(Beta_PA_Veg_soil[[var]]) || is.character(Beta_PA_Veg_soil[[var]])
  
  emm_df <- NULL
  sig_labels <- NULL
  
  # ----------------------------
  # Post-hoc comparisons for categorical variables
  # ----------------------------
  if (is_factor) {
    emm <- emmeans(model, specs = as.formula(paste("~", var)))
    emm_df <- as.data.frame(emm)
    
    # Run pairwise contrasts only if ANOVA significant
    if (anova_tbl$P < 0.05) {
      contrasts <- contrast(emm, method = "pairwise", adjust = "sidak") %>%
        as.data.frame()
      
      sig_labels <- contrasts %>%
        filter(p.value < 0.1) %>%
        mutate(
          p.label = case_when(
            p.value <= 0.001 ~ "***",
            p.value <= 0.01 ~ "**",
            p.value <= 0.05 ~ "*",
            p.value <= 0.06 ~ "+"
          ),
          xmin = str_trim(str_extract(contrast, "(^[^-]+)")),
          xmax = str_trim(ifelse(
            str_detect(contrast, "\\(.*\\)"),
            str_extract(contrast, "(?<=\\().*?(?=\\))"),
            str_extract(contrast, "(?<=-)\\s*.*$")
          ))
        )
      
      # Position brackets above the data
      group_max <- Beta_PA_Veg_soil %>%
        group_by(.data[[var]]) %>%
        summarise(ymax = max(Beta, na.rm = TRUE) + 0.3, .groups = "drop")
      
      sig_labels <- sig_labels %>%
        left_join(group_max, by = c("xmin" = var)) %>%
        rename(y.position = ymax)
    }
  }
  
  # ----------------------------
  # ANOVA p-value label
  # ----------------------------
  model_label <- paste0("p=", signif(anova_tbl$P, 3))
  model_R <- paste0("R^2:", signif(smry$r.squared, 2))
  
  # ----------------------------
  # Create plots
  # ----------------------------
  if (is_factor) {
    plot <- Beta_PA_Veg_soil %>%
      ggplot(aes_string(x = var, y = "Beta", fill = var)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 2, alpha = 0.6) +
      geom_violin(alpha = 0.7, width = 0.4, linewidth = 1.5) +
      geom_jitter(alpha=.5,width=.1,size=3,aes(color=Site))+
      geom_point(data = emm_df, aes_string(x = var, y = "emmean"),
                 shape = 21, size = 6, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
      geom_errorbar(data = emm_df, aes_string(x = var, ymin = "lower.CL", ymax = "upper.CL"),
                    width = 0.15, linewidth = 1.5, inherit.aes = FALSE) +
      annotate("text", x = 1, y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE) + 0.5,
               label = model_label, size = 8, fontface = "bold", hjust = 0) +
      {if (!is.null(sig_labels)) geom_bracket(
        data = sig_labels,
        aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
        inherit.aes = FALSE,
        tip.length = 0.01,
        label.size = 6,
        size = 1.1
      )} +
      labs(y = "Effect size (β_fire)", x = var) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 30, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")
      )
  } else {
    # Scatter + regression line for numeric predictors
    plot <- Beta_PA_Veg_soil %>%
      ggplot(aes_string(x = var, y = "Beta")) +
      geom_point(aes(color = Site), size = 3, alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      annotate("text", x = min(Beta_PA_Veg_soil[[var]], na.rm = TRUE),
               y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE) + 0.5,
               label = model_label, size = 8, fontface = "bold", hjust = 0) +
      annotate("text", x = min(Beta_PA_Veg_soil[[var]], na.rm = TRUE),
               y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE)+.25,
               label = model_R, size = 8, fontface = "bold", hjust = 0) +
      labs(y = "Effect size (β_fire)", x = var) +
      theme_classic() +
      theme(
        axis.text.x = element_text(hjust = 0.5, size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 16, face = "bold")
      )
  }
  
  # Save results & plots
  saveRDS(anova_tbl, paste0("HMSC_MER/Output/Processed_Data/", var, "_anova_summary.rds"))
  ggsave(paste0(output_dir, "/Effect_Size_", var, "_Plot.png"), plot, width = 10, height = 6, dpi = 300)
  
  return(list(anova = anova_tbl, plot = plot))
}

# ----------------------------
# Run analysis for all variables
# ----------------------------
results_list <- lapply(vars, analyze_var)

# Combine all ANOVA results into one summary table
all_anova <- bind_rows(lapply(results_list, `[[`, "anova"))

all_anova
# Save combined ANOVA table
#write_csv(all_anova, "HMSC_MER/Output/Processed_Data/All_ANOVA_Summary.csv")

# Show all plots in R
for (res in results_list) {
  print(res$plot)
}
