library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(ggpubr)
library(multcomp)

# read in soil comm data
myco_tax_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds")
myco_dat_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_RA_soil.rds")

# read in hyph comm data
myco_dat_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds")
myco_tax_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds")

# All tax load
myco_tax_all <- bind_rows(myco_tax_soil, myco_tax_hyph) %>% distinct()

# load most up to date PA model
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

# First lets just use PA HMSC
Beta_estimates_PA <- Beta1$mean

Beta_estimates_PA <- Beta_estimates_PA %>%
    dplyr::select(OTU = Species, Beta = Severity)

# ================================================================
# Select OTUs associated with Ecm_lineage (soil)
# ================================================================
soil_OTU_lineage <- myco_tax_soil %>%
    filter(!is.na(Ecm_lineage)) %>%
    distinct(OTU, Ecm_lineage)

# Join with Beta_estimates
Beta_PA_lineage_soil <- soil_OTU_lineage %>%
    inner_join(Beta_estimates_PA, by = "OTU")

# Quick look at lineage categories present
unique(Beta_PA_lineage_soil$Ecm_lineage)

# ================================================================
# PA MODEL: lm(Beta ~ Ecm_lineage)
# ================================================================
model_lineage <- lm(Beta ~ Ecm_lineage, data = Beta_PA_lineage_soil)

# Model performance summary
model_performance(model_lineage)
check_model(model_lineage)
summary(model_lineage)
Anova(model_lineage)

smry <- summary(model_lineage)
anov <- Anova(model_lineage, type = "II")

anova_tbl_lineage <- tibble(
    Term = rownames(anov)[1],
    sum_sq = anov$`Sum Sq`[1],
    sum_sq_resid = anov$`Sum Sq`[2],
    F = anov$`F value`[1],
    DF_num = anov$Df[1],
    DF_denom = anov$Df[2],
    R2 = smry$r.squared,
    adj_R2 = smry$adj.r.squared,
    P = anov$`Pr(>F)`[1]
)

anova_tbl_lineage

# ================================================================
# Post-hoc pairwise contrasts for Ecm_lineage
# ================================================================
lineage_max <- Beta_PA_lineage_soil %>%
    group_by(Ecm_lineage) %>%
    summarise(y.position = max(Beta, na.rm = TRUE) + 0.05)

anova_p <- Anova(model_lineage)$`Pr(>F)`[1]
p_sig <- case_when(
    anova_p <= 0.0001 ~ "p <0.0001",
    anova_p <= 0.001 ~ "p <0.001",
    anova_p <= 0.1 ~ paste0("p=", as.character(round(anova_p, 3))),
    TRUE ~ "ns"
)

# emmeans + compact letter display
emm_lineage <- emmeans(model_lineage, ~Ecm_lineage, adjust = "sidak")
cld_lineage <- cld(emm_lineage, adjust = "sidak", Letters = letters) %>%
    as.data.frame() %>%
    mutate(.group = str_trim(.group))
cld_lineage

# Combine emmeans letters with y-position
label_df_lineage <- left_join(cld_lineage, lineage_max, by = "Ecm_lineage")
label_df_lineage

# Pairwise contrasts
contrasts_lineage <- emmeans(model_lineage, pairwise ~ Ecm_lineage, adjust = "sidak")
post_hoc_lineage <- contrasts_lineage$contrasts %>%
    as.data.frame() %>%
    filter(p.value < 0.10) %>%
    mutate(
        p.label = case_when(
            p.value <= 0.001 ~ "***",
            p.value <= 0.01 ~ "**",
            p.value <= 0.05 ~ "*",
            p.value <= 0.1 ~ "+"
        )
    )
post_hoc_lineage

# ================================================================
# emmeans with 95% CIs for plotting
# ================================================================
emm_df_lineage <- emmeans(model_lineage, ~Ecm_lineage) %>%
    as.data.frame()

# Create overall model p-value annotation
model_label_lineage <- paste0(p_sig)

# ================================================================
# Final violin plot
# ================================================================
lineage_plot <- Beta_PA_lineage_soil %>%
    ggplot(aes(x = Ecm_lineage, y = Beta, fill = Ecm_lineage)) +
    geom_hline(
        yintercept = 0, linetype = "dashed", color = "grey60",
        linewidth = 1.5, alpha = 0.6
    ) +
    geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
    # Add estimated means and CIs
    geom_point(
        data = emm_df_lineage,
        aes(x = Ecm_lineage, y = emmean),
        shape = 21, size = 5, fill = "white", color = "black",
        stroke = 1.5, inherit.aes = FALSE
    ) +
    geom_errorbar(
        data = emm_df_lineage,
        aes(x = Ecm_lineage, ymin = lower.CL, ymax = upper.CL),
        width = 0.1, linewidth = 1.3, inherit.aes = FALSE
    ) +
    labs(
        y    = expression(paste("Likelihood of occurrence post-fire (", beta[fire], ")")),
        x    = "EcM lineage",
        tag  = "e)"
    ) +
    theme_classic() +
    theme(
        legend.position  = "none",
        axis.text.x      = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"),
        axis.text.y      = element_text(size = 14, face = "bold"),
        axis.title.x     = element_text(size = 18, vjust = .3, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        axis.line        = element_line(linewidth = 1.5, colour = "black")
    ) +
    scale_fill_grey(start = 0.9, end = 0.1) +
    # Compact letter display
    geom_text(
        data = label_df_lineage,
        aes(x = Ecm_lineage, y = y.position, label = .group),
        inherit.aes = FALSE,
        parse = TRUE,
        size = 8,
        fontface = "bold"
    ) +
    annotate(
        "text",
        x = 0.2,
        y = max(Beta_PA_lineage_soil$Beta + 0.15, na.rm = TRUE),
        label = model_label_lineage,
        hjust = -0.3, vjust = 1, size = 8
    )

lineage_plot

ggsave("plots/Effect_Size_Ecm_lineage_Plot.png", lineage_plot,
    width = 14, height = 8, dpi = 300
)

save(anova_tbl_lineage, post_hoc_lineage,
    file = "HMSC_MER/Output/Processed_Data/Ecm_lineage_summary.RDS"
)
