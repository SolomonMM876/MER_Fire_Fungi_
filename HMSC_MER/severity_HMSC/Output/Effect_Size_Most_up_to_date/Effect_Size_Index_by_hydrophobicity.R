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
# Select OTUs associated with hydrophobicity (soil)
# Hydrophobicity is derived from exploration type per Agerer et al. 2012:
#   Hydrophilic: contact, short-distance_coarse, short-distance_delicate,
#                medium-distance_smooth
#   Hydrophobic: medium-distance_fringe, long-distance, mat
# ================================================================
soil_OTU_hydro <- myco_tax_soil %>%
    filter(!is.na(hydrophobicity)) %>%
    distinct(OTU, hydrophobicity)

# Join with Beta_estimates
Beta_PA_hydro_soil <- soil_OTU_hydro %>%
    inner_join(Beta_estimates_PA, by = "OTU") %>%
    mutate(hydrophobicity = factor(hydrophobicity,
        levels = c("hydrophilic", "hydrophobic")
    ))

# Check counts per group
Beta_PA_hydro_soil %>% count(hydrophobicity)

# ================================================================
# PA MODEL: lm(Beta ~ hydrophobicity)
# Binary predictor — equivalent to Welch t-test via lm
# ================================================================
model_hydro <- lm(Beta ~ hydrophobicity, data = Beta_PA_hydro_soil)

# Model performance summary
model_performance(model_hydro)
check_model(model_hydro)
summary(model_hydro)
Anova(model_hydro)

smry <- summary(model_hydro)
anov <- Anova(model_hydro, type = "II")

anova_tbl_hydro <- tibble(
    Term         = rownames(anov)[1],
    sum_sq       = anov$`Sum Sq`[1],
    sum_sq_resid = anov$`Sum Sq`[2],
    F            = anov$`F value`[1],
    DF_num       = anov$Df[1],
    DF_denom     = anov$Df[2],
    R2           = smry$r.squared,
    adj_R2       = smry$adj.r.squared,
    P            = anov$`Pr(>F)`[1]
)

anova_tbl_hydro

# ================================================================
# emmeans & pairwise contrast (hydrophilic vs hydrophobic)
# ================================================================
hydro_max <- Beta_PA_hydro_soil %>%
    group_by(hydrophobicity) %>%
    summarise(y.position = max(Beta, na.rm = TRUE) + 0.05)

anova_p <- Anova(model_hydro)$`Pr(>F)`[1]
p_sig <- case_when(
    anova_p <= 0.0001 ~ "p <0.0001",
    anova_p <= 0.001 ~ "p <0.001",
    anova_p <= 0.1 ~ paste0("p=", as.character(round(anova_p, 3))),
    TRUE ~ "ns"
)

# emmeans + compact letter display
emm_hydro <- emmeans(model_hydro, ~hydrophobicity, adjust = "sidak")
cld_hydro <- cld(emm_hydro, adjust = "sidak", Letters = letters) %>%
    as.data.frame() %>%
    mutate(.group = str_trim(.group))
cld_hydro

# Combine emmeans letters with y-position
label_df_hydro <- left_join(cld_hydro, hydro_max, by = "hydrophobicity")
label_df_hydro

# Pairwise contrast (only one pair here)
contrasts_hydro <- emmeans(model_hydro, pairwise ~ hydrophobicity, adjust = "sidak")
post_hoc_hydro <- contrasts_hydro$contrasts %>%
    as.data.frame() %>%
    mutate(
        p.label = case_when(
            p.value <= 0.001 ~ "***",
            p.value <= 0.01 ~ "**",
            p.value <= 0.05 ~ "*",
            p.value <= 0.1 ~ "+",
            TRUE ~ "ns"
        )
    )
post_hoc_hydro

# ================================================================
# emmeans with 95% CIs for plotting
# ================================================================
emm_df_hydro <- emmeans(model_hydro, ~hydrophobicity) %>%
    as.data.frame() %>%
    mutate(hydrophobicity = factor(hydrophobicity,
        levels = c("hydrophilic", "hydrophobic")
    ))

model_label_hydro <- paste0(p_sig)

# ================================================================
# Final violin plot
# (2-group binary comparison, cleaner layout)
# ================================================================
# Readable axis labels
hydro_labels <- c(
    "hydrophilic" = "Hydrophilic\n(contact, short-distance,\nmed-smooth)",
    "hydrophobic" = "Hydrophobic\n(med-fringe,\nlong-distance, mat)"
)

hydro_plot <- Beta_PA_hydro_soil %>%
    ggplot(aes(x = hydrophobicity, y = Beta, fill = hydrophobicity)) +
    geom_hline(
        yintercept = 0, linetype = "dashed", color = "grey60",
        linewidth = 1.5, alpha = 0.6
    ) +
    geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
    # Estimated means and 95% CIs
    geom_point(
        data = emm_df_hydro,
        aes(x = hydrophobicity, y = emmean),
        shape = 21, size = 5, fill = "white", color = "black",
        stroke = 1.5, inherit.aes = FALSE
    ) +
    geom_errorbar(
        data = emm_df_hydro,
        aes(x = hydrophobicity, ymin = lower.CL, ymax = upper.CL),
        width = 0.08, linewidth = 1.3, inherit.aes = FALSE
    ) +
    labs(
        y   = expression(paste("Likelihood of occurrence post-fire (", beta[fire], ")")),
        x   = "Hydrophobicity (Agerer et al. 2012)",
        tag = "f)"
    ) +
    theme_classic() +
    theme(
        legend.position = "none",
        axis.text.x     = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"),
        axis.text.y     = element_text(size = 14, face = "bold"),
        axis.title.x    = element_text(size = 18, vjust = .3, face = "bold"),
        axis.title.y    = element_text(size = 18, face = "bold"),
        axis.line       = element_line(linewidth = 1.5, colour = "black")
    ) +
    scale_x_discrete(labels = hydro_labels) +
    scale_fill_grey(start = 0.85, end = 0.2) +
    # Compact letter display
    geom_text(
        data = label_df_hydro,
        aes(x = hydrophobicity, y = y.position, label = .group),
        inherit.aes = FALSE,
        parse = TRUE,
        size = 8,
        fontface = "bold"
    ) +
    annotate(
        "text",
        x = 0.55,
        y = max(Beta_PA_hydro_soil$Beta + 0.15, na.rm = TRUE),
        label = model_label_hydro,
        hjust = 0, vjust = 1, size = 8
    )

hydro_plot

ggsave("plots/Effect_Size_hydrophobicity_Plot.png", hydro_plot,
    width = 10, height = 8, dpi = 300
)

save(anova_tbl_hydro, post_hoc_hydro,
    file = "HMSC_MER/Output/Processed_Data/hydrophobicity_summary.RDS"
)
