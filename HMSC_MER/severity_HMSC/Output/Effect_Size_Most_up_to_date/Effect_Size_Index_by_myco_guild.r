library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(multcomp)

# ─────────────────────────────────────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────────────────────────────────────

# Read in soil community data
myco_tax_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds")
myco_dat_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_RA_soil.rds")

# Read in hyphal community data
myco_dat_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds")
myco_tax_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds")

# Merge taxonomic tables
myco_tax_all <- bind_rows(myco_tax_soil, myco_tax_hyph) %>% distinct()

# Load most up-to-date presence-absence HMSC model output
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

# ─────────────────────────────────────────────────────────────────────────────
# Extract beta-fire estimates
# ─────────────────────────────────────────────────────────────────────────────

Beta_estimates_PA <- Beta1$mean

Beta_estimates_PA <- Beta_estimates_PA %>%
    dplyr::select(OTU = Species, Beta = Severity)

# ─────────────────────────────────────────────────────────────────────────────
# Prepare guild-level data
# NOTE: adjust the column name 'guild' below if your taxonomy table uses a
# different variable name (e.g., 'mycorrhizal_type', 'myco_guild', etc.)
# Expected levels: "EcM", "AM", "ErM", "OM" (or similar)
# ─────────────────────────────────────────────────────────────────────────────

soil_OTU_guild <- myco_tax_soil %>%
    distinct(OTU, guild2, phylum, family, genus)

# Join beta estimates with guild assignments
Beta_PA_guild_soil <- soil_OTU_guild %>%
    left_join(Beta_estimates_PA, by = "OTU") %>%
    filter(!is.na(guild2), !is.na(Beta)) %>%
    mutate(guild = as.factor(guild2))

# ─────────────────────────────────────────────────────────────────────────────
# Remove guilds with fewer than 2 OTUs
# ─────────────────────────────────────────────────────────────────────────────

guild_otu_counts <- Beta_PA_guild_soil %>%
    group_by(guild2) %>%
    summarise(n_OTU = n_distinct(OTU), .groups = "drop")

guilds_keep <- guild_otu_counts %>%
    filter(n_OTU >= 2) %>%
    pull(guild2)

message(sprintf(
    "Retaining %d of %d guilds (≥2 OTUs): %s",
    length(guilds_keep),
    nrow(guild_otu_counts),
    paste(guilds_keep, collapse = ", ")
))

Beta_PA_guild_soil <- Beta_PA_guild_soil %>%
    filter(guild2 %in% guilds_keep) %>%
    mutate(guild2 = droplevels(as.factor(guild2)))

# ─────────────────────────────────────────────────────────────────────────────
# Linear model: beta-fire ~ guild
# ─────────────────────────────────────────────────────────────────────────────

model_guild <- lm(Beta ~ guild2, data = Beta_PA_guild_soil)

# Check model assumptions
check_model(model_guild)

# Summarise model
smry <- summary(model_guild)
anov <- Anova(model_guild, type = "II")

# ─────────────────────────────────────────────────────────────────────────────
# ANOVA summary table
# ─────────────────────────────────────────────────────────────────────────────

anova_tbl_guild <- tibble(
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

anova_tbl_guild

# ─────────────────────────────────────────────────────────────────────────────
# Coefficient summary table
# ─────────────────────────────────────────────────────────────────────────────

smry_tbl_guild <- smry$coefficients %>%
    as.data.frame() %>%
    rename(p.value = `Pr(>|t|)`) %>%
    rownames_to_column(var = "Coefficient") %>%
    mutate(
        Coefficient = str_remove(Coefficient, "guild"),
        Term = "Guild",
        p.label = case_when(
            p.value <= 0.001 ~ "***",
            p.value <= 0.01 ~ "**",
            p.value <= 0.05 ~ "*",
            p.value <= 0.1 ~ "+",
            TRUE ~ ""
        )
    )

smry_tbl_guild

# ─────────────────────────────────────────────────────────────────────────────
# Post-hoc pairwise comparisons (Sidak-adjusted)
# ─────────────────────────────────────────────────────────────────────────────

emm_guild <- emmeans(model_guild, ~guild2, adjust = "sidak")

cld_guild <- cld(emm_guild, adjust = "sidak", Letters = letters) %>%
    as.data.frame() %>%
    mutate(.group = str_trim(.group))

contrasts_guild <- emmeans(model_guild, pairwise ~ guild2, adjust = "sidak")

post_hoc_guild <- contrasts_guild$contrasts %>%
    as.data.frame() %>%
    mutate(Variable = "Guild")

post_hoc_guild

# ─────────────────────────────────────────────────────────────────────────────
# Label positions for plot
# ─────────────────────────────────────────────────────────────────────────────

guild_max <- Beta_PA_guild_soil %>%
    group_by(guild2) %>%
    summarise(y.position = max(Beta, na.rm = TRUE) + 0.05, .groups = "drop")

label_df_guild <- left_join(cld_guild, guild_max, by = "guild2")

# ANOVA p-value label
P_guild <- Anova(model_guild)$`Pr(>F)`[1]
model_label_guild <- case_when(
    P_guild <= 0.0001 ~ "p <0.0001",
    P_guild <= 0.001 ~ "p <0.001",
    P_guild <= 0.01 ~ "p <0.01",
    P_guild <= 0.05 ~ "p <0.05",
    P_guild <= 0.1 ~ "+",
    TRUE ~ "ns"
)

# Estimated marginal means with 95% CIs
emm_df_guild <- emmeans(model_guild, ~guild2) %>%
    as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────
# Plot: beta-fire by mycorrhizal guild
# ─────────────────────────────────────────────────────────────────────────────

guild_plot <- Beta_PA_guild_soil %>%
    ggplot(aes(x = guild2, y = Beta, fill = guild2)) +
    geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "grey60", linewidth = 1.5, alpha = 0.6
    ) +
    geom_violin(
        alpha = 0.7, width = 0.5,
        linewidth = 1, color = "black"
    ) +
    geom_point(
        data = emm_df_guild,
        aes(x = guild2, y = emmean),
        shape = 21, size = 5, fill = "white",
        color = "black", stroke = 2,
        inherit.aes = FALSE
    ) +
    geom_errorbar(
        data = emm_df_guild,
        aes(x = guild2, ymin = lower.CL, ymax = upper.CL),
        width = 0.1, linewidth = 1.3,
        inherit.aes = FALSE
    ) +
    geom_text(
        data = label_df_guild,
        aes(x = guild2, y = y.position, label = .group),
        inherit.aes = FALSE,
        size = 8, fontface = "bold"
    ) +
    annotate(
        "text",
        x = 0.2,
        y = max(Beta_PA_guild_soil$Beta, na.rm = TRUE) + 0.2,
        label = model_label_guild,
        hjust = -0.3, vjust = 1, size = 8
    ) +
    labs(
        y   = expression(paste("Likelihood of occurrence post-fire (", beta[fire], ")")),
        x   = "Mycorrhizal Guild",
        tag = "a)"
    ) +
    scale_fill_grey(start = 0.9, end = 0.1) +
    theme_classic() +
    theme(
        legend.position  = "none",
        axis.text.x      = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y      = element_text(size = 14, face = "bold"),
        axis.title.x     = element_text(size = 18, face = "bold"),
        axis.title.y     = element_text(size = 18, face = "bold"),
        axis.line        = element_line(linewidth = 1.5, colour = "black")
    )

guild_plot

# ─────────────────────────────────────────────────────────────────────────────
# Save outputs
# ─────────────────────────────────────────────────────────────────────────────

save(
    post_hoc_guild,
    anova_tbl_guild,
    smry_tbl_guild,
    file = "HMSC_MER/Output/Processed_Data/Guild_summary.RDS"
)

ggsave(
    "plots/Effect_Size_guild_Plot.png",
    guild_plot,
    width = 10, height = 8, dpi = 300
)
