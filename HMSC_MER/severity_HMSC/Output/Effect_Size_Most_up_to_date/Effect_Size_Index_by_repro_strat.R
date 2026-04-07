library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(ggpubr)
library(multcomp)
# load Beta1 and 2 from soil data HMSC
# load('HMSC_MER/results/mcmc_output.RData')

# read in soil comm data
myco_tax_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds")
myco_dat_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_RA_soil.rds")

# read in hyph comm data
myco_dat_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds")
myco_tax_hyph <- readRDS("Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds")

# All tax load
myco_tax_all <- bind_rows(myco_tax_soil, myco_tax_hyph) %>% distinct()


# load most up to date PA model
# load('HMSC_MER/results/Beta1.RData')
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")


# read in Carlos AMF data
# read in hyph comm data
# AMF_repo_type<-read.csv('Processed_data/Seq_dat/Hyph/AMF_Spore_Database_Volume.csv')

# First lets just use PA HMSC
Beta_estimates_PA <- Beta1$mean


Beta_estimates_PA <- Beta_estimates_PA %>%
  dplyr::select(OTU = Species, Beta = Severity)


# read in df extracted from Brundrett 2008 website
final_repo_strat <- readRDS("Processed_data/Seq_dat/datasets_external/repro_strat_df.Rds")

# Select OTUs associated with repo_strategy
soil_OTU_repo_strat <- myco_tax_soil %>%
  left_join(final_repo_strat) %>%
  # mutate(
  #   repo_strategy = case_when(
  #     phylum == "Glomeromycota" ~ "Glomeromycota",  # or "microscopic", "none", etc.
  #     TRUE ~ repo_strategy
  #   )
  # ) %>%
  filter(!is.na(repo_strategy)) %>%
  distinct(OTU, repo_strategy)

# select OTUs associated with phylym
soil_OTU_phy <- myco_tax_all %>%
  distinct(OTU, phylum)
###############################################
### PA MODEL###############
#########################

# Join and merge with Beta_estimates
Beta_PA_repo_strat_soil <- soil_OTU_repo_strat %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>%
  left_join(soil_OTU_phy)



# Plot
ggplot(Beta_PA_repo_strat_soil, aes(x = repo_strategy, y = Beta, fill = repo_strategy)) +
  geom_boxplot(
    alpha = 0.7, outlier.shape = NA, width = .3, position = position_dodge(width = 0.6) # wider dodge increases spacing
  ) +
  geom_jitter(
    alpha = 0.9, size = 1, position = position_dodge(width = 0.6) # must match boxplot dodge
  ) +
  labs(y = "Beta estimates", x = "repo_strategy", title = "Beta values by repo_strategy for hyph taxa PA model") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set2")

hist(Beta_PA_repo_strat_soil$Beta)

# PA model
model <- lm(Beta ~ repo_strategy, data = Beta_PA_repo_strat_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model)

smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_repo <- tibble(
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

anova_tbl_repo

anova_p <- Anova(model)$`Pr(>F)`[1] # First term in the model
# Get significance level
p_sig <- case_when(
  anova_p <= 0.0001 ~ "p <0.0001",
  anova_p <= 0.001 ~ "p <0.001",
  anova_p <= 0.01 ~ "p <0.01",
  anova_p <= 0.05 ~ "p <0.05",
  anova_p <= 0.1 ~ "+",
  TRUE ~ "ns"
)


# Create annotation label
model_label <- paste0(p_sig)
#### create axis labels
repo_strategy_labels_named <- c(
  "Hypogeous or Semi-hypogeous" = "Hypogeous/\nSemi-hypogeous",
  "Resupinate" = "Resupinate\n",
  "Epigeous" = "Epigeous\n",
  "Glomeromycota" = "Glomeromycota\n"
)
repo_strategy_levels <- names(repo_strategy_labels_named)


# Post-hoc pairwise comparisons for interaction
emm <- emmeans(model, ~repo_strategy)

# Tukey pairwise contrasts: Type within each repo_strategy
contrasts_type_within_repo_strategy <- emmeans(model, pairwise ~ repo_strategy, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within_repo_strategy$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.1) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ paste0("p=", round(p.value, 3))
    ),
    xmin = str_trim(str_extract(contrast, "(^[^-]+)")), # before "-"
    xmax = str_trim(ifelse(
      str_detect(contrast, "\\(.*\\)"),
      str_extract(contrast, "(?<=\\().*?(?=\\))"), # extract inside parentheses
      str_extract(contrast, "(?<=-)\\s*.*$") # fallback: extract after "-"
    ))
  )

post_hoc_repo <- contrasts_type_within_repo_strategy$contrasts %>%
  as.data.frame() %>%
  mutate(Variable = "Reproduction strategy")

post_hoc_repo

repo_strategy_max <- Beta_PA_repo_strat_soil %>%
  group_by(repo_strategy) %>%
  summarise(y.position = max(Beta, na.rm = TRUE) + 0.05)

sig_results <- sig_labels %>%
  left_join(repo_strategy_max, by = c("xmax" = "repo_strategy"))

# Get estimated marginal means and compact letter display
cld_repo_strategy <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = str_trim(.group)) # Clean whitespace

# Combine emmeans letters with y-position
label_df <- left_join(cld_repo_strategy, repo_strategy_max, by = "repo_strategy")

label_df



# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~repo_strategy) %>%
  as.data.frame() %>%
  mutate(repo_strategy = factor(repo_strategy, levels = repo_strategy_levels))


# Main boxplot
repo <- Beta_PA_repo_strat_soil %>%
  mutate(repo_strategy = factor(repo_strategy, levels = repo_strategy_levels)) %>%
  ggplot(aes(x = repo_strategy, y = Beta, fill = repo_strategy)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +

  # Violin and jitter
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +

  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x = repo_strategy, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x = repo_strategy, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +

  # Labels and formatting
  labs(
    y = expression(paste("Likelihood of occurrence post-fire (", beta[fire], ")")),
    x = "Fruiting body type", tag = "b)"
  ) +
  theme_classic() +
  scale_fill_grey(start = .9, end = 0.1) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18, vjust = .4, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.5, colour = "black")
  ) +
  scale_x_discrete(labels = repo_strategy_labels_named) +
  geom_text(
    data = label_df,
    aes(x = repo_strategy, y = y.position, label = .group),
    inherit.aes = FALSE, parse = TRUE,
    size = 8,
    fontface = "bold"
  ) +
  # geom_bracket(
  #   data = sig_results,
  #   aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
  #   inherit.aes = FALSE,
  #   tip.length = 0.01,
  #   label.size = 7,
  #   size = 1.1
  # )+
  annotate("text",
    x = 0.2, y = max(Beta_PA_repo_strat_soil$Beta, na.rm = TRUE) + .2,
    label = model_label, hjust = -0.3, vjust = 1, size = 8
  )
repo

ggsave("plots/Effect_Size_repo_strat_Plot.png", repo, width = 12, height = 8, dpi = 300)


save(post_hoc_repo, anova_tbl_repo, file = "HMSC_MER/Output/Processed_Data/repo_strat_summary.RDS")
