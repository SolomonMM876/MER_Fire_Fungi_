library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
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
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

# First lets just use PA HMSC
Beta_estimates_PA <- Beta1$mean


Beta_estimates_PA <- Beta_estimates_PA %>%
  dplyr::select(OTU = Species, Beta = Severity)

# Select OTUs associated with glom
soil_OTU_AM <- myco_tax_soil %>%
  distinct(OTU, phylum, family, genus) %>%
  filter(phylum == "Glomeromycota")

library(tibble)

# table from Weber 2019
amf_traits <- tribble(
  ~Guild, ~Intraradical_hyphae, ~Extraradical_hyphae, ~Families, ~Citations,
  "Rhizophilic", "High", "Low", "Glomeraceae", "1,2,3,4,5",
  "Rhizophilic", "High", "Low", "Claroideoglomeraceae", "2",
  "Rhizophilic", "High", "Low", "Paraglomeraceae", "b",
  "Edaphophilic", "Low", "High", "Gigasporaceae", "1,2,5",
  "Edaphophilic", "Low", "High", "Diversisporaceae", "2,5",
  "Ancestral", "Low", "Low", "Archaeosporaceae", "b",
  "Ancestral", "Low", "Low", "Ambisporaceae", "b",
  "Ancestral", "Low", "Low", "Pacisporaceae", "5",
  "Ancestral", "Low", "Low", "Acaulosporaceae", "1,2"
)

Soil_OTU_rhizo <- soil_OTU_AM %>%
  left_join(amf_traits, by = c("family" = "Families")) %>%
  filter(!is.na(family))

###############################################
### PA MODEL###############
#########################

# Join and merge with Beta_estimates
Soil_OTU_rhizo_B <- Soil_OTU_rhizo %>%
  left_join(Beta_estimates_PA, by = "OTU") %>%
  mutate(family = as.factor(family))


### Guild####




ggplot(Soil_OTU_rhizo_B, aes(x = Guild, y = Beta, fill = Guild)) +
  geom_boxplot(
    alpha = 0.7, outlier.shape = NA, width = .3, # wider dodge increases spacing
  ) +
  geom_jitter(
    alpha = 0.9, size = 1, width = 0.05 # must match boxplot dodge
  ) +
  labs(y = "Effect size (β_fire)", x = "Guild") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


# PA model
model <- lm(Beta ~ Guild, data = Soil_OTU_rhizo_B)

# Check model assumptions visually
check_model(model)
summary(model)
Anova(model, test = "F")

smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_Guild <- tibble(
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

anova_tbl_Guild

smry_tbl_Guild <- smry$coefficients %>%
  as.data.frame() %>%
  rename(p.value = `Pr(>|t|)`) %>%
  rownames_to_column(var = "Coefficent") %>%
  mutate(
    Coefficent = str_remove(Coefficent, "Guild"),
    Term = "Guild ",
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.1 ~ "+",
      TRUE ~ ""
    )
  )
smry_tbl_Guild

# Get label positions (max y for each Guild )
Guild_max <- Soil_OTU_rhizo_B %>%
  group_by(Guild) %>%
  summarise(y.position = max(Beta, na.rm = TRUE) + 0.02, .groups = "drop")

# Combine with significance
smry_intrarad <- left_join(Guild_max, smry_tbl_Guild, by = c("Guild" = "Coefficent")) %>%
  filter(p.label != "") # Keep only significant labels
smry_intrarad

# Get estimated marginal means and compact letter display
emm <- emmeans(model, ~Guild, adjust = "sidak")


cld_Guild <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = str_trim(.group)) # Clean whitespace

# Combine emmeans letters with y-position
label_df <- left_join(cld_Guild, Guild_max, by = "Guild")

label_df


# Get ANOVA p-value for title
P <- Anova(model)$`Pr(>F)`[1]
model_label <- case_when(
  P <= 0.0001 ~ "p <0.0001",
  P <= 0.001 ~ "p <0.001",
  P <= 0.01 ~ "p <0.01",
  P <= 0.05 ~ "p <0.05",
  P <= 0.1 ~ paste0("p=", as.character(round(P, 2))),
  TRUE ~ "ns"
)
# Create annotation label
model_label <- paste0(model_label)

# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~Guild) %>%
  as.data.frame()

# pairwise contrasts: Type within each Guild
contrasts_type_within_Guild <- emmeans(model, pairwise ~ Guild, adjust = "sidak")

post_hoc_Guild <- contrasts_type_within_Guild$contrasts %>%
  as.data.frame() %>%
  mutate(Variable = "Guild")

post_hoc_Guild

# Extract significant results
sig_labels <- contrasts_type_within_Guild$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.11) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.06 ~ paste0("p=", as.character(round(p.value, 3))),
    ),
    xmin = str_trim(str_extract(contrast, "(^[^-]+)")), # before "-"
    xmax = str_trim(ifelse(
      str_detect(contrast, "\\(.*\\)"),
      str_extract(contrast, "(?<=\\().*?(?=\\))"), # extract inside parentheses
      str_extract(contrast, "(?<=-)\\s*.*$") # fallback: extract after "-"
    ))
  ) %>%
  left_join(Guild_max, by = c("xmin" = "Guild"))



library(ggpubr)


# Main plot
AM_rhiz_guild <- Soil_OTU_rhizo_B %>%
  ggplot(aes(x = Guild, y = Beta, fill = Guild)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 2, alpha = 0.6) +

  # Violin and jitter
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1, color = "black") +
  # geom_jitter(aes(shape=family, fill=genus),alpha=.7,fill='black',stroke=1.8, width=.18, size=5)+
  # scale_shape_manual(values = c(1,21, 22, 23, 24,0,25), name='Family') +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x = Guild, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x = Guild, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +

  # Labels and formatting
  labs(y = expression(paste("Likelihood of occurrence post-fire (", beta[fire], ")")), x = "AM Guild", tag = "c)") +
  theme_classic() +
  scale_fill_grey(start = .9, end = 0.1) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.5, colour = "black")
  ) +
  # #Add brackets
  # geom_bracket(
  #   data = sig_labels,
  #   aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
  #   inherit.aes = FALSE,
  #   tip.length = 0.01,
  #   label.size = 7,
  #   size = 1.1
  # )+
  # Add text labels
  geom_text(
    data = label_df,
    aes(x = Guild, y = y.position + 0.01, label = .group),
    inherit.aes = FALSE, parse = TRUE,
    size = 8,
    fontface = "bold"
  ) +
  annotate("text",
    x = 0.4, y = max(Soil_OTU_rhizo_B$Beta, na.rm = TRUE),
    label = model_label, hjust = -0.3, vjust = 1, size = 8
  )

AM_rhiz_guild


plotly::ggplotly(AM_rhiz_guild)

ggsave("plots/Effect_Size_AM_Guild_Plot.png", AM_rhiz_guild, width = 10, height = 8, dpi = 300)

# library(patchwork)
#
# AM_rhiz<-(AM_rhiz_intra+
#             labs(tag = "a)") +
#             theme(plot.tag = element_text(size = 65, face = "bold")))+
#   (AM_rhiz_extra+
#      labs(tag = "b)") +
#      theme(plot.tag = element_text(size = 65, face = "bold")))+
#   plot_layout( widths = c(1, 1.1))
#
#
# AM_rhiz
#
# ggsave("plots/Effect_Size_Intra_Guild_Plot.png", AM_rhiz, width = 35, height = 25, dpi = 300)


save(post_hoc_Guild, anova_tbl_Guild, smry_tbl_Guild, file = "HMSC_MER/Output/Processed_Data/Intra_Extra_radical_hyph_summary.RDS")
