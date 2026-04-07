# ==============================================================================
# SCRIPT 02: Beta Coefficient Analysis (Self-Contained)
# DESCRIPTION: Analyzes influence of grouping variables on fire response (Beta).
# ==============================================================================

# 1. SETUP & LIBRARIES
library(tidyverse)
library(lme4)
library(car)
library(performance)
library(ggpubr)
library(emmeans)
library(multcomp)
library(ggtext)
library(ggeffects)
library(tibble)
library(broom)
library(patchwork)

# SET WORKING DIRECTORY

# Create output folders inside project_script
dir.create("plots", showWarnings = FALSE)
dir.create("processed_data", showWarnings = FALSE)

# 2. LOAD DATA (FROM CLEAN CSVs)
# Since we are inside 'project_script', we just look in 'data/'
Beta_estimates_PA <- read_csv("data/beta_estimates.csv", show_col_types = FALSE) # this is derived data
# ^^^^Replace with "HMSC_Output/models_fitted.RDS" from script 1 if desired
Mer_veg           <- read_csv("data/site_metadata.csv", show_col_types = FALSE)
myco_tax_all      <- read_csv("data/taxonomy.csv", show_col_types = FALSE)
myco_dat_soil     <- read_csv("data/soil_abundance.csv", show_col_types = FALSE)
myco_dat_hyph     <- read_csv("data/hyph_abundance.csv", show_col_types = FALSE)
repo_strat_df     <- read_csv("data/traits_repo_strat.csv", show_col_types = FALSE)
amf_traits        <- read_csv("data/traits_amf_guilds.csv", show_col_types = FALSE)

# 3. PREPARE & JOIN DATA FOR ANALYSES

# --- A. Site Analysis Data ---
soil_Site_OTU <- myco_dat_soil %>% distinct(Site, Plot, OTU)

Beta_PA_Site_soil <- soil_Site_OTU %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  left_join(Mer_veg, by = "Site")

# --- B. Reproductive Strategy Data  ---
soil_OTU_repo_strat <- myco_tax_all %>%
  left_join(repo_strat_df, by = join_by(family, genus)) %>%
  mutate(repo_strategy = case_when(phylum == "Glomeromycota" ~ "Glomeromycota", TRUE ~ repo_strategy)) %>% 
  filter(!is.na(repo_strategy)) %>% 
  distinct(OTU, repo_strategy) 

Beta_PA_repo_strat_soil <- soil_OTU_repo_strat %>%
  inner_join(Beta_estimates_PA, by = "OTU")

# --- C. AMF Guild Data ---
soil_OTU_AM <- myco_tax_all %>%
  distinct(OTU, phylum, family, genus) %>% 
  filter(phylum == 'Glomeromycota')

Soil_OTU_rhizo_B <- soil_OTU_AM %>% 
  left_join(amf_traits, by = c('family' = 'Families')) %>% 
  filter(!is.na(family)) %>%
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(family = as.factor(family)) %>%
  filter(!is.na(Beta))

# --- D. Exploration Type Data ---
soil_OTU_explo <- myco_tax_all %>%
  filter(!is.na(exploration_type)) %>% 
  distinct(OTU, exploration_type)

Beta_PA_explo_soil <- soil_OTU_explo %>%
  inner_join(Beta_estimates_PA, by = "OTU")

# --- E. Phylum Data ---
soil_OTU_phy_only <- myco_tax_all %>%
  distinct(OTU, phylum)

Beta_PA_phy_soil <- soil_OTU_phy_only %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>%
  mutate(phylum = as.factor(phylum))

# ==============================================================================
# HELPER FUNCTION
# ==============================================================================
format_pval <- function(p) {
  case_when(
    p <= 0.0001 ~ "p <0.0001",
    p <= 0.001 ~ "p <0.001",
    p <= 0.01 ~ "p <0.01",
    p <= 0.05 ~ "p <0.05",
    p <= 0.1 ~ paste0("p=", round(p, 3)),
    TRUE ~ "ns"
  )
}

# ==============================================================================
# ANALYSIS SECTIONS
# ==============================================================================

# --- A: SITE EFFECT ---
model_site <- lm(Beta ~ Site, data = Beta_PA_Site_soil)
anova_site <- Anova(model_site, type = "II")
smry_site <- summary(model_site)
a_tbl_site <- tibble(Term = "Site", F = anova_site$`F value`[1], P = anova_site$`Pr(>F)`[1], R2 = smry_site$r.squared)

emm_site <- emmeans(model_site, pairwise ~ Site, adjust = "sidak")
cld_site <- cld(emm_site, adjust = "sidak", sort = FALSE, alpha = 0.05, Letters = letters)
cld_df_site <- as.data.frame(cld_site) %>%
  mutate(letter = stringr::str_trim(stringr::str_remove_all(.group, "[\\(\\)\\[\\]\\s]+")))

site_levels_aridity <- Beta_PA_Site_soil %>% distinct(Site, aridity) %>% arrange(aridity) %>% pull(Site)
site_max <- Beta_PA_Site_soil %>% group_by(Site) %>% summarise(y.position = max(Beta, na.rm=T) + 0.1)
label_df_site <- left_join(cld_df_site, site_max, by = "Site")

site_labels_named <- c(
  "NSMSEQCasWoo01" = "Forest<br>*Melaleuca/Casuarina*", "QDMSEQRainfo01" = "Subtropical-Rainforest",
  "SAMMDDMallee01" = "Triodia-Mallee<br>*Eucalyptus-dumosa*", "SAMMDDMallee02" = "Chenopod-Mallee<br>*Eucalyptus-oleosa*",
  "VCMSECEucFor01" = "Forest<br>*Eucalyptus-muelleriana*", "VCMSECRainfo01" = "Temperate-Rainforest",
  "WAMDALEucWoo01" = "Savannah<br>*Corymbia-paractia*", "WAMDALRainfo01" = "VineThicket",
  "WAMDALEucWoo02" = "Savannah<br>*Triodia-schinzii*", "WAMCOOShrubl01" = "Shrubland<br>*Allocasuarina-acutivalvis*",
  "WAMCOOEucWoo01" = "Woodland<br>*Eucalyptus-salmonophloia*", "WAMCOOEucWoo02" = "Chenopod-Mallee<br>*Eucalyptus-aff.-oleosa*",
  "WAMCOOEucWoo03" = "Woodland<br>*Eucalyptus-salubris*", "WAMESPShrubl01" = "Shrubland<br>*Kingia-australis*"
)

plot_site <- Beta_PA_Site_soil %>%
  mutate(Site = factor(Site, levels = site_levels_aridity)) %>%
  ggplot(aes(x = Site, y = Beta, fill = Site)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.6, linewidth = 1) +
  geom_point(data = as.data.frame(emm_site$emmeans), aes(x = Site, y = emmean), shape = 21, size = 10, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = as.data.frame(emm_site$emmeans), aes(x = Site, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_site, aes(x = Site, y = y.position, label = letter), inherit.aes = FALSE, size = 8, vjust = 0) +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + scale_x_discrete(labels = site_labels_named) +
  theme(legend.position = "none", axis.text.x = element_markdown(angle = 35, hjust = 1, size = 14)) +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Site", tag='a)') +
  annotate("text", x = 1, y = max(Beta_PA_Site_soil$Beta)+0.2, label = format_pval(a_tbl_site$P), hjust = 0, size = 8)

plot_site

ggsave("plots/Effect_Size_site_Plot.png", plot_site, width = 15, height = 8, dpi = 300)

# --- B: FIRE INTERVAL ---
model_fi <- lmer(Beta ~ Fire_interval + (1|Site), data = Beta_PA_Site_soil)
anova_fi <- Anova(model_fi, type = 2, test = 'F')
a_tbl_fi <- tibble(Term = "Fire_interval", F = anova_fi$F[1], P = anova_fi$`Pr(>F)`[1], R2=NA)

emm_fi <- emmeans(model_fi, ~ Fire_interval)
cld_fi <- cld(emm_fi, adjust = "sidak", Letters = letters) %>% as.data.frame() %>% mutate(.group = str_trim(.group))
fi_max <- Beta_PA_Site_soil %>% group_by(Fire_interval) %>% summarise(y.position = max(Beta, na.rm=T) + 0.05)
label_df_fi <- left_join(cld_fi, fi_max, by = "Fire_interval") %>% mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'), labels = c('Short','Medium','Long')))
emm_df_fi <- as.data.frame(emm_fi) %>% mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'), labels = c('Short','Medium','Long')))

plot_fi <- Beta_PA_Site_soil %>%
  mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'), labels = c('Short','Medium','Long'))) %>%
  ggplot(aes(x = Fire_interval, y = Beta, fill = Fire_interval)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  geom_point(data = emm_df_fi, aes(x= Fire_interval,y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df_fi, aes(x= Fire_interval, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_fi, aes(x = Fire_interval, y = y.position, label = .group), inherit.aes = FALSE, size = 8, fontface = "bold") +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + theme(legend.position = "none") +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Expected Fire Interval", tag='b)') +
  annotate("text", x = 0.6, y = max(Beta_PA_Site_soil$Beta)+0.15, label = format_pval(a_tbl_fi$P), hjust = 0, size = 8)

plot_fi

ggsave("plots/Effect_Size_Fire_interval.png", plot_fi, width = 12, height = 8, dpi = 300)

# --- C: ARIDITY ---
model_aridity <- lmer(Beta ~ aridity + (1|Site), data = Beta_PA_Site_soil)
anova_aridity <- Anova(model_aridity, type = 2, test = 'F')
a_tbl_aridity <- tibble(Term = "Aridity", F = anova_aridity$F[1], P = anova_aridity$`Pr(>F)`[1], R2=NA)
pred_aridity <- ggpredict(model_aridity, terms = "aridity")

plot_aridity <- Beta_PA_Site_soil %>%
  ggplot(aes(x = aridity, y = Beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_line(data = pred_aridity, aes(x = x, y = predicted, group = 1), linewidth = 2, inherit.aes = FALSE) +
  geom_ribbon(data = pred_aridity, aes(x = x, ymin = conf.low, ymax = conf.high), inherit.aes = FALSE, alpha = 0.2) +
  geom_point(size = 2, shape = 21, fill = "black", alpha = .7) +
  theme_classic() +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Aridity Index", tag='c)') +
  annotate("text", x = min(Beta_PA_Site_soil$aridity), y = max(Beta_PA_Site_soil$Beta)+0.15, label = format_pval(a_tbl_aridity$P), hjust = 0, size = 8)

plot_aridity

ggsave("plots/Effect_Size_aridity.png", plot_aridity, width = 12, height = 8, dpi = 300)

# --- D: REPRODUCTIVE STRATEGY ---
model_repo <- lm(Beta ~ repo_strategy, data = Beta_PA_repo_strat_soil)
anova_repo <- Anova(model_repo, type = "II")
smry_repo <- summary(model_repo)
a_tbl_repo <- tibble(Term = "Repo_Strategy", F = anova_repo$`F value`[1], P = anova_repo$`Pr(>F)`[1], R2 = smry_repo$r.squared)

emm_repo <- emmeans(model_repo, ~ repo_strategy)
cld_repo <- cld(emm_repo, adjust = "sidak", Letters = letters) %>% as.data.frame() %>% mutate(.group = str_trim(.group))
repo_max <- Beta_PA_repo_strat_soil %>% group_by(repo_strategy) %>% summarise(y.position = max(Beta, na.rm=T) + 0.05)
label_df_repo <- left_join(cld_repo, repo_max, by = "repo_strategy")
emm_df_repo <- as.data.frame(emm_repo)
repo_labels <- c("Hypogeous or Semi-hypogeous" = "Hypogeous/\nSemi-hypogeous", "Resupinate" = "Resupinate\n", "Epigeous" = "Epigeous\n", "Glomeromycota" = "Glomeromycota\n")

plot_repo <- Beta_PA_repo_strat_soil %>%
  mutate(repo_strategy = factor(repo_strategy, levels = names(repo_labels))) %>%
  ggplot(aes(x = repo_strategy, y = Beta, fill = repo_strategy)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  geom_point(data = emm_df_repo, aes(x = repo_strategy, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df_repo, aes(x = repo_strategy, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_repo, aes(x = repo_strategy, y = y.position, label = .group), inherit.aes = FALSE, size = 8, fontface = "bold") +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + theme(legend.position = "none") +
  scale_x_discrete(labels = repo_labels) +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Fruiting body type", tag='d)') +
  annotate("text", x = 0.6, y = max(Beta_PA_repo_strat_soil$Beta)+0.2, label = format_pval(a_tbl_repo$P), hjust = 0, size = 8)

plot_repo

ggsave("plots/Effect_Size_repo_strat_Plot.png", plot_repo, width = 12, height = 8, dpi = 300)

# --- E: AMF GUILD ---
model_guild <- lm(Beta ~ Guild, data = Soil_OTU_rhizo_B)
anova_guild <- Anova(model_guild, type = "II")
smry_guild <- summary(model_guild)
a_tbl_guild <- tibble(Term = "Guild", F = anova_guild$`F value`[1], P = anova_guild$`Pr(>F)`[1], R2 = smry_guild$r.squared)

emm_guild <- emmeans(model_guild, ~ Guild, adjust = "sidak")
cld_guild <- cld(emm_guild, adjust = "sidak", Letters = letters) %>% as.data.frame() %>% mutate(.group = str_trim(.group))
guild_max <- Soil_OTU_rhizo_B %>% group_by(Guild) %>% summarise(y.position = max(Beta, na.rm=T) + 0.02)
label_df_guild <- left_join(cld_guild, guild_max, by = "Guild")
emm_df_guild <- as.data.frame(emm_guild)

plot_guild <- Soil_OTU_rhizo_B %>%
  ggplot(aes(x = Guild, y = Beta, fill = Guild)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 2, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1, color = 'black') +
  geom_point(data = emm_df_guild, aes(x = Guild, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df_guild, aes(x = Guild, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_guild, aes(x = Guild, y = y.position + 0.01, label = .group), inherit.aes = FALSE, size = 8, fontface = "bold") +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + theme(legend.position = "none") +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Guild", tag='e)') +
  annotate("text", x = 0.6, y = max(Soil_OTU_rhizo_B$Beta)+0.2, label = format_pval(a_tbl_guild$P), hjust = 0, size = 8)

plot_guild

ggsave("plots/Effect_Size_Guild_Plot.png", plot_guild, width = 10, height = 8, dpi = 300)

# --- F: EXPLORATION TYPE ---
model_explo <- lm(Beta ~ exploration_type, data = Beta_PA_explo_soil)
anova_explo <- Anova(model_explo, type = "II")
smry_explo <- summary(model_explo)
a_tbl_explo <- tibble(Term = "Exploration_Type", F = anova_explo$`F value`[1], P = anova_explo$`Pr(>F)`[1], R2 = smry_explo$r.squared)

emm_explo <- emmeans(model_explo, ~ exploration_type, adjust = "sidak")
cld_explo <- cld(emm_explo, adjust = "sidak", Letters = letters) %>% as.data.frame() %>% mutate(.group = str_trim(.group))
explo_max <- Beta_PA_explo_soil %>% group_by(exploration_type) %>% summarise(y.position = max(Beta, na.rm=T) + 0.05)
label_df_explo <- left_join(cld_explo, explo_max, by = "exploration_type")
emm_df_explo <- as.data.frame(emm_explo)
exploration_labels_named <- c("contact"="contact\n", "short-distance_coarse"="short\ndistance\ncoarse", "short-distance_delicate"="short\ndistance\ndelicate", "mat"="mat\n", "medium-distance_fringe"="medium\ndistance\nfringe", "medium-distance_smooth"="medium\ndistance\nsmooth", "long-distance"="long-\ndistance\n")
explo_levels <- names(exploration_labels_named)

plot_explo <- Beta_PA_explo_soil %>%
  mutate(exploration_type = factor(exploration_type, levels = explo_levels)) %>%
  ggplot(aes(x = exploration_type, y = Beta, fill = exploration_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  geom_point(data = emm_df_explo %>% mutate(exploration_type = factor(exploration_type, levels = explo_levels)), aes(x = exploration_type, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df_explo %>% mutate(exploration_type = factor(exploration_type, levels = explo_levels)), aes(x = exploration_type, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_explo, aes(x = exploration_type, y = y.position, label = .group), inherit.aes = FALSE, size = 8, fontface = "bold") +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + theme(legend.position = "none") +
  scale_x_discrete(labels = exploration_labels_named) +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Exploration type", tag='f)') +
  annotate("text", x = 0.6, y = max(Beta_PA_explo_soil$Beta)+0.15, label = format_pval(a_tbl_explo$P), hjust = 0, size = 8)

plot_explo

ggsave("plots/Effect_Size_explo_Plot.png", plot_explo, width = 14, height = 8, dpi = 300)

# --- G: PHYLUM ---
model_phy <- lm(Beta ~ phylum, data = Beta_PA_phy_soil)
anova_phy <- Anova(model_phy, type = "II")
smry_phy <- summary(model_phy)
a_tbl_phy <- tibble(Term = "Phylum", F = anova_phy$`F value`[1], P = anova_phy$`Pr(>F)`[1], R2 = smry_phy$r.squared)

emm_phy <- emmeans(model_phy, ~ phylum, adjust = "sidak")
cld_phy <- cld(emm_phy, adjust = "sidak", Letters = letters) %>% as.data.frame() %>% mutate(.group = str_trim(.group))
phylum_max <- Beta_PA_phy_soil %>% group_by(phylum) %>% summarise(y.position = max(Beta, na.rm=T) + 0.05)
label_df_phy <- left_join(cld_phy, phylum_max, by = "phylum")
emm_df_phy   <- as.data.frame(emm_phy)

plot_phy <- Beta_PA_phy_soil %>%
  ggplot(aes(x = phylum, y = Beta, fill = phylum)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.5, linewidth = 1, color = "black") +
  geom_point(data = emm_df_phy, aes(x = phylum, y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df_phy, aes(x = phylum, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  geom_text(data = label_df_phy, aes(x = phylum, y = y.position, label = .group), inherit.aes = FALSE, size = 8, fontface = "bold") +
  scale_fill_grey(start = .9, end = 0.1) + theme_classic() + theme(legend.position = "none") +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')), x = "Phylum", tag='g)') +
  annotate("text", x = 0.6, y = max(Beta_PA_phy_soil$Beta)+0.2, label = format_pval(a_tbl_phy$P), hjust = 0, size = 8)

plot_phy

ggsave("plots/Effect_Size_phylum_Plot.png", plot_phy, width = 10, height = 8, dpi = 300)

# --- H: RELATIVE ABUNDANCE RATIO ---
hyph_count <- myco_dat_hyph %>% summarise(total = sum(count, na.rm = TRUE)) %>% pull()
soil_count <- myco_dat_soil %>% summarise(total = sum(count, na.rm = TRUE)) %>% pull()
hyph_RA <- myco_dat_hyph %>% group_by(OTU) %>% summarise(RA_hyph = sum(count,na.rm=T)/hyph_count, .groups="drop")
soil_RA <- myco_dat_soil %>% group_by(OTU) %>% summarise(RA_soil = sum(count,na.rm=T)/soil_count, .groups="drop")
RA_joined <- full_join(hyph_RA, soil_RA, by = "OTU") %>% mutate(RA_hyph = replace_na(RA_hyph, 0))
min_ratio <- RA_joined %>% filter(RA_hyph>0 & RA_soil>0) %>% mutate(ratio=RA_hyph/RA_soil) %>% summarise(min=min(ratio,na.rm=T)) %>% mutate(pc=round(min/2,4)) %>% pull(pc)
RA_ratio_PA <- RA_joined %>% mutate(ratio = RA_hyph/RA_soil, log_ratio = log10(ratio + min_ratio)) %>% inner_join(Beta_estimates_PA, by = "OTU") %>% mutate(Abs_Beta = abs(Beta))

analyze_beta_type <- function(data, y_var, beta_type_label, plot_color) {
  model <- lm(as.formula(paste(y_var, "~ log_ratio")), data = data)
  smry <- summary(model)
  anov <- Anova(model, test.statistic = "F")
  anova_tbl <- tibble(Type = beta_type_label, F = anov$`F value`[1], P = anov$`Pr(>F)`[1], R2 = smry$r.squared)
  
  label_text_p <- if(anova_tbl$P > 0.10) "ns" else paste0("italic(p):", round(anova_tbl$P, 2))
  label_text_R <- if(anova_tbl$P > 0.10) "" else paste0("R^2: ", round(anova_tbl$R2, 2))
  parse_flag <- anova_tbl$P <= 0.10
  
  plot <- ggplot(data, aes(x = log_ratio, y = !!sym(y_var))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1, alpha = 0.6) +
    geom_point(size = 4, alpha = 0.6, color = "black") +
    geom_smooth(method = "lm", se = FALSE, color = plot_color, linewidth = 2) +
    theme_classic() +
    labs(y = eval(parse(text = beta_type_label)), x = expression(log[10]~"(RA Hyphae / RA Soil)")) +
    annotate("text", x = min(data$log_ratio, na.rm = TRUE)+0.4, y = max(data[[y_var]], na.rm = TRUE)-0.03, label = label_text_p, parse = parse_flag, hjust = 0, size = 8) +  
    annotate("text", x = min(data$log_ratio, na.rm = TRUE)+0.4, y = max(data[[y_var]], na.rm = TRUE)-0.055, label = label_text_R, parse = parse_flag, hjust = 0, size = 8) +
    theme(axis.text = element_text(size = 16, face = "bold"), axis.title = element_text(size = 18, face = "bold"), axis.line = element_line(linewidth = 2))
  return(list(anova = anova_tbl, plot = plot))
}

res_raw <- analyze_beta_type(RA_ratio_PA, "Beta", "expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')'))", "purple4")
res_abs <- analyze_beta_type(RA_ratio_PA, "Abs_Beta", "expression(paste('Absolute likelihood of occurrence post-fire (', beta[fire], ')'))", "green")
res_pos <- analyze_beta_type(RA_ratio_PA %>% filter(Beta > 0), "Beta", "expression(paste('Positive likelihood of occurrence post-fires (', beta[fire], '>0)'))", "blue")
res_neg <- analyze_beta_type(RA_ratio_PA %>% filter(Beta < 0), "Beta", "expression(paste('Negative likelihood of occurrence post-fire (', beta[fire], '<0)'))", "red")

anova_table_all_RA <- bind_rows(res_abs$anova, res_raw$anova, res_pos$anova, res_neg$anova) %>% mutate(Term="RA_Ratio")
plot_a <- res_raw$plot + labs(tag = "a)") + theme(plot.tag = element_text(size = 16, face = "bold"))
plot_b <- res_abs$plot + labs(tag = "b)") + theme(plot.tag = element_text(size = 16, face = "bold"))
plot_c <- res_pos$plot + labs(tag = "c)") + theme(plot.tag = element_text(size = 16, face = "bold"))
plot_d <- res_neg$plot + labs(tag = "d)") + theme(plot.tag = element_text(size = 16, face = "bold"))
final_RA_plot <- (plot_a + plot_b) / (plot_c + plot_d)

final_RA_plot

ggsave("plots/RA_OTU_Effect_Size_Plot.png", final_RA_plot, width = 14, height = 14, dpi = 300)

# ==============================================================================
# COMBINE STATS TABLES
# ==============================================================================
master_anova <- bind_rows(
  a_tbl_site, a_tbl_fi, a_tbl_aridity, a_tbl_repo, a_tbl_guild, a_tbl_explo, a_tbl_phy,
  anova_table_all_RA
)

master_anova

# Extract contrasts for master post-hoc table
post_hoc_list <- list(
  site = as.data.frame(emmeans(emm_site, pairwise ~ Site, adjust="sidak")$contrasts) %>% mutate(Variable="Site"),
  fi = as.data.frame(emmeans(emm_fi, pairwise ~ Fire_interval, adjust="sidak")$contrasts) %>% mutate(Variable="Fire_Interval"),
  repo = as.data.frame(emmeans(emm_repo, pairwise ~ repo_strategy, adjust="sidak")$contrasts) %>% mutate(Variable="Repo_Strategy"),
  guild = as.data.frame(emmeans(emm_guild, pairwise ~ Guild, adjust="sidak")$contrasts) %>% mutate(Variable="Guild"),
  explo = as.data.frame(emmeans(emm_explo, pairwise ~ exploration_type, adjust="sidak")$contrasts) %>% mutate(Variable="Exploration"),
  phy = as.data.frame(emmeans(emm_phy, pairwise ~ phylum, adjust="sidak")$contrasts) %>% mutate(Variable="Phylum")
)
master_posthoc <- bind_rows(post_hoc_list)

master_posthoc

write_csv(master_anova, "processed_data/Master_ANOVA_Results.csv")
write_csv(master_posthoc, "processed_data/Master_PostHoc_Results.csv")

print("Analysis Complete. All outputs saved in 'project_script/'.")