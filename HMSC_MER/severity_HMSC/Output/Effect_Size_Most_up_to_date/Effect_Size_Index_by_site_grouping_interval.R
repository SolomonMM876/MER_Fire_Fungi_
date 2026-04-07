library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(ggpubr)
library(readxl)

# ----------------------------
# Load data
# ----------------------------
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

TSF <- readRDS('Processed_data/metadata/fire/Fire_regime_MER.Rdata') %>%
  distinct(Site, severity, TSF_years)

MER_Site_description<-read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx",
                                 sheet = "SP suggestions (highlighted)") #%>% 
  #mutate(Fire_interval= if_else(Site=='WAMCOOEucWoo02','moderate',Fire_interval))

MER_aridity<-readRDS('Processed_data/metadata/Aridity_MER.Rdata')

colnames(MER_Site_description)
Mer_veg <- MER_Site_description %>%
  left_join(TSF) %>% 
  dplyr::select(Site, Site_nickname, Last_fire_severity, severity,TSF_years,
                Veg_type,Veg_type_alt,Veg_type_alt_2, Fire_interval, Min_TFI = `Min TFI`,
                Simple_fire_response_dominants, Fire_response) %>% 
  mutate(Min_TFI=as.numeric(Min_TFI)) %>% 
  left_join(MER_aridity)

load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

test_vars<-Mer_veg %>% dplyr::select(Site_nickname, aridity, severity,TSF=TSF_years, Veg_type_alt,Fire_interval,Simple_fire_response_dominants,Fire_response)

write.csv(test_vars, file ='Processed_data/metadata/veg/Site_test_var.csv', row.names = FALSE )

Beta_estimates_PA <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Severity)

soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site, Plot, OTU,guild2) %>%
  left_join(Mer_veg)

Beta_PA_Veg_soil <- soil_Site_OTU %>%
  left_join(Beta_estimates_PA, by = "OTU")

Beta_PA_Veg_soil %>% distinct(Simple_fire_response_dominants)
# -----------------------------------------------
# Variables to test
# -----------------------------------------------
vars <- c('aridity', 'TSF_years', 'Veg_type',# 'Veg_type_alt', 'Veg_type_alt_2',
          'Fire_interval', 'Simple_fire_response_dominants', 'Fire_response')

# -----------------------------------------------
# Function to run model + extract ANOVA
# -----------------------------------------------
run_anova <- function(var) {
  # Define model formula dynamically
  f <- as.formula(paste0("Beta ~ ", var, " + (1|Site)"))
  
  # Fit LMM
  model <- tryCatch(
    lmer(f, data = Beta_PA_Veg_soil),
    error = function(e) return(NULL)
  )
  
  if (is.null(model)) return(NULL)
  
  # Compute ANOVA
  anov <- Anova(model, type = 2, test = "F")
  r2_vals <- performance::r2(model)
  
  # Build tidy table
  tibble(
    Variable = var,
    sum_sq = anov$`Sum Sq`[1],
    sum_sq_resid = anov$`Sum Sq`[2],
    F_value = anov$F,
    DF_num = anov$Df,
    DF_denom = anov$Df.res,
    P = anov$`Pr(>F)`,
    R2_marginal = r2_vals$R2_marginal,
    R2_conditional = r2_vals$R2_conditional
  )
}

# -----------------------------------------------
# Run across all variables & bind together
# -----------------------------------------------
anova_results_site_vars <- vars %>%
  map_dfr(run_anova)

# -----------------------------------------------
# Preview results
# -----------------------------------------------
print(anova_results_site_vars)







# -----------------------------------------------
# Function to run model + extract ANOVA
# Returns one row per term (main effects + interaction)
# -----------------------------------------------
run_anova <- function(var) {
  
  f <- as.formula(paste0("Beta ~ ", var, " * guild2 + (1|Site)"))
  
  model <- tryCatch(
    lmer(f, data = Beta_PA_Veg_soil),
    error = function(e) {
      message("Model failed for variable: ", var, " — ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(model)) return(NULL)
  
  anov     <- Anova(model, type = 2, test = "F")
  r2_vals  <- performance::r2(model)
  
  # anov has one row per term — bind all rows, labelled
  tibble(
    Variable        = var,
    Term            = rownames(anov),          # e.g. "var", "guild2", "var:guild2"
    F_value         = anov$F,
    DF_num          = anov$Df,
    DF_denom        = anov$Df.res,
    P               = anov$`Pr(>F)`,
    sig             = case_when(               # quick significance flag
      P <= 0.001 ~ "***",
      P <= 0.01  ~ "**",
      P <= 0.05  ~ "*",
      P <= 0.1   ~ ".",
      TRUE       ~ "ns"
    ),
    R2_marginal     = r2_vals$R2_marginal,     # fixed effects only
    R2_conditional  = r2_vals$R2_conditional   # fixed + random
  )
}

# -----------------------------------------------
# Run across all variables & bind together
# -----------------------------------------------
anova_results_site_vars <- vars %>%
  map_dfr(run_anova)

# -----------------------------------------------
# Preview — interaction rows highlighted
# -----------------------------------------------
print(anova_results_site_vars)

# Isolate just the interaction terms for easy scanning
interaction_results <- anova_results_site_vars %>%
  filter(str_detect(Term, ":"))

cat("\n--- Interaction terms only ---\n")
print(interaction_results)














###PA MODEL###############
#########################

#PA model
model <- lmer(Beta ~ Fire_interval + (1|Site), data = Beta_PA_Veg_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model, test='F')

smry <- summary(model)
anov <- Anova(model, type = 2,test='F')
r2_vals <- performance::r2(model)  # returns tibble with R2_marginal and R2_conditional

anova_tbl_Fire_interval <- tibble(
  Term = rownames(anov),
  F = anov$F,
  DF_num = anov$Df,
  DF_denom = anov$Df.res,
  R2_marginal = r2_vals$R2_marginal,
  R2_conditional = r2_vals$R2_conditional,
  P = anov$`Pr(>F)`[1]
)

anova_tbl_Fire_interval

anova_p <- anova_tbl_Fire_interval$P  # First term in the model
# Get significance level
model_label <- case_when(
  anova_p <= 0.0001 ~ "p <0.0001",
  anova_p <= 0.001 ~ "p <0.001",
  anova_p <= 0.01 ~ "p <0.01",
  anova_p <= 0.1 ~ paste0('p=',round(anova_p,3)),
  TRUE ~ "ns"
)

# Tukey pairwise contrasts: Type within each Fire_interval
contrasts_type_within <- emmeans(model, pairwise ~ Fire_interval, adjust = "sidak")

# Extract significant results
sig_labels <- contrasts_type_within$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.1) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.1 ~ paste0('p=',round(p.value,3))
    )
    )

post_hoc_veg<-contrasts_type_within$contrasts %>% 
  as.data.frame() %>% 
  mutate(Variable='Fire_interval')

post_hoc_veg

Fire_interval_max <- Beta_PA_Veg_soil %>%
  group_by(Fire_interval) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

# Post-hoc pairwise comparisons for interaction
emm <- emmeans(model, ~ Fire_interval)

library(multcomp)
# Get estimated marginal means and compact letter display
cld_Fire_interval  <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>% 
  mutate(.group = str_trim(.group))  # Clean whitespace

# Combine emmeans letters with y-position
label_df <- left_join(cld_Fire_interval ,Fire_interval_max, by = "Fire_interval")%>% 
  mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'),
                                labels = c('Short','Medium','Long')))

label_df



# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~ Fire_interval) %>%
  as.data.frame()%>% 
  mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'),
                                labels = c('Short','Medium','Long')))

Beta_PA_Veg_soil %>%
  group_by(Fire_interval) %>%
  summarise(n_OTUs = n_distinct(OTU))

Beta_PA_Veg_soil %>%
  distinct(Fire_interval, Site) %>% arrange(Fire_interval)

# Main boxplot
Fire_interval<-Beta_PA_Veg_soil%>% 
  mutate(Fire_interval = factor(Fire_interval, levels = c('short','moderate','long'),
                                labels = c('Short','Medium','Long'))) %>% 
  ggplot( aes(x = Fire_interval, y = Beta, fill=Fire_interval)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  # Violin and jitter
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  #geom_jitter(alpha=.5,aes(color=Site), width = .2)+
  
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x= Fire_interval,y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x= Fire_interval, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  # Labels and formatting
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')),
       x = "Expected Fire Interval", tag='b)') +
  theme_classic() +
  scale_fill_grey(start = .9, end = 0.1) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18,vjust=.4, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.5, colour = "black")
  ) +
  #scale_x_discrete(labels = Fire_interval_labels_named)+
  geom_text(
    data = label_df,
    aes(x = Fire_interval , y = y.position, label = .group),
    inherit.aes = FALSE,parse = TRUE,
    size = 8,
    fontface = "bold"
  )+
  annotate("text", x = 0.2, y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE),
           label = model_label, hjust = -0.3, vjust = 1, size = 8)
Fire_interval

ggsave("plots/Effect_Size_Fire_interval.png", Fire_interval, width = 12, height = 8, dpi = 300)


save(post_hoc_veg,anova_tbl_Fire_interval,anova_results_site_vars, file = 'HMSC_MER/Output/Processed_Data/Fire_interval_summary.RDS')


#####################
###########
###aridity###############
#########################

#PA model
model <- lmer(Beta ~ aridity + (1|Site), data = Beta_PA_Veg_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model, test='F')

smry <- summary(model)
anov <- Anova(model, type = 2,test='F')
r2_vals <- performance::r2(model)  # returns tibble with R2_marginal and R2_conditional

anova_tbl_aridity <- tibble(
  Term = rownames(anov),
  F = anov$F,
  DF_num = anov$Df,
  DF_denom = anov$Df.res,
  R2_marginal = r2_vals$R2_marginal,
  R2_conditional = r2_vals$R2_conditional,
  P = anov$`Pr(>F)`[1]
)

anova_tbl_aridity

anova_p <- anova_tbl_aridity$P  # First term in the model
# Get significance level
model_label <- case_when(
  anova_p <= 0.0001 ~ "p <0.0001",
  anova_p <= 0.001 ~ "p <0.001",
  anova_p <= 0.01 ~ "p <0.01",
  anova_p <= 0.1 ~ paste0('p=',round(anova_p,3)),
  TRUE ~ "ns"
)

library(ggeffects)

# Generate predicted line from model
pred <- ggpredict(model, terms = "aridity") %>% 
  as.data.frame() %>% 
  mutate(aridity = factor(x))   # ensure same factor structure

aridity_max <- Beta_PA_Veg_soil %>%
  group_by(aridity) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

# Post-hoc pairwise comparisons for interaction
emm <- emmeans(model, ~ aridity)

library(multcomp)
# Get estimated marginal means and compact letter display
cld_aridity  <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>% 
  mutate(.group = str_trim(.group))  # Clean whitespace

# Combine emmeans letters with y-position
label_df <- left_join(cld_aridity ,aridity_max, by = "aridity")

label_df



# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~ aridity) %>%
  as.data.frame()%>% 
  mutate(aridity = factor(aridity))

Beta_PA_Veg_soil %>%
  group_by(aridity) %>%
  summarise(n_OTUs = n_distinct(OTU))

Beta_PA_Veg_soil %>%
  distinct(aridity, Site) %>% arrange(aridity)

# Main boxplot
aridity<-Beta_PA_Veg_soil%>% 
  ggplot( aes(x = aridity, y = Beta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_line(data = pred, aes(x = x, y = predicted, group = 1),
            linewidth = 2, color = "black", inherit.aes = FALSE) +
  geom_ribbon(data = pred,aes(x = x, ymin = conf.low, ymax = conf.high),
              inherit.aes = FALSE, alpha = 0.2, color = NA)+
  geom_point( aes(x = aridity, y = Beta),
             size = 2, shape = 21, fill = "black", alpha=.7,
             inherit.aes = FALSE) +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')),
       x = "Aridity Index", tag='c)') +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18,vjust=.4, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.5, colour = "black")
  ) +
  geom_text(
    data = label_df,
    aes(x = aridity , y = y.position, label = .group),
    inherit.aes = FALSE,parse = TRUE,
    size = 8,
    fontface = "bold"
  )+
  annotate("text", x = -Inf, y = max(Beta_PA_Veg_soil$Beta, na.rm = TRUE)+.15,
           label = model_label, hjust = -0.3, vjust = 1, size = 8)
aridity

ggsave("plots/Effect_Size_aridity.png", aridity, width = 12, height = 8, dpi = 300)


save(post_hoc_veg,anova_tbl_aridity,anova_results_site_vars, file = 'HMSC_MER/Output/Processed_Data/aridity_summary.RDS')
