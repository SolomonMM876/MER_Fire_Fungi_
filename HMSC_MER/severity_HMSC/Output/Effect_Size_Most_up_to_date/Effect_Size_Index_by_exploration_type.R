library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(ggpubr)
library(multcomp)

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#All tax load
myco_tax_all<-bind_rows(myco_tax_soil,myco_tax_hyph) %>% distinct()


#load most up to date PA model
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean


Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)


# Select OTUs associated with exploration_type

soil_OTU_explo <- myco_tax_soil %>%
  filter(!is.na(exploration_type)) %>% 
  distinct(OTU,exploration_type) 

###PA MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_PA_explo_soil<- soil_OTU_explo %>%
  inner_join(Beta_estimates_PA, by = "OTU") 



# Plot
ggplot(Beta_PA_explo_soil, aes(x = exploration_type, y = Beta, fill=exploration_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta estimates", x = "exploration_type", title = "Beta values by exploration_type for hyph taxa PA model") +
  theme_minimal() + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(palette = "Set2")


#PA model
model <- lm(Beta ~ exploration_type, data = Beta_PA_explo_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model)


smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_explo <- tibble(
  Term = rownames(anov)[1],
  sum_sq= anov$`Sum Sq`[1],
  sum_sq_resid= anov$`Sum Sq`[2],
  F = anov$`F value`[1],
  DF_num = anov$Df[1],
  DF_denom = anov$Df[2],
  R2 = smry$r.squared,
  adj_R2 = smry$adj.r.squared,
  P = anov$`Pr(>F)`[1]
)

anova_tbl_explo


exploration_type_max <- Beta_PA_explo_soil %>%
  group_by(exploration_type) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.05)

anova_p <- Anova(model)$`Pr(>F)`[1]  # First term in the model
# Get significance level
p_sig <- case_when(
  anova_p <= 0.0001 ~ "p <0.0001",
  anova_p <= 0.001 ~ "p <0.001",
  anova_p <= 0.1    ~ paste0('p=',as.character(round(anova_p, 3))),
  TRUE ~ "ns"
)


# ================================================================
# Post-hoc pairwise contrasts for exploration_type
# ================================================================
# Get estimated marginal means and compact letter display
emm <- emmeans(model, ~ exploration_type ,adjust = "sidak")
cld_explo  <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>% 
  mutate(.group = str_trim(.group))  # Clean whitespace
cld_explo

# Combine emmeans letters with y-position
label_df <- left_join(cld_explo , exploration_type_max, by = "exploration_type")

label_df




# Tukey pairwise contrasts: Type within each exploration_type
contrasts_type_within_exploration_type <- emmeans(model, pairwise ~ exploration_type, adjust = "sidak")

# Extract significant results
post_hoc_explo <- contrasts_type_within_exploration_type$contrasts %>%
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

post_hoc_explo



# contrasts_type_within_exploration_type <- emmeans(model, pairwise ~ exploration_type, adjust = "sidak")
# 
# sig_labels_explo <- contrasts_type_within_exploration_type$contrasts %>%
#   as.data.frame() %>%
#   filter(p.value < 0.1) %>%                    # keep interesting contrasts
#   mutate(
#     p.label = case_when(
#       p.value <= 0.001 ~ "***",
#       p.value <= 0.01  ~ "**",
#       p.value <= 0.05  ~ "*",
#       p.value <= 0.1   ~ paste0("p=", round(p.value, 3)),
#       TRUE             ~ ""
#     ),
#     # split into left / right around the " - " separator (keeps whatever is inside parentheses or not)
#     left_part  = str_trim(str_split_fixed(contrast, " - ", 2)[,1]),
#     right_part = str_trim(str_split_fixed(contrast, " - ", 2)[,2]),
#     
#     # if a part contains parentheses, extract inside them; otherwise use the whole part
#     xmin = if_else(
#       str_detect(left_part, "\\("),
#       str_extract(left_part, "(?<=\\().*?(?=\\))"),
#       left_part
#     ),
#     xmax = if_else(
#       str_detect(right_part, "\\("),
#       str_extract(right_part, "(?<=\\().*?(?=\\))"),
#       right_part
#     ),
#     
#     # final cleanup: remove any stray parentheses or surrounding spaces
#     xmin = str_remove_all(xmin, "[\\(\\)]") %>% str_trim(),
#     xmax = str_remove_all(xmax, "[\\(\\)]") %>% str_trim()
#   ) %>%
#   select(contrast, estimate, SE, df, t.ratio, p.value, p.label, xmin, xmax)
# 
# # Inspect
# sig_labels_explo
# # Join y.position (align xmax factor to dataset)
# sig_results_explo <- sig_labels_explo %>%
#   left_join(exploration_type_max, by = c("xmax" = "exploration_type"))

# ================================================================
# Get emmeans with 95% CIs for plotting
# ================================================================


# Create annotation label
model_label <- paste0(p_sig)
#### create axis labels
exploration_labels_named <- c(
  "contact"                  = "contact\n",
  "short-distance_coarse"   = "short\ndistance\ncoarse",
  "short-distance_delicate" = "short\ndistance\ndelicate",
  "mat"                     = "mat\n",
  "medium-distance_fringe"  = "medium\ndistance\nfringe",
  "medium-distance_smooth"  = "medium\ndistance\nsmooth",
  "long-distance"           = "long-\ndistance\n"
)

exploration_levels <- names(exploration_labels_named)


emm_df <- emmeans(model, ~ exploration_type) %>%
  as.data.frame() %>%
  mutate(exploration_type = factor(exploration_type, levels = exploration_levels))
# Main boxplot
explo<-Beta_PA_explo_soil %>% 
  mutate(exploration_type = factor(exploration_type, levels = exploration_levels)) %>% 
  ggplot( aes(x = exploration_type, y = Beta, fill=exploration_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x= exploration_type,y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x= exploration_type, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')),
       x = "Exploration type", tag='d)') +
  theme_classic() +
  theme(    legend.position = "none",
            axis.text.x = element_text(hjust = 0.5,vjust=1, size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 18,vjust=.3, face = "bold"),
            axis.title.y = element_text(size = 18, face = "bold"),
            axis.line = element_line(linewidth = 1.5, colour = "black")
        )+
  scale_x_discrete(labels = exploration_labels_named)+
  scale_fill_grey(start = .9, end = 0.1) +
  #Add text labels
  geom_text(
    data = label_df,
    aes(x = exploration_type , y = y.position, label = .group),
    inherit.aes = FALSE,parse = TRUE,
    size = 8,
    fontface = "bold"
  )+
  annotate("text", x = 0.2, y = max(Beta_PA_explo_soil$Beta+0.15, na.rm = TRUE),
           label = model_label, hjust = -.3, vjust = 1, size = 8)

explo

ggsave("plots/Effect_Size_explo_Plot.png", explo, width = 14, height = 8, dpi = 300)

save(anova_tbl_explo,post_hoc_explo, file = 'HMSC_MER/Output/Processed_Data/explo_summary.RDS')


# ================================================================
# Ruhlandiella-specific plot
# ================================================================

ruhl_OTUs <- myco_tax_soil %>%
  filter(genus == "Ruhlandiella") %>%
  distinct(OTU, exploration_type)

# Subset of Beta estimates for Ruhlandiella OTUs only
Beta_PA_ruhl <- Beta_PA_explo_soil %>%
  filter(OTU %in% ruhl_OTUs$OTU) %>%
  mutate(exploration_type = factor(exploration_type, levels = exploration_levels))

ruhl <- Beta_PA_explo_soil %>%
  mutate(exploration_type = factor(exploration_type, levels = exploration_levels)) %>%
  ggplot(aes(x = exploration_type, y = Beta, fill = exploration_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 1) +
  geom_point(data = emm_df, aes(x = exploration_type, y = emmean),
             shape = 21, size = 5, fill = "white", color = "black", stroke = 1.5,
             inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x = exploration_type, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  # Ruhlandiella raw OTUs
  geom_point(data = Beta_PA_ruhl, aes(x = exploration_type, y = Beta),
             shape = 21, size = 3, fill = "red", color = "black", stroke = 1,
             position = position_jitter(width = 0.05, seed = 42),
             inherit.aes = FALSE) +
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')),
       x = "Exploration type", tag = 'e)') +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(hjust = 0.5, vjust = 1, size = 14, face = "bold"),
    axis.text.y  = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 18, vjust = .3, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line    = element_line(linewidth = 1.5, colour = "black")
  ) +
  scale_x_discrete(labels = exploration_labels_named) +
  scale_fill_grey(start = .9, end = 0.1) +
  geom_text(data = label_df, aes(x = exploration_type, y = y.position, label = .group),
            inherit.aes = FALSE, parse = TRUE, size = 8, fontface = "bold") +
  annotate("text", x = 0.2, y = max(Beta_PA_explo_soil$Beta + 0.15, na.rm = TRUE),
           label = model_label, hjust = -.3, vjust = 1, size = 8)

ruhl

ggsave("plots/Effect_Size_Ruhlandiella_Plot.png", ruhl, width = 14, height = 8, dpi = 300)

