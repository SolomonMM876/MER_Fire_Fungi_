library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)
library(multcomp)


#load Beta1 and 2 from soil data HMSC
#load('HMSC_MER/results/mcmc_output.RData')

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
Beta_estimates_PA <- Beta1$mean


Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)

# Select OTUs associated with phylum
soil_OTU_phy <- myco_tax_soil %>%
  distinct(OTU,phylum,family,genus,exploration_type) 

###############################################
###PA MODEL###############
#########################

#Join and merge with Beta_estimates
Beta_PA_phy_soil<- soil_OTU_phy %>%
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(phylum=as.factor(phylum))

# Plot
ggplot(Beta_PA_phy_soil, aes(x = phylum, y = Beta, fill=phylum)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,  # wider dodge increases spacing
) +
  geom_jitter( alpha = 0.9, size = 1, width = 0.05  # must match boxplot dodge
) +
  labs(y = "Effect size (β_fire)", x = "Phylum") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


#PA model
model <- lm(Beta ~phylum, data = Beta_PA_phy_soil)

# Check model assumptions visually
check_model(model)
summary(model)
Anova(model,test= 'F')

smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_phy <- tibble(
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

anova_tbl_phy

smry_tbl_phy<-  smry$coefficients%>% 
  as.data.frame() %>% 
  rename(p.value = `Pr(>|t|)`)%>%
  rownames_to_column(var='Coefficent') %>% 
  mutate(Coefficent= str_remove(Coefficent,'phylum'),
         Term= 'Phylum',
             p.label = case_when(
               p.value <= 0.001 ~ "***",
               p.value <= 0.01  ~ "**",
               p.value <= 0.05  ~ "*",
               p.value <= 0.1   ~ '+',
               TRUE             ~ ""
             )
           )
smry_tbl_phy

# Get label positions (max y for each phylum)
phylum_max <- Beta_PA_phy_soil %>%
  group_by(phylum) %>%
  summarise(y.position = max(Beta, na.rm = TRUE) + 0.05, .groups = "drop")

# Combine with significance
smry_phy <- left_join( phylum_max,smry_tbl_phy, by = c('phylum'="Coefficent")) %>%
  filter(p.label != "")  # Keep only significant labels
smry_phy

# Get estimated marginal means and compact letter display
emm <- emmeans(model, ~ phylum,adjust = "sidak")


cld_phylum <- cld(emm, adjust = "sidak", Letters = letters) %>%
  as.data.frame() %>% 
  mutate(.group = str_trim(.group))  # Clean whitespace


# Combine emmeans letters with y-position
label_df <- left_join(cld_phylum, phylum_max, by = "phylum")

label_df


# Get ANOVA p-value for title
P <- Anova(model)$`Pr(>F)`[1]
model_label <- case_when(
  P <= 0.0001 ~ "p <0.0001",
  P <= 0.001 ~ "p <0.001",
  P <= 0.01  ~ "p <0.01",
  P <= 0.05  ~ "p <0.05",
  P <= 0.1 ~ "+",
  TRUE ~ "ns"
)
# Create annotation label
model_label <- paste0(model_label)

# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~ phylum) %>%
  as.data.frame()

# pairwise contrasts: Type within each phylum
contrasts_type_within_phylum <- emmeans(model, pairwise ~ phylum, adjust = "sidak")

post_hoc_phy<-contrasts_type_within_phylum$contrasts %>% 
  as.data.frame() %>% 
  mutate(Variable='Phylum')

post_hoc_phy


# Extract significant results
sig_labels <- contrasts_type_within_phylum$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.1) %>%
  mutate(
    p.label = case_when(
      p.value <= 0.001 ~ "***",
      p.value <= 0.01 ~ "**",
      p.value <= 0.05 ~ "*",
      p.value <= 0.06 ~ '+'
    ),
    xmin = str_trim(str_extract(contrast, "(^[^-]+)")),         # before "-"
    xmax = str_trim(ifelse(
      str_detect(contrast, "\\(.*\\)"),
      str_extract(contrast, "(?<=\\().*?(?=\\))"),   # extract inside parentheses
      str_extract(contrast, "(?<=-)\\s*.*$")         # fallback: extract after "-"
    ))) %>% 
  left_join(phylum_max, by=c('xmin'= 'phylum'))






# Main plot
phy <- Beta_PA_phy_soil %>% 
  ggplot(aes(x = phylum, y = Beta, fill = phylum)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  # Violin and jitter
  geom_violin(alpha = 0.7, width = 0.5, linewidth = 1, color = "black") +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x= phylum,y = emmean), shape = 21, size = 5, fill = "white", color = "black", stroke = 2, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x= phylum, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  # Labels and formatting
  labs(y = expression(paste('Likelihood of occurrence post-fire (', beta[fire], ')')),
       x = "Phylum", tag='a)')+
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
  #Add text labels
  geom_text(
    data = label_df,
    aes(x = phylum, y = y.position, label = .group),
    inherit.aes = FALSE,parse = TRUE,
    size = 8,
    fontface = "bold"
  )+
  annotate("text", x = 0.2, y = max(Beta_PA_phy_soil$Beta, na.rm = TRUE)+.2,
           label = model_label, hjust = -0.3, vjust = 1, size = 8)

phy


save(post_hoc_phy,anova_tbl_phy,smry_tbl_phy, file = 'HMSC_MER/Output/Processed_Data/Phylum_summary.RDS')

ggsave("plots/Effect_Size_phylum_Plot.png", phy, width = 10, height = 8, dpi = 300)








