library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(ggpubr)
library(emmeans)

#load Beta  soil data HMSC

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#load most up to date PA model
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean

Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)

MER_aridity<-readRDS('Processed_data/metadata/Aridity_MER.Rdata')


# Calculate Beta across entire dataset

soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site,Plot,OTU) 

hyph_Site_OTU <- myco_dat_hyph %>%
  distinct(Site,OTU) 

##############################
#######PA MODEL###############
##########################


# Join and merge with Beta_estimates
Beta_PA_Site_soil<- soil_Site_OTU %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  left_join(MER_aridity)

# Plot for both
ggplot(Beta_PA_Site_soil, aes(x = Site, y = Beta)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  geom_jitter( alpha = 0.9, size = 1,     position = position_dodge(width = 0.6)  # must match boxplot dodge
  ) +
  labs(y = "Beta values per OTU", x = "Site", title = "Beta values by Site") +
  theme_minimal() + 
  theme(     axis.text.x = element_text(angle = 45, hjust = 1))


#PA model
model <- lm(Beta ~ Site, data = Beta_PA_Site_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model)

smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_site <- tibble(
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

anova_tbl_site

smry_tbl_site<-  smry$coefficients%>% 
  as.data.frame() %>% 
  rownames_to_column(var='Coefficent') %>% 
  mutate(Coefficent= str_remove(Coefficent,'site'),
         Term= 'Site')
smry_tbl_site

site_max <- Beta_PA_Site_soil %>%
  group_by(Site) %>%
  summarise(y.position = max(Beta, na.rm = TRUE)+0.1)

# Tukey pairwise contrasts: Type within each site
contrasts_type_within_site <- emmeans(model, pairwise ~ Site, adjust = "sidak")

post_hoc_site<-contrasts_type_within_site$contrasts %>% 
  as.data.frame() %>%
  filter(p.value < 0.10) 







# Extract significant results
bracket_df<- contrasts_type_within_site$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.11) %>%
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
   left_join(site_max, by=c('xmin'= 'Site'))

 
 # Extract significance info from model summary
 site_summary <- summary(model)$coefficients %>%
   as.data.frame() %>%
   rownames_to_column("Site") %>%
   rename(p.value = `Pr(>|t|)`)%>%
   mutate(
     Site=str_remove(Site, "Site"),
     p.label = case_when(
       p.value <= 0.001 ~ "***",
       p.value <= 0.01  ~ "**",
       p.value <= 0.05  ~ "*",
       p.value <= 0.1   ~ paste0("p=", round(p.value, 3)),
       TRUE             ~ ""
     )
   )
site_summary 




sig_results <- site_summary %>%
  left_join(site_max, by = c("Site"="Site") )

sig_results


anova_p <- Anova(model)$`Pr(>F)`[1]  # First term in the model
# Get significance level
p_sig <- case_when(
  anova_p <= 0.0001 ~ "p <0.0001",
  anova_p <= 0.001 ~ "p <0.001",
  anova_p <= 0.002 ~ "p <0.002",
  anova_p <= 0.05 ~ "p <0.05",
  anova_p <= 0.1 ~ "+",
  TRUE ~ "ns"
)

# Create annotation label
model_label <- paste0( p_sig)
#### Order sites by aridity
site_levels_aridity <- Beta_PA_Site_soil %>% 
  distinct(Site, aridity) %>% 
  arrange(aridity) %>%         # change to arrange(desc(aridity)) if needed
  pull(Site)

#### Create axis labels
site_labels_named <- c(
  "NSMSEQCasWoo01" = "Forest<br>*Melaleuca/Casuarina*",
  "QDMSEQRainfo01" = "Subtropical-Rainforest",
  "SAMMDDMallee01" = "Triodia-Mallee<br>*Eucalyptus-dumosa*",
  "SAMMDDMallee02" = "Chenopod-Mallee<br>*Eucalyptus-oleosa*",
  "VCMSECEucFor01" = "Forest<br>*Eucalyptus-muelleriana*",
  "VCMSECRainfo01" = "Temperate-Rainforest",
  "WAMDALEucWoo01" = "Savannah<br>*Corymbia-paractia*",
  "WAMDALRainfo01" = "VineThicket",
  "WAMDALEucWoo02" = "Savannah<br>*Triodia-schinzii*",
  "WAMCOOShrubl01" = "Shrubland<br>*Allocasuarina-acutivalvis*",
  "WAMCOOEucWoo01" = "Woodland<br>*Eucalyptus-salmonophloia*",
  "WAMCOOEucWoo02" = "Chenopod-Mallee<br>*Eucalyptus-aff.-oleosa*",
  "WAMCOOEucWoo03" = "Woodland<br>*Eucalyptus-salubris*",
  "WAMESPShrubl01" = "Shrubland<br>*Kingia-australis*"
)

site_levels <- names(site_labels_named)


### Emmeans ordered by aridity
emm_df <- emmeans(model, ~ Site) %>%
  as.data.frame() %>%
  mutate(Site = factor(Site, levels = site_levels_aridity))

# ---------------------------
# Add CLD letters (compact letter display)
# ---------------------------

# compute emmeans object (re-usable)
emm_obj <- emmeans(model, ~ Site)
library(multcomp)
# get CLD letters (adjust method to match your pairwise adjust; here we use 'sidak' to match contrasts)
# cld() sometimes returns columns named '.group', '.cld', 'Letters', or '.letters' depending on versions/options.
cld_raw <- cld(emm_obj, adjust = "sidak", sort = FALSE, alpha = 0.05, Letters = letters)

cld_df <- as.data.frame(cld_raw) %>%
  mutate(
    letter = coalesce(as.character(.group) %>% na_if("")),
    letter = stringr::str_trim(stringr::str_remove_all(letter, "[\\(\\)\\[\\]\\s]+")),
    Site = factor(Site, levels = site_levels_aridity)
  )

emm_with_letters <- site_max %>%
  mutate(Site = factor(Site, levels = site_levels_aridity)) %>%
  left_join(cld_df %>% dplyr::select(Site, letter), by = "Site") %>%
  mutate(y.position = y.position * 0.90)



# Main boxplot
site<-Beta_PA_Site_soil %>% 
  mutate(Site = factor(Site, levels = site_levels)) %>% 
  ggplot( aes(x = Site, y =Beta, fill=Site)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  
  geom_violin(alpha = 0.7, width = 0.6, linewidth = 1) +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x= Site,y = emmean), shape = 21, size = 10, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x= Site, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  
  labs(y = expression(paste('Effect size (', beta[fire], ')')),
       x = "Site") +
  theme_classic() +
  theme(    legend.position = "none",
            axis.text.x = element_text(hjust = 0.5,vjust=.2, size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 18,vjust=0.5, face = "bold"),
            axis.title.y = element_text(size = 18, face = "bold"),
            axis.line = element_line(linewidth = 1.5, colour = "black"),
            plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
  )+
  scale_x_discrete(labels = site_labels_named)+
  # geom_bracket(
  #   data = bracket_df,
  #   aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
  #   inherit.aes = FALSE,
  #   tip.length = 0.01,
  #   label.size = 5,
  #   size = 1.1
  # )+
  # geom_text(data = sig_explo,
  #           aes(x = site_type, y = y.position, label = p.label),
  #           inherit.aes = FALSE,
  #           fontface= 'bold',
  #           size = 6)+
  annotate("text", x = 0.2, y = max(Beta_PA_Site_soil$Beta+0.05, na.rm = TRUE),
           label = model_label, hjust = -.3, vjust = 1, size = 7)


library(ggtext)


# Main boxplot with CLD letters
site <- Beta_PA_Site_soil %>%
  mutate(Site = factor(Site, levels = site_levels_aridity)) %>%
  ggplot(aes(x = Site, y = Beta, fill = Site)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 1.5, alpha = 0.6) +
  geom_violin(alpha = 0.7, width = 0.6, linewidth = 1) +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x = Site, y = emmean),
             shape = 21, size = 10, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x = Site, ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 1.3, inherit.aes = FALSE) +
  # CLD letters (from emm_with_letters)
  geom_text(data = emm_with_letters, aes(x = Site, y = y.position, label = letter),
            inherit.aes = FALSE, size = 8, vjust = 0) +
  labs(y = expression(paste('Likelihood of occurence post-fire (', beta[fire], ')')), x = "Site", tag='a)') +
  theme_classic() +
  scale_fill_grey(start = .9, end = 0.1) +
  theme(legend.position = "none",
        axis.text.x = element_markdown(angle = 35, hjust = 1, vjust = 1, face = "bold", size = 14),
        axis.title.x = element_text(size = 18, vjust = 0.5, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.line = element_line(linewidth = 1.5, colour = "black"),
        plot.margin = margin(t = 10, r = 10, b = 40, l = 10)) +
  scale_x_discrete(labels = site_labels_named, guide = guide_axis(n.dodge=1)) +
  annotate("text", x = 0.2, y = max(Beta_PA_Site_soil$Beta + 0.2, na.rm = TRUE),
           label = model_label, hjust = -.3, vjust = 1, size = 8)

site
save(anova_tbl_site,smry_tbl_site,post_hoc_site, file = 'HMSC_MER/Output/Processed_Data/site_summary.RDS')

ggsave("plots/Effect_Size_site_Plot.png", site, width = 15, height = 8, dpi = 300)
