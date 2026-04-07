library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(ggpubr)
library(emmeans)

#load Beta1 and 2 from soil data HMSC
#load('HMSC_MER/results/mcmc_output.RData')

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')


#read in veg type assignments
library(readxl)
MER_Site_description <- read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx",
                                   sheet = "SP suggestions (highlighted)")

Mer_veg<-MER_Site_description %>% 
  select(Site,Veg_type=Veg_type_alt_2)

#load most up to date PA model
load('HMSC_MER/results/Beta1.RData')
load('HMSC_MER/results/Beta1_host_freq.RData')

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean

Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)


# Calculate Beta across entire dataset

soil_Site_OTU <- myco_dat_soil %>%
  distinct(Site,Plot,OTU) %>% 
  left_join(Mer_veg)

# hyph_Site_OTU <- myco_dat_hyph %>%
#   distinct(Site,OTU) 

##############################
#######PA MODEL###############
##########################


# Join and merge with Beta_estimates
Beta_PA_Veg_type_soil<- soil_Site_OTU %>%
  #filter(OTU %in% myco_tax_hyph$OTU) %>% 
  left_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(    Type= "Soil")

# Plot for both
ggplot(Beta_PA_Veg_type_soil, aes(x = Veg_type, y = Beta)) +

  geom_jitter( aes(color=Site), alpha = 0.9, size = 3,     width = 0.3 # must match boxplot dodge
  ) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width=.3,     position = position_dodge(width = 0.6)  # wider dodge increases spacing
  ) +
  labs(y = "Beta values per OTU", x = "Veg_type", title = "Beta values by Veg_type") +
  theme_minimal() + 
  theme(     axis.text.x = element_text(angle = 45, hjust = 1))


#PA model#PA modelSite
model <- lm(Beta ~ Veg_type, data = Beta_PA_Veg_type_soil)

# Model performance summary
model_performance(model)
# Check model assumptions visually
check_model(model)
summary(model)
Anova(model)

smry <- summary(model)
anov <- Anova(model, type = "II")

anova_tbl_Veg_type <- tibble(
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

anova_tbl_Veg_type

smry_tbl_Veg_type<-  smry$coefficients%>% 
  as.data.frame() %>% 
  rownames_to_column(var='Coefficent') %>% 
  mutate(Coefficent= str_remove(Coefficent,'Veg_type'),
         Term= 'Veg_type')
smry_tbl_Veg_type


# Tukey pairwise contrasts: Type within each Veg_type
contrasts_type_within_Veg_type <- emmeans(model, pairwise ~ Veg_type, adjust = "sidak")

# Extract significant results
sig_labels<- contrasts_type_within_Veg_type$contrasts %>%
  as.data.frame() %>%
  filter(p.value < 0.11)%>%
   mutate(
     p.label = case_when(
       p.value <= 0.001 ~ "***",
       p.value <= 0.01 ~ "**",
       p.value <= 0.05 ~ "*",
       p.value <= 0.1 ~ paste0('p=',round(p.value,3))
     ),
     xmin = str_trim(str_extract(contrast, "(^[^-]+)")),         # before "-"
     xmax = str_trim(ifelse(
       str_detect(contrast, "\\(.*\\)"),
       str_extract(contrast, "(?<=\\().*?(?=\\))"),   # extract inside parentheses
       str_extract(contrast, "(?<=-)\\s*.*$")         # fallback: extract after "-"
     )))
 
 post_hoc_veg<-contrasts_type_within_Veg_type$contrasts %>% 
   as.data.frame() %>% 
   mutate(Variable='Veg_type')
 
 post_hoc_veg
 
 veg_type_max <- Beta_PA_Veg_type_soil %>%
   group_by(Veg_type) %>%
   summarise(y.position = max(Beta, na.rm = TRUE)+0.05)
 
 sig_results <- sig_labels %>%
   left_join(veg_type_max, by = c("xmax"="Veg_type") )
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
#### create axis labels

# Get emmeans with 95% CIs
emm_df <- emmeans(model, ~ Veg_type) %>%
  as.data.frame()


# Main boxplot
Veg_type<-Beta_PA_Veg_type_soil %>% 
  ggplot( aes(x = Veg_type, y =Beta, fill=Veg_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 2, alpha = 0.6) +
  
  geom_violin(alpha = 0.7, width = 0.4, linewidth = 2) +
  # Add estimated means and CIs
  geom_point(data = emm_df, aes(x= Veg_type,y = emmean), shape = 21, size = 10, fill = "white", color = "black", stroke = 1.5, inherit.aes = FALSE) +
  geom_errorbar(data = emm_df, aes(x= Veg_type, ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 2, inherit.aes = FALSE) +
  
  labs(y = "Effect size (β_fire)", x = "Veg_type") +
  theme_classic() +
  theme(    legend.position = "none",
            axis.text.x = element_text(hjust = 0.5,vjust=.2, size = 14, face = "bold"),
            axis.text.y = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(size = 14,vjust=0.5, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.line = element_line(linewidth = 2, colour = "black"),
            plot.margin = margin(t = 10, r = 10, b = 40, l = 10)
  )+
  geom_bracket(
    data = sig_results,
    aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.label),
    inherit.aes = FALSE,
    tip.length = 0.01,
    label.size = 10,
    size = 1.1
  )+
  annotate("text", x = 0.2, y = max(Beta_PA_Veg_type_soil$Beta+0.05, na.rm = TRUE),
           label = model_label, hjust = -.3, vjust = 1, size = 10)

Veg_type
save(anova_tbl_Veg_type,smry_tbl_Veg_type, file = 'HMSC_MER/Output/Processed_Data/Veg_type_summary.RDS')

ggsave("plots/Effect_Size_Veg_type_Plot.png", Veg_type, width = 40, height = 20, dpi = 300)
