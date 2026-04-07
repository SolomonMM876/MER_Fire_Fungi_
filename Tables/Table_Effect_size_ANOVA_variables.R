library(tidyverse)
library(flextable)
library(officer)


#load phylum summary stats
load( 'HMSC_MER/Output/Processed_Data/Phylum_summary.RDS')

#load repo strat summary stats
load('HMSC_MER/Output/Processed_Data/repo_strat_summary.RDS')

#load exploration summary stats
load('HMSC_MER/Output/Processed_Data/explo_summary.RDS')

#load AM hyphal production stats
load('HMSC_MER/Output/Processed_Data/Intra_Extra_radical_hyph_summary.RDS')

#load site summary stats
load('HMSC_MER/Output/Processed_Data/site_summary.RDS')

#load fire interval stats
load('HMSC_MER/Output/Processed_Data/Fire_interval_summary.RDS')

anova_results_site_vars<-anova_results_site_vars %>% 
  filter(Variable %in% c('severity','TSF_years','Veg_type','Fire_interval','Simple_fire_response_dominants'))
  



ANOVA_terms <- bind_rows(anova_tbl_phy, anova_tbl_explo,anova_tbl_Guild, anova_tbl_repo,anova_tbl_site) %>% 
  mutate(
    across(sum_sq:adj_R2, ~ round(.x, 2)),
    P = case_when(
      P <= 0.0001 ~ "<0.0001",
      P <= 0.001 ~ "<0.001",
      TRUE       ~ as.character(round(P, 3))),
      `Sum of squares` = paste0(sum_sq, ", ", sum_sq_resid),
      Df = paste0(DF_num, ", ", DF_denom)
    ) %>% 
  dplyr::select(Term,`Sum of squares`,Df,F,R2,adj_R2,P)





# Build flextable ANOVA
ft_ANOVA <- ANOVA_terms %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of OTU grouping variable on fire effect size (Beta)") %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_ANOVA


anova_results_site_vars

ANOVA_site_terms <- anova_results_site_vars %>% 
  mutate(
    across(`F_value`:R2_conditional, ~ round(.x, 2)),
    P = case_when(
      P <= 0.0001 ~ "<0.0001",
      P <= 0.001 ~ "<0.001",
      TRUE       ~ as.character(round(P, 3))),
    Df = paste0(DF_num, ", ", DF_denom)
  ) %>% 
  dplyr::select(Variable,Df,`F_value`,`R^2 marginal`=R2_marginal,`R^2 conditional`=R2_conditional,P)





# Build flextable ANOVA
ft_ANOVA_site <- ANOVA_site_terms %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of Site grouping charecteristics on fire effect size (Beta)") %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_ANOVA_site












#### create site labels
site_labels_named <- c(
  "NSMSEQCasWoo01"  = "NSW-Forest-Flood",
  "QDMSEQRainfo01"  = "QLD-Rainfo-Subtrop",
  "SAMMDDMallee01"  = "SA-Mallee-Spinifex",
  "SAMMDDMallee02"  = "SA-Mallee-Chenopod",
  "VCMSECEucFor01"  = "VIC-Forest-Euc",
  "VCMSECRainfo01"  = "VIC-Rainfo-Temperate",
  "WAMDALEucWoo01"  = "WA-Savanna-Corymbia",
  "WAMDALRainfo01"  = "WA-Vine-Thicket",
  "WAMDALEucWoo02"  = "WA-Savanna-Minyjuru",
  "WAMCOOShrubl01"  = "WA-Sandplain-AlloCas",
  "WAMCOOEucWoo01"  = "WA-Woodland-SalmGum",
  "WAMCOOEucWoo02"  = "WA-Woodland-Oleosa",
  "WAMCOOEucWoo03"  = "WA-Woodland-Gimlet",
  "WAMESPShrubl01"  = "WA-Shrub-Kingia"
)
t<-site_labels_named %>% as.data.frame()

# summary_stat_terms<-bind_rows(smry_tbl_phy %>% mutate(Term='Phylum') %>% rename(`Pr(>|t|)`=p.value),
#                               smry_tbl_explo%>% mutate(Term='Exploration type'),
#                               smry_tbl_Guild%>% mutate(Term='Guild') %>% rename(`Pr(>|t|)`=p.value),
#                               smry_tbl_repo%>% mutate(Term='Fruiting body type'),
#                               smry_tbl_site%>% mutate(Term='Site')) %>% 
#   dplyr::select(Term,Coefficient =Coefficent, Estimate, `Std. Error`,`t value`,P=`Pr(>|t|)`) %>% 
#   mutate(
#     across(Estimate:P, ~ round(.x, 2)),
#     P = case_when(
#       P <= 0.001 ~ "<0.001",
#       P <= 0.01  ~ "<0.01",
#       P <= 0.05  ~ "<0.05",
#       TRUE       ~ as.character(round(P, 3))
#     ),
#     Coefficient=str_replace(Coefficient,'_'," "),
#     Coefficient=str_remove(Coefficient, 'Site'),
#     Coefficient = str_replace_all(Coefficient, site_labels_named)
#   ) 
# 
# 
# 
# 
# # Build flextable ANOVA
# ft_summary <- summary_stat_terms %>%
#   flextable() %>%
#   theme_booktabs() %>%
#   add_header_lines("Supp Table X. The influence of X on fire effect size (Beta) ") %>%
#   align(align = "left", part = "all") %>%
#   hline(i=c(3,10,14,16,18)) %>% 
#   set_table_properties(layout = "autofit")
# 
# # View in RStudio Viewer
# ft_summary



post_hoc_stat<-bind_rows(post_hoc_phy,post_hoc_repo,post_hoc_explo %>% mutate(Variable='Exploration type'),post_hoc_Guild,post_hoc_site %>% mutate(Variable='Site'),post_hoc_veg) %>% 
  dplyr::select(Term=Variable,contrast, Estimate=estimate, `Std. Error`=SE,df,`t value`=t.ratio,P=p.value) %>% 
  mutate(
    across(Estimate:P, ~ round(.x, 2)),
    P = case_when(
      P <= 0.001 ~ "<0.001",
      P <= 0.01  ~ "<0.01",
      P <= 0.05  ~ as.character(round(P, 3)),
      TRUE       ~ as.character(round(P, 3))
    ),
    Term=str_replace(Term,'Reproduction strategy','Fruiting body type'),
    Term=str_replace(Term,'Extraradical_hyphae','Extraradical hyphal levels')
    
  )




# Build flextable ANOVA
ft_post_hoc <- post_hoc_stat %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Post hoc differences between groups with Sidak correction ") %>%
  align(align = "left", part = "all") %>%
  hline(i=c(3,9,12,15,60)) %>% 
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_post_hoc





# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_ANOVA) %>% 
  body_add_break() %>% 
  body_add_flextable(ft_ANOVA_site) %>% 
  body_add_break() %>% 
  body_add_flextable(ft_post_hoc)

print(doc, target = "Tables/Output/Table_ANOVA_Beta_Index_groupings.docx")
