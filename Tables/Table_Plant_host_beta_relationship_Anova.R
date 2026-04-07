library(tidyverse)
library(flextable)
library(officer)

#load Anova table and interaction smmry
load('HMSC_MER/Output/Processed_Data/Myco_veg_host_ANOVA_output.RDS')


host_beta_ANOVA<-anova_table%>% 
  mutate(across(Estimate:R2_marginal, ~ round(.x, 2)),
         Df_res=round(Df_res,0),
         P = case_when(
           P <= 0.0001 ~ "<0.0001",
           P <= 0.001 ~ "<0.001",
           P <= 0.01  ~ "<0.01",
           P <= 0.05  ~ "<0.05",
           TRUE       ~ as.character(round(P, 3))),
          Df = paste0(Df, ", ", Df_res))%>% 
  dplyr::select(Factor,Estimate,Std_Error, Df,F,P,`R² marginal`=R2_marginal)



# Build flextable
ft_host_ANOVA <- host_beta_ANOVA %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of relationship between host frequency and Beta") %>%
  align(align = "left", part = "all") %>%
 # hline(i=c(4,8)) %>% 
  set_table_properties(layout = "autofit") 

# View in RStudio Viewer
ft_host_ANOVA

# 
# host_beta_smry_int<-smry_table%>% 
#   mutate(across(Estimate:P, ~ round(.x, 3)),
#          P = case_when(
#            P <= 0.0001 ~ "<0.0001",
#            P <= 0.001 ~ "<0.001",
#            P <= 0.01  ~ "<0.01",
#            P <= 0.05  ~ "<0.05",
#            TRUE       ~ as.character(round(P, 3))
#          ) )%>% 
#   select(Metric,Term,Estimate,`Std. Error`,`t value`,P)
# VarCorr()
# 
# 
# # Build flextable
# ft_host_int_smry <- host_beta_smry_int %>%
#   flextable() %>%
#   theme_booktabs() %>%
#   add_header_lines("Supp Table X. Interaction summary of relationship between host frequency and Beta") %>%
#   align(align = "left", part = "all") %>%
#   hline(i=c(11)) %>% 
#   set_table_properties(layout = "autofit") 
# 
# # View in RStudio Viewer
# ft_host_int_smry
# 

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_host_ANOVA) %>% 
  body_add_break() 
  #body_add_flextable(ft_host_int_smry)

print(doc, target = "Tables/Output/Table_host_freq_vs_Beta_ANOVA.docx")
