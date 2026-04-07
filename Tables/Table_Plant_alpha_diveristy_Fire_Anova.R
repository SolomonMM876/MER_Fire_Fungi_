library(tidyverse)
library(flextable)
library(officer)

#load alpha diversity for all plant spp
load('HMSC_MER/Output/Processed_Data/Vegitation_alpha_metrics.RDS')

#load for just mycorrhizal species
load('HMSC_MER/Output/Processed_Data/Vegitation_myco_alpha_metrics.RDS')
anova_table_all_myco
anova_table_veg<-bind_rows(anova_table_all_veg %>% mutate(Group='All vegetation')
                             ,anova_table_all_myco)

alpha_ANOVA<-anova_table_veg%>% 
  mutate(across(Estimate:P, ~ round(.x, 2)),
         DF_denom=round(DF_denom,0)) %>% 
  select(Group,Metric,Estimate,Std_Error,DF_num,DF_denom,F,P)
  





# Build flextable
ft_summary <- alpha_ANOVA %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of alpha metrics response to fire") %>%
  align(align = "left", part = "all") %>%
  hline(i=c(4,8,12,16)) %>% 
  set_table_properties(layout = "autofit") %>% 
  compose(j = "P", 
          value = as_paragraph(
            as_chunk(sprintf("%.2f", P), props = fp_text(bold = TRUE, color = "black"))
          ),
          i = ~ P < 0.05)  # Only bold significant P values

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_alpha_ANOVA_veg_fire.docx")
