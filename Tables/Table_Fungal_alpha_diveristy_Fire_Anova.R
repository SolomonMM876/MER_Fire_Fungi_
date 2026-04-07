library(tidyverse)
library(flextable)
library(officer)

#load alpha diversity for all soil data
load('HMSC_MER/Output/Processed_Data/soil_myco_alpha_metrics.RDS')

#load for bag data
load('HMSC_MER/Output/Processed_Data/bag_myco_alpha_metrics.RDS')



anova_table_myco<-bind_rows(anova_table_formatted_soil ,anova_table_formatted_bag)

alpha_ANOVA_myco<-anova_table_myco%>% 
  mutate(across(Estimate:P, ~ round(.x, 2)),
         DF_denom=round(DF_denom,0),
         #`Sum of squares` = paste0(sum_sq, ", ", sum_sq_resid),
         Df = paste0(DF_num, ", ", DF_denom)) %>% 
  dplyr::select(Type,Metric,Estimate,Std_Error,Df,F,P)
  





# Build flextable
ft_summary <- alpha_ANOVA_myco %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of mycorrhizal fungal alpha metrics response to fire") %>%
  align(align = "left", part = "all") %>%
  hline(i=c(4,8)) %>% 
  set_table_properties(layout = "autofit") %>% 
  compose(j = "P", 
          value = as_paragraph(
            as_chunk(sprintf("%.2f", P), props = fp_text(bold = TRUE))
          ),
          i = ~ P < 0.05)  # Only bold significant P values

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_alpha_ANOVA_myco_fungi_fire.docx")
