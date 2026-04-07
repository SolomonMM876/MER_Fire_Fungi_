library(tidyverse)
library(flextable)
library(officer)



load('HMSC_MER/Output/Processed_Data/Effect_Size_RA_Beta.RDS')


RA_Beta_ANOVA <- anova_table_all %>% 
  mutate(
    across(sum_sq_resid:adj_R2, ~ round(.x, 2)),
    across(Estimate:sum_sq, ~ round(.x, 3)),
    `Sum of squares` = paste0(sum_sq, ", ", sum_sq_resid),
    Df = paste0(DF_num, ", ", DF_denom)
  ) %>%
  dplyr::select(
    Type,
    Estimate,
    `Std error` = Std_Error,
    `Sum of squares`,
    Df,
    F,
    P,
    R2
  )
  



# Build flextable
ft_summary <- RA_Beta_ANOVA %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. ANOVA of log10 ratio of OTU enrichment related to effect size (β_fire )") %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_Effect_Size_RA.docx")
