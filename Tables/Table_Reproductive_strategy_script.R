library(tidyverse)
library(flextable)
library(officer)



final_repo_strat<-readRDS('Processed_data/Seq_dat/datasets_external/repro_strat_df.Rds') %>% 
  select(family,genus,repo_strategy,notes) %>% 
  mutate(across(family:repo_strategy, ~ ifelse(is.na(.), "NA", .))
  )






# Build flextable
ft_summary <- final_repo_strat %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Fungal families and genera with reproductive strategy, with justification.") %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "autofit") 

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_reproductive_strategy.docx")
