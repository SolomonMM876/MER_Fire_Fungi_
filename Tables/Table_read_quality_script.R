library(tidyverse)
library(flextable)
library(officer)
library(scales)  # for comma()
# Format numbers with commas inside "X (Y-Z)" strings or just plain numbers
add_commas_to_string_numbers <- function(x) {
  str_replace_all(
    x,
    "(\\d{1,3})(?=(\\d{3})+(?!\\d))",  # Regex for thousands separator
    "\\1,"
  )
}



summary_column_2nd_myc<-read_rds('Processed_data/seq_summary/2nd_myc_summary.rds')


summary_column_1st_myc<-read_rds('Processed_data/seq_summary/1st_myc_summary.rds')

summary_combined <- summary_column_1st_myc %>%
  rename(`1st_hyph` = First_myc) %>%
  left_join(summary_column_2nd_myc %>% rename(`2nd_hyph` = Second_myc), by = "metric") %>%
  #left_join(summary_column_3rd_myc %>% rename(`soil_1st` = value), by = "metric") %>%
  rename(Metric = metric)

summary_combined_fmt <- summary_combined %>%
  mutate(across(c(`1st_hyph`, `2nd_hyph`), ~ add_commas_to_string_numbers(.)))


# Build flextable
ft_summary <- summary_combined_fmt %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Table X. An overview of the sequencing results, including fungal reads per sample and proportion of reads assigned to mycorrhizal taxa. Values are the total number of reads retained after each step.") %>%
  align(align = "left", part = "all") %>%
  autofit()

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_seq_summary.docx")

