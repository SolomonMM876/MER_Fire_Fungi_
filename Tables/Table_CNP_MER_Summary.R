library(flextable)
library(officer)
library(dplyr)
library(tidyr)
library(stringr)

CNP_MER<-readRDS('Processed_data/Stoich/CNP_MER.RDS')



# Function to calculate mean ± SE and format nicely
mean_se <- function(x) {
  m <- mean(x, na.rm = TRUE)
  se <- sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
  return(sprintf("%.1f ± %.1f", m, se))
}

# Create summary table
CNP_summary <- CNP_MER %>%
  group_by(Fire_Treatment) %>%
  summarise(
    `% C` = mean_se(C),
    `% N` = mean_se(N),
    `% P` = mean_se(P),
    `C:N` = mean_se(CN),
    `C:P` = mean_se(CP),
    `N:P` = mean_se(NP),
    .groups = "drop"
  ) %>%
  # Add overall mean ± SE
  bind_rows(
    CNP_MER %>%
      summarise(
        `% C` = mean_se(C),
        `% N` = mean_se(N),
        `% P` = mean_se(P),
        `C:N` = mean_se(CN),
        `C:P` = mean_se(CP),
        `N:P` = mean_se(NP),
      ) %>%
      mutate(Fire_Treatment = "Overall")
  ) %>%
  # Reorder rows: Overall first, then B, then U
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("Overall", "B", "U"))) %>%
  arrange(Fire_Treatment)

# View the table
CNP_summary

# Rename column for nicer display
CNP_summary <- CNP_summary %>%
  rename("Treatment" = Fire_Treatment)

# Build flextable
ft_summary <- CNP_summary %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Summary of soil nutrients (mean ± SE) for burnt, unburnt, and overall samples") %>%
  align(align = "center", part = "all") %>%
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>% 
  body_add_flextable(ft_summary)

# Export table to Word
print(doc, target = "Tables/Output/Table_Myc_CNP_Summary.docx")
