library(tidyverse)
library(flextable)
library(officer)

# Load data
myco_tax_final <- readRDS("processed_data/metadata/veg/myco_host_df.rds")

# Create cleaned host table
myco_host <- myco_tax_final %>%
  select(Genus = genus,
         Family = family,
         Myco_type,
         Notes) %>%
  distinct() %>%
  mutate(
    # Add superscript "¹" when Notes match
    Myco_type = if_else(
      Notes == "Classification based on other genera in family",
      paste0(Myco_type, "\u00B9"),  # Unicode superscript 1
      Myco_type
    ),
    Notes=str_remove(Notes,"Classification based on other genera in family"),
    Myco_type = if_else(
      Notes == "Classification based Soudzilovskaia et al 2020 FungalRoot databse",
      paste0(Myco_type, "\u00B2"),  # Unicode superscript 2?
      Myco_type
    ),
    Notes=str_remove(Notes,"Classification based Soudzilovskaia et al 2020 FungalRoot databse")
  )

# Create flextable
ft_myco <- myco_host %>%
  arrange(Family,Genus) %>% 
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines(
    "Table S8. Assignment of mycorrhizal type to plant species based on genus-level classifications from Brundrett (2009). When genus-level data were unavailable, family-level classifications were used and noted."
  ) %>%
  compose(j = "Genus", value = as_paragraph(as_i(Genus))) %>%  # Italicize genus
  compose(j = "Family", value = as_paragraph(as_i(Family))) %>%  # Italicize family
  compose(j = "Myco_type", value = as_paragraph(as_chunk(Myco_type))) %>%
  autofit() %>%
  align(align = "left", part = "all")

ft_myco

# Export to Word
doc <- read_docx() %>%
  body_add_flextable(ft_myco)

print(doc, target = "Tables/Output/AM_ECM_host_ID_Table.docx")
