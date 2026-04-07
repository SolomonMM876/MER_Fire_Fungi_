library(tidyverse)
library(flextable)
library(officer)

#create table 
#this is all the familys and genera from hyphal and soil samples with probable mycorrhizal status 
Ecto_included <- tribble(
  ~family,            ~genera,       ~Included_Excluded, ~Reason,
  "Thelephoraceae",   NA,            "Included",         "Many are likely to be ectomycorrhizal (Cannon & Kirk, 2007)",
  "Thelephoraceae",   "Thelephora",  "Included",         "Known ectomycorrhizal genus (https://www.cabi.org/isc/datasheet/53570)",
  "Myxotrichaceae",   "Myxotrichum", "Included",         "Dalpé Y. 1989 shown to form ericoid mycorrhizas",
  "Myxotrichaceae",   NA,            "Excluded",         "Some are ectomycorrhizal, but most are not (Thiem et al., 2018)",
  "Tricholomataceae", NA,            "Included",         "FunGuild notes suggest most species are ectomycorrhizal",
  "Tricholomataceae", "Bonomyces",   "Included",         "Most species are mycorrhizal in this family (https://mycorrhizas.info/ecmf.html)",
  "Amanitaceae",      "Amanita",     "Included",         "Most are ectomycorrhizal (Yang et al., 1999)",
  "Leotiaceae",       "Pezoloma",    "Included",         "(Midgley et al. 2017) ericaceous fungi found on roots",
  "Hymenochaetaceae", NA,            "Excluded",         "Lots are not ectomycorrhizal (Salvador-Montoya et al., 2022)",
  "Pluteaceae",       NA,            "Included",         "Forming ectomycorrhizas with roots of broadleaved trees, or saprobic on rotten wood, plant remains or humus (Cannon & Kirk 2007)",
  "Helotiaceae",      "Neocrinula",  "Excluded",         "Usually saprobic on herbaceous or woody tissues, some species fungicolous. A few species are known to form mycorrhizae (Cannon & Kirk 2007)",
  "Helotiaceae",      NA,            "Excluded",         "Lack of information, FunGuild notes suggest unlikely ectomycorrhizal",
  "Hysterangiaceae",  NA,            "Included",         "Known ectomycorrhizal family (Weber et al. 2013, Nuske et al. 2019)",
  "Boletaceae",       NA,            "Included",         "FunGuild notes suggest most species are ectomycorrhizal"
) %>% arrange(family,genera)


saveRDS(Ecto_included,'Processed_data/Seq_dat/Ecto_Included.rds')

# Clean and format
ecto_clean <- Ecto_included %>%
  rename(genus = genera) %>%
  mutate(genus = if_else(is.na(genus), "NA", genus)) %>%
  arrange(family, desc(genus == "NA"), genus) %>%  # Sort within family: known genus first, then "NA"
  group_by(family) %>%
  mutate(family = ifelse(row_number() == 1, family, "")) %>%
  ungroup()

# Create flextable
ft_ecto <- ecto_clean %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Fungal families and genera included or excluded as ectomycorrhizal taxa, with justification.") %>%
  compose(j = "genus", value = as_paragraph(as_i(genus))) %>%  # Italicize genus
  autofit() %>%
  align(align = "left", part = "all")
ft_ecto

# Export to Word
doc <- read_docx() %>%
  body_add_flextable(ft_ecto)

print(doc, target = "Tables/Output/Ecto_Included_Table.docx")
