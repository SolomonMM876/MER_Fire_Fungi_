library(tibble)
library(dplyr)

# First table: sites
sites <- tribble(
  ~Site,                ~Site_nickname,        ~Veg_type,
  "NSMSEQCasWoo01",     "NSW-Forest-Flood",    "Casuarina Swamp Oak Forest.",
  "QDMSEQRainfo01",     "QLD-Rainfo-Subtrop",  "Tropical/subtropical rainforest",
  "SAMMDDMallee01",     "SA-Mallee-Spinifex",  "Shrubby/Triodia mallee",
  "SAMMDDMallee02",     "SA-Mallee-Chenopod",  "Tussock grass/chenopod mallee",
  "VCMSECEucFor01",     "VIC-Forest-Euc",      "Resprouter eucalypt forest",
  "VCMSECRainfo01",     "VIC-Rainfo-Temperate","Warm temperate rainforest",
  "WAMDALEucWoo01",     "WA-Savanna-Corymbia", "Wet-dry tropical eucalypt woodlands",
  "WAMDALRainfo01",     "WA-Vine-Thicket",     "Clittoral rainforest and coastal vine thickets",
  "WAMDALEucWoo02",     "WA-Savanna-Minyjuru", "Wet-dry tropical eucalypt woodlands",
  "WAMCOOShrubl01",     "WA-Sandplain-AlloCas","under development",
  "WAMCOOEucWoo01",     "WA-Woodland-SalmGum", "obligate-seeder eucalypt woodland",
  "WAMCOOEucWoo02",     "WA-Woodland-Oleosa",  "resprouter eucalypt woodland/chenopod mallee",
  "WAMCOOEucWoo03",     "WA-Woodland-Gimlet",  "obligate-seeder eucalypt woodland",
  "WAMESPShrubl01",     "WA-Shrub-Kingia",     "under development"
)

# Second table: floristic types
floristics <- tribble(
  ~Floristic_type,       ~Understorey_characteristics,                                         ~Rainfall_mm, ~Soils,                                             ~Fire,
  "Maalok",              "obligate-seeder sparse or absent with dense litter",                ">300",      "heavy clays or calcareous sands",                  "infrequent, stand-replacing",
  "Mesic heathy mallee", "dense sclerophyll shrubs from temperate genera",                     ">300",      "duplex with sandy topsoil, coarse sands",          "frequent, fire-age dependent",
  "Shrubby mallee",      "dominated by shrubs from arid plant genera",                         "<350",      "loams",                                             "moderate, ephemeral driven",
  "Triodia mallee",      "Triodia dominant with shrubs from arid plant genera",                "<350",      "light sands to sandy loams",                       "moderate, ephemeral driven",
  "Tussock grass mallee","tussock grasses dominant, with semi-succulent chenopods",            ">300",      "heavier soils",                                    "infrequent crown fire, poorly known",
  "Chenopod mallee",     "semi-succulent chenopods dominant, with tussock grasses",            "<350",      "heavier, alkaline/saline soils",                   "infrequent, ephemeral driven"
)

# Manual match table: link Veg_type to Floristic_type
veg_to_floristic <- tribble(
  ~Veg_type,                               ~Floristic_type,
  "Shrubby/Triodia mallee",                "Triodia mallee",
  "Tussock grass/chenopod mallee",         "Chenopod mallee",
  "obligate-seeder eucalypt woodland",     "Maalok",
  "resprouter eucalypt woodland/chenopod mallee", "Chenopod mallee"
  # all others will be NA for now
)

# Join to get Floristic_type
sites_with_flor <- sites %>%
  left_join(veg_to_floristic, by = "Veg_type") %>%
  left_join(floristics, by = "Floristic_type") %>%
  mutate(source = if_else(!is.na(Floristic_type), "(Prober et al. 2023)", NA_character_))

sites_with_flor

library(writexl)

write_xlsx(sites_with_flor,'Processed_data/metadata/veg/MER_Site_description.xlsx')
library(readxl)

site_flor<-read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx")

library(tidyverse)
library(flextable)
library(officer)



# Build flextable ANOVA
ft_AusEcoModel<- site_flor %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Table X. Description of sites, vegitation type and expected fire frequency based on... Site_nickname is based on Aus state aberviation- broad ecosystem type- specific vegitation found at the site. 
                   Expected fire interval is based on relevant archetype model") %>%
  align(align = "left", part = "all") %>%
  set_table_properties(layout = "autofit")

# View in RStudio Viewer
ft_AusEcoModel





# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_AusEcoModel) %>% 
  body_add_break() 

print(doc, target = "Tables/Output/Table_AusEcoModels_updated.docx")
