library(flextable)
library(tidyverse)
library(officer)
library(readxl)

#Read in Suzzanne Prober  AusEcoModel spreadsheet
MER_Site_description<-read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx",
                                 sheet = "SP suggestions (highlighted)")
#read in site precip data
precip_site<-readRDS("Processed_data/metadata/precip_MER.RDS")

#read in site pH data
MER_pH<-readRDS("processed_data/nutrients/pH_MER.rds") %>% 
  group_by(Site) %>% 
  summarise(pH=round(mean(Soil_pH, na.rm=TRUE),1))

#Read in soil details
CSBP_Lab_MER_Soil <- read_excel("raw_data/Soil_data/CSBP_Lab_250130_Solomon McMahan_001.xls", 
                                               sheet = "Results")
# make a named vector as a lookup
texture_lookup <- c(
  "1"   = "Sand",
  "1.5" = "Sand/Loam",
  "2"   = "Loam",
  "2.5" = "Loam/Clay",
  "3"   = "Clay",
  "3.5" = "Heavy Clay"
)

Soil_text <- CSBP_Lab_MER_Soil %>%
  mutate(
    Texture_num = as.numeric(as.character(Texture)),
    Texture_cat = recode(as.character(Texture_num), !!!texture_lookup)
  ) %>% 
  dplyr::select(Site=`Sample Name 2`,Soil_type=Texture_cat)
  



AusEcoModel<-MER_Site_description%>% 
  left_join(precip_site) %>% #precip data
  left_join(MER_pH) %>% 
  left_join(Soil_text) %>% 
  dplyr::select(Site=Site_nickname,`Floristic type`=Floristic_type,`Broad veg type`=Veg_type_alt,`Mean Rainfall per year`=mean_annual_precip,
                pH,Soil_type,`Most recent fire severity`= Last_fire_severity,
                `Fire response of dominant veg`=Fire_response,Source=source_freq) 
# Build flextable ANOVA
ft_AusEcoModel<- AusEcoModel %>%
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

print(doc, target = "Tables/Output/Table_AusEcoModels.docx")
