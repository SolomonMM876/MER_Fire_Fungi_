library(tidyverse)
library(flextable)
library(officer)


#from myco_host_script
Plots_not_sampled<-readRDS('Processed_data/metadata/veg/Plots_w_no_veg_data.rds') %>% 
  mutate(Type='Vegetation') %>% 
  rename(Burnt=B, Unburnt=U)

#from missing samples script
load('Processed_data/Seq_dat/Plots_w_missing_data.rds')

#### create site labels
site_labels_named <- c(
  "NSMSEQCasWoo01"  = "NSW-Forest-Flood",
  "QDMSEQRainfo01"  = "QLD-Rainfo-Subtrop",
  "SAMMDDMallee01"  = "SA-Mallee-Spinifex",
  "SAMMDDMallee02"  = "SA-Mallee-Chenopod",
  "VCMSECEucFor01"  = "VIC-Forest-Euc",
  "VCMSECRainfo01"  = "VIC-Rainfo-Temperate",
  "WAMDALEucWoo01"  = "WA-Savanna-Corymbia",
  "WAMDALRainfo01"  = "WA-Vine-Thicket",
  "WAMDALEucWoo02"  = "WA-Savanna-Minyjuru",
  "WAMCOOShrubl01"  = "WA-Sandplain-AlloCas",
  "WAMCOOEucWoo01"  = "WA-Wood-SalmGum",
  "WAMCOOEucWoo02"  = "WA-Wood-Oleosa",
  "WAMCOOEucWoo03"  = "WA-Wood-Gimlet",
  "WAMESPShrubl01"  = "WA-Shrub-Kingia"
)

summary_table <- bind_rows(Plots_not_sampled, table_soil_missing, table_hyph_missing) %>%
  mutate(Type = factor(Type, levels = c("Vegetation", "Soil", "In-growth bags"))) %>%
  select(Type, everything()) %>%
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, sort(levels(Site)))  # reorder alphabetically
  ) %>% 
  arrange(Type, Site)
  

# View result
print(summary_table)



# Build flextable
ft_summary <- summary_table %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Table SX.Number of OTUs (operational taxonomic units) at each posterior probability") %>%
  hline(i=c(4,10)) %>% 
  align(align = "center", part = "all") %>%
  set_table_properties(layout = "autofit") 

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_missing_plots.docx")
