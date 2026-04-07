#### Site locations map ####

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(maps)
library(plotly)
library(readxl)
library(oz)

#--- Load site metadata ---#
MER_Site_description <- read_excel(
  "Processed_data/metadata/veg/MER_Site_description.xlsx",
  sheet = "SP suggestions (highlighted)"
)

All_Sites <- readRDS("raw_data/MER_Site_data/All_Sites.RDS")


#### Create axis labels
site_labels_named <- c(
  "NSMSEQCasWoo01" = "Forest\n(Melaleuca/Casuarina*)",
  "QDMSEQRainfo01" = "Subtropical-Rainforest",
  "SAMMDDMallee01" = "Triodia-Mallee\n(Eucalyptus-dumosa)",
  "SAMMDDMallee02" = "Chenopod-Mallee\n(Eucalyptus-oleosa)",
  "VCMSECEucFor01" = "Forest\n(Eucalyptus-muelleriana)",
  "VCMSECRainfo01" = "Temperate-Rainforest",
  "WAMDALEucWoo01" = "Savannah\n(Corymbia-paractia)",
  "WAMDALRainfo01" = "VineThicket",
  "WAMDALEucWoo02" = "Savannah\n(Sersalisia-sericea)",
  "WAMCOOShrubl01" = "Shrubland\n(Allocasuarina-acutivalvis)",
  "WAMCOOEucWoo01" = "Woodland\n(Eucalyptus-salmonophloia)",
  "WAMCOOEucWoo02" = "Chenopod-Mallee\n(Eucalyptus-affinis-oleosa)",
  "WAMCOOEucWoo03" = "Woodland\n(Eucalyptus-salubris)",
  "WAMESPShrubl01" = "Shrubland\n(Kingia australis)"
)

#--- Manual Nudge Assignments ---#
manual_nudges <- tribble(
  ~Site,             ~nudge_x, ~nudge_y,
  "NSMSEQCasWoo01",   13.0,     -3.0,     # Far East 
  "QDMSEQRainfo01",   12.0,     4.0,     # Far North-East
  "SAMMDDMallee01",   -3.0,      -14,     # Far North
  "SAMMDDMallee02",  -12.0,     -6.0,     # Far South
  "VCMSECEucFor01",   17.0,     -8.0,    # Far East
  "VCMSECRainfo01",   12.0,     -3.0,     # Far South-East
  "WAMDALEucWoo01",  -16.0,     10.0,     # Far North-West
  "WAMDALRainfo01",   4.0,      8.0,     # Far North
  "WAMDALEucWoo02",  -17.0,    0.0,     # Far South-West
  "WAMCOOShrubl01",  -18.0,     0.0,     # Far West
  "WAMCOOEucWoo01",  -21.0,      8.0,     # Inland NW
  "WAMCOOEucWoo02",   16.0,      3.0,     # Inland East
  "WAMCOOEucWoo03",   1.0,      10.0,     # Inland North
  "WAMESPShrubl01",   -10,     -7.0      # Far South
)




# 2. Create Map_Info with Correct Join Sequence
Map_Info <- All_Sites %>%
  group_by(Site) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(Site, Lat, Long) %>%
  # JOIN NUDGES FIRST (while Site matches the short codes)
  left_join(manual_nudges, by = "Site") %>%
  mutate(
    nudge_x = replace_na(nudge_x, 0),
    nudge_y = replace_na(nudge_y, 0)
  ) %>%
  # JOIN METADATA
  left_join(MER_Site_description %>% 
              rename(Broad_veg = Veg_type_alt), 
            by = "Site") %>% 
  # RENAME SITES LAST
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, sort(levels(Site))))

#--- Get Australia basemap ---#
australia_map <- map_data("world", region = "Australia")

#--- Static ggplot map with shapes ---#
p_map <- ggplot() +
  geom_polygon(data = australia_map, aes(x = long, y = lat, group = group),
               fill = "grey90", color = "black", linewidth = 1.2) +
  geom_point(data = All_Sites, aes(x = Long, y = Lat), 
             size = 6, shape=21, stroke=3, fill="white") +
  geom_text_repel(
    data = Map_Info, 
    aes(x = Long, y = Lat, label = Site),
    # Nudges are now active and valid
    nudge_x = Map_Info$nudge_x,
    nudge_y = Map_Info$nudge_y,
    size = 7,
    lineheight = 0.75,       # <--- Added: Reduces space between lines (Default is ~1.2)
    min.segment.length = 0,
    segment.size = 1.5,
    point.padding = 0.5,
    force = 0.5,
    max.overlaps = Inf, # Essential for manual placement
    color = "black",   
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15
  ) +
  coord_fixed(ratio = 1.1, xlim = c(95, 175), ylim = c(-50, -5)) + 
  theme_minimal(base_line_size = 1.5) +
  labs(x = "Latitude", y = "Longitude") +
  theme(
    axis.text = element_text(size=16, face='bold'),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.line = element_line(linewidth = 1.5, colour = "black")
  )
p_map
ggsave("plots/Site_map.png", p_map, width = 15, height = 12, dpi = 300)
