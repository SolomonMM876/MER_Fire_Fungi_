#############################################
# Script: Site-level Annual Precipitation
# Purpose: Calculate mean annual precipitation across plots per site
# Author: [Your Name]
# Date: [Date]
#############################################

library(tidyverse)
library(terra)

#--------------------------------------------
# 1. Import site data
#--------------------------------------------
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  dplyr::select(Site, Plot, lat = Lat, lon = Long)

#--------------------------------------------
# 2. Load WorldClim monthly precipitation rasters
#--------------------------------------------
precip_files <- list.files(
  "C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/precip",
  pattern = "\\.tif$",
  full.names = TRUE
)

# Ensure files are ordered Jan–Dec
precip_files <- precip_files[order(precip_files)]

# Load into SpatRaster stack
precip_stack <- rast(precip_files)
names(precip_stack) <- month.abb  # Jan, Feb, ..., Dec

#--------------------------------------------
# 3. Extract precipitation values for each plot
#--------------------------------------------
site_vect <- vect(All_Sites, geom = c("lon", "lat"), crs = "EPSG:4326")
precip_values <- extract(precip_stack, site_vect)

# Combine with site metadata
precip_result <- bind_cols(All_Sites, precip_values)

#--------------------------------------------
# 4. Calculate annual precipitation per plot
#--------------------------------------------
precip_plot <- precip_result %>%
  mutate(precip_annual = rowSums(across(all_of(month.abb)), na.rm = TRUE)) %>%
 dplyr::select(Site, Plot, precip_annual)

#--------------------------------------------
# 5. Average annual precipitation across plots per site
#--------------------------------------------
precip_site <- precip_plot %>%
  group_by(Site) %>%
  summarise(mean_annual_precip = mean(precip_annual, na.rm = TRUE), .groups = "drop") %>% 
  mutate(mean_annual_precip=round(mean_annual_precip,-1))

#--------------------------------------------
# 6. Save outputs
#--------------------------------------------
saveRDS(precip_site, "Processed_data/metadata/precip_MER.RDS")
write.csv(precip_site, "Processed_data/metadata/precip_MER.csv", row.names = FALSE)

