library(tidyverse)
library(sf)
library(lubridate)


#IMPORT and Clean########
#set path
#gdb_path<-'raw_data/fire_history_aus_2_7_25/Bushfire_Boundaries_Historical_2024_V3.gdb'
#read in all data
#fire_data <- st_read(gdb_path, layer = "Bushfire_Boundaries_Historical_V3")
#subset for data from 1970 onward and no burn date
# fire_data_1970on <- fire_data %>%
#   filter(!is.na(ignition_date)) %>%
#   filter(ignition_date >= ymd("1970-01-01"))
#st_write(fire_data_1970on, "raw_data/fire_history_aus_2_7_25/fire_data_1970on.gpkg")

#ID problematic polygons
# # Create a logical vector indicating validity
# validity_vec <- st_is_valid(fire_data_1970on)
# 
# # Add this as a new column to your data
# fire_data_1970on <- fire_data_1970on %>%
#   mutate(is_valid = validity_vec)
# 
# invalid_geoms <- fire_data_1970on %>%
#   filter(!is_valid)

# # Plot invalid geometries (red) and valid ones (grey)
# library(ggplot2)
# 
# ggplot() +
#   geom_sf(data = fire_data_1970on %>% filter(is_valid), fill = "grey80", color = NA) +
#   geom_sf(data = invalid_geoms, fill = "red", alpha = 0.6) +
#   labs(title = "Invalid Fire Polygons in Red") +
#   theme_minimal()

#Fix invalid geometries (if you want to clean them)
#fire_data_1970on <- st_make_valid(fire_data_1970on)

##############

#load cleaned data 
fire_data_1970on <- st_read("raw_data/fire_history_aus_2_7_25/fire_data_1970on.gpkg")

#read in all site locations
All_Sites_sf<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% 
  select(Site,Plot,Fire_Treatment,Lat,Long) %>% 
  st_as_sf(coords = c("Long", "Lat"), crs = 4283)  # GDA94 (EPSG:4283)


sf::sf_use_s2(FALSE)

# Spatial join: get all fire events for each plot
fires_per_plot <- st_join(All_Sites_sf, fire_data_1970on, join = st_intersects, left = FALSE)


previous_fires<-fires_per_plot %>% 
  mutate(
    ignition_year = year(ymd(ignition_date))
  ) %>%
  distinct(Site,ignition_year) %>% 
  as.data.frame() %>% 
  #filter(ignition_year < 2018) %>% 
  group_by(Site) %>%
  summarise(previous_fire_year = max(ignition_year, na.rm = TRUE),
            next_previous_fire_year = {
              sorted_years <- sort(ignition_year)
              n <- length(sorted_years)
              if (n %% 2 == 1) {
                # Odd number: take the middle
                sorted_years[(n + 1) / 2]
              } else {
                # Even number: return NA (no true middle value)
                NA_integer_
              }
            },
            oldest_fire_year = min(ignition_year, na.rm = TRUE),
            .groups = "drop")%>%
  mutate(
    next_previous_fire_year = if_else(next_previous_fire_year == previous_fire_year, NA_integer_, next_previous_fire_year),
    oldest_fire_year = if_else(oldest_fire_year == previous_fire_year, NA_integer_, oldest_fire_year)
  )

#save previous fire data
saveRDS(previous_fires, 'Processed_data/metadata/fire/previous_fires.Rdata')




#############extract fire history for Broome sites############


library(sf)
library(dplyr)

# # 1. Read in your site data (assuming it's in GDA94 / EPSG:4283)
# All_Sites_sf <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
#   select(Site, Plot, Fire_Treatment, Lat, Long) %>%
#   filter(str_detect(Site, 'WAMDAL')) %>%
#   st_as_sf(coords = c("Long", "Lat"), crs = 4283)  # GDA94
# 
# # 2. Read the NAFI fire frequency shapefile
# # Update the path to your shapefile directory
# nafi_fire <- st_read("raw_data/fire_history_aus_2_7_25/North_Aus_Fire_hist/FF_2000_24_gda94.shp")  # or similar file name
# 
# # 3. Ensure both layers use the same CRS (assume NAFI is in GDA94)
# nafi_fire <- st_transform(nafi_fire, 4283)  # EPSG:4283
# 
# # 4. Spatial join: extract the GRIDCODE (fire frequency) for each site
# #The data are vector polygons derived from an image with a resolution of 250m per pixel and each polygon is tagged with an attribute (gridcode) that corresponds to the number of years that pixel has been detected as being burnt in the 25 year period 2000 to 2024. The data is based on monthly fire scar mapping and extends across the WA rangelands, across the entire NT (down to 26 degrees S) and across northern Qld down to approx. 20 degrees S. The monthly fire scar mapping of landscapes north of 20 degrees S has been validated by aerial transects and checked using high resolution imagery across northern Australia 
# site_fire_history <- st_join(All_Sites_sf, nafi_fire, left = TRUE)
# 
# # 5. View or save the result
# print(site_fire_history)
# 
# # Optional: Export to CSV
# site_fire_history_df <- site_fire_history %>%
#   st_drop_geometry()  # Remove geometry to save as flat table
# 
# saveRDS(site_fire_history_df, 'Processed_data/metadata/fire/Broom_site_fire_history_df.Rdata')


