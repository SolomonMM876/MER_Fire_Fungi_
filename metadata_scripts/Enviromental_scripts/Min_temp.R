library(tidyverse)


# Read in your site data
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, lat = Lat, lon = Long) 

# Load the Global Aridity Index (AI) monthly raster stack
#this is from a folder outside of my wd
tmin_files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/tmin", 
                         pattern = "\\.tif$", 
                         full.names = TRUE)
# Step 3: Sort files by month — assumes file names are in order like tmin_01.tif ... tmin_12.tif
tmin_files <- tmin_files[order(tmin_files)]

# Step 4: Load into SpatRaster stack
tmin_stack <- rast(tmin_files)

# Step 5: Name layers by month
names(tmin_stack) <- month.abb  # Jan, Feb, ..., Dec

# Step 6: Create spatial object from sites
site_vect <- vect(All_Sites, geom = c("lon", "lat"), crs = "EPSG:4326")

# Step 7: Extract tmin values
tmin_values <- extract(tmin_stack, site_vect)

# Step 8: Combine with site info
tmin_result <- bind_cols(All_Sites, tmin_values)

# View output
head(tmin_result)

tmin_result <- tmin_result %>%
  mutate(tmin_annual_min = pmin(!!!syms(month.abb), na.rm = TRUE)) %>%
  select(Site, Plot, tmin_annual_min) %>% 
  group_by(Site) %>% 
  summarise(tmin=min(tmin_annual_min, na.rm = TRUE))

saveRDS(tmin_result,'Processed_data/metadata/tmin_MER.Rdata')
