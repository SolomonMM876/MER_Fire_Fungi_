library(tidyverse)


# Read in your site data
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, lat = Lat, lon = Long) 

# Load the Global Aridity Index (AI) monthly raster stack
#this is from a folder outside of my wd
tmax_files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/tmax", 
                         pattern = "\\.tif$", 
                         full.names = TRUE)
# Step 3: Sort files by month — assumes file names are in order like tmax_01.tif ... tmax_12.tif
tmax_files <- tmax_files[order(tmax_files)]

# Step 4: Load into SpatRaster stack
tmax_stack <- rast(tmax_files)

# Step 5: Name layers by month
names(tmax_stack) <- month.abb  # Jan, Feb, ..., Dec

# Step 6: Create spatial object from sites
site_vect <- vect(All_Sites, geom = c("lon", "lat"), crs = "EPSG:4326")

# Step 7: Extract Tmax values
tmax_values <- extract(tmax_stack, site_vect)

# Step 8: Combine with site info
tmax_result <- bind_cols(All_Sites, tmax_values)

# View output
head(tmax_result)

tmax_result <- tmax_result %>%
  mutate(tmax_annual_max = pmax(!!!syms(month.abb), na.rm = TRUE)) %>%
  select(Site, Plot, tmax_annual_max) %>% 
  group_by(Site) %>% 
  summarise(tmax=max(tmax_annual_max, na.rm = TRUE))

saveRDS(tmax_result,'Processed_data/metadata/Tmax_MER.Rdata')
