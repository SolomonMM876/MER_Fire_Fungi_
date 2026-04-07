# Load necessary libraries
library(terra)
library(dplyr)

# Note: In the Global-Aridity dataset, which uses this formulation, Aridity Index values increase
# for more humid conditions, and decrease with more arid conditions.


# Read in your site data
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, lat = Lat, lon = Long) 

# Load the Global Aridity Index (AI) monthly raster stack
#this is from a folder outside of my wd
ai_rasters <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/aridity/Global-AI_v3_monthly", 
                         pattern = "\\.tif$", 
                         full.names = TRUE)

# Sort files to make sure months are in order from 01 to 12
ai_rasters <- ai_rasters[order(ai_rasters)]

# Stack the rasters (each layer represents a month)
ai_stack <- rast(ai_rasters)

# Give the layers clear names: Jan to Dec
month_names <- month.abb  # "Jan", "Feb", ..., "Dec"
names(ai_stack) <- month_names

# Create spatial points from your site data
sites_vect <- vect(All_Sites, geom = c("lon", "lat"), crs = "EPSG:4326")

# Extract AI values for each site
aridity_values <- terra::extract(ai_stack, sites_vect)

# Combine with site identifiers
result <- bind_cols(All_Sites, aridity_values)

##########################################################################################
#some plots have missing values
missing_val<-result %>% 
  filter(Jan==0) %>% select(Site:lon)

# Snap each lat/lon to the nearest raster cell center
snap_to_grid <- function(coord, res = 1/120) {
  round(coord / res) * res  # 1/120 = 0.0083333° (30 arc-sec)
}

missing_val_snapped <- missing_val %>%
  mutate(
    lon_snap = snap_to_grid(lon),
    lat_snap = snap_to_grid(lat)
  )

# Create spatial points
snapped_vect <- vect(missing_val_snapped, geom = c("lon_snap", "lat_snap"), crs = "EPSG:4326")
# Extract using snapped coordinates
snapped_extract <- terra::extract(ai_stack, snapped_vect)
# Join results
missing_val_snapped <- bind_cols(missing_val_snapped, snapped_extract) %>% mutate(snap='snapped')
missing_val_snapped
# Step 1: Identify rows in `result` where all months are zero
result_fil <- result %>%
  filter(!if_all(Jan:Dec, ~ .x == 0)) 

# Step 2: Join updated values to the original result
result_updated <- result_fil %>%
  bind_rows(missing_val_snapped)


###############find values for empty plot data############

# Step 1: Extract the latitude of VCRF15 that has 0 values
vcrf15_lat <- result_updated %>%
  filter(Plot == "VCRF15") %>%
  pull(lat)

# Step 2: Find the plot with the closest latitude (excluding VCRF15)
closest_by_lat <- result_updated %>%
  filter(Plot != "VCRF15") %>%
  mutate(lat_diff = abs(lat - vcrf15_lat)) %>%
  arrange(lat_diff) %>%
  slice(1)

# Step 3: Replace VCRF15's aridity values with values from closest_by_lat
month_cols <- month.abb

result_updated <- result_updated %>%
  mutate(across(all_of(month_cols),
                ~ if_else(Plot == "VCRF15",
                          closest_by_lat[[cur_column()]],
                          .x)))


##################finally################



# View result
head(result_updated)
#calc annnual averages
result_updated$AI_annual <- rowMeans(result_updated[, month_names], na.rm = TRUE)

####Final Df
# Aridity Index values increase for more humid conditions, and decrease with more arid conditions.
MER_aridity<-result_updated %>% dplyr::select(Site,Plot,AI_annual) %>% 
  group_by(Site) %>% 
  summarise(aridity=mean(AI_annual)*0.0001)

saveRDS(MER_aridity,'Processed_data/metadata/Aridity_MER.Rdata')
