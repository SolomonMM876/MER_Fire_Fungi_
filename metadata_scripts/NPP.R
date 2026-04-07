
library(MODISTools)
library(dplyr)
library(purrr)

#site coordinates and sampling dates
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, lat=Lat,lon=Long)

All_biomass<-readRDS('processed_data/Bag_data/MER_biomass.RDS') %>% 
  select(Site,end=Harvest_date,start=Deployment_Date) %>% distinct() %>% 
  #select distinct time for sampling each plot
  group_by(Site) %>%
  summarise(
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    .groups = "drop"
  )

#testing to see if NPP from full period since burn explained more variation
# #import df from Fire_regime_MER
# date_recent_fire<-MER_fire_regime %>% distinct(Site, year_of_most_recent_fire, month_of_most_recent_fire) %>% 
#   mutate(
#     recent_fire = ymd(paste(year_of_most_recent_fire, month_of_most_recent_fire, "01"))
#   ) %>% select(Site,recent_fire)
# 


sites<-All_Sites %>% left_join(All_biomass) %>%
  #left_join(date_recent_fire) %>% 
  rename(site_name=Plot) %>% select(-Site) 
# sites <- sites %>% mutate(site_id = row_number())



#select product
#mt_products() %>% filter(grepl("MOD13", product))
#mt_bands("MOD13Q1")


######################
#this download for sampling period
# Efficient, vectorized MODIS NDVI download with pmap
ndvi_data <- pmap_dfr(
  sites,
  function(site_name, lat, lon, start, end) {
    mt_subset(
      product = "MOD13Q1",
      band = "250m_16_days_NDVI",
      lat = lat,
      lon = lon,
      start = start,
      end = end,
      site_name = site_name,
      internal = TRUE,
      progress = FALSE
    )
  }
)
###################################
# #this downloads since last burn
# # Efficient, vectorized MODIS NDVI download with pmap
# ndvi_data <- pmap_dfr(
#   sites,
#   function(site_name, lat, lon, recent_fire, end) {
#     mt_subset(
#       product = "MOD13Q1",
#       band = "250m_16_days_NDVI",
#       lat = lat,
#       lon = lon,
#       start = recent_fire,
#       end = end,
#       site_name = site_name,
#       internal = TRUE,
#       progress = FALSE
#     )
#   }
# )

# Rescale and format
ndvi_MER <- ndvi_data %>%
  mutate(
    NDVI = value * 0.0001,
    date = as.Date(calendar_date)
  ) %>% 
  filter(NDVI>0) %>% 
  group_by(site,latitude,longitude) %>%
  summarise(
    mean_NDVI = mean(NDVI, na.rm = TRUE),
    sd_NDVI = sd(NDVI, na.rm = TRUE),
    n_observations = n(),
    se_NDVI=sd_NDVI/sqrt(n_observations)
  )

#######################
#Adjust for sampling period
#I have no idea why I had to do this stupid rounding step to get the joins to work, but it did now
ndvi_plot_post_burn<-left_join(ndvi_MER %>%   mutate(lat_round = round(latitude, 6),
                                          lon_round = round(longitude, 6)),
                    sites %>% 
                      mutate(lat_round = round(lat, 6),
                             lon_round = round(lon, 6)), by = c("lat_round", "lon_round")) %>% 
  ungroup() %>% 
  select(Plot=site_name,mean_NDVI,se_NDVI,n_observations)

ndvi_plot_post_burn<-ndvi_plot_post_burn %>% 
  left_join(All_Sites %>% select(Site,Plot))

ndvi_site_post_burn<-ndvi_plot_post_burn %>% 
  left_join(All_Sites %>% select(Site,Plot))%>% 
  group_by(Site) %>% 
  summarise(NPP_since_fire=mean(mean_NDVI))

head(ndvi_plot_post_burn)

saveRDS(ndvi_site_post_burn,'Processed_data/metadata/NPP_MER_by_site.Rdata')
saveRDS(ndvi_plot_post_burn,'Processed_data/metadata/NPP_MER_by_plot.Rdata')
