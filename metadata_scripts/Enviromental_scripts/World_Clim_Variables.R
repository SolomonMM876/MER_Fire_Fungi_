#load libraries and import df
library(raster)
library(dplyr)
library(tibble)
library(readxl)

#import data
#source('MER_All_Site_Data.R')
# Read in your site data
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  dplyr::select(Site, Plot, lat = Lat, lon = Long) 


sites <- All_Sites[,c('Site', 'lat', 'lon')] %>%
  distinct(Site, .keep_all = T)%>%
  column_to_rownames(var='Site')




#extract monthly precipitation averages for all locations
#not sure why I cant change file path here, but this is all that works
precip.files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/precip", ".tif", full.names=TRUE)
precip <- stack(precip.files)

month <- c("precip_Jan", "precip_Feb", "precip_Mar", "precip_Apr", "precip_May", "precip_Jun", "precip_Jul", "precip_Aug",
           "precip_Sep", "precip_Oct", "precip_Nov", "precip_Dec")
names(precip) <- month

precip.data<-raster::extract(precip,sites)

###################bioclimatic var###
bio.files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/bio", ".tif", full.names=TRUE)
bio.var <- stack(bio.files)
#Bio 1 and Bio12 are mean annual temperature and annual precipitation
Temp__Precip_Annual_ <- bio.var[[c(1,12)]]
names(Temp__Precip_Annual_) <- c("Annual_Temp","Annual_Prec")

Temp__Precip_Annual_data<-raster::extract(Temp__Precip_Annual_,sites)


########Tavg######
# Tavg.files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/tavg", ".tif", full.names=TRUE)
# Tavg.var <- stack(Tavg.files)
# 
# month <- c("Tavg_Jan", "Tavg_Feb", "Tavg_Mar", "Tavg_Apr", "Tavg_May", "Tavg_Jun", "Tavg_Jul", "Tavg_Aug",
#            "Tavg_Sep", "Tavg_Oct", "Tavg_Nov", "Tavg_Dec")
# names(Tavg.var) <- month
# 
# Tavg_data<-raster::extract(Tavg.var,sites)

###elevation####
elev.files <- list.files("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/world clim/elevation", ".tif", full.names=TRUE)
elev.var <- stack(elev.files)

elev_data<-raster::extract(elev.var,sites)

###Cbind all extracted data in one df


Site_Env_data<-cbind(sites,precip.data,Temp__Precip_Annual_data,elev_data)
Site_Env_data<-rownames_to_column(Site_Env_data, var='site')

#remove df's used to make
rm(sites,precip.data,Temp__Precip_Annual_data,elev_data,bio.var,elev.var,Temp__Precip_Annual_)

