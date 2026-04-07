# Install and load necessary packages
library(raster)
library(sf)
library("chirps")
library(tidyverse)

#import data
All_Sites<-read.csv("MER_Site_Data/All_Sites.csv")

MER_unique_sites <- All_Sites[,c('Site', 'Long', 'Lat')] %>%
  distinct(Site, .keep_all = T)%>%
  group_by(Site)%>%
  mutate(ID = cur_group_id())


#extract precipitation data

#extract lon and lat
location <-as.data.frame(MER_unique_sites[,c('Long','Lat')])

###this step takes a while
Precip <- get_chirps(location, dates = c("2021-01-01", "2023-11-24"), server = "ClimateSERV")

#merge df's based on ID
colnames(Precip)[1]<-'ID'
Precip$ID<-as.character(Precip$ID)
MER_unique_sites$ID<-as.character(MER_unique_sites$ID)

Sites_Precip<-full_join(MER_unique_sites,Precip, by=c('ID'))

