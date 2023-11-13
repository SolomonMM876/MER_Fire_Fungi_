#libraries
library(readxl)
library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
library(splitstackshape)
library(tidyr)
library(car)

#import and clean data so it is proper format
#NSMSEQCasWoo01
NSMSEQCasWoo01<- read_excel("Site_Labels_BaileyKing.xlsx", sheet = 2)
NSMSEQCasWoo01<-NSMSEQCasWoo01[-1,]
NSMSEQCasWoo01$Site<-'NSMSEQCasWoo01'
colnames(NSMSEQCasWoo01)[1]<-'Plot'
colnames(NSMSEQCasWoo01)[2]<-'Plot_num'
colnames(NSMSEQCasWoo01)[5]<-'Lat'
colnames(NSMSEQCasWoo01)[6]<-'Long'
colnames(NSMSEQCasWoo01)[7]<-'Deployment_Date'
colnames(NSMSEQCasWoo01)[9]<-'Notes'
NSMSEQCasWoo01$Lat <- signif(as.numeric(NSMSEQCasWoo01$Lat), 6)
NSMSEQCasWoo01$Long <- signif(as.numeric(NSMSEQCasWoo01$Long), 7)
NSMSEQCasWoo01$Deployment_Date<-as.Date(as.numeric(NSMSEQCasWoo01$Deployment_Date),origin = "1899-12-30")

#WAMCOOEucWoo03
WAMCOOEucWoo03<- read_excel("Great_Western_Woodlands_Completed.xlsx", sheet = 4)
WAMCOOEucWoo03<-WAMCOOEucWoo03[-1,-1]
WAMCOOEucWoo03$Site<-'WAMCOOEucWoo03'
colnames(WAMCOOEucWoo03)[1]<-'Plot'
colnames(WAMCOOEucWoo03)[2]<-'Plot_num'
colnames(WAMCOOEucWoo03)[5]<-'Lat'
colnames(WAMCOOEucWoo03)[6]<-'Long'
colnames(WAMCOOEucWoo03)[7]<-'Deployment_Date'
WAMCOOEucWoo03$Lat <- signif(as.numeric(WAMCOOEucWoo03$Lat), 6)
WAMCOOEucWoo03$Long <- signif(as.numeric(WAMCOOEucWoo03$Long), 7)
WAMCOOEucWoo03$Deployment_Date<-as.Date(as.numeric(WAMCOOEucWoo03$Deployment_Date),origin = "1899-12-30")

#WAMCOOShrubl01 
WAMCOOShrubl01 <- read_excel("Great_Western_Woodlands_Completed.xlsx", sheet = 2)
WAMCOOShrubl01 <-WAMCOOShrubl01 [-1,-1]
WAMCOOShrubl01 $Site<-'WAMCOOShrubl01 '
colnames(WAMCOOShrubl01 )[1]<-'Plot'
colnames(WAMCOOShrubl01)[2]<-'Plot_num'
colnames(WAMCOOShrubl01 )[5]<-'Lat'
colnames(WAMCOOShrubl01 )[6]<-'Long'
colnames(WAMCOOShrubl01 )[7]<-'Deployment_Date'
WAMCOOShrubl01 $Lat <- signif(as.numeric(WAMCOOShrubl01$Lat), 6)
WAMCOOShrubl01 $Long <- signif(as.numeric(WAMCOOShrubl01$Long), 7)
#WAMCOOShrubl01 $Deployment_Date<-as.Date(as.numeric(WAMCOOShrubl01 $Deployment_Date),origin = "1899-12-30") #not needed here

#WAMESPShrubl01
WAMESPShrubl01 <- read_excel("Site_Labels_Waychinicup.xlsx", sheet = 2)

#conver utm to lon lat
utm2lonlat(easting, northing, zone = 1, hemisphere = "N", km = FALSE)

#combine df's to master sheet
combined_df_MER<-rbind(NSMSEQCasWoo01,WAMCOOShrubl01,WAMCOOShrubl01)

