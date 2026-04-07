#libraries
library(readxl)
library(tidyverse)
library(patchwork)
library(oce) #change format of GPS cord
library(lubridate)


#import and clean data so it is proper format
#NSMSEQCasWoo01
NSMSEQCasWoo01<- read_excel("raw_data/MER_Site_Data/Site_Labels_BaileyKing.xlsx", sheet = 2)
NSMSEQCasWoo01<-NSMSEQCasWoo01[-1,]
NSMSEQCasWoo01$Site<-'NSMSEQCasWoo01'
colnames(NSMSEQCasWoo01)[1]<-'Plot'
colnames(NSMSEQCasWoo01)[2]<-'Plot_num'
colnames(NSMSEQCasWoo01)[5]<-'Lat'
colnames(NSMSEQCasWoo01)[6]<-'Long'
colnames(NSMSEQCasWoo01)[7]<-'Deployment_Date'
colnames(NSMSEQCasWoo01)[9]<-'Notes'
NSMSEQCasWoo01$Lat <- signif(as.numeric(NSMSEQCasWoo01$Lat), 10)
NSMSEQCasWoo01$Long <- signif(as.numeric(NSMSEQCasWoo01$Long), 10)
NSMSEQCasWoo01$Deployment_Date<-as.Date(as.numeric(NSMSEQCasWoo01$Deployment_Date),origin = "1899-12-30")

#WAMCOOEucWoo01
WAMCOOEucWoo01<- read_excel("raw_data/MER_Site_Data/Great_Western_Woodlands_Completed.xlsx", sheet = 9)
WAMCOOEucWoo01<-WAMCOOEucWoo01[-1,-1]
WAMCOOEucWoo01$Site<-'WAMCOOEucWoo01'
colnames(WAMCOOEucWoo01)[1]<-'Plot'
colnames(WAMCOOEucWoo01)[2]<-'Plot_num'
colnames(WAMCOOEucWoo01)[5]<-'Lat'
colnames(WAMCOOEucWoo01)[6]<-'Long'
colnames(WAMCOOEucWoo01)[7]<-'Deployment_Date'
WAMCOOEucWoo01$Lat <- signif(as.numeric(WAMCOOEucWoo01$Lat), 10)
WAMCOOEucWoo01$Long <- signif(as.numeric(WAMCOOEucWoo01$Long), 10)
WAMCOOEucWoo01$Deployment_Date<-as.Date(as.numeric(WAMCOOEucWoo01$Deployment_Date),origin = "1899-12-30")

#WAMCOOEucWoo02
WAMCOOEucWoo02<- read_excel("raw_data/MER_Site_Data/Great_Western_Woodlands_Completed.xlsx", sheet = 7)
WAMCOOEucWoo02<-WAMCOOEucWoo02[-1,-1]
WAMCOOEucWoo02$Site<-'WAMCOOEucWoo02'
colnames(WAMCOOEucWoo02)[1]<-'Plot'
colnames(WAMCOOEucWoo02)[2]<-'Plot_num'
colnames(WAMCOOEucWoo02)[5]<-'Lat'
colnames(WAMCOOEucWoo02)[6]<-'Long'
colnames(WAMCOOEucWoo02)[7]<-'Deployment_Date'
WAMCOOEucWoo02$Lat <- signif(as.numeric(WAMCOOEucWoo02$Lat), 10)
WAMCOOEucWoo02$Long <- signif(as.numeric(WAMCOOEucWoo02$Long), 10)
WAMCOOEucWoo02$Deployment_Date<-as.Date(as.numeric(WAMCOOEucWoo02$Deployment_Date),origin = "1899-12-30")

#WAMCOOEucWoo03
WAMCOOEucWoo03<- read_excel("raw_data/MER_Site_Data/Great_Western_Woodlands_Completed.xlsx", sheet = 4)
WAMCOOEucWoo03<-WAMCOOEucWoo03[-1,-1]
WAMCOOEucWoo03$Site<-'WAMCOOEucWoo03'
colnames(WAMCOOEucWoo03)[1]<-'Plot'
colnames(WAMCOOEucWoo03)[2]<-'Plot_num'
colnames(WAMCOOEucWoo03)[5]<-'Lat'
colnames(WAMCOOEucWoo03)[6]<-'Long'
colnames(WAMCOOEucWoo03)[7]<-'Deployment_Date'
WAMCOOEucWoo03$Lat <- signif(as.numeric(WAMCOOEucWoo03$Lat), 10)
WAMCOOEucWoo03$Long <- signif(as.numeric(WAMCOOEucWoo03$Long), 10)
WAMCOOEucWoo03$Deployment_Date<-as.Date(as.numeric(WAMCOOEucWoo03$Deployment_Date),origin = "1899-12-30")

#WAMCOOShrubl01 
WAMCOOShrubl01 <- read_excel("raw_data/MER_Site_Data/Great_Western_Woodlands_Completed.xlsx", sheet = 2)
WAMCOOShrubl01 <-WAMCOOShrubl01 [-1,-1]
WAMCOOShrubl01 $Site<-'WAMCOOShrubl01'
colnames(WAMCOOShrubl01 )[1]<-'Plot'
colnames(WAMCOOShrubl01)[2]<-'Plot_num'
colnames(WAMCOOShrubl01 )[5]<-'Lat'
colnames(WAMCOOShrubl01 )[6]<-'Long'
colnames(WAMCOOShrubl01 )[7]<-'Deployment_Date'
WAMCOOShrubl01 $Lat <- signif(as.numeric(WAMCOOShrubl01$Lat), 10)
WAMCOOShrubl01 $Long <- signif(as.numeric(WAMCOOShrubl01$Long), 10)
#WAMCOOShrubl01 $Deployment_Date<-as.Date(as.numeric(WAMCOOShrubl01 $Deployment_Date),origin = "1899-12-30") #not needed here

#VCMSECEucFor01
VCMSECEucFor01 <- read_excel("raw_data/MER_Site_Data/Site_Labels_East Gipplands_Robyn.xlsx", sheet = 2)
VCMSECEucFor01 <-VCMSECEucFor01 [-1,-1]
VCMSECEucFor01 $Site<-'VCMSECEucFor01'
colnames(VCMSECEucFor01 )[1]<-'Plot'
colnames(VCMSECEucFor01)[2]<-'Plot_num'
colnames(VCMSECEucFor01 )[7]<-'Deployment_Date'
VCMSECEucFor01$Deployment_Date<-as.Date(as.numeric(VCMSECEucFor01 $Deployment_Date),origin = "1899-12-30") 
#load GPS points in dif format
VCMSECEucFor01_GPS <- read_excel("raw_data/MER_Site_Data/Site_Labels_East Gipplands_Robyn.xlsx", sheet = 4)
result <- utm2lonlat(VCMSECEucFor01_GPS$Easting, VCMSECEucFor01_GPS$Northing,55,hemisphere = "S")
VCMSECEucFor01<-cbind(VCMSECEucFor01,result)
VCMSECEucFor01<-VCMSECEucFor01[,-c(5,6,10)]
colnames(VCMSECEucFor01)[9]<-'Long'
colnames(VCMSECEucFor01 )[10]<-'Lat'
VCMSECEucFor01$Lat <- gsub("-", "", VCMSECEucFor01$Lat)
rm(result,VCMSECEucFor01_GPS)


#WAMESPShrubl01
WAMESPShrubl01 <- read_excel("raw_data/MER_Site_Data/Site_Labels_Waychinicup.xlsx", sheet = 2)
WAMESPShrubl01 <-WAMESPShrubl01 [-1,]
WAMESPShrubl01 $Site<-'WAMESPShrubl01'
colnames(WAMESPShrubl01 )[1]<-'Plot'
colnames(WAMESPShrubl01)[2]<-'Plot_num'
colnames(WAMESPShrubl01 )[5]<-'Lat'
colnames(WAMESPShrubl01 )[6]<-'Long'
colnames(WAMESPShrubl01 )[7]<-'Deployment_Date'
WAMESPShrubl01 $Lat <- as.numeric(WAMESPShrubl01$Lat)
WAMESPShrubl01 $Long <- as.numeric(WAMESPShrubl01$Long)
WAMESPShrubl01$Deployment_Date<-as.Date(as.numeric(WAMESPShrubl01 $Deployment_Date),origin = "1899-12-30")
#change format GPS points
result <- utm2lonlat(WAMESPShrubl01$Lat, WAMESPShrubl01$Long,50,hemisphere = "S") #timezone 50
WAMESPShrubl01<-cbind(WAMESPShrubl01,result)
WAMESPShrubl01<-WAMESPShrubl01[,-c(5,6)]
colnames(WAMESPShrubl01)[9]<-'Long'
colnames(WAMESPShrubl01 )[10]<-'Lat'
WAMESPShrubl01$Lat <- gsub("-", "", WAMESPShrubl01$Lat)
WAMESPShrubl01 $Lat <- signif(as.numeric(WAMESPShrubl01$Lat), 10)
WAMESPShrubl01 $Long <- signif(as.numeric(WAMESPShrubl01$Long), 10)
rm(result)


#SAMMDDMallee01
SAMMDDMallee01 <- read_excel("raw_data/MER_Site_Data/Site Labels_Murraylands.xlsx", sheet = 1)
SAMMDDMallee01 <-SAMMDDMallee01 [-1,-1]
SAMMDDMallee01 $Site<-'SAMMDDMallee01'
colnames(SAMMDDMallee01 )[1]<-'Plot'
colnames(SAMMDDMallee01)[2]<-'Plot_num'
colnames(SAMMDDMallee01)[5] <- 'Easting'   # Was 'Lat'
colnames(SAMMDDMallee01)[6] <- 'Northing'  # Was 'Long'
colnames(SAMMDDMallee01 )[7]<-'Deployment_Date'
SAMMDDMallee01$Easting <- as.numeric(SAMMDDMallee01$Easting)
SAMMDDMallee01$Northing <- as.numeric(SAMMDDMallee01$Northing)
SAMMDDMallee01$Deployment_Date<-as.Date(as.numeric(SAMMDDMallee01 $Deployment_Date),origin = "1899-12-30")
#change format GPS points
result <- utm2lonlat(easting = SAMMDDMallee01$Easting, 
                     northing = SAMMDDMallee01$Northing, 
                     zone = 54,        # Changed from 50 to 54
                     hemisphere = "S")

SAMMDDMallee01<-cbind(SAMMDDMallee01,result)
SAMMDDMallee01<-SAMMDDMallee01[,-c(5,6)]
colnames(SAMMDDMallee01)[9]<-'Long'
colnames(SAMMDDMallee01 )[10]<-'Lat'
SAMMDDMallee01$Lat <- gsub("-", "", SAMMDDMallee01$Lat)
SAMMDDMallee01 $Lat <- signif(as.numeric(SAMMDDMallee01$Lat), 10)
SAMMDDMallee01 $Long <- signif(as.numeric(SAMMDDMallee01$Long), 10)
rm(result)



#SAMMDDMallee02
SAMMDDMallee02 <- read_excel("raw_data/MER_Site_Data/Site Labels_Murraylands.xlsx", sheet = 3)
SAMMDDMallee02 <-SAMMDDMallee02 [-1,-1]
SAMMDDMallee02 $Site<-'SAMMDDMallee02'
colnames(SAMMDDMallee02 )[1]<-'Plot'
colnames(SAMMDDMallee02)[2]<-'Plot_num'
colnames(SAMMDDMallee02 )[5]<-'Easting'   # Was 'Lat'
colnames(SAMMDDMallee02 )[6]<-'Northing'  # Was 'Long'
colnames(SAMMDDMallee02 )[7]<-'Deployment_Date'
SAMMDDMallee02$Easting <- as.numeric(SAMMDDMallee02$Easting)
SAMMDDMallee02$Northing <- as.numeric(SAMMDDMallee02$Northing)
SAMMDDMallee02$Deployment_Date<-as.Date(as.numeric(SAMMDDMallee02 $Deployment_Date),origin = "1899-12-30")
#change format GPS points
result <- utm2lonlat(easting = SAMMDDMallee02$Easting, 
                     northing = SAMMDDMallee02$Northing, 
                     zone = 54,        # Changed from 50 to 54
                     hemisphere = "S")
SAMMDDMallee02<-cbind(SAMMDDMallee02,result)
SAMMDDMallee02<-SAMMDDMallee02[,-c(5,6)]
colnames(SAMMDDMallee02)[9]<-'Long'
colnames(SAMMDDMallee02 )[10]<-'Lat'
SAMMDDMallee02$Lat <- gsub("-", "", SAMMDDMallee02$Lat)
SAMMDDMallee02 $Lat <- signif(as.numeric(SAMMDDMallee02$Lat), 10)
SAMMDDMallee02 $Long <- signif(as.numeric(SAMMDDMallee02$Long), 10)
rm(result)













#VCMSECRainfo01
VCMSECRainfo01<- read_excel("raw_data/MER_Site_Data/Site_Labels_East Gipplands_SP.xlsx", sheet = 2)
VCMSECRainfo01<-VCMSECRainfo01[-1,]
VCMSECRainfo01$Site<-'VCMSECRainfo01'
colnames(VCMSECRainfo01)[1]<-'Plot'
colnames(VCMSECRainfo01)[2]<-'Plot_num'
colnames(VCMSECRainfo01)[5]<-'Lat'
colnames(VCMSECRainfo01)[6]<-'Long'
colnames(VCMSECRainfo01)[7]<-'Deployment_Date'
VCMSECRainfo01$Lat <- signif(as.numeric(VCMSECRainfo01$Lat), 10)
VCMSECRainfo01$Long <- signif(as.numeric(VCMSECRainfo01$Long), 10)
VCMSECRainfo01$Deployment_Date<-as.Date(as.numeric(VCMSECRainfo01$Deployment_Date),origin = "1899-12-30")

#WAMDALRainfo01
WAMDALRainfo01<- read_excel("raw_data/MER_Site_Data/Site_Labels_Broome.xlsx", sheet = 'Bag Deployment WAMDALRainfo01')
WAMDALRainfo01<-WAMDALRainfo01[-c(1),-1]
WAMDALRainfo01$Site<-'WAMDALRainfo01'
colnames(WAMDALRainfo01)[1]<-'Plot'
colnames(WAMDALRainfo01)[2]<-'Plot_num'
WAMDALRainfo01[3] <- substr(WAMDALRainfo01$'Fire Treatment', 1,1)
colnames(WAMDALRainfo01)[5]<-'Lat'
WAMDALRainfo01$Lat <- abs(as.numeric(WAMDALRainfo01$Lat))
colnames(WAMDALRainfo01)[6]<-'Long'
colnames(WAMDALRainfo01)[7]<-'Deployment_Date'
WAMDALRainfo01$Lat <- signif(as.numeric(WAMDALRainfo01$Lat), 10)
WAMDALRainfo01$Long <- signif(as.numeric(WAMDALRainfo01$Long), 10)
WAMDALRainfo01$Deployment_Date<-as.Date(as.numeric(WAMDALRainfo01$Deployment_Date),origin = "1899-12-30")


#WAMDALEucW0001
WAMDALEucW0001<- read_excel("raw_data/MER_Site_Data/Site_Labels_Broome.xlsx", sheet = 2)
WAMDALEucW0001<-WAMDALEucW0001[-c(1),-1]
WAMDALEucW0001$Site<-'WAMDALEucWoo01'
colnames(WAMDALEucW0001)[1]<-'Plot'
colnames(WAMDALEucW0001)[2]<-'Plot_num'
WAMDALEucW0001[3] <- substr(WAMDALEucW0001$'Fire Treatment', 1,1)
colnames(WAMDALEucW0001)[5]<-'Lat'
WAMDALEucW0001$Lat <- abs(as.numeric(WAMDALEucW0001$Lat))
colnames(WAMDALEucW0001)[6]<-'Long'
colnames(WAMDALEucW0001)[7]<-'Deployment_Date'
WAMDALEucW0001$Lat <- signif(as.numeric(WAMDALEucW0001$Lat), 10)
WAMDALEucW0001$Long <- signif(as.numeric(WAMDALEucW0001$Long), 10)
WAMDALEucW0001$Deployment_Date<-as.Date(as.numeric(WAMDALEucW0001$Deployment_Date),origin = "1899-12-30")



#WAMDALEucW0002
WAMDALEucW0002<- read_excel("raw_data/MER_Site_Data/Site_Labels_Broome.xlsx", sheet = 4)
WAMDALEucW0002<-WAMDALEucW0002[-c(1),-1]
WAMDALEucW0002$Site<-'WAMDALEucW0002'
colnames(WAMDALEucW0002)[1]<-'Plot'
colnames(WAMDALEucW0002)[2]<-'Plot_num'
WAMDALEucW0002[3] <- substr(WAMDALEucW0002$'Fire Treatment', 1,1)
colnames(WAMDALEucW0002)[5]<-'Lat'
WAMDALEucW0002$Lat <- abs(as.numeric(WAMDALEucW0002$Lat))
colnames(WAMDALEucW0002)[6]<-'Long'
colnames(WAMDALEucW0002)[7]<-'Deployment_Date'
WAMDALEucW0002$Lat <- signif(as.numeric(WAMDALEucW0002$Lat), 10)
WAMDALEucW0002$Long <- signif(as.numeric(WAMDALEucW0002$Long), 10)
WAMDALEucW0002$Deployment_Date<-as.Date(as.numeric(WAMDALEucW0002$Deployment_Date),origin = "1899-12-30")

#QDMSEQRainfo01
QDMSEQRainfo01<- read_excel("raw_data/MER_Site_Data/Site_Labels_Lamington.xlsx", sheet = 2)
QDMSEQRainfo01<-QDMSEQRainfo01[-c(1),]
QDMSEQRainfo01$Site<-'QDMSEQRainfo01'
colnames(QDMSEQRainfo01)[1]<-'Plot'
colnames(QDMSEQRainfo01)[2]<-'Plot_num'
QDMSEQRainfo01[3] <- substr(QDMSEQRainfo01$'Fire Treatment', 1,1)
colnames(QDMSEQRainfo01)[5]<-'Lat'
QDMSEQRainfo01$Lat <- abs(as.numeric(QDMSEQRainfo01$Lat))
colnames(QDMSEQRainfo01)[6]<-'Long'
colnames(QDMSEQRainfo01)[7]<-'Deployment_Date'
QDMSEQRainfo01$Lat <- signif(as.numeric(QDMSEQRainfo01$Lat), 10)
QDMSEQRainfo01$Long <- signif(as.numeric(QDMSEQRainfo01$Long), 10)
QDMSEQRainfo01$Deployment_Date<-as.Date(as.numeric(QDMSEQRainfo01$Deployment_Date),origin = "1899-12-30")





#combine df's to master sheet
All_Sites<-rbind(NSMSEQCasWoo01,VCMSECEucFor01,VCMSECRainfo01, WAMCOOEucWoo01,WAMCOOEucWoo02,
                       WAMCOOEucWoo03,WAMCOOShrubl01,WAMESPShrubl01,WAMDALRainfo01,WAMDALEucW0001,
                       WAMDALEucW0002,SAMMDDMallee01,SAMMDDMallee02,QDMSEQRainfo01) %>%
  mutate(Lat=as.numeric(paste0('-',Lat)))#add '-' to all latitudes 

colnames(All_Sites)[3]<-'Fire_Treatment'
colnames(All_Sites)[4]<-'Operator'

#remove all other df's
rm(NSMSEQCasWoo01,VCMSECEucFor01,VCMSECRainfo01, WAMCOOEucWoo01,WAMCOOEucWoo02,
WAMCOOEucWoo03,WAMCOOShrubl01,WAMESPShrubl01,WAMDALRainfo01,WAMDALEucW0001,WAMDALEucW0002,
SAMMDDMallee01,SAMMDDMallee02,QDMSEQRainfo01)


All_Sites<-All_Sites %>% 
  mutate(  Plot = str_remove(Plot, "[BU]$"),            # Remove burn assignment
           Plot = str_replace(Plot, "^NSS00", "NSS0"),
           Plot = str_replace(Plot, "^NSSO", "NSS"),
           Site = str_replace(Site, "WAMDALEucW00","WAMDALEucWoo"),
             Deployment_Date = as.Date(Deployment_Date),
           Site=as.factor(Site),
           Plot=as.factor(Plot),
           Fire_Treatment=as.factor(Fire_Treatment)
  ) %>% 
  rename(Notes_install=Notes,
         Fire_Treatment=Fire_Treatment)

saveRDS(All_Sites, file = "raw_data/MER_Site_Data/All_Sites.RDS") 

temp<-All_Sites %>% 
  filter(Site %in% c('SAMMDDMallee02','SAMMDDMallee01' )) %>% 
  select(Site,Plot,Lat,Long,Fire_Treatment)

#write.csv(temp, 'SA_sites.csv', row.names = FALSE)
