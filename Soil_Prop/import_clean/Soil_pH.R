library(readxl)
library(tidyverse)


MER_Soil_ph<- read_excel("raw_data/Soil_data/Soil_Ph_Kumari_24.xlsx")

#read in nute sample info
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS")





# Modify the Plot column in both data frames, removing the last character if it's U or B
MER_Soil_ph <- MER_Soil_ph %>%
  mutate(Plot = if_else(str_ends(Plot, "[UB]"), str_sub(Plot, 1, nchar(Plot) - 1), Plot)) %>% 
  select(-Site) %>% 
  mutate(Plot = str_replace(Plot, "NSSO", "NSS")) #my dyslexia left me to right VCFR instead of the correct VCRF
    


# Perform the anti join

anti_join(MER_Soil_ph,All_Sites, by = "Plot")

duplicates <- MER_Soil_ph %>%
  group_by(Plot) %>%
  filter(n() > 1) %>%
  ungroup()

MER_pH<-All_Sites %>% 
  left_join(MER_Soil_ph) %>% 
  select(Site,Plot,Soil_pH=pH)

saveRDS(MER_pH, "processed_data/nutrients/pH_MER.rds")
