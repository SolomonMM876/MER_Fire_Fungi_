library(tidyverse)


PO4<- read_csv("raw_data/Soil_data/Nutrients/Resins_rerun/O-PHOS-SOL-RR-7S-20250606.csv")


temp<-PO4%>%
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate( DIL = as.numeric(DIL)) %>%   # Convert DIL to numeric
  filter(!str_detect(Sample_Number,'STANDARD|C C')) %>% 
  mutate(Adj_result=DIL*Result)%>% select(-c(`Test Name`, Time, Units))

#Run Ortho P script
t<-inner_join(All_ortho_P_filtered %>% select(-c(source_file,`Test Name`,DIL, Time, Units)),temp, by = join_by(Sample_Number))

plot(t$Result.x,t$Adj_result) #not a good relationship


NH4_NO3<- read_csv("raw_data/Soil_data/Nutrients/Resins_rerun/NH4-NO3-SOLO-RR-7S.csv")

Faizun_NH4<-NH4_NO3%>%
  mutate(`Sample ID`=if_else(`Sample ID`=='SOLO-1014X3 (NH)','1014X3',`Sample ID`)) %>% 
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate( DIL = as.numeric(DIL)) %>%   # Convert DIL to numeric
  filter(!str_detect(Sample_Number,'STANDARD|CC|C C')) %>% 
  mutate(Adj_result=DIL*Result)%>% select(-c(`Sample Details`)) %>% 
  filter(Result>0)
saveRDS(Faizun_NH4,'raw_data/Soil_data/Nutrients/Resins_rerun/NH4_new_dat.Rdata')
