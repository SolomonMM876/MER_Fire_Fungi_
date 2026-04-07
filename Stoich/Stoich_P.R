library(tidyverse)
library(readxl)
library(ggplot2)

#Data from kazi for Total P
Total_P <- read_excel("raw_data/Stoich/Tube_biomass_MER.xlsx")

#IDs for Tube data from Kazi
Tube_label <- read_excel("raw_data/Tube_biomass_MER_ID.xlsx") %>% 
  select(Tube_ID,Total_P_ID) %>% 
  mutate(Tube_ID=as.factor(Tube_ID),
         Total_P_ID=as.factor(Total_P_ID))


Standards_P<-Total_P %>% 
  filter(str_detect(`Sample ID_main`,'Apple')) %>% 
  rename(Tube_ID=`Tube no`,weight_mg=`Weight (mg)`,percent_P=`P (%)`) %>% 
  mutate(standard='Apple_leaf',
         Tube_ID=as.factor(Tube_ID))

Standards_P %>% 
  ggplot(aes(x=`Sample ID_main`, y= percent_P)) +
           geom_point()

Standards_P %>% 
  ggplot(aes(x=standard, y= percent_P)) +
  geom_violin(draw_quantiles = TRUE)+
  geom_jitter()

MER_TP<-Total_P%>% 
  filter(`Project Name`=='MER') %>% 
  rename(Total_P_ID=`Sample ID_main`,weight_mg=`Weight (mg)`,percent_P=`P (%)`) %>% 
  select(Total_P_ID,weight_mg,percent_P) %>% 
  mutate(Total_P_ID=as.factor(Total_P_ID)) %>% 
  left_join(Tube_label) %>% 
  filter(!Total_P_ID=='Blank')


Site_Tube_info<-readRDS("raw_data/Soil_data/Nute_plot_info.RDS") 

MER_P<-Site_Tube_info %>% 
  distinct(Site,Plot,Tube_ID) %>% 
  filter(!is.na(Tube_ID)) %>% 
  left_join(MER_TP) %>% 
  filter(!is.na(percent_P))


MER_P %>% 
  ggplot(aes(x=weight_mg , y= percent_P)) +
  geom_point()

MER_P %>% 
  ggplot(aes(x=Site, y= percent_P)) +
  geom_violin(drop=FALSE)+
  geom_jitter(width=.1)

saveRDS(MER_P, 'Processed_data/Stoich/Total_P_MER.Rdata')

