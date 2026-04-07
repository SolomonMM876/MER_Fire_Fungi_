library(tidyverse)

#PO4
MER_Ortho<-readRDS('Processed_data/nutrients/Otho_P_dat.RDS')
MER_Ortho_Sum<-readRDS('Processed_data/nutrients/Otho_P_summary.RDS') %>% 
  select(-n)

#NH4
MER_NH4<-readRDS("processed_data/nutrients/Ammonia_MER.rds")%>% 
  select(-n)

#NO3
MER_NO3<-readRDS("processed_data/nutrients/NO3_MER.rds")%>% 
  select(-n)

#pH
MER_pH<-readRDS("processed_data/nutrients/pH_MER.rds")


MER_nute<-MER_pH %>% 
  left_join(MER_NH4 )%>% 
  left_join(MER_Ortho_Sum) %>% 
  left_join(MER_NO3) %>% 
  rename(PO4=PO4_mg_kg,NH4=Amm_mg_kg,NO3=NO3_mg_kg) %>% 
  select(-se_ortho)

# Step 1: Identify sites with any NA in nutrient values
sites_with_na <- MER_nute %>%
  group_by(Site) %>%
  filter(any(is.na(PO4), is.na(NO3), is.na(NH4))) %>%
  pull(Site) %>%
  unique()

PO4<-MER_nute %>%
  filter(Site %in% sites_with_na) %>%
  ggplot() +
  geom_col(aes(x = Plot, y = PO4)) +
  theme_bw() +
  labs(x = "Plot (Grouped by Site)", y = "Orthophos (mg/kg)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

NO3<-MER_nute %>%
  filter(Site %in% sites_with_na) %>%
  ggplot() +
  geom_col(aes(x = Plot, y = NO3)) +
  theme_bw() +
  labs(x = "Plot (Grouped by Site)", y = "Nitrate (mg/kg)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


NH4<-MER_nute %>%
  filter(Site %in% sites_with_na) %>%
  ggplot() +
  geom_col(aes(x = Plot, y = NH4)) +
  theme_bw() +
  labs(x = "Plot (Grouped by Site)", y = "Ammonia (mg/kg)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(patchwork)

PO4/NO3/NH4

saveRDS(MER_nute, 'Processed_data/nutrients/MER_All_Soil.rds')





####################
auscribe_plots<-readRDS('Processed_data/metadata/veg/auscribe_veg_both.Rdata') %>% 
  distinct(Site,plot,plot_join) %>% rename(Plot=plot_join)


MER_nute_plot<-auscribe_plots %>% 
  left_join(MER_nute) %>% 
  rename(Sol_plot=Plot,auscribe_plot=plot)

write.csv(MER_nute_plot, file='Processed_data/nutrients/Soil_properties_MER.csv',row.names = FALSE)
