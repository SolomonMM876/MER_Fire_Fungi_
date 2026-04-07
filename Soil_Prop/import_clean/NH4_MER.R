library(tidyverse)
library(readxl)
library(ggplot2)
library(purrr)

#Read in ortho_P data###########
# Define the folder path
folder_path <- "raw_data/Soil_data/Nutrients/NH4"

# Get a list of all file paths in the folder
file_list_NH <- list.files(folder_path, full.names = TRUE)
#file_list<-file_list_NH[3:15]#select the files we want

# Read and combine all files
All_NH4<- file_list_NH %>%
  set_names() %>%  # Preserve filenames for later
  map_dfr(read_csv, .id = "source_file") # Read files and add filename column

# Extract only the filename (not full path)
All_NH4_filtered <- All_NH4 %>%
  mutate(source_file = basename(source_file))%>%
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate(DIL = str_remove(DIL, "DIL|dil"),  # Remove 'DIL' from the dilution column
         DIL = as.numeric(DIL),  # Convert DIL to numeric
         `Manual Dil` = if_else(!is.na(DIL), `Manual Dil` * DIL, `Manual Dil`),
         Sample_Number = case_when( #replace values incorrectly input by Faizun below
           str_detect(Sample_Number, "655") ~ "655",
           str_detect(Sample_Number, "652") ~ "652",
           str_detect(Sample_Number, "354") ~ "354",
           TRUE ~ Sample_Number
         ))

#read in samples Faizun reran
# Faizun_NH4<-readRDS('raw_data/Soil_data/Nutrients/Resins_rerun/NH4_new_dat.Rdata')
# 
# All_NH4_filtered <- All_NH4_filtered %>%
#   filter(!Sample_Number %in% Faizun_NH4$Sample_Number) %>%
#   bind_rows(Faizun_NH4)



All_NH4_filtered %>% 
  filter(Result>2)->temp

temp
#These are the samples I pooled to have enough HCL
All_names<-read_excel("raw_data/Soil_data/Nutrients/Ortho_phos/Resin_Nutrient_IDs.xlsx",col_names = FALSE)%>%
  rename(ID=1, Names=2,Manual_Dil=4, Combined=3) %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia left me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT'),
         Names= str_replace(Names, 'WABEW15','WABEW18'))#there is no wabew15 only 18, must have been misswritten

MER_names_pooled<-read_excel("raw_data/Soil_data/Nutrients/MER_Combined.xlsx") %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia left me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT'),
         Names= str_replace(Names, 'WABEW15','WABEW18'))#there is no wabew15 only 18, must have been misswritten


#join by name col
unjoined_values<-MER_names_pooled%>% mutate(ID=as.factor(ID))%>%
  anti_join(All_NH4_filtered, by= c("ID"="Sample_Number"))


joined_values<-MER_names_pooled%>% mutate(ID=as.factor(ID))%>%
  left_join(All_NH4_filtered, by= c("ID"="Sample_Number")) %>%
  mutate( Total_Dil= `Manual Dil`) 


############


#read in nute sample info
Site_Plot_rep<-readRDS("raw_data/Soil_data/Nute_plot_info.RDS") 

Site_Plot<-Site_Plot_rep %>% 
  filter(!is.na(Nutrient_sub) &! Nutrient_sub=='0') %>% 
  group_by(Site,Plot) %>% 
  summarise(n=n(),
            Nutrient_sub=sum(Nutrient_sub, na.rm = TRUE),
            HCL= n*7.5,
            t=Nutrient_sub/n)

#join plot ID by sample ID sheet
MER_NH4<-Site_Plot%>%
  left_join(joined_values, by=c("Plot"="Names")) %>% 
  filter(!( !is.na(ID) & is.na(Result) ))

temp<-All_names %>% mutate(ID=as.factor(ID))%>%
  right_join(All_NH4_filtered, by= c("ID"="Sample_Number")) %>% 
  filter(ID %in% c(130:460)) %>% 
  filter(!str_detect(Names, 'BLANK')) %>% 
  mutate(Names=str_remove(Names,'.$')) %>% 
  rename(Plot=Names) %>% 
  left_join(Site_Plot, by=c('Plot'))%>%
  mutate( Total_Dil= `Manual Dil`*`Auto Dil`) 

combined_NH4 <- bind_rows(MER_NH4,temp) %>% 
  #remove one plot rep that is negative, use pooled sample instead
  filter(!(Plot == "WABVT13" &Result < 0))

# Step 1: Create a version of joined_values with duplicates averaged
summarised_dupes <- combined_NH4 %>%
  group_by(Site,Plot) %>%
  filter(n() > 1) %>%            # Only duplicated Plots
  filter(Result>0) %>%  #remove na values and remove negative results, which are weird anyway
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  left_join(distinct(combined_NH4 %>% select(-where(is.numeric)), Plot), by = "Plot") %>%
  relocate(Plot)  # Optional: keep ID as first column


# Step 2: Remove all duplicated rows (by ID) from the original dataframe
non_dupes <- combined_NH4 %>%
  filter(!(Plot %in% summarised_dupes$Plot))

# Step 3: Combine the non-duplicates and the averaged duplicates
joined_values_cleaned <- bind_rows(non_dupes, summarised_dupes) %>% 
  left_join(Site_Plot %>% select(Site_from_ref = Site, Plot), by = "Plot") %>%
  mutate(Site = coalesce(Site, Site_from_ref)) %>%
  select(-Site_from_ref,-Site.x,-Site.y)


#ID samples that are missing
missing<-Site_Plot%>%
# filter(str_detect(Names, "WAEW"))%>%
  anti_join(joined_values_cleaned,by=c("Plot"))%>%
  filter((!is.na(Nutrient_sub)|Nutrient_sub==0.00))

missing
####################################################

#calc mg/kg Ammon

LOQ_NH4NO3<-(.0757+.0587)/2

Final_NH4_MER<-joined_values_cleaned%>%  
  mutate(Dil_result = Result*Total_Dil,
         Amm_mg_kg= Dil_result *(HCL/Nutrient_sub)) %>%  #7.5mL used for 1 g of resin
  select(Site,Plot,n,Amm_mg_kg)


dupes <- Final_NH4_MER %>%
  group_by(Site,Plot) %>%
  filter(n() > 1)

temp<-Site_Plot %>% 
  anti_join(Final_NH4_MER)

library(ggplot2)
Final_NH4_summary <- Final_NH4_MER %>%
  group_by(Site, Plot) %>%
  summarise(mean = mean(as.numeric(Amm_mg_kg), na.rm = TRUE),
            se = (sd(as.numeric(Amm_mg_kg), na.rm = TRUE)) / sqrt(n())) %>%
  ungroup()

NH4_MER<-Final_NH4_summary %>% 
  ggplot() +
  geom_col(aes(x = Plot, y = mean))+
  geom_errorbar(aes(x=Plot,ymin = mean - se, ymax = mean + se), width = 0.2) +
  theme_bw() +
  labs(x = "Plot (Grouped by Site)", y = "Ammonia (mg/kg)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")  # Group by Site



saveRDS(Final_NH4_MER, "processed_data/nutrients/Ammonia_MER.rds")


library(patchwork)

#PO4_MER/NO3_MER/NH4_MER
