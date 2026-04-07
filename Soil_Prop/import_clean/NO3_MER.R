library(tidyverse)
library(readxl)
library(ggplot2)
library(purrr)

#Read in ortho_P data###########
# Define the folder path
folder_path <- "raw_data/Soil_data/Nutrients/NO3"

# Get a list of all file paths in the folder
file_list_NO3 <- list.files(folder_path, full.names = TRUE)
file_to_exclude <- file.path(folder_path, c("NO3_rerun","NO3-SOLOMON-L8-RR FOR HIGH DIL.csv"))
file_list <- setdiff(file_list_NO3, file_to_exclude)


# Read and combine all files
All_NO3<- file_list %>%
  set_names() %>%  # Preserve filenames for later
  map_dfr(read_csv, .id = "source_file") # Read files and add filename column

# Extract only the filename (not full path)
All_NO3_filtered <- All_NO3 %>%
  mutate(source_file = basename(source_file))%>%
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate(DIL = str_remove(DIL, "DIL|dil|D|\\s*\\(A\\)"),  # Remove 'DIL' from the dilution column
          DIL = str_remove(DIL, "\\s*\\(A\\)"),  # Remove '(A)' from the dilution column
         DIL = as.numeric(DIL),  # Convert DIL to numeric
         `Manual Dil` = if_else(!is.na(DIL), `Manual Dil` * DIL, `Manual Dil`),
         Sample_Number = case_when( #replace values incorrectly input by Faizun below
           str_detect(Sample_Number, "655") ~ "655",
           str_detect(Sample_Number, "652") ~ "652",
           str_detect(Sample_Number, "354") ~ "354",
           TRUE ~ Sample_Number
         ))


#select ID file
NO3_names <- read_excel("raw_data/Soil_data/Nutrients/Ortho_phos/Resin_Nutrient_IDs.xlsx",col_names = FALSE)%>%
  rename(ID=1, Names=2,Manual_Dil=4, Combined=3) %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia left me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT'),
         Names= str_replace(Names, 'WABEW15','WABEW18'))#there is no wabew15 only 18, must have been misswritten

#quick check for neg values

temp<-All_NO3_filtered %>% 
  filter(Result<0) %>% 
  left_join(NO3_names %>% mutate(ID=as.factor(ID)), by= c("Sample_Number"="ID"))


MER_names_pooled<-read_excel("raw_data/Soil_data/Nutrients/MER_Combined.xlsx") %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia left me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT'))

#join by name col
unjoined_values<-NO3_names%>% mutate(ID=as.factor(ID))%>%
  anti_join(All_NO3_filtered, by= c("ID"="Sample_Number"))


joined_values<-MER_names_pooled%>% mutate(ID=as.factor(ID))%>%
  left_join(All_NO3_filtered, by= c("ID"="Sample_Number")) 



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
MER_NO3<-Site_Plot%>%
  left_join(joined_values, by=c("Plot"="Names")) %>% 
  filter(!( !is.na(ID) & is.na(Result) ))

temp<-NO3_names %>% mutate(ID=as.factor(ID))%>%
  right_join(All_NO3_filtered, by= c("ID"="Sample_Number")) %>% 
  filter(ID %in% c(130:460)) %>% 
  filter(!str_detect(Names, 'BLANK')) %>% 
  mutate(Names=str_remove(Names,'.$')) %>% 
  rename(Plot=Names) %>% 
  left_join(Site_Plot)

combined_NO3 <- bind_rows(MER_NO3,temp) %>% 
  filter( !is.na(Result) ) #remove empty rows
  
# Step 1: Create a version of joined_values with duplicates averaged
summarised_dupes <- combined_NO3 %>%
  group_by(ID) %>%
  filter(n() > 1) %>% 
  #remove blanks
  filter(!is.na(source_file)) %>% 
  #remove diluted samples, values do not allign well and better to have similar runs instead of averaging
  filter(is.na(DIL))%>%            # Only duplicated IDs
  summarise(Result=mean(Result)) %>% 
 left_join(distinct(combined_NO3 %>% select(Site,Plot,ID)), by = "ID") %>%
  relocate(ID)  # Optional: keep ID as first column

# Step 2: Remove all duplicated rows (by ID) from the original dataframe
non_dupes <- combined_NO3 %>%
  filter(!(ID %in% summarised_dupes$ID))

# Step 3: Combine the non-duplicates and the averaged duplicates
joined_values_cleaned <- bind_rows(non_dupes, summarised_dupes) %>% 
  left_join(Site_Plot %>% select(Site_from_ref = Site, Plot), by = "Plot") %>%
  mutate(Site = coalesce(Site, Site_from_ref)) 

#Repeat this step but at plot level
dupes <- joined_values_cleaned %>%
  group_by(Plot) %>%
  filter(n() > 1) %>% 
  filter(!is.na((ID))) %>% 
  filter(!( !is.na(DIL) & Result < 0 )) %>% #remove samples that were dillutted and and results became too low
  filter(!is.na(Note)) %>% 
  group_by(Plot) %>%
  summarise(Result=mean(Result)) %>% 
  left_join(distinct(combined_NO3 %>% select(Site,Plot)), by = "Plot") %>%
  ungroup()

# Step 2: Remove all duplicated rows (by ID) from the original dataframe
non_dupes_2 <- joined_values_cleaned %>%
  filter(!(Plot %in% dupes$Plot))

# Step 3: Combine the non-duplicates and the averaged duplicates
joined_values_final <- bind_rows(non_dupes_2, dupes) %>% 
  filter(!is.na(Result))

dupes <- joined_values_final %>%
  group_by(Plot) %>%
  filter(n() > 1) 

#add empty nutrient values
joined_values_final <- joined_values_final %>%
  left_join(Site_Plot %>% select(Plot, Nutrient_sub_site = Nutrient_sub,HcL_Site=HCL), by = c("Site","Plot"))%>%
  mutate(Nutrient_sub = coalesce(Nutrient_sub, Nutrient_sub_site),
         HCL = coalesce(HCL, HcL_Site)) %>%
  select(-Nutrient_sub_site,-HcL_Site) 
############

#ID samples that are missing
missing<-Site_Plot%>%
# filter(str_detect(Names, "WAEW"))%>%
  anti_join(joined_values_final,by=c("Plot"))%>%
  filter((!is.na(Nutrient_sub)|Nutrient_sub==0.00))

missing
#Also missing SAM105
####################################################

#calc mg/kg NO3
#select the blank values from raw data file
Blank_value<-All_NO3_filtered%>%
  #Adjust as needed
  filter( str_detect(`Sample ID`, "BLANK"))%>%
  rename(Plot= `Sample ID`)%>%
  summarise(Blank_avg = mean(Result, na.rm = TRUE))%>%
  pull()

LOQ_NH4NO3<-(0.00070+0.002)/2

Final_NO3<-joined_values_final%>%  
  mutate(NO3_dil = Result,#Manual Dil is 1 for all cols
         NO3_dil = ifelse(NO3_dil < 0, LOQ_NH4NO3, NO3_dil),
         NO3_mg_kg= NO3_dil *(HCL/Nutrient_sub)) #7.5mL used for 1 g of resin

Final_NO3 %>% 
  select(Site,Plot,n,NO3_mg_kg)->Final_NO3_MER



library(ggplot2)
Final_NO3_summary <- Final_NO3_MER %>%
  group_by(Site, Plot) %>%
  summarise(mean = mean(as.numeric(NO3_mg_kg), na.rm = TRUE),
            se = (sd(as.numeric(NO3_mg_kg), na.rm = TRUE)) / sqrt(n())) %>%
  ungroup()

NO3_MER<-Final_NO3_summary %>% 
  ggplot() +
  geom_col(aes(x = Plot, y = mean))+
  geom_errorbar(aes(x=Plot,ymin = mean - se, ymax = mean + se), width = 0.2) +
  theme_bw() +
  labs(x = "Plot (Grouped by Site)", y = "Nitrate (mg/kg)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")  # Group by Site

NO3_MER

saveRDS(Final_NO3_MER, "processed_data/nutrients/NO3_MER.rds")



