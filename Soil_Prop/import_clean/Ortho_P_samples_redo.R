library(tidyverse)
library(readxl)
library(ggplot2)
library(purrr)

# Define the folder path
folder_path <- "raw_data/Nutrients/Ortho_phos"

# Get a list of all file paths in the folder
file_list_o <- list.files(folder_path, full.names = TRUE)
# Exclude the specific file
file_to_exclude <- file.path(folder_path, c("rerun","Resin_Nutrient_IDs.xlsx"))
file_list <- setdiff(file_list_o, file_to_exclude)

# Read and combine all files
All_ortho_P <- file_list %>%
  set_names() %>%  # Preserve filenames for later
  map_dfr(read_csv, .id = "source_file") # Read files and add filename column

# Extract only the filename (not full path)
All_ortho_P_filtered <- All_ortho_P %>%
  mutate(source_file = basename(source_file)
  )%>%
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate(DIL = str_remove(DIL, "DIL|dil"),  # Remove 'DIL' from the dilution column
         DIL = as.numeric(DIL),  # Convert DIL to numeric
         `Manual Dil` = if_else(!is.na(DIL), `Manual Dil` * DIL, `Manual Dil`),
         Total_Dil= `Manual Dil`*`Auto Dil`)


temp<-All_ortho_P_filtered %>% 
  filter(Result>1)
All_ortho_P_filtered %>% 
  filter(Result>2)->temp


#Some samples had to be rerun because they were missing
# Define the folder path
folder_path_rerun <- "raw_data/Nutrients/Ortho_phos/rerun"

# Get a list of all file paths in the folder
file_list_redo <- list.files(folder_path_rerun, full.names = TRUE)

# Read and combine all files
redo_ortho_P <- file_list_redo %>%
  set_names() %>%  # Preserve filenames for later
  map_dfr(read_csv, .id = "source_file") %>% # Read files and add filename column
  mutate(source_file = basename(source_file)
  )%>%
  separate(`Sample ID`, into = c("Sample_Number", "DIL"), sep = "X|x", remove = FALSE) %>%
  mutate(DIL = str_remove(DIL, "DIL|dil"),  # Remove 'DIL' from the dilution column
         DIL = as.numeric(DIL),  # Convert DIL to numeric
         `Manual Dil` = if_else(!is.na(DIL), `Manual Dil` * DIL, `Manual Dil`),
         Total_Dil= `Manual Dil`*`Auto Dil`)



#select ID file
Ortho_names <- read_excel("raw_data/Nutrients/Ortho_phos/Resin_Nutrient_IDs.xlsx",col_names = FALSE)%>%
  rename(ID=1, Names=2,Manual_Dil=4, Combined=3) %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia --left me to write VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) 

#join by name col
unjoined_values<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  anti_join(All_ortho_P_filtered, by= c("ID"="Sample_Number"))%>% 
  filter(!Combined=='combined_MER') %>% #combined samples and created duplicates, but only relevant for avail N
  filter(!Manual_Dil=='10.9XManual_DIL')#dilluted samples for N

All_ortho_P<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  left_join(All_ortho_P_filtered, by= c("ID"="Sample_Number")) %>%#%>%select(`Names`,`Test Name`)%>%distinct()
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia --left me to write VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) %>% 
  group_by(ID) %>%
  slice_min(Result) %>% 
  ungroup()


redone_values<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  right_join(redo_ortho_P, by= c("ID"="Sample_Number")) %>%#%>%select(`Names`,`Test Name`)%>%distinct()
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia lef me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) %>% 
  filter(!is.na(Names))


#read in nute sample info
#Sample_ID
Nute_info<-readRDS('raw_data/Soil_dat/Nute_plot_info.RDS')

#join with redone samples
temp<-redone_values %>% 
  filter(Result<1) %>% 
  left_join(Nute_info %>% 
              select(Site,Plot) %>% 
              distinct(),  by=c("Names"="Plot")) 


#join with done samples
samps_redone<-All_ortho_P %>% 
  filter(Result>1) %>% 
  left_join(Nute_info,  by=c("Names"="Plot_Rep")) %>% 
  anti_join(temp, by= c('Plot'='Names')) %>% 
  mutate(Label=if_else(is.na(Plot),Names,Plot)) %>% 
  distinct() %>% 
  arrange(desc(Label)) %>% 
  # left_join(Ortho_names, by=c('Label'='Names'))
  
  
  MER_Redo<-samps_redone %>% 
  filter(!is.na(Site)) %>% 
  group_by(Label) %>% 
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

The_rest <-samps_redone %>% 
  filter(is.na(Site)) %>% 
  select(colnames(MER_Redo))

samps<-rbind(The_rest,MER_Redo) %>% 
  left_join(Ortho_names, by=c('Label'='Names'))

#sampes that need to be redone
samps_redone

write.csv(samps_redone, 'raw_data/Initial_Concentrated_Samps.csv',row.names = FALSE)