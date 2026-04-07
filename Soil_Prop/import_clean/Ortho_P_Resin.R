library(tidyverse)
library(readxl)
library(ggplot2)
library(purrr)

# Define the folder path
folder_path <- "raw_data/Soil_data/Nutrients/Ortho_phos"

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
         Total_Dil= `Manual Dil`)
           

temp<-All_ortho_P_filtered %>% 
  filter(Result>1)
All_ortho_P_filtered %>% 
  filter(Result>2)->temp


#Some samples had to be rerun because they were missing
# Define the folder path
folder_path_rerun <- "raw_data/Soil_data/Nutrients/Ortho_phos/rerun"

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
         Total_Dil= `Manual Dil`)

#select ID file
Ortho_names <- read_excel("raw_data/Soil_data/Nutrients/Ortho_phos/Resin_Nutrient_IDs.xlsx",col_names = FALSE)%>%
  rename(ID=1, Names=2,Manual_Dil=4, Combined=3) %>% 
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia --left me to write VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) 

#join by name col see wehat is not joined = 0 
unjoined_values<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  anti_join(All_ortho_P_filtered, by= c("ID"="Sample_Number"))%>% 
  filter(!Combined=='combined_MER') %>% #combined samples and created duplicates, but only relevant for avail N
  filter(!Manual_Dil=='10.9XManual_DIL')#dilluted samples for N

#combine names and redout for all samps
All_ortho_P<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  left_join(All_ortho_P_filtered, by= c("ID"="Sample_Number")) %>%#%>%select(`Names`,`Test Name`)%>%distinct()
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia --left me to write VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) %>% 
  group_by(ID) %>%
  slice_min(Result) %>% 
  ungroup()

#samples that were redone to check machine accuracy
redone_values<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  inner_join(redo_ortho_P, by= c("ID"="Sample_Number")) %>%#%>%select(`Names`,`Test Name`)%>%distinct()
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia lef me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) 


#read in nute sample info
#Sample_ID
Nute_info<-readRDS('raw_data/Soil_data/Nute_plot_info.RDS')

#join with redone samples
redone_nute<-redone_values %>% 
  filter(Result<1) %>% 
  left_join(Nute_info %>% 
              select(Site,Plot) %>% 
              distinct(),  by=c("Names"="Plot")) 


#Samples I thought had to be redone
# samps_redone<-All_ortho_P %>% 
#   filter(Result>1) %>% 
#   left_join(Nute_info,  by=c("Names"="Plot_Rep")) %>% 
#   anti_join(redone_nute, by= c('Plot'='Names')) %>% 
#   mutate(Label=if_else(is.na(Plot),Names,Plot)) %>% 
#   distinct() %>% 
#   arrange(desc(Label))

################################################



#join by name col
Ortho_labeled<-Ortho_names%>% mutate(ID=as.factor(ID))%>%
  mutate(Names = str_replace(Names, "VCFR", "VCRF"), #my dyslexia lef me to right VCFR instead of the correct VCRF
         Names= str_replace(Names, 'WABUT','WABVT')) %>% 
  left_join(All_ortho_P)

Missing_values<-Ortho_labeled%>%
  filter(is.na(Result)&!is.na(Names)) %>% 
  filter(!Combined=='combined_MER'|is.na(Combined))#combined samples and created duplicates
  

#######################################################################################

#join by name col
MER_Ortho<-Nute_info%>%
  left_join(Ortho_labeled, by= c("Plot_Rep"="Names")) 

#select unmatching values
unmatched<-MER_Ortho%>%
  filter(is.na(Result)&!is.na(Nutrient_sub)) %>% 
  #remove QDR samples nute_sub=0
  filter(!Plot=='QDRF06')

temp<-unmatched %>% 
  left_join(redone_nute)


#Now add a count of the number of bags collected from each plot, this is to calculate how much HCl we used
Ortho_P_MER<-MER_Ortho%>%
  #remove empty reps
  filter(!Nutrient_sub=='0') %>% 
  filter(!is.na(source_file)) %>% 
  left_join(MER_Ortho%>%
              #remove empty reps
              filter(!Nutrient_sub=='0') %>% 
              filter(!is.na(source_file)) %>% 
              group_by(Plot) %>%
              summarise(n= n()))



#This is your final df, note that the duplicated row now only have one value for Result and the rest of the values
#are na, this shouldnt be an issue for the rest of your calculations
Ortho_P_MER

#select the blank values from raw data file
Blank_value<-All_ortho_P%>%
  #Adjust as needed
  filter( str_detect(Names, "BLANK|Blank")) %>% 
  rename(Plot= `Sample ID`)%>%
  summarise(Blank_avg = mean(Result, na.rm = TRUE))%>%
  pull()

Final_Ortho_numbers<-Ortho_P_MER%>%  
  mutate( Dil_result= Result* Total_Dil,
          Ortho_blanked = (Dil_result-Blank_value),
          Ortho_blanked = ifelse(Ortho_blanked < 0, 0.0288/2, Ortho_blanked),
          Ortho_resin= Ortho_blanked *(7.5/Nutrient_sub )) %>% #7.5mL used for 1 g of resin
  select(Site:Rep,ID, n,Ortho_resin,Nutrient_sub)

#some samples did not have resins or biomass so these samples have no biomass as well
Ortho_all_dat<-Final_Ortho_numbers %>% 
  filter(!Nutrient_sub<0.1)

# Step 1: Summarize the data
Ortho_summary <- Final_Ortho_numbers %>%
  group_by(Plot,Site,n) %>%
  summarise(
    PO4_mg_kg = mean(Ortho_resin, na.rm = TRUE),
    se_ortho = sd(Ortho_resin, na.rm = TRUE) / sqrt(n())
  )


temp<-Nute_info %>% select(-Rep,-Nutrient_sub) %>% distinct() %>% 
  anti_join(Ortho_summary)
# Step 2: Plot
PO4_MER<-Ortho_summary %>% 
  filter(PO4_mg_kg<1000) %>% 
  ggplot(aes(x = Plot, y = PO4_mg_kg)) +
  geom_col() +
  geom_errorbar(aes(ymin = PO4_mg_kg - se_ortho, ymax = PO4_mg_kg + se_ortho), 
                width = 0.2) +
  labs(x = "Plot (Grouped by Site)", y = "PO4(mg/kg)") +
  facet_grid(. ~ Site, scales = "free_x", space = "free_x")+  # Group by Site
  theme(axis.text.x = element_text(angle = 90, hjust = 0))


PO4_MER


saveRDS(Final_Ortho_numbers,'Processed_data/nutrients/Otho_P_dat.RDS')
saveRDS(Ortho_summary,'Processed_data/nutrients/Otho_P_summary.RDS')

