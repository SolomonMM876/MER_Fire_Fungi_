library(tidyverse)


#read in Auscribe year 1 data already extract from veg_import_clean
auscribe_veg1<-readRDS('Processed_data/metadata/veg/auscribe_veg_yr1.Rdata')


#read in combined auscribe datasets extract from veg_import_clean
auscribe_veg<-readRDS('Processed_data/metadata/veg/auscribe_veg_both.Rdata')


temp<-anti_join(auscribe_veg,auscribe_veg1)

library(APCalign)
#load for austraits
tax_resources <- load_taxonomic_resources()

# Extract distinct species names as a character vector
spp_names <- auscribe_veg %>%
  distinct(herbarium_determination) %>%
  filter(!is.na(herbarium_determination)) %>% 
  pull(herbarium_determination) %>%  # extract as a character vector
  iconv(from = "", to = "UTF-8", sub = "byte")  # Convert to valid UTF-8


# Align the species names using APC resources
aligned_spp <- align_taxa(spp_names, resources = tax_resources)

upd_spp_list<-update_taxonomy(aligned_spp, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
  select(original_name,suggested_name,aligned_name,genus:taxon_rank)%>%
  filter(!is.na(suggested_name))#remove indeterminate IDs

#read table of myco IDs
brundrett_2017 <- read_csv("Processed_data/metadata/veg/datasets_published/brundrett_2017_myco.csv")%>%
  #correct family names so joins work
  mutate(Family = str_replace(Family, "^(Euphorbaceae|Fabaceae)\\b.*", "\\1"),
         Family = str_replace(Family, "Euphorbaceae", "Euphorbiaceae"))

#read in North american classification dataset
Brundrett_Tedersoo_2020.csv<- read_csv("Processed_data/metadata/veg/datasets_published/Brundrett_Tedersoo_2020.csv")

library(readxl)
#read in global fungalroot db for myco association Nadejda A. Soudzilovskaia et al 2020
Soudzilovskaia2020<-read_excel('Processed_data/metadata/veg/datasets_published/Soudzilovskaia_etal_2020.xlsx', sheet = 'Table S2',skip = 1)


# First join: Match trait data based on Genus, some taxa are dual, but we just care if they are or are not
genus_myco <- upd_spp_list %>%
  left_join(
    brundrett_2017 %>%
      rename(Myco_type = Primary.Obs) %>%
      group_by(Genus) %>%
      summarise(Myco_type = paste(unique(Myco_type), collapse = " | "), .groups = "drop"),  
    by = c("genus" = "Genus")
  ) %>%
  mutate(
    Myco_type = case_when(
      genus %in% c("Thysanotus", "Cassytha", "Olax",'Anarthria', 'Crassula') ~ "NM",
      genus %in% c('Doryanthes',"Brunoniella",'Caesia') ~ "Unkown",
      TRUE ~ Myco_type),
    Notes = case_when(
      genus == "Thysanotus" ~ "Unique sub-epidermal association - McGee 1988, Brundrett & Abbott 1991.",
      genus == "Cassytha" ~ "Can have members of family that are mycorrhizal, but all species observed are hemi-parasitic growing on stems and are NM as such.",
      genus == "Brunoniella" ~ "Type not determined due to lack of records, removed from analysis.",
      #taxon_name == "Doryanthes excelsa" ~ "Type not determined due to lack of records, classified as NM for analysis.",
      #taxon_name == "Caesia sp." ~ "Type not determined due to lack of records, classified as NM for analysis.",
      genus == "Olax" ~ "Non-mycorrhizal Bellgard et al. 1994",
      genus == "Anarthria" ~ "Non-myorrhizal Brundrett 2009 (https://mycorrhizas.info/nmplants.html)",
      genus == "Crassula" ~ "Non-myorrhizal Brundrett 2009 Most samples have NM roots, some few examples of VAM *More sampling required.(https://mycorrhizas.info/nmplants.html)",
      is.na(Myco_type) ~ "Classification based on other genera in family",
      TRUE ~ ""))

#join remaining based on family
family_myco <- genus_myco %>%
  # Join to Brundrett data again, this time by Family
  left_join(
    brundrett_2017 %>%
      rename(Myco_type_fam = Primary.Obs) %>%
      group_by(Family) %>%
      summarise(
        Myco_type_fam = paste(unique(Myco_type_fam), collapse = " | "),
        .groups = "drop"
      ),
    by = c("family" = "Family")  # Join on family names
  ) %>%
  # If genus-level Myco_type is still NA, use the family-level classification
  mutate(Myco_type = coalesce(Myco_type, Myco_type_fam)) %>%
  # Drop the intermediate family trait column
  select(-Myco_type_fam) 

#combine classifications
myco_tax<-family_myco%>%
  left_join(genus_myco) %>% 
  rename(herbarium_determination=original_name)

#ID unknown plant taxa
unk_myco_plants<-myco_tax%>%
  filter(is.na(Myco_type))%>%
  distinct()

#ID unknown myco types by Soudzilovskaia et al 2020
Soudzilovskaia_genus_myco<-unk_myco_plants %>% 
  left_join(Soudzilovskaia2020, by=c('genus'='Genus')) %>% 
  filter(!is.na(`Mycorrhizal type`)) %>% 
  mutate(Myco_type = coalesce(Myco_type, `Mycorrhizal type`),
         Myco_type=if_else(Myco_type=='uncertain','Unkown',Myco_type),
         Myco_type=if_else(Myco_type=='AM','VAM',Myco_type),
         Notes= 'Classification based Soudzilovskaia et al 2020 FungalRoot databse') %>% 
  select(-`Mycorrhizal type`)


# Remove original unknown rows and add updated identifications
myco_tax_final <- myco_tax %>%
  filter(!(is.na(Myco_type)) ) %>%
  bind_rows(Soudzilovskaia_genus_myco)

#quick check
myco_tax_final%>%filter(is.na(Myco_type))%>%distinct()

myco_tax_final%>%filter(Myco_type=='Unkown')%>%distinct()

saveRDS(myco_tax_final,'processed_data/metadata/veg/myco_host_df.rds')
#######################count Frequency myco hosts##############################################################

# Step 1: Join vegetation with mycorrhizal classifications, remove 'Unknown' types
auscribe_myco <- auscribe_veg %>% 
  left_join(myco_tax_final) %>% 
  filter(!str_detect(Myco_type, "Unkown"))

# Step 2: Summarize per plot
myco_freq_summary <- auscribe_myco %>% 
  rename(Plot = plot_join) %>% 
  group_by(Site, Plot) %>%
  summarise(
    total_species = sum(!is.na(aligned_name)),  # total occurrences
    total_myco_hosts = sum(str_detect(Myco_type, "VAM|Ericoid|ECM"), na.rm = TRUE),
    AM_hosts = sum(str_detect(Myco_type, "VAM"), na.rm = TRUE),
    ECM_hosts = sum(str_detect(Myco_type, "ECM"), na.rm = TRUE),
    non_myco_host = sum(!str_detect(Myco_type, "VAM|Ericoid|ECM"), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    freq_total_myco = total_myco_hosts / total_species,
    freq_non_host= non_myco_host/total_species,
    freq_AM = AM_hosts / total_species,
    freq_ECM = ECM_hosts / total_species
  )

All_Sites<-readRDS("raw_data/MER_Site_data/All_Sites.RDS")

#Two sites in broome were not samples WAMDALRainfo01 and WAMDALEucWoo02
missing_plots<-anti_join(All_Sites,myco_freq_summary, by='Plot') %>% distinct(Site,Plot)

missing_plots %>% 
  left_join(All_Sites)

# THIS CREATES ERRORS AND IS NOT WORTHWHILE
# #to fill in gaps in first year data collection I have taken data from second data collection, which exists for some plots not sampled in first yr
# #load year 2 data
# 
# #load Auscribe_JoinedTables df's YEAR 2 already extract from veg_import_clean
# auscribe_veg_yr2<-readRDS('Processed_data/metadata/veg/auscribe_veg_yr2.Rdata')
# 
# #extract missing data from plots
# missing_plots %>% inner_join(auscribe_veg_yr2, by='Plot') %>% distinct(Site.x,Plot)
# 
# auscribe_veg_yr2_miss_plot<-missing_plots %>% inner_join(auscribe_veg_yr2, by='Plot') %>% 
#   select(-Site.y) %>% rename(Site=Site.x)
# 
# # Extract distinct species names as a character vector again
# spp_names_y2 <- auscribe_veg_yr2_miss_plot %>%
#   distinct(herbarium_determination) %>%
#   filter(!is.na(herbarium_determination)) %>% 
#   pull(herbarium_determination) %>%  # extract as a character vector
#   iconv(from = "", to = "UTF-8", sub = "byte")  # Convert to valid UTF-8
# 
# 
# # Align the species names using APC resources
# aligned_spp_y2 <- align_taxa(spp_names_y2, resources = tax_resources)
# 
# upd_spp_list_y2<-update_taxonomy(aligned_spp_y2, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
#   select(original_name,suggested_name,aligned_name,genus:taxon_rank)%>%
#   filter(!is.na(suggested_name))#remove indeterminate IDs
# 
# anti_join(upd_spp_list_y2,upd_spp_list)
# 
# myco_tax_y2<-left_join(upd_spp_list_y2,myco_tax_final %>%
#                          distinct(genus,Myco_type)) %>% 
#   rename(herbarium_determination=original_name)
# 
# 
# #REPEAT Step 1: Join vegetation with mycorrhizal classifications, remove 'Unknown' types
# auscribe_myco_y2 <- auscribe_veg_yr2 %>% 
#   left_join(myco_tax_y2) %>% 
#   filter(!str_detect(Myco_type, "Unkown"))
# 
# # Step 2: Summarize per plot
# myco_freq_summary_y2 <- auscribe_myco_y2 %>%
#   group_by(Site,Plot) %>% # add plot here to get by plot level
#   summarise(
#     total_species = n_distinct(aligned_name),
#     total_myco_hosts = n_distinct(aligned_name[str_detect(Myco_type, "VAM|Ericoid|ECM")]),
#     AM_hosts = n_distinct(aligned_name[str_detect(Myco_type, "VAM")]),
#     ECM_hosts = n_distinct(aligned_name[str_detect(Myco_type, "ECM")]),
#     non_myco_host = sum(!str_detect(Myco_type, "VAM|Ericoid|ECM"), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     freq_total_myco = total_myco_hosts / total_species,
#     freq_non_host= non_myco_host/total_species,
#     freq_AM = AM_hosts / total_species,
#     freq_ECM = ECM_hosts / total_species
#   ) %>% 
#   filter(Plot=='WASP02')

####year 2 data is from WAMESPShrubl01 
###############
#####Combine y1 and y2 data
#myco_freq_final_sum_plot<-bind_rows(myco_freq_summary,myco_freq_summary_y2)

myco_freq_final_sum_plot<-myco_freq_summary

#  Summarize by site
myco_freq_final_site<- myco_freq_final_sum_plot %>%
  group_by(Site) %>% # add plot here to get by plot level
  summarise(
      total_species = sum(total_species, na.rm = TRUE),
      total_myco_hosts = sum(total_myco_hosts, na.rm = TRUE),
      total_non_myco_hosts = sum(non_myco_host, na.rm = TRUE),
      AM_hosts = sum(AM_hosts, na.rm = TRUE),
      ECM_hosts = sum(ECM_hosts, na.rm = TRUE)
  ) %>% 
  mutate(
    freq_total_myco = total_myco_hosts / total_species,
    freq_non_host= total_non_myco_hosts/total_species,
    freq_AM = AM_hosts / total_species,
    freq_ECM = ECM_hosts / total_species
  ) 


missing_plots<-anti_join(All_Sites,myco_freq_final_sum_plot, by='Plot') %>% distinct(Site,Plot,Fire_Treatment)

missing_plots_summary <- missing_plots %>% 
  group_by(Site, Fire_Treatment) %>% 
  count() %>% 
  pivot_wider(
    names_from = Fire_Treatment,
    values_from = n,
    values_fill = 0   # fill missing combos with 0
  )
missing_plots_summary

saveRDS(missing_plots_summary,'Processed_data/metadata/veg/Plots_w_no_veg_data.rds')


saveRDS(myco_freq_final_sum_plot,'Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata')
saveRDS(myco_freq_final_site,'Processed_data/metadata/veg/myco_freq_Site_level_df.Rdata')



###############################################################
#Calculate difference in frquency of myco hosts between burnt and unburnt plots

myco_host_dif<-myco_freq_final_sum_plot %>% 
  left_join(All_Sites %>% select(Site,Plot,Fire_Treatment)) %>% 
  group_by(Fire_Treatment,Site) %>% 
  summarise(
    total_species = sum(total_species, na.rm = TRUE),
    total_myco_hosts = sum(total_myco_hosts, na.rm = TRUE),
    total_non_myco_hosts = sum(non_myco_host, na.rm = TRUE),
    AM_hosts = sum(AM_hosts, na.rm = TRUE),
    ECM_hosts = sum(ECM_hosts, na.rm = TRUE)
  ) %>%
  mutate(
    freq_total_myco = total_myco_hosts / total_species,
    freq_non_host=total_non_myco_hosts/total_species,
    freq_AM = AM_hosts / total_species,
    freq_ECM = ECM_hosts / total_species
  )

myco_host_dif %>% 
  ggplot(aes(x=Fire_Treatment,y=freq_total_myco))+
  geom_col()+
  facet_wrap(~Site)

myco_host_dif %>% 
  ggplot(aes(x=Fire_Treatment,y=freq_AM))+
  geom_col()+
  facet_wrap(~Site)

myco_host_dif %>% 
  ggplot(aes(x=Fire_Treatment,y=freq_ECM))+
  geom_col()+
  facet_wrap(~Site)

myco_host_dif %>% 
  ggplot(aes(x=Fire_Treatment,y=freq_non_host))+
  geom_col()+
  facet_wrap(~Site)


#Pivot wider to get both treatments in same row
myco_host_diff_wide <- myco_host_dif %>%
  select(Site, Fire_Treatment, freq_total_myco, freq_non_host, freq_AM, freq_ECM) %>%
  pivot_wider(
    names_from = Fire_Treatment,
    values_from = c(freq_total_myco, freq_non_host, freq_AM, freq_ECM),
    names_sep = "_"
  )

#select host freq in unburnt plots
host_freq_unburnt<-myco_host_diff_wide %>% 
  select(Site,freq_AM_U,freq_ECM_U,freq_total_myco_U,freq_non_host_U)


saveRDS(host_freq_unburnt,'Processed_data/metadata/veg/myco_freq_Site_unburnt_plots.Rdata')


# Step 4: Calculate differences (Burnt - Unburnt)
myco_host_diff_final <- myco_host_diff_wide %>%
  mutate(
    diff_total_myco = freq_total_myco_U- freq_total_myco_B,
    diff_non_host   = freq_non_host_U - freq_non_host_B,
    diff_AM         = freq_AM_U - freq_AM_B,
    diff_ECM        = freq_ECM_U - freq_ECM_B
  ) %>% 
  select(Site,diff_total_myco,diff_non_host,diff_AM,diff_ECM)

saveRDS(myco_host_diff_final,'Processed_data/metadata/veg/myco_freq_Site_differences.Rdata')

