#libraries
library(tidyverse)
library(lubridate)

# #extract files
#only do once
# library(archive)
# 
# # List contents
# archive::archive("raw_data/site_veg_data/Auscribe_JoinedTables.7z")
# 
# # Extract all
# archive::archive_extract("raw_data/site_veg_data/Auscribe_JoinedTables.7z", 
# dir = "raw_data/site_veg_data/Auscribe_JoinedTables_update")



#load Auscribe_JoinedTables df's

#contains the point intercept records from Auscribe_allTables in a more accessible format.
#Records have been joined to show their site, plot and herbarium ID
point_intercept<- read_csv("raw_data/site_veg_data/Auscribe_JoinedTables_update/point_intercept_records.csv")

#contains the site and location data for the visits contained in the other tables in this folder. 
site_location_visit_name<-read_csv("raw_data/site_veg_data/Auscribe_JoinedTables_update/site_location_visit_name.csv")

#contains the vegetation records records from Auscribe_allTables in a more accessible format
vegetation_records <- read_csv("raw_data/site_veg_data/Auscribe_JoinedTables_update/vegetation_records.csv")




#make an understandable df

auscribe_dat<-site_location_visit_name %>% 
  left_join(point_intercept) %>% 
  left_join(vegetation_records) %>% 
  relocate(Site,site_location_name,mer_plot_type) %>% 
  select(-c(Plot.size,method,field_name,hasBeenSighted,site_location_visit_id,
            site_location_id,id, visit_start_date)) %>% 
  rename(plot=site_location_name)

#######read in site data for sites that I collected data from
All_Sites<-readRDS("raw_data/MER_Site_data/All_Sites.RDS")

auscribe_dat %>% distinct(Site)


#redo labels for auscribe data so I can join by plot
auscribe_sites <- auscribe_dat %>%
  mutate(
    plot_join = plot %>%
      str_replace('WAMCOO0009', 'WAEW09') %>%
      str_replace('WAMCOO000', 'WASS0') %>%
      str_replace('WAMCOO00', 'WAEW') %>%
      str_replace('WAMESP00', 'WASP') %>%
      str_replace('SAMMDD0009', 'SAM209') %>%
      str_replace('SAMMDD000', 'SAM10') %>%
      str_replace('SAMMDD001', 'SAM21') %>%
      str_replace('QDMSEQ00', 'QDRF') %>%
      str_replace('VCMSEC000', 'VCEF0') %>%
      str_replace('VCMSEC00', 'VCRF') %>%
      str_replace('NSMSEQ00', 'NSS') %>% 
      str_replace('WAMDAL0009','WABVT09') %>% 
      str_replace("WAMDAL001([0-6])", "WABVT1\\1") %>% 
      str_replace("WAMDAL00", 'WABEW')
  ) %>% 
  filter(!Site %in% c('QDMBBSEucWoo01','QDMWETEucWoo01','NSMSECEucFor01','NSMNNCEucFor01', 'WAMJAFEucFor01')) %>% #sites not part of my study
  filter(!plot %in% paste0("NSMSEQ00", stringr::str_pad(9:12, width = 2, pad = "0"))) %>%  #removes extra plots at NSMSEQCasWoo01
  filter(!plot %in% paste0("VCMSEC00", stringr::str_pad(9:12, width = 2, pad = "0"))) %>%  #removes extra plots at VCMSECRainfo01
  select(Site,plot,plot_join) %>% distinct()%>% 
  group_by(Site) %>% 
  mutate(n=n()) %>% 
  left_join(All_Sites, by=c('Site','plot_join'='Plot')) %>% 
  ungroup()

auscribe_dat %>% distinct(Site)
auscribe_sites %>% distinct(Site)  

#join by plot names that exist in auscribe site df
auscribe_veg<-auscribe_dat %>% 
  inner_join(auscribe_sites %>% select(plot,plot_join)) %>% 
  filter(!is.na(substrate))#remove rows with all na values

missing_plots_y1<-anti_join(All_Sites,auscribe_veg, by=c('Plot'='plot_join'))

#quality check
anti_join(auscribe_veg,auscribe_sites, by=c('plot_join'='plot_join')) %>% 
  distinct(Site,plot)

saveRDS(auscribe_veg, 'Processed_data/metadata/veg/auscribe_veg_yr1.Rdata')


########fulcrum data#######
#LAST time since fire!!!

#load Fulcrum_data df's

#load plot level data
mer_pilot_plot_monitoring_module <- readr::read_csv("raw_data/site_veg_data/Fulcrum_data/mer_pilot_plot_monitoring_module.csv")%>%
  #extract the plot_name into the correct col
  mutate(
    plot_name_other = if_else(
      is.na(plot_name_other) & stringr::str_count(plot_name, ",") >= 2,
      sapply(strsplit(plot_name, ","), function(x) x[3]),
      plot_name_other)) %>% 
  #select the first row input
  mutate(date = ymd(date)) %>%  
  group_by(plot_name_other) %>%
  slice_min(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup()

t<-mer_pilot_plot_monitoring_module %>% select(plot_name_other,starts_with('percentage'))

#select cols I want and relevant plots
MER_plot_severity<-mer_pilot_plot_monitoring_module %>%
  select(plot_name_other,gps_altitude,estimated_fire_intensity_at_plot) %>% 
  mutate(plot = str_remove(plot_name_other, "-.*"))


#join and select focal sites 
plot_fire_severity<-MER_plot_severity %>% 
  inner_join(auscribe_sites %>% 
              select(Site,plot,plot_join,n,Fire_Treatment))

saveRDS(plot_fire_severity, 'Processed_data/metadata/fire/MER_plot_fire_severity.Rdata')



mer_veg_burn_plot<-mer_pilot_plot_monitoring_module %>%
  select(plot_name_other,percentage_of_all_individuals_that_are_dead_or_defoliated_e:percentage_of_individuals_with_intact_canopy_shrubs)%>% 
  mutate(plot = str_remove(plot_name_other, "-.*")) %>% 
  inner_join(auscribe_sites %>% 
               select(Site,plot,plot_join,n,Fire_Treatment))

# Count missing values per Site for each column
missing_summary <- mer_veg_burn_plot %>%
  group_by(Site) %>%
  summarise(across(
    percentage_of_all_individuals_that_are_dead_or_defoliated_e:
      percentage_of_individuals_with_intact_canopy_shrubs,
    ~ sum(is.na(.)),
    .names = "missing_{.col}"
  ))

# View the result
print(missing_summary)

write.csv(missing_summary,file = 'missing_summary_fulcrum.csv')


#load site level data
mer_pilot_site_description_module <- read_csv("raw_data/site_veg_data/Fulcrum_data/mer_pilot_site_description_module.csv")


#extract most recent burn data
MER_Site_fire_info<-mer_pilot_site_description_module %>%
  rename(Site=site_name_auscribe) %>% 
  filter(!Site %in% c('QDMBBSEucWoo01','QDMWETEucWoo01','NSMSECEucFor01','NSMNNCEucFor01', 'WAMJAFEucFor01')) %>% #sites not part of my study
  select(Site,year_of_most_recent_fire,month_of_most_recent_fire,
         has_there_been_another_significant_fire_at_the_site_in_the_last_10_years_that_should_be_recorded
         ) %>% 
  # Step 1: Update year and month for site VCMSECRainfo01
  mutate(
    year_of_most_recent_fire = if_else(Site == "VCMSECRainfo01", 2020L, year_of_most_recent_fire),
    month_of_most_recent_fire = if_else(Site == "VCMSECRainfo01", "Jan", month_of_most_recent_fire)
  )%>%
  # Step 2: Remove row where Site == WAMDALRainfo01 and year_of_most_recent_fire == 2020
  filter(!(Site == "WAMDALRainfo01" & year_of_most_recent_fire == 2020) # "Removed: confirmed by national fire history dataset"
         ) %>% 
  filter(!is.na(year_of_most_recent_fire))

saveRDS(MER_Site_fire_info, 'Processed_data/metadata/fire/MER_Site_fire_info.Rdata')




##############
#load year 2 data

#load Auscribe_JoinedTables df's

#contains the point intercept records from Auscribe_allTables in a more accessible format.
#Records have been joined to show their site, plot and herbarium ID
point_intercept_y2<- read_csv("raw_data/site_veg_data/MER_year_2_dat/pointIntercept_records_y2_cleaned.csv")

#contains the site and location data for the visits contained in the other tables in this folder.
site_location_visit_name_y2<-read_csv("raw_data/site_veg_data/MER_year_2_dat/site_location_visit_name_y2.csv")

#contains the vegetation records records from Auscribe_allTables in a more accessible format
vegetation_records_y2 <- read_csv("raw_data/site_veg_data/MER_year_2_dat/vegetation_records_y2.csv")



auscribe_dat_y2<-point_intercept_y2 %>%
  relocate(Site,name,mer_plot_type) %>%
  select(-c(in_canopy_branch,in_canopy_sky,transectId,plotId,plot_size,method))


#make an understandable df

auscribe_dat_y2<-point_intercept_y2%>% 
  left_join(vegetation_records_y2) %>% 
  relocate(Site,plotName,mer_plot_type) %>% 
  select(-c(method,field_name,
            id, visit_start_date)) %>% 
  rename(plot=plotName)


t<-auscribe_dat_y2%>% 
  filter(Site %in% c('WAMDALEucWoo01','WAMDALEucWoo02','WAMDALRainfo01'))


#redo labels for auscribe data so I can join by plot
auscribe_sites_y2 <- auscribe_dat_y2 %>%
  filter(!mer_plot_type=='BT') %>% 
  mutate(
    plot_join = plot %>%
      str_replace('WAMCOO0009', 'WAEW09') %>%
      str_replace('WAMCOO000', 'WASS0') %>%
      str_replace('WAMCOO00', 'WAEW') %>%
      str_replace('WAMESP00', 'WASP') %>%
      str_replace('SAMMDD0009', 'SAM209') %>%
      str_replace('SAMMDD000', 'SAM10') %>%
      str_replace('SAMMDD001', 'SAM21') %>%
      str_replace('QDMSEQ00', 'QDRF') %>%
      str_replace('VCMSEC000', 'VCEF0') %>%
      str_replace('VCMSEC00', 'VCRF') %>%
      str_replace('NSMSEQ00', 'NSS') %>% 
      str_replace('WAMDAL0009','WABVT09') %>% 
      str_replace("WAMDAL001([0-6])", "WABVT1\\1") %>% 
      str_replace("WAMDAL00", 'WABEW')
  ) %>% 
  filter(!Site %in% c('QDMBBSEucWoo01','QDMWETEucWoo01','NSMSECEucFor01','NSMNNCEucFor01', 'WAMJAFEucFor01')) %>% #sites not part of my study
  filter(!plot %in% paste0("NSMSEQ00", stringr::str_pad(9:12, width = 2, pad = "0"))) %>%  #removes extra plots at NSMSEQCasWoo01
  filter(!plot %in% paste0("VCMSEC00", stringr::str_pad(9:12, width = 2, pad = "0")) #removes extra plots at VCMSECRainfo01
         ) %>% 
  filter(!is.na(plot_join)) %>% 
  left_join(All_Sites, by=c('Site','plot_join'='Plot'))
 
# Count number of plots per site
plot_counts <- auscribe_sites_y2 %>%
  group_by(Site, mer_plot_type) %>%
  summarise(n_plots = n_distinct(plot), .groups = "drop")
  
auscribe_veg_y2<-auscribe_sites_y2 %>% 
  inner_join(missing_plots_y1)
  
temp<-auscribe_veg_y2 %>% distinct(Site,plot_join)

setdiff(names(auscribe_veg_y2), names(auscribe_veg))

colnames(auscribe_veg)
colnames(auscribe_veg_y2)
auscribe_veg_y2 %>% select(Site,herbarium_determination, herbariumId )

#labeling is stupid, this is what I did to fix it
auscribe_veg_y2 <- auscribe_veg_y2 %>%
  mutate(herbarium_determination = coalesce(herbarium_determination, herbariumId)) %>%
  select(-herbariumId) 

auscribe_veg_y2 <- auscribe_veg_y2 %>%
  select(any_of(names(auscribe_veg)))

library(lubridate)

# Combine datasets and label source
auscribe_veg_both <- bind_rows(
  auscribe_veg %>%
    select(-dead) %>%
    mutate(year = "first"),
  
  auscribe_veg_y2 %>%
    mutate(
      visit_end_date = as.POSIXct(dmy(visit_end_date), tz = "UTC"),
      year = "second"
    )
)

# Create a new column with labeled dates (Unknown if NA)
auscribe_veg_both <- auscribe_veg_both %>%
  mutate(visit_end_date_label = if_else(
    is.na(visit_end_date), 
    as.POSIXct(NA),  # Keep NA for proper histogram binning
    visit_end_date
  ))

# Count missing visit_end_date rows
unknown_count <- auscribe_veg_both %>%
  filter(is.na(visit_end_date))  %>%
  count(year, name = "n")

# Plot with an extra "Unknown" bar
ggplot(auscribe_veg_both %>% filter(!is.na(visit_end_date)), 
       aes(x = visit_end_date, fill = interaction(year))) +
  geom_histogram(
    binwidth = 7 * 24 * 60 * 60,
    color = 'black',
    alpha = 0.7
  ) +
  geom_bar(
    data = unknown_count,
    aes(x = as.POSIXct(max(auscribe_veg_both$visit_end_date, na.rm = TRUE) + 7*24*60*80),
        y = n,
        fill = year,
        color ='black'),
    stat = "identity",
    inherit.aes = FALSE,
    width = 7 * 24 * 60 * 60
  ) +
  scale_fill_manual(values = c("first" = "#1f78b4", "second" = "#33a02c")) +
  scale_x_datetime(date_breaks = "1 month", date_labels = "%b-%Y") +
  labs(
    title = "Number of Observations per Visit Date",
    x = "Month (Unknown grouped on far right)",
    y = "Number of Observations",
    fill = "Dataset"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

missing_plots_y1<-anti_join(All_Sites,auscribe_veg_both, by=c('Plot'='plot_join')) %>% distinct(Site,Plot)


saveRDS(auscribe_veg_both, 'Processed_data/metadata/veg/auscribe_veg_both.Rdata')


########fulcrum data y2#######
#LAST time since fire!!!

#load Fulcrum_data df's

#load plot level data
mer_pilot_plot_monitoring_module_2 <- readr::read_csv("raw_data/site_veg_data/fulcrum_y2_aug24/mer_pilot_plot_monitoring_module_aug24/mer_pilot_plot_monitoring_module.csv")%>%
  #extract the plot_name into the correct col
  mutate(
    plot_name_other = if_else(
      is.na(plot_name_other) & stringr::str_count(plot_name, ",") >= 2,
      sapply(strsplit(plot_name, ","), function(x) x[3]),
      plot_name_other)) %>% 
  #select the first row input
  mutate(date = ymd(date)) %>%  
  group_by(plot_name_other) %>%
  slice_min(order_by = date, n = 1, with_ties = FALSE) %>%
  ungroup()

t<-mer_pilot_plot_monitoring_module_2 %>% select(plot_name_other,starts_with('percentage'))

#select cols I want and relevant plots
MER_plot_severity_2<-mer_pilot_plot_monitoring_module_2 %>%
  select(plot_name_other,gps_altitude,estimated_fire_intensity_at_plot) %>% 
  mutate(plot = str_remove(plot_name_other, "-.*"))


#join and select focal sites 
plot_fire_severity_2<-MER_plot_severity_2 %>% 
  inner_join(auscribe_sites %>% 
               select(Site,plot,plot_join,n,Fire_Treatment))

#saveRDS(plot_fire_severity, 'Processed_data/metadata/fire/MER_plot_fire_severity.Rdata')



mer_veg_burn_plot_2<-mer_pilot_plot_monitoring_module_2 %>%
  select(plot_name_other,percentage_of_all_individuals_that_are_dead_or_defoliated_e:percentage_of_individuals_with_intact_canopy_shrubs)%>% 
  mutate(plot = str_remove(plot_name_other, "-.*")) %>% 
  inner_join(auscribe_sites %>% 
               select(Site,plot,plot_join,n,Fire_Treatment))

# Count missing values per Site for each column
missing_summary_2 <- mer_veg_burn_plot_2 %>%
  group_by(Site) %>%
  summarise(across(
    percentage_of_all_individuals_that_are_dead_or_defoliated_e:
      percentage_of_individuals_with_intact_canopy_shrubs,
    ~ sum(is.na(.)),
    .names = "missing_{.col}"
  ))

# View the result
print(missing_summary_2)

#write.csv(missing_summary,file = 'missing_summary_fulcrum_2.csv')
