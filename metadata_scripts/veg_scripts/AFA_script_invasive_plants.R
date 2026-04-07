library(tidyverse)
library(APCalign)


#read in national dataset from MArtin-Flores et al. 2023
AFA_national<-read_csv('Processed_data/metadata/veg/datasets_published/AFA_Martin_Fores_etal_2023/AFA_national_dataset.csv')


#read in Auscribe year 1 data already extract from veg_import_clean
auscribe_veg<-readRDS('Processed_data/metadata/veg/auscribe_veg_yr1.Rdata')


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


##############################################################################

upd_spp_list <- update_taxonomy(
  aligned_spp,
  taxonomic_splits = "most_likely_species",
  resources = tax_resources
) %>%
  select(original_name, suggested_name, aligned_name, genus:taxon_rank) %>%
  filter(!is.na(suggested_name))  # Remove indeterminate names

# Join updated taxonomy back into auscribe_veg
veg_tax <- auscribe_veg %>%
  left_join(upd_spp_list, by = c("herbarium_determination" = "original_name")) %>% 
# Create a state abbreviation column
  mutate(State_code = str_sub(Site, 1, 3))

# Define a function to subset by state and join unified status
get_state_df <- function(state_abbrev, state_colname) {
  veg_tax %>%
    filter(State_code == state_abbrev) %>%
    left_join(
      AFA_national %>%
        select(APC_canonical_name, unified_status = all_of(state_colname)),
      by = c("suggested_name" = "APC_canonical_name")
    )
}

# Apply function to each state
NSW <- get_state_df("NSM", "NSW_unified_status")
QLD <- get_state_df("QDM", "QLD_unified_status")
SAM <- get_state_df("SAM", "SA_unified_status")
VCM <- get_state_df("VCM", "VIC_unified_status")
WAM <- get_state_df("WAM", "WA_unified_status")

#combine sites back together
veg_combined <- bind_rows(NSW, QLD, SAM, VCM, WAM)


temp<-veg_combined %>% 
  filter(is.na(unified_status)) %>% 
#exclude plants that were unable to be identified or only identified to genus or family status
  filter(taxon_rank=='species',!is.na(herbarium_determination)) %>% 
  distinct(aligned_name,Site,State_code)

unique(veg_combined$unified_status)


# Step 1: Rename and filter valid unified status
veg_plot_summary <- veg_combined %>%
  filter(!is.na(unified_status)) %>%
  rename(Plot = plot_join) 

veg_species<-veg_plot_summary%>% 
# Step 2: Get total species per plot
  group_by(Site, Plot) %>%
  summarise(
    total_species = sum(!is.na(aligned_name)),
    .groups = "drop"
  )

# Step 3: Count frequency of each unified_status per plot
status_freq_per_plot <- veg_plot_summary %>%
  group_by(Site, Plot, unified_status) %>%
  summarise(
    count = n(),
    .groups = "drop"
  )

# Step 4: Join total species to compute relative frequency
status_freq_summary <- status_freq_per_plot %>%
  left_join(veg_species, by = c("Site", "Plot")) %>%
  mutate(
    rel_freq = count / total_species  # relative frequency per plot
  )

saveRDS(status_freq_summary, 'processed_data/metadata/veg/AFA_native_invasive_freq_MER.Rds')

All_Sites<-readRDS("raw_data/MER_Site_data/All_Sites.RDS")