# ==============================================================================
# SCRIPT 00: MASTER DATA PREPARATION
# DESCRIPTION: Reads raw project files and creates a self-contained 'data' folder
#              with all CSVs required for HMSC fitting and Beta analysis.
# ==============================================================================

library(tidyverse)
library(readxl)

# 1. SETUP DIRECTORY STRUCTURE
# ------------------------------------------------------------------------------
main_dir <- "project_script"
data_dir <- file.path(main_dir, "data")

dir.create(main_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

print(paste("Created data directory at:", data_dir))


# 2. EXPORT INPUTS FOR HMSC MODEL FITTING (For Script 01)
# ------------------------------------------------------------------------------

# --- A. Site Metadata ---
# Input: All_Sites.RDS
# Output: input_site_metadata.csv
print("Processing Site Metadata...")
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS")

site_meta_clean <- All_Sites %>% 
  dplyr::select(Site, Plot, Lat, Long, Fire_Treatment)

write_csv(site_meta_clean, file.path(data_dir, "input_site_metadata.csv"))

# --- B. Host Frequency Covariates ---
# Input: myco_freq_plot_level_df.Rdata
# Output: input_host_frequency.csv
print("Processing Host Frequency...")
# Note: readRDS is safer if the file is truly RDS. If it's an environment RData, we load it.
# Based on your previous code using readRDS for this file:
myco_freq <- readRDS('Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata') %>% 
  dplyr::select(Plot, freq_AM, freq_ECM)

write_csv(myco_freq, file.path(data_dir, "input_host_frequency.csv"))

# --- C. Soil Community Matrix (Raw Counts) ---
# Input: wide_myco.rds
# Output: input_soil_community.csv
print("Processing Soil Community Matrix...")
soil_comm <- readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')

# Ensure Plot column exists and is first
# (Assuming the raw file has Plot as a column based on your join logic)
write_csv(soil_comm, file.path(data_dir, "input_soil_community.csv"))


# 3. EXPORT DATA FOR BETA ANALYSIS (For Script 02)
# ------------------------------------------------------------------------------

# --- D. Beta Estimates ---
# Input: Beta1_host_freq.RData (HMSC output)
# Output: beta_estimates.csv
print("Processing Beta Estimates...")
load('HMSC_MER/results/Beta1_host_freq.RData') # Loads 'Beta1' object

beta_clean <- Beta1$mean %>%
  dplyr::select(OTU = Species, Beta = Fire_TreatmentB)

write_csv(beta_clean, file.path(data_dir, "beta_estimates.csv"))

# --- E. Detailed Site Metadata (Analysis) ---
# Inputs: Aridity, Fire Regime, Site Description Excel
# Output: site_metadata_analysis.csv
print("Processing Analysis Metadata...")
MER_aridity <- readRDS('Processed_data/metadata/Aridity_MER.Rdata')
TSF <- readRDS('Processed_data/metadata/fire/Fire_regime_MER.Rdata') %>% 
  distinct(Site, severity, TSF_years)
MER_Desc <- read_excel("Processed_data/metadata/veg/MER_Site_description.xlsx", 
                       sheet = "SP suggestions (highlighted)")

site_metadata_analysis <- MER_Desc %>%
  left_join(TSF, by = "Site") %>%
  dplyr::select(Site, Site_nickname, Last_fire_severity, severity, TSF_years,
                Veg_type, Veg_type_alt, Veg_type_alt_2, Fire_interval, 
                Min_TFI = `Min TFI`, Simple_fire_response_dominants, Fire_response) %>%
  mutate(Min_TFI = as.numeric(Min_TFI)) %>%
  left_join(MER_aridity, by = "Site")

write_csv(site_metadata_analysis, file.path(data_dir, "site_metadata_analysis.csv"))

# --- F. Taxonomy (Combined) ---
# Inputs: Soil and Hyph Taxonomy RDS
# Output: taxonomy.csv
print("Processing Taxonomy...")
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_tax_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')
taxonomy_combined <- bind_rows(myco_tax_soil, myco_tax_hyph) %>% distinct()

write_csv(taxonomy_combined, file.path(data_dir, "taxonomy.csv"))

# --- G. Abundance Data (For RA Calculations) ---
# Inputs: RA RDS files
# Output: soil_abundance.csv, hyph_abundance.csv
print("Processing Abundance Data...")
myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
write_csv(myco_dat_soil, file.path(data_dir, "soil_abundance.csv"))

myco_dat_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
write_csv(myco_dat_hyph, file.path(data_dir, "hyph_abundance.csv"))

# --- H. Functional Traits ---
# Inputs: Reproductive strategy RDS
# Output: traits_repo_strat.csv, traits_amf_guilds.csv
print("Processing Traits...")
repo_strat <- readRDS('Processed_data/Seq_dat/datasets_external/repro_strat_df.Rds')
write_csv(repo_strat, file.path(data_dir, "traits_repo_strat.csv"))

amf_traits <- tribble(
  ~Guild,         ~Intraradical_hyphae, ~Extraradical_hyphae, ~Families,              ~Citations,
  "Rhizophilic",  "High",               "Low",                "Glomeraceae",          "1,2,3,4,5",
  "Rhizophilic",  "High",               "Low",                "Claroideoglomeraceae", "2",
  "Rhizophilic",  "High",               "Low",                "Paraglomeraceae",      "b",
  "Edaphophilic", "Low",                "High",               "Gigasporaceae",        "1,2,5",
  "Edaphophilic", "Low",                "High",               "Diversisporaceae",     "2,5",
  "Ancestral",    "Low",                "Low",                "Archaeosporaceae",     "b",
  "Ancestral",    "Low",                "Low",                "Ambisporaceae",        "b",
  "Ancestral",    "Low",                "Low",                "Pacisporaceae",        "5",
  "Ancestral",    "Low",                "Low",                "Acaulosporaceae",      "1,2"
)
write_csv(amf_traits, file.path(data_dir, "traits_amf_guilds.csv"))

print("SUCCESS: All data files exported to 'project_script/data/'")
print("You may now run Script 01 or Script 02.")