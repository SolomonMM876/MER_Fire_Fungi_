library(tidyverse)
library(readxl)

# read in all site locations
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% dplyr::select(Site, Plot, Fire_Treatment)

# read in fire history from Fire history_import
previous_fires <- readRDS("Processed_data/metadata/fire/previous_fires.Rdata")

# read in fire freq for Broome plots from Fire history_import from 2000-2024
# site_fire_history_df<-readRDS('Processed_data/metadata/fire/Broom_site_fire_history_df.Rdata')

# read in fulcrum last fire data at site level from veg import_clean
MER_Site_fire_info <- readRDS("Processed_data/metadata/fire/MER_Site_fire_info.Rdata")

# read in plot level fire severity data from veg import_clean
plot_fire_severity <- readRDS("Processed_data/metadata/fire/MER_plot_fire_severity.Rdata")

# Read in Suzzanne Prober  AusEcoModel spreadsheet
AusEcoModel <- read_excel("Processed_data/metadata/veg/MER eco_type fire details 3Jul25_SP.xlsx")

MER_fire_regime <- All_Sites %>%
  left_join(plot_fire_severity %>% dplyr::select(Site, Plot = plot_join, severity = estimated_fire_intensity_at_plot)) %>%
  left_join(MER_Site_fire_info %>% dplyr::select(-has_there_been_another_significant_fire_at_the_site_in_the_last_10_years_that_should_be_recorded)) %>%
  left_join(previous_fires %>% dplyr::select(Site, next_previous_fire_year, oldest_fire_year)) %>%
  mutate(
    severity = if_else(Fire_Treatment == "U" & is.na(severity), "Unburnt", severity),
    note = case_when(
      Site == "WAMDALRainfo01" ~ "potential fire in 2024 via fire dataset",
      Site == "SAMMDDMallee02" ~ "potential sig fire in last 10 years fulcrum",
      Site == "WAMDALEucWoo02" ~ "potential sig fire in last 10 years fulcrum",
      Site == "VCMSECRainfo01" ~ "data from Forest Information Portal (https://maps.ffm.vic.gov.au/fip/index.html?viewer=fip)",
      Site == "VCMSECEucFor01" ~ "data from Forest Information Portal (https://maps.ffm.vic.gov.au/fip/index.html?viewer=fip) and digital atlas of Australia (https://digital.atlas.gov.au/pages/national-bushfire-data)",
      Site == "WAMCOOShrubl01" ~ "Allocasuarina is an obligate-seeder, i.e. killed by fire. Fire frequency in this location should be >40 years",
      Site %in% c("WAMCOOEucWoo01", "WAMCOOEucWoo03") ~ "See eucalypt woodland booklet, infrequent (>300 year) stand-replacing (intense) fires. Large mature trees in this site likely unburnt for 250 or more yearsprior to the MER fire",
      Site == "WAMCOOEucWoo02" ~ "See eucalypt woodland and mallee booklets, based on large size of burnt stems Id say these sites were at least 80-100 years since prior fire"
    ),
    # confirmed via email with site managers
    severity = case_when(
      Plot %in% c("VCEF07", "VCEF08") ~ "moderate",
      Plot %in% c("WASP05") ~ "very high",
      TRUE ~ severity
    ),
    next_previous_fire_year = case_when(
      Site == "VCMSECRainfo01" ~ 1992,
      Site == "VCMSECEucFor01" ~ 1981,
      TRUE ~ next_previous_fire_year
    ),
    oldest_fire_year = case_when(
      Site == "VCMSECRainfo01" ~ 1981,
      TRUE ~ oldest_fire_year
    )
  ) %>%
  mutate(
    severity_num = case_when(
      str_detect(severity, regex("extreme", ignore_case = TRUE)) ~ 5,
      str_detect(severity, regex("very high", ignore_case = TRUE)) ~ 4,
      str_detect(severity, regex("high", ignore_case = TRUE)) ~ 3,
      str_detect(severity, regex("moderate", ignore_case = TRUE)) ~ 2,
      str_detect(severity, regex("low", ignore_case = TRUE)) ~ 1,
      TRUE ~ NA_real_ # Optional: for cases that don't match anything
    )
  ) %>%
  left_join(AusEcoModel %>% dplyr::select(Site, fire_interval))

MER_fire_regime



MER_regime <- MER_fire_regime %>%
  filter(Fire_Treatment == "B") %>%
  group_by(Site) %>%
  summarise(
    severity = mean(severity_num),
    year_of_most_recent_fire,
    fire_interval,
    .groups = "drop"
  ) %>%
  distinct()
MER_regime

# read in collection data
All_biomass <- readRDS("processed_data/Bag_data/MER_biomass.RDS") %>% distinct(Site, Harvest_date)

TSF <- MER_Site_fire_info %>%
  left_join(All_biomass) %>%
  # Convert month names to numbers and create a proper fire date
  mutate(
    fire_month_num = match(month_of_most_recent_fire, month.abb), # Convert Jan -> 1, Feb -> 2, etc.
    fire_date = make_date(year_of_most_recent_fire, fire_month_num, 1),
    # Calculate difference in months, rounding to nearest whole month
    TSF_months = round(interval(fire_date, Harvest_date) / months(1)),
    TSF_years = round(TSF_months / 12, 1)
  ) %>%
  select(-fire_month_num, -fire_date) # Clean up helper columns

MER_regime <- MER_regime %>%
  left_join(TSF %>% select(Site, TSF_years))

saveRDS(MER_regime, "Processed_data/metadata/fire/Fire_regime_MER.rds")

# Plot-level fire severity: unburnt plots assigned fire_severity_num = 0
# severity_num already computed above via str_detect:
#   1 = Low, 2 = Moderate, 3 = High, 4 = Very High, 5 = Extreme
# Here we add 0 for Unburnt so all plots sit on the same continuous scale.
MER_regime_plot <- MER_fire_regime %>%
  mutate(
    fire_severity_num = if_else(Fire_Treatment == "U", 0, severity_num)
  )

saveRDS(MER_regime_plot, "Processed_data/metadata/fire/Fire_regime_plot.rds")
