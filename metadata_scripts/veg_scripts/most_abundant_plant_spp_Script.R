library(tidyverse)

#read in combined auscribe datasets extract from veg_import_clean
auscribe_veg<-readRDS('Processed_data/metadata/veg/auscribe_veg_both.Rdata')


#library(APCalign)
#load for austraits
# tax_resources <- load_taxonomic_resources()
# 
# # Extract distinct species names as a character vector
# spp_names <- auscribe_veg %>%
#   distinct(herbarium_determination) %>%
#   filter(!is.na(herbarium_determination)) %>% 
#   pull(herbarium_determination) %>%  # extract as a character vector
#   iconv(from = "", to = "UTF-8", sub = "byte")  # Convert to valid UTF-8
#
# Align the species names using APC resources
# aligned_spp <- align_taxa(spp_names, resources = tax_resources)
# 
# upd_spp_list<-update_taxonomy(aligned_spp, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
#   select(original_name,suggested_name,aligned_name,genus:taxon_rank)%>%
#   filter(!is.na(suggested_name))#remove indeterminate IDs
# 
# All_Sites<-readRDS("raw_data/MER_Site_data/All_Sites.RDS")
# 
# auscribe_veg
# # 1. Join taxonomy (Run this once)
# auscribe_combined <- auscribe_veg %>%
#   left_join(upd_spp_list, by = c("herbarium_determination" = "original_name")) %>%
#   mutate(final_species = coalesce(suggested_name, herbarium_determination))

# 2. Filter, Calculate Site Totals, and Rank
top_per_site <- auscribe_veg %>%
  # Remove unwanted growth forms and NAs
  filter(!growth_form %in% c("Vine", "Fern", "Tussock grass","Forb", "Sedge","Shrub"),
         !is.na(growth_form), !is.na(herbarium_determination)) %>%
  # Groups by species to count how many unique Sites they appear in
  group_by(Site, herbarium_determination) %>%
  mutate(n_plots_present = n_distinct(plot)) %>% 
  ungroup() %>% 
# Calculate TOTAL hits per Site (Denominator)
add_count(Site, name = "site_total_hits") %>%
  # Count frequency of each species within that Site and Growth Form
  # CRITICAL: Include n_plots_present in count() to preserve the column
  count(Site,mer_plot_type, growth_form, herbarium_determination, site_total_hits, n_plots_present, name = "spp_hits") %>%
  # Calculate percentage (Species Hits / Total Site Hits)
  mutate(percent_site_cover = (spp_hits / site_total_hits) * 100) %>% 
  #  Re-group by Site and Plot Type before slicing
  group_by(Site, mer_plot_type) %>%
  slice_max(order_by = percent_site_cover, n = 2, with_ties = TRUE)

# View results
head(top_per_site)

# Filter for specific sites as requested
temp <- top_per_site %>% 
  filter(Site %in% c('NSMSEQCasWoo01'))

site_spp_dist <- auscribe_veg %>%
  filter(Site == 'NSMSEQCasWoo01',
         herbarium_determination %in% c('Casuarina glauca', 'Melaleuca quinquenervia')) %>%
  # Count hits for each species in each plot
  count(plot,mer_plot_type, herbarium_determination, name = "hits_per_plot") %>%
  # Reshape for easier comparison
  pivot_wider(names_from = herbarium_determination, 
              values_from = hits_per_plot, 
              values_fill = 0)

print(site_spp_dist)


top_trees_wide <- auscribe_veg %>%
  # 1. Filter: Keep only trees (by excluding others) and valid data
  filter(!growth_form %in% c("Vine", "Fern", "Tussock grass", "Forb", "Sedge", "Shrub"),
         !is.na(growth_form), 
         !is.na(herbarium_determination)) %>%
  
  # 2. Summarise counts per Species within Site AND Plot Type
  group_by(Site, mer_plot_type, herbarium_determination) %>%
  summarise(
    spp_hits = n(), 
    n_plots_present = n_distinct(plot), 
    .groups = "drop"
  ) %>%
  
  # 3. Calculate Totals and Percentages within the Site/Plot Type group
  group_by(Site, mer_plot_type) %>%
  mutate(
    total_hits = sum(spp_hits),
    percent_cover = (spp_hits / total_hits) * 100
  ) %>%
  
  # 4. Filter for presence in >3 plots
  #filter(n_plots_present > 2) %>%
  
  # 5. Rank: Select top 2 and assign a rank number (1 or 2)
  slice_max(order_by = percent_cover, n = 2, with_ties = FALSE) %>%
  mutate(rank = row_number()) %>% 
  ungroup() %>%
  
  # 6. Reshape: Select relevant columns and pivot to wide format
  select(Site, mer_plot_type, herbarium_determination, rank) %>%
  pivot_wider(
    names_from = rank,
    values_from = herbarium_determination,
    names_prefix = "rank_"
  ) %>%
  
  # 7. Rename columns to final requirement
  rename(
    most_abundant_tree = rank_1, 
    second_most_abundant_tree = rank_2
  )

# View Result
print(top_trees_wide, n=30)
