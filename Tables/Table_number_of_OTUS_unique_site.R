library(tidyverse)

dat_myco_RA_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')




# Count how many unique sites each OTU is present in
otu_site_counts <- dat_myco_RA_soil %>%
  filter(resampled_count > 0) %>%  # keep only OTUs actually detected
  distinct(Site, OTU) %>%         # ensure unique combinations
  count(OTU, name = "n_sites")    # count how many sites each OTU occurs in

head(otu_site_counts)

otu_shared_summary <- otu_site_counts %>%
  count(n_sites, name = "n_OTUs") %>%   # count number of OTUs per site-occurrence level
  arrange(desc(n_sites))                # start from 14 sites down to 1

otu_shared_summary

#number of unique OTUs at each site
otus_per_site <- dat_myco_RA_soil %>%
  filter(resampled_count > 0) %>%
  distinct(Site, OTU) %>%
  count(Site, name = "n_OTUs") %>%
  arrange(desc(n_OTUs))

otus_per_site


list(
  OTUs_shared_across_sites = otu_shared_summary,
  OTUs_per_site = otus_per_site
)
