library(tidyverse)

dat_myco_RA_hyph<-readRDS('Processed_data/Seq_dat/hyph/myco_RA_hyph.rds')

#### create site labels
site_labels_named <- c(
  "NSMSEQCasWoo01"  = "NSW-Forest-Flood",
  "QDMSEQRainfo01"  = "QLD-Rainfo-Subtrop",
  "SAMMDDMallee01"  = "SA-Mallee-Spinifex",
  "SAMMDDMallee02"  = "SA-Mallee-Chenopod",
  "VCMSECEucFor01"  = "VIC-Forest-Euc",
  "VCMSECRainfo01"  = "VIC-Rainfo-Temperate",
  "WAMDALEucWoo01"  = "WA-Savanna-Corymbia",
  "WAMDALRainfo01"  = "WA-Vine-Thicket",
  "WAMDALEucWoo02"  = "WA-Savanna-Minyjuru",
  "WAMCOOShrubl01"  = "WA-Sandplain-AlloCas",
  "WAMCOOEucWoo01"  = "WA-Woodland-SalmGum",
  "WAMCOOEucWoo02"  = "WA-Woodland-Oleosa",
  "WAMCOOEucWoo03"  = "WA-Woodland-Gimlet",
  "WAMESPShrubl01"  = "WA-Shrub-Kingia"
)

theme_pub <- theme_classic(base_size = 18) +  # Bigger base text size
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold")
  )

# Count how many unique sites each OTU is present in, keeping guild2
otu_site_counts <- dat_myco_RA_hyph %>%
  filter(resampled_count > 0) %>% 
  distinct(Site, OTU, guild2) %>%      # keep guild2 for later grouping
  count(OTU, guild2, name = "n_sites") # count unique sites per OTU per guild

# Summarize how many OTUs are shared across sites, grouped by guild2
otu_shared_summary <- otu_site_counts %>%
  count(n_sites, guild2, name = "n_OTUs") %>%
  arrange(desc(n_sites))

# Plot
overlap_OTU_sites_hyph<-otu_shared_summary %>% 
ggplot(aes(x = n_sites, y = n_OTUs, fill = guild2)) +
  geom_col(position = "stack", alpha = 0.9) +
  geom_text(aes(label = n_OTUs), 
            position = position_stack(vjust = 0.5), size = 5, color = "white") +
  scale_x_continuous(breaks = seq(1, max(otu_shared_summary$n_sites), 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    tag='a)',
    x = "Number of Sites OTU Occurs In",
    y = "Number of OTUs",
    fill = "Guild"
  ) +
  theme_pub+ 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )






# Count number of unique OTUs per Site and guild2
otus_per_site_guild <- dat_myco_RA_hyph %>%
  filter(resampled_count > 0) %>%
  distinct(Site, OTU, guild2) %>%  # keep guild info per OTU per Site
  count(Site, guild2, name = "n_OTUs") %>%
  arrange(desc(n_OTUs))

# Count number of unique plots per site
plots_per_site <- dat_myco_RA_hyph %>%
  distinct(Site, Plot) %>%  # <- assuming column 'plot' exists in your data
  count(Site, name = "n_plots")

# Join plots per site to total OTUs for positioning labels
label_positions <- otus_per_site_guild %>%
  group_by(Site) %>%
  summarise(max_OTUs = sum(n_OTUs), .groups = "drop") %>%
  left_join(plots_per_site, by = "Site") %>% 
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, sort(levels(Site)))  # reorder alphabetically
  )

# Create stacked bar plot of OTU richness per site by guild
OTUs_site_hyph<-otus_per_site_guild %>% 
mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
       Site = fct_relevel(Site, sort(levels(Site)))  # reorder alphabetically
) %>% 
ggplot( aes(x = Site, y = n_OTUs, fill = guild2)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = n_OTUs),
            position = position_stack(vjust = 0.5),
            size = 3, color = "white") +
  # Add number of plots above each bar
  geom_text(data = label_positions,
            aes(x = Site, y = max_OTUs + 3, label = paste0("n=", n_plots)),
            inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    tag = 'a)',
    x = "Site",
    y = "Number of OTUs",
    fill = "Guild"
  ) +
  theme_pub+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )
OTUs_site_hyph
