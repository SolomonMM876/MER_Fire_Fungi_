library(tidyverse)

dat_myco_RA_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

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
  "WAMCOOEucWoo01"  = "WA-Wood-SalmGum",
  "WAMCOOEucWoo02"  = "WA-Wood-Oleosa",
  "WAMCOOEucWoo03"  = "WA-Wood-Gimlet",
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
otu_site_counts <- dat_myco_RA_soil %>%
  filter(count > 0) %>% 
  distinct(Site, OTU, guild2) %>%      # keep guild2 for later grouping
  count(OTU, guild2, name = "n_sites") # count unique sites per OTU per guild

# Summarize how many OTUs are shared across sites, grouped by guild2
otu_shared_summary <- otu_site_counts %>%
  count(n_sites, guild2, name = "n_OTUs") %>%
  arrange(desc(n_sites)) %>% 
  mutate(label_x = case_when(
    n_sites == 4 & guild2 == "EctoMycorrhizal" ~ n_sites - 0.1,
    n_sites == 4 & guild2 == "Ericoid Mycorrhizal" ~ n_sites + 0.1,
    n_sites == 6 & guild2 == "EctoMycorrhizal" ~ n_sites - 0.2,
    n_sites == 6 & guild2 == "Arbsucular Mycorrhizal" ~ n_sites + 0.3,
    TRUE ~ as.numeric(n_sites)),
    label_y = case_when(
      n_sites == 6 & guild2 == "EctoMycorrhizal" ~ n_sites - 5,
      n_sites == 6 & guild2 == "Arbsucular Mycorrhizal" ~ n_sites + 5,
      TRUE ~ as.numeric(n_sites)
  ))

# If you also want the overall percent (ignoring guild2):
otu_shared_overall <- otu_site_counts %>%
  count(n_sites, name = "n_OTUs") %>%
  mutate(
    total_OTUs = sum(n_OTUs),
    percent = 100 * n_OTUs / total_OTUs
  )

# Plot
overlap_OTU_sites_soil<-otu_shared_summary %>% 
  mutate(guild2= str_replace(guild2,'EctoMycorrhizal','Ectomycorrhizal')) %>% 
ggplot(aes(x = n_sites, y = n_OTUs, fill = guild2)) +
  geom_col(position = "stack", alpha = 0.9) +
  geom_text(aes(x = label_x, label = n_OTUs), 
            position = position_stack(vjust = 0.5), size = 6, color = "black") +
  geom_text(
    data = otu_shared_overall,
    aes(x = n_sites, y = n_OTUs, label = paste0(round(percent, 1), "%")),
    inherit.aes = FALSE,  # avoid conflicting aesthetics
    vjust = -0.5,         # nudge above bar
    size = 8,
    color = "black"
  ) +
  scale_x_continuous(breaks = seq(1, max(otu_shared_summary$n_sites), 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
   # tag='a)',
    x = "Number of Sites OTU Occurs In",
    y = "Number of OTUs",
    fill = "Guild"
  ) +
  theme_pub+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

overlap_OTU_sites_soil
ggsave(filename = 'Plots/OTUs_overlap_site_soil.png', overlap_OTU_sites_soil, width = 15, height = 11, dpi = 300)

#repeat above, but dont group by guild############

# Count how many unique sites each OTU is present in, keeping guild2
otu_site_counts <- dat_myco_RA_soil %>%
  filter(count > 0) %>% 
  distinct(Site, OTU) %>%      # keep guild2 for later grouping
  count(OTU, name = "n_sites") # count unique sites per OTU per guild

# Summarize how many OTUs are shared across sites, grouped by guild2
otu_shared_summary <- otu_site_counts %>%
  count(n_sites, name = "n_OTUs")

# If you also want the overall percent (ignoring guild2):
otu_shared_overall <- otu_site_counts %>%
  count(n_sites, name = "n_OTUs") %>%
  mutate(
    total_OTUs = sum(n_OTUs),
    percent = 100 * n_OTUs / total_OTUs
  )

# Plot
overlap_OTU_sites_soil_no_guild<-otu_shared_summary %>% 
  ggplot(aes(x = n_sites, y = n_OTUs)) +
  geom_col(position = "stack", alpha = 0.9, fill='black') +
  geom_text(
    data = otu_shared_overall,
    aes(x = n_sites, y = n_OTUs, label = paste0(round(percent, 1), "%")),
    inherit.aes = FALSE,  # avoid conflicting aesthetics
    vjust = -0.5,         # nudge above bar
    size = 8,
    color = "black"
  ) +
  scale_x_continuous(breaks = seq(1, max(otu_shared_summary$n_sites), 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    # tag='a)',
    x = "Number of Sites OTU Occurs In",
    y = "Number of OTUs") +
  theme_pub+
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "right"
  )

overlap_OTU_sites_soil_no_guild



























# Count number of unique OTUs per Site and guild2
otus_per_site_guild <- dat_myco_RA_soil %>%
  filter(count > 0) %>%
  distinct(Site, OTU, guild2) %>%  # keep guild info per OTU per Site
  count(Site, guild2, name = "n_OTUs") %>%
  arrange(desc(n_OTUs))

# Count number of unique plots per site
plots_per_site <- dat_myco_RA_soil %>%
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
OTUs_site_soil<-otus_per_site_guild %>% 
mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
       Site = fct_relevel(Site, sort(levels(Site)))  # reorder alphabetically
) %>% 
ggplot( aes(x = Site, y = n_OTUs, fill = guild2)) +
  geom_col(alpha = 0.9) +
  geom_text(aes(label = n_OTUs),
            position = position_stack(vjust = 0.5),
            size = 4, color = "black") +
  # Add number of plots above each bar
  geom_text(data = label_positions,
            aes(x = Site, y = max_OTUs + 3, label = paste0("n=", n_plots)),
            inherit.aes = FALSE,
            size = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    tag = 'b)',
    x = "Site",
    y = "Number of OTUs",
    fill = "Guild"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,),
    legend.position = "right"
  )
library(patchwork)

#overlap_OTU_sites<-overlap_OTU_sites_hyph+overlap_OTU_sites_soil
overlap_OTU_sites_soil
Soil_OTUS<-overlap_OTU_sites_soil+OTUs_site_soil

Soil_OTUS
ggsave(filename = 'Plots/OTU_overlap_soil.png', overlap_OTU_sites_soil, width = 15, height = 11, dpi = 300)

#OTUs_site_soil<-OTUs_site_hyph+OTUs_site_soil

OTUs_site_soil
ggsave(filename = 'Plots/OTUs_site_soil.png', OTUs_site_soil, width = 15, height = 11, dpi = 300)

