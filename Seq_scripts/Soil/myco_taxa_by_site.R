library(tidyverse)

library(ggplot2)
library(scales)

dat_myco_RA_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#### create site labels
# site_labels_named <- c(
#   "NSMSEQCasWoo01"  = "NSW-Forest-Casuarina",
#   "QDMSEQRainfo01"  = "QLD-Rainfo-Subtrop",
#   "SAMMDDMallee01"  = "SA-Mallee-Spinifex",
#   "SAMMDDMallee02"  = "SA-Mallee-Chenopod",
#   "VCMSECEucFor01"  = "VIC-Forest-Euc",
#   "VCMSECRainfo01"  = "VIC-Rainfo-Temperate",
#   "WAMDALEucWoo01"  = "WA-Savanna-Corymbia",
#   "WAMDALRainfo01"  = "WA-Vine-Thicket",
#   "WAMDALEucWoo02"  = "WA-Savanna-Minyjuru",
#   "WAMCOOShrubl01"  = "WA-Sandplain-AlloCas",
#   "WAMCOOEucWoo01"  = "WA-Woodland-SalmGum",
#   "WAMCOOEucWoo02"  = "WA-Woodland-Oleosa",
#   "WAMCOOEucWoo03"  = "WA-Woodland-Gimlet",
#   "WAMESPShrubl01"  = "WA-Shrub-Kingia"
# )

site_labels_named <- c(
  "NSMSEQCasWoo01" = "Forest-Melaleuca/Casuarina",
  "QDMSEQRainfo01" = "Subtropical-Rainforest",
  "SAMMDDMallee01" = "Triodia-Mallee-Eucalyptus-dumosa",
  "SAMMDDMallee02" = "Chenopod-Mallee-Eucalyptus-oleosa",
  "VCMSECEucFor01" = "Forest-Eucalyptus-muelleriana",
  "VCMSECRainfo01" = "Temperate-Rainforest",
  "WAMDALEucWoo01" = "Savannah-Corymbia-paractia",
  "WAMDALRainfo01" = "VineThicket",
  "WAMDALEucWoo02" = "Savannah-Sersalisia-sericea",
  "WAMCOOShrubl01" = "Shrubland-Allocasuarina-acutivalvis",
  "WAMCOOEucWoo01" = "Woodland-Eucalyptus-salmonophloia",
  "WAMCOOEucWoo02" = "Chenopod-Mallee-Eucalyptus-affinis-oleosa",
  "WAMCOOEucWoo03" = "Woodland-Eucalyptus-salubris",
  "WAMESPShrubl01" = "Shrubland-Kingia australis"
)
# CUSTOM THEME FOR READABILITY
# =========================
theme_pub <- theme_classic(base_size = 18) +  # Bigger base text size
  theme(
    axis.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold"),
    plot.tag = element_text(size = 20, face = "bold")
  )


# -----------------------------
# Filter for Glomeromycota
# -----------------------------
dat_glom <- dat_myco_RA_soil %>%
  filter(phylum == "Glomeromycota")

dat_glom %>%
  distinct(OTU, family) %>%            # Ensure unique OTU-family pairs
  count(family, name = "n_OTUs") %>%   # Count OTUs per family
  mutate(total_OTUs = sum(n_OTUs),
         percentage = n_OTUs / total_OTUs * 100) %>%
  arrange(desc(percentage))


# -----------------------------
# PLOT 1: All Glomeromycota OTUs
# -----------------------------

# Create a named color vector for all families
families <- dat_glom %>% 
  filter(!is.na(family)) %>% 
  pull(family) %>% 
  factor() %>% 
  levels()

family_colors <- setNames(hue_pal()(length(families)), families)

family_colors

#create summary
dat_glom_summary <- dat_glom %>% 
  group_by(Site, Plot, family) %>% 
  summarise(RA_fam_samp = sum(count) / unique(reads_samp), .groups = "drop") %>% 
  filter(!is.na(family)) %>% 
  mutate(
    Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
    Site = fct_relevel(Site, rev(levels(Site)))
  )

n_sites <- nlevels(dat_glom_summary$Site)

p_glom_all <- dat_glom %>% 
  filter(!is.na(family)) %>% 
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, rev(levels(Site)))  # reorder alphabetically
) %>% 
  ggplot( aes(x = RA_samp, 
              y = Site, # set order & labels
              color = family)) +
  geom_point(alpha = 0.7, position = position_dodge(width = .8), size=4) +
  geom_hline(yintercept = seq(1.5, n_sites - 0.5, 1),
             color = "grey70", linetype = "dashed", linewidth = 0.4) +
  labs( tag='a)', x = "Relative Abundance (% log10)", y = "") +
  scale_x_log10()+
  #scale_size(guide = 'none')+
  scale_color_manual(values = family_colors) +  # <-- apply colors here
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 5)))
p_glom_all

#by family per plot
# p_glom_all <- dat_glom_summary %>%
#   ggplot(aes(x = RA_fam_samp, y = Site, color = family)) +
#   geom_point(alpha = 0.7, position = position_dodge(width = .8), size = 4) +
#   # horizontal separators:
#   geom_hline(yintercept = seq(1.5, n_sites - 0.5, 1),
#              color = "grey70", linetype = "dashed", linewidth = 0.4) +
#   labs(tag = "a)", x = "Relative Abundance (% log10)", y = "") +
#   scale_x_log10() +
#   scale_color_manual(values = family_colors) +
#   theme_pub +
#   guides(color = guide_legend(override.aes = list(size = 5)))
# p_glom_all

# -----------------------------
# PLOT 2: Glomeromycota by Family 
### -----------------------------
# Filter for summary stats: keep only groups with >=3 observations
dat_glom_summary <- dat_glom %>%
  filter(!is.na(family)) %>%
  group_by(family, Fire_Treatment) %>%
  filter(n() >= 2) %>%
  ungroup()

n_families <- dat_glom %>% filter(!is.na(family)) %>% distinct(family) %>%  pull(family) %>%   factor() %>%nlevels()

p_glom_family_soil <- dat_glom %>%
  filter(!is.na(family)) %>%
  mutate(family = fct_reorder(family, RA_samp, .fun = sum, .desc = FALSE)) %>%  # Reorder genera
  ggplot(aes(x = RA_samp, y = family, color = Fire_Treatment)) +
  geom_point(aes(color=Fire_Treatment),alpha = 0.7,size=4, position = position_dodge(width = 0.8)) +
  stat_summary( data=dat_glom_summary, aes(x = RA_samp, y = family, group = Fire_Treatment), # <-- don't map color here!
               fun = mean,
               fun.min = function(z) mean(z) - sd(z)/sqrt(length(z)),
               fun.max = function(z) mean(z) + sd(z)/sqrt(length(z)),
               geom = "pointrange",
               size = 0.5,
               position = position_dodge(width = 0.8),
               color='black') +
  geom_hline(yintercept = seq(1.5, n_families - 0.5, 1),
             color = "grey70", linetype = "dashed", linewidth = 0.4) +
  scale_x_log10() +
  scale_color_manual(
    values = c("B" = "#E64B35", "U" = "#4DBBD5"),  # <-- Choose your colors
    labels = c("B" = "Burnt", "U" = "Unburnt")  ) +
  labs( tag= 'a)',x = "Relative Abundance (%log10)", y = "", color= 'Fire Treatment') +
  theme_pub 


p_glom_family_soil
# -----------------------------
# Repeat for EctoMycorrhizal
# -----------------------------
dat_ecto <- dat_myco_RA_soil %>%
  filter(guild2 == "EctoMycorrhizal")

# Create a named color vector for all families
genera <- dat_ecto %>% 
  filter(!is.na(genus)) %>% 
  pull(genus) %>% 
  factor() %>% 
  levels()

genus_colors <- setNames(hue_pal()(length(genera)), genera)

genus_colors

n_Sites_ecto <- dat_ecto %>% filter(!is.na(Site)) %>% distinct(Site) %>%  pull(Site) %>%   factor() %>%nlevels()


p_ecto_all <- dat_ecto %>% 
  filter(!is.na(genus)) %>% 
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, rev(levels(Site)))  # reorder alphabetically
  ) %>% 
  ggplot(aes(x = RA_samp, y = Site, color = genus)) +
  geom_point(alpha = 0.7, position = position_dodge(width = .8), size=4) +
  scale_size(guide = 'none')+
  labs( tag= 'b)',x = "Relative Abundance (% log10)", y = "") +
  scale_x_log10()+
  geom_hline(yintercept = seq(1.5, n_Sites_ecto - 0.5, 1),
             color = "grey70", linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(values = genus_colors)+
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 5)))
p_ecto_all



# Filter for summary stats: keep only groups with >=3 observations
dat_ecto_summary <- dat_ecto %>%
  filter(!is.na(genus)) %>%
  group_by(genus, Fire_Treatment) %>%
  filter(n() >= 3) %>%
  ungroup()

n_genera_ecto<- dat_ecto %>% filter(!is.na(genus)) %>% distinct(genus) %>%  pull(genus) %>%   factor() %>%nlevels()


p_ecto_genus_soil <- dat_ecto %>%
  filter(!is.na(genus)) %>%
  mutate(genus = fct_reorder(genus, RA_samp, .fun = sum, .desc = FALSE)) %>%  # Reorder genera
  ggplot(aes(x = RA_samp, y = genus, color = Fire_Treatment)) +
  geom_point(aes(color=Fire_Treatment),alpha = 0.7, size=4, position = position_dodge(width = 0.8)) +
  stat_summary(data = dat_ecto_summary,
               aes(x = RA_samp, y = genus, group = Fire_Treatment), # <-- don't map color here!
               fun = mean,
               fun.min = function(z) mean(z) - sd(z)/sqrt(length(z)),
               fun.max = function(z) mean(z) + sd(z)/sqrt(length(z)),
               geom = "pointrange",
               size = 0.5,
               position = position_dodge(width = 0.8),
               color='black') +
  scale_x_log10()+
  geom_hline(yintercept = seq(1.5, n_genera_ecto - 0.5, 1),
             color = "grey70", linetype = "dashed", linewidth = 0.4) +
  scale_color_manual(
    values = c("B" = "#E64B35", "U" = "#4DBBD5"),  # <-- Choose your colors
    labels = c("B" = "Burnt", "U" = "Unburnt")
  ) +
  labs( tag= 'b)',x = "Relative Abundance (% log10)", y = "", color= 'Fire Treatment') +
  theme_pub 

  p_ecto_genus_soil
# -----------------------------
# Show Plots
# -----------------------------
p_glom_all
p_glom_family_soil
p_ecto_all
p_ecto_genus_soil
####################################
library(patchwork)
Comm_by_Site<-p_glom_all/p_ecto_all

Comm_by_Site
Fire_response_by_phylo<-p_glom_family_soil/p_ecto_genus_soil
Fire_response_by_phylo

ggsave(filename='Plots/Fire_response_by_phylo.png', Fire_response_by_phylo, width = 15, height = 20, dpi = 300)
ggsave(filename = 'Plots/Comm_by_Site.png', Comm_by_Site, width = 15, height = 20, dpi = 300)

