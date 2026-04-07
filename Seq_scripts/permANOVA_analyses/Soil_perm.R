library(tidyverse)
library(vegan)
library(ggrepel)

# Load data
wide_myco_soil <- readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)

# Merge metadata
wide_myco_soil <- All_Sites %>%
  left_join(wide_myco_soil) %>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

site_list <- unique(wide_myco_soil$Site)

# Color palette for plotting
burn_colors <- c("B" = "darkred", "U" = "orange")

# Store results
permanova_summary <- list()
permanova_full <- list()
cap_scores_list <- list()
cap_var_explained <- list() 

# Loop over sites
for (site_name in site_list) {
  site_data <- wide_myco_soil %>% filter(Site == site_name)
  mat_myco <- site_data %>% select(starts_with("ITSall"))
  
  if (nrow(site_data) < 3 || length(unique(site_data$Fire_Treatment)) < 2) {
    permanova_summary[[site_name]] <- tibble(Site = site_name, R2 = NA_real_, p_value = NA_real_, N = nrow(site_data))
    next
  }
  
  # PERMANOVA
  perm <- adonis2(mat_myco ~ Fire_Treatment, data = site_data, distance = 'robust.aitchison', by = 'margin')
  fire_row <- perm["Fire_Treatment", ]
  
  permanova_summary[[site_name]] <- tibble(
    Site = site_name,
    R2 = fire_row$R2,
    p_value = fire_row$`Pr(>F)`,
    N = nrow(site_data)
  )
  
  permanova_full[[site_name]] <- perm %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(Site = site_name, .before = 1)
  
  # CAP analysis
  cap_mod <- capscale(mat_myco ~ Fire_Treatment, data = site_data, distance = 'robust.aitchison', add = TRUE)
  
  
  # Variance explained by CAP1
  proportions<-round(cap_mod$CCA$eig/cap_mod$tot.chi *100, 1) # proportion of variation associated with each axis

  cap_var_explained[[site_name]] <- tibble(Site = site_name, var_CAP1 = proportions)
  
  #Site scores
  scrs <- scores(cap_mod, tidy = TRUE) %>%
    filter(score == "sites") %>%
    bind_cols(site_data %>% select(Site, Plot, Fire_Treatment)) %>%
    mutate(Site = site_name)
  
  cap_scores_list[[site_name]] <- scrs
}

# Combine results
permanova_summary_df <- bind_rows(permanova_summary)
permanova_full_df <- bind_rows(permanova_full)
cap_var_df <- bind_rows(cap_var_explained)
cap_all_sites_df <- bind_rows(cap_scores_list) %>% 
  left_join(cap_var_df)%>%
  mutate(label = paste0( " CAP1 (",round(var_CAP1, 1), "%)")) 



# View summary stats
print(permanova_summary_df)

saveRDS(permanova_summary_df,'Processed_data/Seq_dat/permANOVA_data/soil_perm.Rdata')


facet_labels <- cap_all_sites_df %>%
  filter(Site %in% veg_long$Site) %>% 
  distinct(Site, var_CAP1) %>%
  mutate(Site = as.character(Site),
         label = paste0(Site, " (CAP1 ", round(var_CAP1, 1), "%)")) %>% 
  select(Site,label) %>% 
  deframe()


# Plot facetted CAP by Site
cap_all_sites_df %>% 
  filter(Site %in% veg_long$Site) %>% 
ggplot( aes(x = CAP1, y = MDS1, colour = Fire_Treatment)) +
  geom_point(size = 3) +
  facet_wrap(~ Site, scales = "free",
             labeller= as_labeller(facet_labels),
             strip.position = "bottom") +
  scale_colour_manual(values = burn_colors) +
  theme_bw(base_size = 14) +
  labs(
         x = "CAP1",
         y = "MDS1", color = "Burn") +
  theme(strip.text = element_text(size = 14))->cap_plot

print(cap_plot)

# ggsave("plots/facetted_CAP_by_site.png", cap_plot, width = 15, height = 10, dpi = 300)


