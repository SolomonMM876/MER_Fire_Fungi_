library(tidyverse)
library(vegan)
library(ggrepel)

# Load data
wide_myco_hyph<-readRDS('Processed_data/Seq_dat/Hyph/wide_myco.rds')
myco_tax_Hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)

# Merge metadata
wide_myco_hyph <- All_Sites %>%
  left_join(wide_myco_hyph) %>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

site_list <- unique(wide_myco_hyph$Site)

# Color palette for plotting
burn_colors <- c("B" = "darkred", "U" = "orange")

# Store results
permanova_summary <- list()
permanova_full <- list()
cap_scores_list <- list()

# Loop over sites
for (site_name in site_list) {
  site_data <- wide_myco_hyph %>% filter(Site == site_name)
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
  
  # CAP Analysis: catch overfitted models
  cap_mod <- try(capscale(mat_myco ~ Fire_Treatment, data = site_data, distance = 'robust.aitchison', add = TRUE), silent = TRUE)
  
  if (inherits(cap_mod, "try-error")) {
    warning(paste("CAP failed for site:", site_name))
    next
  }
  
  # Check if residual variance exists
  if (cap_mod$CA$rank == 0) {
    warning(paste("CAP overfitted (no residual variance) for site:", site_name))
    next
  }
  
  scrs <- try(scores(cap_mod, tidy = TRUE), silent = TRUE)
  
  if (!inherits(scrs, "try-error")) {
    scrs_site <- scrs %>%
      filter(score == "sites") %>%
      bind_cols(site_data %>% select(Site, Plot, Fire_Treatment)) %>%
      mutate(Site = site_name)
    
    cap_scores_list[[site_name]] <- scrs_site
  }
}


# Combine results
permanova_summary_df <- bind_rows(permanova_summary)
permanova_full_df <- bind_rows(permanova_full)
cap_all_sites_df <- bind_rows(cap_scores_list)

# View summary stats
print(permanova_summary_df)

saveRDS(permanova_summary_df,'Processed_data/Seq_dat/permANOVA_data/hyph_perm.Rdata')


# Plot facetted CAP by Site
cap_plot <- ggplot(cap_all_sites_df, aes(x = CAP1, y = MDS1, colour = Fire_Treatment)) +
  geom_point(size = 3) +
  facet_wrap(~ Site, scales = "free") +
  scale_colour_manual(values = burn_colors) +
  theme_bw(base_size = 14) +
  labs(x = "CAP1", y = "MDS1", color = "Burn") +
  theme(strip.text = element_text(size = 14))

print(cap_plot)

# ggsave("plots/facetted_CAP_by_site.png", cap_plot, width = 15, height = 10, dpi = 300)


