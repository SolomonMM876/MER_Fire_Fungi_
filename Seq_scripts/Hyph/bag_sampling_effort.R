library(tidyverse)
library(vegan)


#import Hyph seq data
MER_Hyph<-readRDS('Processed_data/Seq_dat/MER_Hyph.rds' )
#redo plot names
MER_Hyph<-MER_Hyph %>% 
  mutate(  Plot = str_replace(Plot, "(U|B)$", ""),#Remove Unburnt|Burnt
           Plot = str_replace(Plot, "^NSS0", "NSS"),#Remove extra 0 added by mistake
           Plot= str_replace(Plot, 'WABUT','WABVT'),#fix spelling error
           Plot = str_replace(Plot, "VCFR", "VCRF"),#Fix mixup of letters
           Plot = str_replace(Plot, "WAE3", "WAEW3"), #missed a W add that back in for two samples
           Plot= str_replace(Plot, 'WABEW15','WABEW18')) %>% # no wabew15
  #remove empty cols
  dplyr::select(-sample,qbit_DNA_conc_ng_uL)

#Check to make sure I have all my samples
#first all site data
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") 


# 1. Create presence/absence matrix (rows = Plot/bags, columns = OTUs)
mat_pa <- All_Sites %>%
  distinct(Site, Plot) %>%
  left_join(MER_Hyph) %>%
  mutate(presence = 1) %>%
  distinct(Site, Plot, OTU, presence) %>%
  pivot_wider(names_from = OTU, values_from = presence, values_fill = 0)

# 2. Loop over each Site to generate species accumulation curves using `specaccum()`
sac_by_site <- mat_pa %>%
  group_split(Site) %>%
  map_dfr(function(df_site) {
    site_name <- unique(df_site$Site)

    
    site_matrix <- df_site %>%
      column_to_rownames("Plot") %>%
      dplyr::select(-Site)%>%
      dplyr::select(where(is.numeric))  # Ensure only OTUs are included
    
    sac <- specaccum(site_matrix, method = "random", permutations = 999)
    
    tibble(
      Site = site_name,
      Samples = sac$sites,
      Species = sac$richness,
      SD = sac$sd
    )
  })

# 3. Plot accumulation curves by Site
ggplot(sac_by_site, aes(x = Samples, y = Species, color = Site)) +
  geom_line(size = 1) +
  #geom_ribbon(aes(ymin = Species - SD, ymax = Species + SD, fill = Site), alpha = 0.2, color = NA) +
  labs(x = "Number of Bags", y = "Species Richness") +
  theme_classic(base_size = 16)
