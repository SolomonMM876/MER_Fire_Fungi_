library(tidyverse)

library(ggplot2)

dat_myco_RA_hyph<-readRDS('Processed_data/Seq_dat/hyph/myco_RA_hyph.rds')

#source('Seq_scripts/soil/myco_tax_by_site.R')
# -----------------------------
# Filter for Glomeromycota
# -----------------------------
dat_glom <- dat_myco_RA_hyph %>%
  filter(phylum == "Glomeromycota")

library(scales)  # For label_number()
# -----------------------------
# PLOT 1: All Glomeromycota OTUs
# -----------------------------

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


p_glom_all_hyph <- dat_glom %>% 
  filter(!is.na(family)) %>% 
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, rev(levels(Site)))  # reorder alphabetically
         ) %>% 
  ggplot( aes(x = RA_samp, y = Site, size = RA_samp, color = family)) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.5)) +
  labs(tag= 'a)', x = "Relative Abundance (% log10)", y = "Site") +
  scale_x_log10(  labels = label_number(accuracy = 0.001, scale_cut = cut_short_scale()) # <-- 3 sig figs, no sci notation
  )+
  scale_size(guide = 'none')+
  scale_color_manual(values = family_colors) +  # <-- apply colors here
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme_pub+
  theme(legend.position = 'none')
p_glom_all_hyph

# -----------------------------
# PLOT 2: Glomeromycota by Family 
### -----------------------------
dat_glom_summary <- dat_glom %>%
  filter(!is.na(family)) %>%
  group_by(family, Fire_Treatment) %>%
  filter(n() >= 4) %>%
  ungroup()



p_glom_family_hyph<- dat_glom %>%
  filter(!is.na(family)) %>%
  mutate(family = fct_reorder(family, RA_samp, .fun = sum, .desc = FALSE)) %>%  # Reorder genera
  ggplot(aes(x = RA_samp, y = family, color = Fire_Treatment)) +
  geom_point(aes(color=Fire_Treatment),alpha = 0.6,size=3, position = position_dodge(width = 0.5)) +
  stat_summary(data=dat_glom_summary,
               aes(x = RA_samp, y = family, group = Fire_Treatment), # <-- don't map color here!
                 fun = mean,
                 fun.min = function(z) mean(z) - sd(z)/sqrt(length(z)),
                 fun.max = function(z) mean(z) + sd(z)/sqrt(length(z)),
                 geom = "pointrange",
                 size = 0.5,
                 position = position_dodge(width = 0.5),
                 color='black') +
  scale_x_log10(  labels = label_number(accuracy = 0.001, scale_cut = cut_short_scale()) # <-- 3 sig figs, no sci notation
) +
  labs( tag= 'a)',x = "Relative Abundance (%log10)", y = "AM family", color= 'Fire Treatment') +
  theme_pub +
  theme(legend.position = 'none')
  

p_glom_family_hyph
# -----------------------------
# Repeat for EctoMycorrhizal
# -----------------------------
dat_ecto <- dat_myco_RA_hyph %>%
  filter(guild2 == "Ectomycorrhizal")

p_ecto_all_hyph <- dat_ecto %>% 
  filter(!is.na(genus)) %>% 
  mutate(Site = factor(Site, levels = names(site_labels_named), labels = site_labels_named),
         Site = fct_relevel(Site, rev(levels(Site)))  # reorder alphabetically
  ) %>% 
  ggplot(aes(x = RA_samp, y = Site, size = RA_samp, color = genus)) +
  geom_point(alpha = 0.6, position = position_dodge(width = 0.2)) +
  scale_size(guide = 'none')+
  labs( tag= 'c)',x = "Relative Abundance (% log10)", y = "Site") +
  scale_x_log10()+
  scale_color_manual(values = genus_colors)+
  theme_pub+
  theme(legend.position = 'none')+
  guides(color = guide_legend(override.aes = list(size = 5)))
p_ecto_all_hyph



# Filter for summary stats: keep only groups with >=3 observations
dat_ecto_summary <- dat_ecto %>%
  filter(!is.na(genus)) %>%
  group_by(genus, Fire_Treatment) %>%
  filter(n() >= 3) %>%
  ungroup()


p_ecto_genus_hyph <- dat_ecto %>%
  filter(!is.na(genus)) %>%
  mutate(genus = fct_reorder(genus, RA_samp, .fun = sum, .desc = FALSE)) %>%  # Reorder genera
  ggplot(aes(x = RA_samp, y = genus, color = Fire_Treatment)) +
  geom_point(aes(color=Fire_Treatment),alpha = 0.6, size=3, position = position_dodge(width = 0.8)) +
  stat_summary(data = dat_ecto_summary,
               aes(x = RA_samp, y = genus, group = Fire_Treatment), # <-- don't map color here!
               fun = mean,
               fun.min = function(z) mean(z) - sd(z)/sqrt(length(z)),
               fun.max = function(z) mean(z) + sd(z)/sqrt(length(z)),
               geom = "pointrange",
               size = 0.5,
               position = position_dodge(width = 0.5),
               color='black') +
  scale_x_log10()+
  labs( tag= 'c)', x = "Relative Abundance (% log10)", y = "Ectomycorrhizal genus", color= 'Fire Treatment') +
  theme_pub+
  theme(legend.position = 'none')

p_ecto_genus_hyph
# -----------------------------
# Show Plots
# -----------------------------
p_glom_all_hyph
p_glom_family_hyph
p_ecto_all_hyph
p_ecto_genus_hyph
####################################


library(patchwork)

glom_all<-p_glom_all_hyph+p_glom_all
glom_all
ecto_all<-p_ecto_all_hyph + p_ecto_all
ecto_all

myco<-glom_all/ecto_all

myco

ggsave(filename = 'Plots/myco_comm.png', myco, width = 20, height = 11, dpi = 300)


glom_family<-p_glom_family_hyph+p_glom_family_soil

ecto_genus<-p_ecto_genus_hyph+p_ecto_genus_soil

phylo_fire<-glom_family/ecto_genus

phylo_fire
ggsave(filename = 'Plots/phylo_fire.png', phylo_fire, width = 18, height = 15, dpi = 300)
