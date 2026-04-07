library(tidyverse)

library(ggplot2)

dat_myco_RA_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')



# Identify the top 10 most abundant genera
top_genera_soil <- dat_myco_RA_soil %>%
  filter(!is.na(genus)) %>% 
  group_by(genus) %>%
  summarise(total_abundance = sum(RA_total_Burn, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(genus)
# Compute ordered levels for genus_grouped
genus_order <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus)) %>%
  group_by(genus_grouped) %>%
  summarise(total_abundance = sum(RA_total_Burn, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%  # ascending so the most abundant ends up on top
  pull(genus_grouped)

genus_order <- c(setdiff(genus_order, "Other"), "Other")

# Make genus_grouped a factor with the correct order
dat_myco_RA_soil <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus),
         genus_grouped = factor(genus_grouped, levels = genus_order))

#Make plot

Fire_Genus_myco<-dat_myco_RA_soil %>% 
  ggplot(aes(x=Fire_Treatment, y=RA_total_Burn, fill=genus_grouped, text=species)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  geom_bar(stat = 'identity', position = position_stack(), width = 0.6) +
  scale_fill_viridis_d(option = "H", name = "Genus",direction = -1) +
  theme_classic() +
  labs(
    y = 'Porportion mycorrhizal sequences',
    x = '',
    tag = "(b)"
  ) +
  theme_classic()+
  theme(
    legend.position = 'right',
    axis.text.x = element_text(hjust = 0.5, size = 25,lineheight = 0.7),
    axis.text.y = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 30),
    strip.text = element_text(size = 20),
    plot.tag = element_text(size = 22, hjust = 10),
  )

Fire_Genus_myco
#########################################################################
#FACET by site#############
#calculate relative abundance of reads
dat_myco_RA_soil<-dat_myco_RA_soil%>%
  group_by(Site,Fire_Treatment)%>%
  summarise(burn_reads_site=sum(count))%>%
  ungroup() %>% 
  left_join(dat_myco_RA_soil) %>% 
  mutate( RA_Site_Burn= count/burn_reads_site)


# Identify the top 10 most abundant genera
top_genera_soil <- dat_myco_RA_soil %>%
  filter(!is.na(genus)) %>% 
  group_by(Site,genus) %>%
  summarise(total_abundance = sum(RA_Site_Burn, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = 10) %>%
  pull(genus)
# Compute ordered levels for genus_grouped
genus_order <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus)) %>%
  group_by(Site,genus_grouped) %>%
  summarise(total_abundance = sum(RA_Site_Burn, na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%  # ascending so the most abundant ends up on top
  pull(genus_grouped)

genus_order <- c(setdiff(genus_order, "Other"), "Other")

# Make genus_grouped a factor with the correct order
dat_myco_RA_soil <- dat_myco_RA_soil %>%
  mutate(genus_grouped = if_else(is.na(genus) | !(genus %in% top_genera_soil), "Other", genus),
         genus_grouped = factor(genus_grouped, levels = genus_order))

#Make plot

Fire_Genus_myco_Site<-dat_myco_RA_soil %>% 
  ggplot(aes(x=Fire_Treatment, y=RA_Site_Burn, fill=genus_grouped, text=species)) + # text aesthetic is for the ggplotly visualisation below
  facet_wrap(~Site)+
  geom_bar(stat = 'identity' ,position = position_stack(), width = 0.4) +
  geom_bar(stat = 'identity', position = position_stack(), width = 0.6) +
  scale_fill_viridis_d(option = "H", name = "Genus",direction = -1) +
  theme_classic() +
  labs(
    y = 'Porportion mycorrhizal sequences',
    x = '',
    tag = "(b)"
  ) +
  theme_classic()+
  theme(
    legend.position = 'right',
    axis.text.x = element_text(hjust = 0.5, size = 25,lineheight = 0.7),
    axis.text.y = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 30),
    strip.text = element_text(size = 20),
    plot.tag = element_text(size = 22, hjust = 10),
  )

Fire_Genus_myco_Site

##############################################################################
#Phylum by site#############

library(readxl)

shipping<-read_excel("C:/Users/90957135/OneDrive - Western Sydney University/MER Fire/MER Fire and Fungi Involvement Sheet.xlsx") %>% 
  select(Site=2,ship_time=13) %>% filter(!is.na(Site))

#Make plot

Fire_phylum_myco_Site<-dat_myco_RA_soil %>% left_join(shipping) %>% 
  ggplot(aes(x=Fire_Treatment, y=RA_Site_Burn, fill=phylum, text=species)) + # text aesthetic is for the ggplotly visualisation below
  facet_wrap(~Site)+
  geom_bar(stat = 'identity', position = position_stack(), width = 0.6) +
  scale_fill_viridis_d(option = "H", name = "phylum",direction = -1) +
  geom_text(
    aes(x = Fire_Treatment, y = 1 + 0.02, label = burn_reads_site),
    inherit.aes = FALSE,
    size = 6
  ) +
  geom_text(aes(x=1.5,y = .5, label = ship_time),
              inherit.aes = FALSE,
              size = 6
  ) +
  theme_classic() +
  labs(
    y = 'Porportion mycorrhizal sequences',
    x = '',
    tag = "(b)"
  ) +
  theme_classic()+
  theme(
    legend.position = 'right',
    axis.text.x = element_text(hjust = 0.5, size = 25,lineheight = 0.7),
    axis.text.y = element_text(hjust = 0.5, size = 25),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 30),
    strip.text = element_text(size = 20),
    plot.tag = element_text(size = 22, hjust = 10),
  )

Fire_phylum_myco_Site

