#LOAD LIBRARARIES
library(forcats)
library(vegan)
library(tidyverse) 
library(readxl)
library(ggplot2)



#From ITS prep script
wide_myco_soil<-readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')
myco_tax_Soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_Soil.rds')
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site,Plot,Fire_Treatment)


###next analysis#######

# next analysis - permanova

# first remove samples that have no mycorrhizal reads because they cause errors below
wide_myco_soil<-All_Sites%>%
  left_join(wide_myco_soil) %>% 
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

# extract the community table, save as a new object
mat_myco<- wide_myco_soil %>% select(starts_with("ITSall"))

adonis.Site<-adonis2(mat_myco~ Site , data=wide_myco_soil, distance='robust.aitchison', by= 'margin')


temp<-adonis2(mat_myco~ Fire_Treatment +Site, data=wide_myco_soil, distance='robust.aitchison', by= 'margin')

temp<-temp%>%
  as.data.frame()%>%
  rownames_to_column(var='Factor') %>% 
  mutate(Sample_Type='Soil')
temp
#write_RDS(temp,'Tables/Soil_permanova.csv', row.names = FALSE)
table(wide_myco_soil$Fire_Treatment)

cap.site <- capscale(mat_myco~ Site , data=wide_myco_soil, distance='robust.aitchison', add=TRUE)

cap.all <- capscale(mat_myco~ Fire_Treatment +Site, data=wide_myco_soil, distance='robust.aitchison', add=TRUE)
anova(cap.all)
summary(cap.all)
Cap1_aov<-as.data.frame( anova(cap.all, by = "margin"))%>%
  rownames_to_column()
cap.all
plot(cap.all)
proportions<-round(cap.all$CCA$eig/cap.all$tot.chi *100, 1) # proportion of variation associated with each axis




# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap.all, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site<- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')
scrs_biplot <- scrs %>% filter(score=='biplot')

burn_colors <- c("B" = "darkred", "U" = "orange")

library(ggrepel)

# first plot - site scores along with centroids for each group
cbind(wide_myco_soil, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2 )) + 
 # facet_wrap(~Site, scales='free_x')+
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  geom_point(aes( colour= Fire_Treatment), size=10, stroke = 3)+ 
  #stat_ellipse(aes(color = Fire_Treatment), level = 0.95, size = 2, linetype = 1) + # Add confidence ellipses
  #geom_text(aes( label = Site), color= 'black', size=4)+
  #scale_shape_manual(values = c(19,1))+
  scale_colour_manual(values = burn_colors) +     # Custom colors for Interval
  labs( x=  paste0("CAP1 (", proportions[1], "%)"), y= paste0("CAP2 (", proportions[2], "%)"),
        colour = "Burn",  # Rename legend for color (Fire.Interval)
  )+
  xlim(c(min(scrs_site[, 'CAP1']-1), max(scrs_site[, 'CAP1'])+.2)) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2'])+1)) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 36, color = 'black'),
    axis.text.y = element_text(size = 36, color = 'black'),
    axis.title.y = element_text(size = 36, color = 'black'),
    axis.title.x = element_text(size = 36),
    legend.position = "top",
    legend.text = element_text(size = 30),        # Larger legend text
    legend.title = element_text(size = 30)        # Optional: larger legend title
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 8)),  # Adjust shape size
    colour = guide_legend(override.aes = list(size = 8))  # Adjust color legend size
  )  ->p2

p2

#ggsave(filename = "plots/30_Soil_PCoA.png", plot = p2, dpi=300, device = "png", width = 65, height = 45, units = "cm")

