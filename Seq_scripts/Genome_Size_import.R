
#LOAD LIBRARARIES
library(tidyverse)


#mycorrhizal taxa
#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
#read in hyph comm data
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')


myco_tax<-bind_rows(myco_tax_hyph,myco_tax_soil) %>% 
  distinct(OTU, .keep_all = TRUE)

gs_dat<-read.csv('Processed_data/Seq_dat/datasets_external/gs_dataset_mycorrhizae_Hiyang_3.3.25.csv')

gs_dat_genus<-gs_dat%>%filter(!is.na(genus))%>%
  group_by(genus)%>%summarise(mean_gs=mean(GS))

tax_gs<-myco_tax%>%
  left_join(gs_dat_genus)%>%filter(!is.na(mean_gs))


saveRDS(tax_gs, 'Processed_data/Seq_dat/datasets_external/tax_genome_size_MER.Rds')
