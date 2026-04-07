library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(performance)
library(emmeans)

#load Beta1 and 2 from soil data HMSC
load('HMSC_MER/results/mcmc_output.RData')

#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#load most up to date PA model
load('HMSC_MER/results/Beta1.RData')

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean
Beta_estimates_COP<-Beta2$mean


Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)

Beta_estimates_COP<-Beta_estimates_COP %>% 
  dplyr::select(OTU=Species,Beta=Severity)

# Select OTUs associated with phylum

soil_OTU_phy <- myco_tax_soil %>%
  distinct() 

###############################################
###PA MODEL###############
#########################

#Join and merge with Beta_estimates
data<- soil_OTU_phy %>%
  inner_join(Beta_estimates_PA, by = "OTU") %>% 
  mutate(abs_Beta=abs(Beta))

library(lme4)

model_phylum <- lmer(Beta ~ 1 + (1|phylum), data = data)
model_order  <- lmer(Beta ~ 1 + (1|phylum/order), data = data)
model_family <- lmer(Beta ~ 1 + (1|phylum/order/family), data = data)
model_genus <- lmer(Beta ~ 1 + (1|phylum/order/family/genus), data = data)

#anova(model_phylum, model_order, model_family, model_genus)
VarCorr(model_genus)


data_nested <- data %>%
  filter(!is.na(phylum), !is.na(order), !is.na(family), !is.na(genus))

model_phylum <- lmer(Beta ~ 1 + (1|phylum), data = data_nested)
model_order  <- lmer(Beta ~ 1 + (1|phylum/order), data = data_nested)
model_family <- lmer(Beta ~ 1 + (1|phylum/order/family), data = data_nested)
model_genus  <- lmer(Beta ~ 1 + (1|phylum/order/family/genus), data = data_nested)

anova(model_phylum, model_order, model_family, model_genus)


summary(model_genus)
VarCorr(model_genus)
