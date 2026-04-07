library(tidyverse)
library(readxl)

#results of CN spiking
CN_nonspike<- read_excel("raw_data/Stoich/CN_Faizun/Nonspike_Solomon_Kazi-RESULT.xlsx", sheet = "MER_Nonspike")

CN_nonspike<-CN_nonspike %>% 
  rename(perc_N_Samp = `N%`,
         perc_C_Samp = `C%`) %>% 
  select(-c(Spike,Plate_ID, `Sample no`))

#spiked samples
CN_spike<- read_excel("raw_data/Stoich/CN_Faizun/Spike_Solomon_Kazi-RESULT.xlsx", sheet = "MER_Spike",skip = 1)

#drifts
CN_drift<- read_excel("raw_data/Stoich/CN_Faizun/drifts-FN-CN.xlsx") 

#separate standards
OL_samps<-CN_drift %>% 
  filter(name %in% c( 'OL 1 mg spiked', 'OL'))

#seperate N and C spikes
CN_drift<-CN_drift %>% 
  filter(!name %in% c( 'OL 1 mg spiked', 'OL')) %>% 
  mutate(N_mg_spike= WT*`N%`/100,
         C_mg_spike= WT* `C%`/100)
#Removes OL standards and computes the actual mg amounts of N and C in the spike materials, using weight × (% / 100).

#calc percent of N and C in spiked samples
CO3_mean_N_perc<-CN_drift%>%filter(name=='NaHCO3')%>% summarise(mean(`N%`) )%>%pull()
CO3_mean_C_perc<-CN_drift%>%filter(name=='NaHCO3')%>% summarise(mean(`C%`) )%>%pull()
NH4_mean_N_perc<-CN_drift%>%filter(name=='(NH)4SO4')%>% summarise(mean(`N%`) )%>%pull()
NH4_mean_C_perc<-CN_drift%>%filter(name=='(NH)4SO4')%>% summarise(mean(`C%`) )%>%pull()



#calc Total spike for C and N for OL
CN_standards<-OL_samps%>%
  filter(name %in% c( 'OL 1 mg spiked')) %>% 
  #Calculating total N and C (mg) in the spiked sample.
  mutate(samp_N_mg= WT*`N%`/100,
         samp_C_mg= WT* `C%`/100) %>% 
  mutate(Carb_spike_t_mg= (NaHCO3*CO3_mean_C_perc/100)+(`(NH)4SO4`*NH4_mean_C_perc/100),#total carbon spike
         Nitr_spike_t_mg= (NaHCO3*CO3_mean_N_perc/100)+(`(NH)4SO4`*NH4_mean_N_perc/100))%>%#total nitrogen spike
  #Subtracting out the known spike amounts.
  mutate(Samp_N_mg= samp_N_mg-Nitr_spike_t_mg,
         Samp_C_mg= samp_C_mg-Carb_spike_t_mg,
         #Dividing by the sample weight to get %N and %C after correcting for spike.
         perc_N_Samp= (Samp_N_mg/sample)*100,
         perc_C_Samp= (Samp_C_mg/sample)*100) 

Standards<-CN_standards %>% 
  select(name,WT, sample,perc_N_Samp,perc_C_Samp) %>% 
  bind_rows(OL_samps %>%   filter(name %in% c( 'OL')) %>% 
              select(name,WT,perc_N_Samp= `N%`,perc_C_Samp= `C%`))

#calc Total spike for C and N
CN_All<-CN_spike%>%
  mutate(samp_N_mg= Total*(`N%`/100),
         samp_C_mg= Total*(`C%`/100)) %>% 
  select(-c(Plate_ID,Spike,`...9`)) %>% 
  mutate(Carb_spike_t_mg= (C_Spike*CO3_mean_C_perc/100)+(N_Spike*NH4_mean_C_perc/100),#total carbon spike
         Nitr_spike_t_mg= (C_Spike*CO3_mean_N_perc/100)+(N_Spike*NH4_mean_N_perc/100))%>%#total nitrogen spike
  mutate(Samp_N_mg= samp_N_mg-Nitr_spike_t_mg,
         Samp_C_mg= samp_C_mg-Carb_spike_t_mg,
         perc_N_Samp= (Samp_N_mg/weight_mg)*100,
         perc_C_Samp= (Samp_C_mg/weight_mg)*100)

#bind spiked samples with nonspiked samples
CN_All_MER<-bind_rows(CN_All%>% select(Tube_ID,weight_mg,perc_N_Samp,perc_C_Samp) %>% 
                        mutate(Spiked='Yes'),CN_nonspike %>% mutate(Spiked='NO'))

Site_Tube_info<-readRDS("raw_data/Soil_data/Nute_plot_info.RDS") 

CN_MER<-Site_Tube_info %>% 
  distinct(Site,Plot,Tube_ID) %>% 
  left_join(CN_All_MER) %>% 
  filter(!is.na(perc_N_Samp))%>% 
  mutate(Note = 
           case_when(
             Tube_ID == "WT16"~ " mod,xtreme,minor D, maybe beads in sample",
             Tube_ID == "Wa03"~ " FULL OF WEIRD POWDER, MAYBE TERMITES?",
             Tube_ID == "3"~ " Lots of beads",
             Tube_ID == "Q1"~ " one bag had extreme damage",
             Tube_ID == "WT11"~ " 2 minor D, 1 Mod D"
           ))


saveRDS(CN_MER, 'Processed_data/Stoich/CN_Final_MER.Rdata')


#Tube_biomass <- read_excel("raw_data/Tube_biomass_MER_ID.xlsx")

CN_MER %>% 
  ggplot()+
  geom_boxplot(aes(x=Spiked, y= perc_N_Samp))

CN_MER %>% 
  ggplot()+
  geom_boxplot(aes(x=Spiked, y= perc_C_Samp))

CN_MER %>% 
  ggplot()+
  geom_point(aes(x=weight_mg , y= perc_C_Samp, color=Spiked))+
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) 

CN_MER %>% 
  ggplot()+
  geom_point(aes(x=weight_mg , y= perc_N_Samp, color=Spiked))+
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) 

library(car)
# -------------------------------
# ANOVA: Does weight_mg or Spiked explain more variation?
# -------------------------------

# Nitrogen %
m1_N <- lm(perc_N_Samp ~ weight_mg, data = CN_MER )
m2_N <- lm(perc_N_Samp ~ Spiked, data = CN_MER)

anova_N <- anova(m1_N, m2_N)

# Carbon %
m1_C <- lm(perc_C_Samp ~ weight_mg, data = CN_MER)
m2_C <- lm(perc_C_Samp ~ Spiked, data = CN_MER)

anova_C <- anova(m1_C, m2_C)

# Summaries + R² for each model
summary(m1_N)$adj.r.squared
summary(m2_N)$adj.r.squared
summary(m1_C)$adj.r.squared
summary(m2_C)$adj.r.squared

# Print results
cat("\nANOVA Nitrogen:\n")
print(anova_N)

cat("\nANOVA Carbon:\n")
print(anova_C)

# --------------------------
# Combined models
# --------------------------

# Nitrogen
m_both_N <- lm(perc_N_Samp ~ weight_mg + Spiked, data = CN_MER)
summary(m_both_N)  # model fit
Anova(m_both_N, type = 2)  # Type II ANOVA: unique effect of each predictor

# Carbon
m_both_C <- lm(perc_C_Samp ~ weight_mg + Spiked, data = CN_MER)
summary(m_both_C)
Anova(m_both_C, type = 2)

# Compare R²
summary(m_both_N)$adj.r.squared
summary(m_both_C)$adj.r.squared



