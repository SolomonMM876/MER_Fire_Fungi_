library(tidyverse)
library(Hmsc)

#set wd
localDir = "HMSC_MER/severity_HMSC"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")
resultDir = file.path(localDir, "results")
if (!dir.exists(dataDir)) dir.create(dataDir, recursive = TRUE)
if (!dir.exists(modelDir)) dir.create(modelDir, recursive = TRUE)
if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)


#read in Site Data
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS")

#select data 
Site<-All_Sites %>% 
  select(Site,Plot,Lat,Long)

MER_fire_regime<-readRDS("Processed_data/metadata/fire/Fire_regime_plot.rds") %>% 
  select(Site,Plot,Severity=fire_severity_num)


##################################
#soil fungal communities
soil_comm<-readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')%>% 
  mutate(source='soil') %>% 
  rowwise() %>%
  mutate(total_reads = sum(c_across(where(is.numeric)), na.rm = TRUE))%>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>% 
  relocate(Plot,total_reads)

#read in host frequency data
myco_freq_plot <- readRDS('Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata') %>% 
  select(Plot,freq_AM,freq_ECM)
# myco_freq_plot %>% left_join(Soil) %>% 
#   ggplot(aes(x=Site))+
#   geom_point(aes(y=freq_AM))+
#   geom_point(aes(y=freq_ECM),color='red')

##############################

################
#join data together
All<-Site %>% 
  left_join(MER_fire_regime) %>% 
  left_join(myco_freq_plot) %>% 
  left_join(soil_comm) 

All %>% filter(is.na(freq_AM)) %>% distinct(Plot) %>% pull()
#9 plots have no host frequency data VCEF07,VCEF08,WASP05,WABEW04,WABEW05, WABEW06, WABEW20,WABEW21, WABEW22

All %>% filter(is.na(total_reads)) %>% distinct()
#VCRF20 no biomass and excluded
#no bag data for NSS06, NSS07and WAEW23
#no nute data WAEW15 and WAEW16
#no NO3 for SAM105
#VCEF01,VCEF07,VCEF08,WAEW09,WAEW19,WASS01,WASS06,WABVT11,QDRF08 all no soil data
data <-All%>%
  mutate( total_reads=log10(total_reads)
  )%>% 
  filter(!is.na(total_reads)
         ) %>%  
  filter(!is.na(freq_AM)) %>%
  filter(!is.na(Severity))

#select factors
XData <-data %>% 
  select(Severity,freq_AM, freq_ECM,total_reads)#

#Select Community Data
YData<-data%>%
  select(starts_with('ITSall'))

Y = as.matrix(YData)

#look at relative abundances of species
library(ggplot2)
data.frame(
  Species = names(YData),
  Relative_Abundance = colSums(YData, na.rm = TRUE)/sum(YData, na.rm = TRUE))%>%
  ggplot( aes(x = reorder(Species, Relative_Abundance), y = Relative_Abundance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip to make labels more readable
  labs(title = "Relative Abundance of Fungal Species",
       x = "Species",
       y = "Relative Abundance") +
  theme_minimal() 
###########################
#Spatial variables

# random effect structure -- spatial latent variable
sData=data %>% group_by(Plot) %>% 
      column_to_rownames('Plot') %>% select(Long,Lat) %>% 
  as.matrix()

#plot points to check
library(maps)
# Plot a map of Australia
map("world", regions = "Australia", 
    xlim = c(110, 155), ylim = c(-45, -10),  # Set bounds manually
    col = "lightgray", fill = TRUE, bg = "white", lwd = 0.5, main = "Australian Sites")
# Add your points
points(sData, pch = 21, bg = "red", col = "black")

# STUDY DESIGN

studyDesign = data.frame(
  Plot = data$Plot,
  Site = data$Site,
  Spatial = as.factor(data$Plot)   # link spatial random level to Plot IDs
  
  ) %>%
  mutate(across(everything(), as.factor))

#set random effects of site and plot
rL_Plot= HmscRandomLevel(units = studyDesign$Plot)
rL_Site = HmscRandomLevel(units = studyDesign$Site)

# random effect structure -- spatial latent variable
rL_spatial <- HmscRandomLevel(sData=sData, longlat=TRUE)

# REGRESSION MODEL FOR ENVIRONMENTAL COVARIATES.
XFormula = ~Severity+ freq_AM+ freq_ECM+ total_reads

####Run models

## prepare models
models <- list()

# presence-absence model
models[[1]] <- Hmsc(Y = 1*(Y>0), XData = XData, XFormula = XFormula, 
                    studyDesign = studyDesign, 
                    ranLevels = list(Site=rL_Site,Plot=rL_Plot,Spatial= rL_spatial),
                    distr="probit",
                    YScale = TRUE)

# abundance model
# Y2 <- Y
# Y2[Y2==0] <- NA
# models[[2]] <- Hmsc(Y = log(Y2), XData = XData, XFormula = XFormula, 
#                     studyDesign = studyDesign, 
#                     ranLevels = list(Plot=rL_Plot, Site=rL_Site),
#                     distr="normal",
#                     YScale = TRUE)

#To understand what was done with the above script, let us look at the objects:
# The study design simply lists the id of each plot
  head(studyDesign)

#The right-hand side (rL) describes the structure of 
# random effect. To see the structure, evaluate
rL_Site
#rL_Plot

rL_spatial

# It is always a good idea to look at how the XData and XFormula are translate to the design matrix X. 
# These can be seen as e.g. by evaluating the following

models[[1]]$XFormula
head(models[[1]]$XData)
head(models[[1]]$X)

# Note further that there is "1 trait" even if we did not define any traits. To see why this is the case, evaluate
head(models[[1]]$Tr)

#good idea to look at the model object as well:
#models[[2]]

# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}

# model parameters
thin <- 20
samples <- 2000
nChains <- 2
transient <- 0.3*thin*samples
iterations<-thin*samples
nParallel<-nChains

# COMBINING AND SAVING MODELS (START)
##################################################################################################

names(models) = c("Severity_Soil_MER_PA")
# export workspace
save(models,thin, samples, nChains, transient,nParallel, file = file.path(dataDir, "unfitted_models_severity.RData"))
