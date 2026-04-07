library(tidyverse)

#read in CN
CN_MER<-readRDS('Processed_data/Stoich/CN_Final_MER.Rdata') %>% 
  select(Site,Plot,CN_wt=weight_mg,N=perc_N_Samp,C=perc_C_Samp,CN_spike=Spiked,Note)
#read in TP
MER_P<-readRDS('Processed_data/Stoich/Total_P_MER.Rdata') %>% 
  select(Site,Plot,TP_wt=weight_mg,P=percent_P)
#read in Site Data
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% 
  select(Site,Plot,Fire_Treatment)

#combine
CNP_MER<-CN_MER %>% 
  left_join(MER_P) %>% 
  left_join(All_Sites)


#compute ratios
CNP_MER<-CNP_MER %>% 
  mutate(
    # Set C and N to NA where CN_wt < 0.4
    C = if_else(CN_wt < 0.4, NA_real_, C),
    N = if_else(CN_wt < 0.4, NA_real_, N),
    
    # Set P to NA where TP_wt < 0.4
    P = if_else(TP_wt < 0.4, NA_real_, P)
  ) %>%
  mutate(CN=(C/N),
         CP=(C/P),
         NP=(N/P)
         ) 

CNP_MER %>% 
  ggplot()+
  geom_point(aes(x=CN_wt , y= CN, color=CN_spike))+
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) 

CNP_MER %>% 
  ggplot()+
  geom_point(aes(x=CN_wt , y= CP, color=CN_spike))+
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) 

CNP_MER %>% 
  ggplot()+
  geom_point(aes(x=CN_wt , y= NP, color=CN_spike))+
  scale_color_manual(values = c("Yes" = "red", "No" = "black")) 

library(lme4)
library(ggplot2)
library(performance)
library(lmerTest) # for p-values

#log transform
CNP_MER <- CNP_MER %>%
  mutate(
    log_CN = log(CN),
    log_CP = log(CP),
    log_NP = log(NP)
  )

saveRDS(CNP_MER, 'Processed_data/Stoich/CNP_MER.RDS')


par(mfrow = c(2, 3))
hist(CNP_MER$CN, main = "CN", xlab = "CN")
hist(CNP_MER$CP, main = "CP", xlab = "CP")
hist(CNP_MER$NP, main = "NP", xlab = "NP")
hist(CNP_MER$log_CN, main = "log_CN", xlab = "CN")
hist(CNP_MER$log_CP, main = "log_CP", xlab = "CP")
hist(CNP_MER$log_NP, main = "log_NP", xlab = "NP")
par(mfrow = c(1, 1))


CNP_MER %>% 
  ggplot(aes(x=Fire_Treatment,y=CN, fill=Fire_Treatment))+
  geom_boxplot()+ 
  geom_jitter(width=.2,color='blue')+
  facet_grid(~Site)+ 
  theme_classic()

CNP_MER %>% 
  ggplot(aes(x=Fire_Treatment,y=NP, fill=Fire_Treatment))+
  geom_boxplot()+
  geom_jitter(width=.2,color='blue')+
  facet_grid(~Site)+ 
  theme_classic()

CNP_MER %>% 
  ggplot(aes(x=Fire_Treatment,y=CP, fill=Fire_Treatment))+
  geom_boxplot()+
  geom_jitter(width=.2,color='blue')+
  facet_grid(~Site)+ 
  theme_classic()

library(patchwork)

# Define a bigger custom theme
big_theme <- theme(
  axis.text = element_text(color = "black", size = 16),
  axis.title = element_text(color = "black", size = 18, face = "bold"),
  legend.position = "right",
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16, face = "bold"),
)

# CN plot
CN_MER <- CNP_MER %>% 
  ggplot(aes(x = Fire_Treatment, y = CN, fill = Fire_Treatment)) +
  geom_boxplot() + 
  geom_jitter(width = 0.2, aes(color = Site), size = 4) +
  theme_classic() +
  big_theme+
  theme(legend.position = "none")

# CP plot
CP_MER <- CNP_MER %>% 
  ggplot(aes(x = Fire_Treatment, y = CP, fill = Fire_Treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = Site), size = 4) +
  theme_classic() +
  big_theme+
  theme(legend.position = "none")


# NP plot
NP_MER <- CNP_MER %>% 
  ggplot(aes(x = Fire_Treatment, y = NP, fill = Fire_Treatment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, aes(color = Site), size = 4) +
  theme_classic() +
  big_theme

# Combine using patchwork
CNP_MER_plot <- CN_MER | CP_MER | NP_MER
CNP_MER_plot



#CN
mod_CN <- lmer(log_CN ~ Fire_Treatment + (1 | Site), data = CNP_MER)
summary(mod_CN)
check_model(mod_CN) # from performance package

#CP
mod_CP <- lmer(log_CP ~ Fire_Treatment + (1 | Site), data = CNP_MER)
summary(mod_CP)
check_model(mod_CP)

#NP
mod_NP <- lmer(log_NP ~ Fire_Treatment + (1 | Site), data = CNP_MER)
summary(mod_NP)
check_model(mod_NP)

#CN SPike
mod_CN_spike <- lmer(log_CN ~ Fire_Treatment + CN_spike + (1 | Site), data = CNP_MER)
summary(mod_CN_spike)
check_model(mod_CN_spike)

anova(mod_CN, mod_CN_spike) # if spike is a fixed effect

