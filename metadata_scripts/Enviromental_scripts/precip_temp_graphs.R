#load libraries
library(ggplot2)
library(dplyr)

source('Temp_Precip_import.R')

#precip per day#########

#facet all sites
ggplot(Sites_Precip,aes(x = date, y = chirps)) +
         geom_line() +
  facet_wrap(~Site, scales = "free_y")

#loop for each site individually
plots_list <- list()
for (i in unique(Sites_Precip$Site)) {
  plot_data <- subset(Sites_Precip, i == Site)
  
 
  
  current_plot<-ggplot(plot_data, aes(date, chirps)) +
    geom_line()+
    scale_x_date(date_breaks = "1 year", date_labels = "%b %Y") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab('Precipitation (mm)') +
    ggtitle(i)
    #facet_wrap(~Burn_Treatment) not needed
  
  plots_list[[i]] <- current_plot
}

#print graphs
for (i in unique(Sites_Precip$Site)) {
  ggsave(plots_list[[i]], file=paste0("plot_y_", i,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")
}

####for monthly averages####


#calc monthly precip
Sites_Precip <-
  Sites_Precip %>%
  mutate(MonthYear = format(date, "%Y-%m")) %>%
  group_by(MonthYear) %>%
  mutate( AvgPrecipitation = mean(chirps, na.rm = TRUE))

ggplot(Sites_Precip,aes(x = date, y = AvgPrecipitation)) +
  geom_col() +
  facet_wrap(~Site, scales = "free_y")



m_plots_list <- list()
for (i in unique(Sites_Precip$Site)) {
  plot_data <- subset(Sites_Precip, i == Site)
  
  
  
  current_plot<-ggplot(plot_data, aes(date, AvgPrecipitation)) +
    geom_col()+
    scale_x_date(date_breaks = "1 year", date_labels = "%b %Y") +
    ylab('Precipitation (mm)') +
    ggtitle(i)
  #facet_wrap(~Burn_Treatment) not needed
  
  m_plots_list[[i]] <- current_plot
}
for (i in unique(Sites_Precip$Site)) {
  ggsave(plots_list[[i]], file=paste0("plot_month_avg_", i,".png"),path='precip_graph/', width = 14, height = 10, units = "cm")
}


