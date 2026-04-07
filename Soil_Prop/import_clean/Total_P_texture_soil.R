library(tidyverse)
library(readxl)


CSBP_Lab_Total_P<- read_excel("raw_data/Soil_data/CSBP_Lab_250130_Solomon McMahan_001.xls",
                              sheet= 2)[-(1:4),]#skip empty top row and 4 samples from Ku-ring-gai



