# ==============================================================================
# SCRIPT 01: HMSC Model Fitting (Self-Contained)
# DESCRIPTION: Fits JSDM (Presence-Absence) using inputs from 'data/' folder.
# ==============================================================================

# 1. SETUP & LIBRARIES
library(tidyverse)
library(Hmsc)
library(maps)
library(ggplot2)
library(reshape2)

# Set Working Directory to Project Folder
# setwd("xxx")

# MCMC Configuration
n_chains <- 2
n_samples <- 2300
n_thin <- 50
n_transient <- 0.3 * n_thin * n_samples
n_parallel <- n_chains

# Output Directory
dir.create("HMSC_Output", showWarnings = FALSE)

# 2. DATA IMPORT & WRANGLING
# Load clean CSVs
Site_Meta <- read_csv("data/input_site_metadata.csv", show_col_types = FALSE)
myco_freq_plot <- read_csv("data/input_host_frequency.csv", show_col_types = FALSE)
soil_comm_raw <- read_csv("data/input_soil_community.csv", show_col_types = FALSE)

# Process Community Data
# (Replicating logic: summing numeric reads, replacing NAs)
soil_comm <- soil_comm_raw %>%
  mutate(source = "soil") %>%
  rowwise() %>%
  mutate(total_reads = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(across(starts_with("ITS"), ~ replace_na(.x, 0))) %>%
  relocate(Plot, total_reads)

# Join Datasets
All <- Site_Meta %>%
  left_join(myco_freq_plot, by = "Plot") %>%
  left_join(soil_comm, by = "Plot")

# Filtering & Transformation
data <- All %>%
  mutate(
    Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")),
    total_reads = log10(total_reads)
  ) %>%
  filter(!is.na(total_reads)) %>%
  filter(!is.na(freq_AM))

# 3. DEFINE MODEL MATRICES
# X: Covariates
XData <- data %>%
  dplyr::select(Fire_Treatment, freq_AM, freq_ECM, total_reads)

# Y: Community Matrix (Binary for PA model)
YData <- data %>% select(starts_with("ITSall"))
Y_PA <- 1 * (as.matrix(YData) > 0)

# S: Spatial Coordinates
sData <- data %>%
  group_by(Plot) %>%
  column_to_rownames("Plot") %>%
  select(Long, Lat) %>%
  as.matrix()

# Study Design
studyDesign <- data.frame(
  Plot = as.factor(data$Plot),
  Site = as.factor(data$Site),
  Spatial = as.factor(data$Plot)
)

# 4. HMSC MODEL SPECIFICATION
rL_Site <- HmscRandomLevel(units = studyDesign$Site)
rL_spatial <- HmscRandomLevel(sData = sData, longlat = TRUE)

XFormula <- ~ Fire_Treatment + freq_AM + freq_ECM + total_reads

models <- list()
models[[1]] <- Hmsc(
  Y = Y_PA,
  XData = XData,
  XFormula = XFormula,
  studyDesign = studyDesign,
  ranLevels = list(Site = rL_Site, Spatial = rL_spatial),
  distr = "probit",
  YScale = TRUE
)

# 5. MODEL FITTING (MCMC)
print("Starting MCMC Sampling...")
models[[1]] <- sampleMcmc(
  models[[1]],
  samples = n_samples,
  thin = n_thin,
  transient = n_transient,
  nChains = n_chains,
  nParallel = n_parallel,
  verbose = TRUE
)

# 6. MODEL EVALUATION & SAVING
# Convergence
mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate = FALSE)$psrf
print(summary(psrf.beta1))

# Fit (Explanatory)
predY1 <- computePredictedValues(models[[1]], expected = FALSE, nParallel = n_parallel)
MF1 <- evaluateModelFit(hM = models[[1]], predY = predY1)
print(paste("Mean Tjur R2:", mean(MF1$TjurR2, na.rm = TRUE)))
print(paste("Mean AUC:", mean(MF1$AUC, na.rm = TRUE)))

# Save Model Object
print("Saving model object...")
# Saving into the project_script subfolder
saveRDS(list(model = models[[1]], data = data), "HMSC_Output/models_fitted.RDS")
print("Model saved to project_script/HMSC_Output/models_fitted.RDS")
