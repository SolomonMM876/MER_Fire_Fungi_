library(Hmsc)
library(tidyverse)

localDir = "HMSC_MER"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")

#########
# model parameters
thin <- 50
samples <- 1500
nChains <- 2
filename=file.path(modelDir, paste0("MER_model_thin_",as.character(thin),"_samples_",as.character(samples),"_chains_",as.character(nChains),'.Rdata'))

#save(models, file=filename)


file.exists(filename)  # Should return TRUE if path is correct

load(filename)

n_pred <- nrow(models[[1]]$XData)  # or however many rows you want
# n_pred <- 10


# Create a new covariate data frame where all samples are "burned"
XData_burned <- models[[1]]$XData %>%
  slice(1:n_pred) %>%
  mutate(
    Fire_Treatment = "B",
    total_reads = mean(models[[1]]$XData$total_reads)
  )

# Use 10 replicated rows of the first entry in studyDesign
studyDesign_burned <- models[[1]]$studyDesign %>%
  slice(rep(1, n_pred))  # repeat the first row for each new prediction

# Use the same random level structure as the model
ranLevels_burned <- models[[1]]$ranLevels

# Predict expected presence probability under Fire_Treatment = "B"
pred_burned <- predict(
  models[[1]],
  XData = XData_burned,
  studyDesign = studyDesign_burned,
  ranLevels = ranLevels_burned,
  expected = TRUE,
  type = "response"
)

# pred_burned: array of dim [samples, species, chains]
# To get overall probabilities per species, take the mean over samples and chains:
pred_array_burn <- array(unlist(pred_burned), dim = c(10, 438, 1600))
burned_prob <- apply(pred_array_burn, 2, mean)  # vector of length = number of OTUs
names(burned_prob) <- colnames(pred_burned[[1]])  # add OTU names

XData_unburned <- models[[1]]$XData %>%
  slice(1:n_pred) %>%
  mutate(
    Fire_Treatment = "U",
    total_reads = mean(models[[1]]$XData$total_reads)
  )

studyDesign_unburned <- studyDesign_burned  # same as above

pred_unburned <- predict(
  models[[1]],
  XData = XData_unburned,
  studyDesign = studyDesign_unburned,
  ranLevels = ranLevels_burned,
  expected = TRUE,
  type = "response"
)

pred_array_unburn <- array(unlist(pred_unburned), dim = c(10, 438, 1600))
unburned_prob <- apply(pred_array_unburn, 2, mean)  # vector of length = number of OTUs
names(unburned_prob) <- colnames(pred_unburned[[1]])  # add OTU names


OTU_likelihoods <- data.frame(
  OTU = models[[1]]$spNames,
  prob_burned = burned_prob,
  prob_unburned = unburned_prob,
  ratio = burned_prob / (burned_prob + unburned_prob)
)

saveRDS(OTU_likelihoods,'HMSC_MER/Output/processed_Data/HMSC_OTU_likelyhood.Rdata')
