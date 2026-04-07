library(Hmsc)
library(tidyverse)

# set local wd

localDir <- "HMSC_MER"
dataDir <- file.path(localDir, "data")
modelDir <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")

#### load models-----------

if (!dir.exists(resultDir)) dir.create(resultDir)
set.seed(1)

# model parameters
thin <- 50
samples <- 2300
nChains <- 2

# load models

# filename=file.path(modelDir, paste0("MER_fire_model_thin_",as.character(thin),"_samples_",as.character(samples),
#                                     "_chains_",as.character(nChains),'_1.Rdata'))
# load(filename)


# filename=file.path(modelDir, paste0("MER_model_thin_",as.character(thin),"_samples_",as.character(samples),
#                                     "_chains_",as.character(nChains),'_2.Rdata'))
# load(filename)
#
models <- list()

# models[[1]] <- models1
# models[[2]] <- models2
# rm(models1, models2)

# save as one model
#########
filename <- file.path(modelDir, paste0(
     "MER_host_freq_model_thin_", as.character(thin), "_samples_", as.character(samples),
     "_chains_", as.character(nChains), ".Rdata"
))
# save(models, file=filename)


file.exists(filename) # Should return TRUE if path is correct

load(filename)

## check model convergence####
mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate = FALSE)$psrf
summary(psrf.beta1)
quantile(psrf.beta1[, 1], seq(0, 1, 0.05))

# traceplot(mpost1$Beta, ask = FALSE, col = 1:6, main = "Traceplot for Beta Parameters")


# mpost2 <- convertToCodaObject(models[[2]])
# psrf.beta2 <- gelman.diag(mpost2$Beta, multivariate=FALSE)$psrf
# summary(psrf.beta2)
# quantile(psrf.beta2[, 1], seq(0, 1, 0.05))

# traceplot(mpost2$Beta, ask = FALSE, col = 1:6, main = "Traceplot for Beta Parameters")


predY1 <- computePredictedValues(models[[1]], expected = FALSE, nParallel = 2)
MF1 <- evaluateModelFit(hM = models[[1]], predY = predY1)
mean(MF1$TjurR2, na.rm = T) # 0.109
mean(MF1$AUC, na.rm = T) # 0.929
# WAIC <- computeWAIC(models[[1]])

# predY2 <- computePredictedValues(models[[2]], expected=FALSE, nParallel=2)
# MF2 <- evaluateModelFit(hM=models[[2]], predY=predY2)
# mean(MF2$R2, na.rm=T)  # 0.956


hist(psrf.beta1, main = "Probit", xlab = "psrf (Beta)")
# hist(psrf.beta2, main="COP", xlab="psrf (Beta)")

# # two-fold cross validation.
partition_1 <- createPartition(models[[1]], nfolds = 2)
preds_new_1 <- computePredictedValues(models[[1]], partition = partition_1, nParallel = 2)
MF_pred_1 <- evaluateModelFit(hM = models[[1]], predY = preds_new_1)
# # two-fold cross validation.
# partition_2 <- createPartition(models[[2]], nfolds=2)
# preds_new_2 <- computePredictedValues(models[[2]], partition=partition_2, nParallel=2)
# MF_pred_2 <- evaluateModelFit(hM=models[[2]], predY=preds_new_2)
#
# save(MF_pred_1,MF_pred_2,file='HMSC_MER/results/mcmc_convergence_pred.RData')
# load('HMSC_ABS/ABS_2nd_Rnd/results/mcmc_convergence_pred.RData')

# plot in 2x2 layout
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2)) # optional: adjust margins

# Plot explanatory vs predictive power using Tjur R²
plot(MF1$TjurR2, MF_pred_1$TjurR2,
     xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Explanatory Power (Tjur R²)",
     ylab = "Predictive Power (Tjur R²)",
     main = paste0(
          "Tjur R²\nmean(expl) = ", round(mean(MF1$TjurR2, na.rm = TRUE), 3),
          ", mean(pred) = ", round(mean(MF_pred_1$TjurR2, na.rm = TRUE), 3)
     )
)
abline(0, 1, col = "blue", lty = 2)
abline(h = 0, v = 0, lty = 3)
# • Points near the diagonal = good generalization
# • Points above the line: model predicts better on test data (rare)
# • Points below the line: model performs worse in prediction → overfitting
# • Tjur R² values near or below 0 suggest poor discrimination between presence and absence

# Plot explanatory vs predictive power using AUC
plot(MF1$AUC, MF_pred_1$AUC,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Explanatory Power (AUC)",
     ylab = "Predictive Power (AUC)",
     main = paste0(
          "AUC\nmean(expl) = ", round(mean(MF1$AUC, na.rm = TRUE), 3),
          ", mean(pred) = ", round(mean(MF_pred_1$AUC, na.rm = TRUE), 3)
     )
)
abline(0, 1, col = "blue", lty = 2)
abline(h = 0.5, v = 0.5, lty = 3)
# • AUC values > 0.7 = acceptable; > 0.8 = good; > 0.9 = excellent
# • Points near y = x = good match between explanatory and predictive performance
# • Points below diagonal = model performs worse on test data → possible overfitting
# • Points near 0.5 on y-axis = model not better than random in prediction


# Plot RMSE (lower = better)
plot(MF1$RMSE, MF_pred_1$RMSE,
     xlim = range(MF1$RMSE, MF_pred_1$RMSE, na.rm = TRUE),
     ylim = range(MF1$RMSE, MF_pred_1$RMSE, na.rm = TRUE),
     xlab = "Explanatory Power (RMSE)",
     ylab = "Predictive Power (RMSE)",
     main = paste0(
          "RMSE\nmean(expl) = ", round(mean(MF1$RMSE, na.rm = TRUE), 3),
          ", mean(pred) = ", round(mean(MF_pred_1$RMSE, na.rm = TRUE), 3)
     )
)
# Add reference line (perfect match between training and prediction error)
abline(0, 1, col = "blue", lty = 2)
# --- What to look for ---
# • Lower RMSE is better (closer to 0)
# • Points below the diagonal = lower RMSE in prediction → excellent generalization (rare)
# • Points above diagonal = model fits training data better than test → possible overfitting
# • Large RMSE values = poor fit (could indicate very rare or noisy OTUs)


# save.image(file='HMSC_MER/results/mcmc_convergence_plots.RData')

########
# Evaluate plot for models[[2]]

# 2x2 layout, adjust as needed
# par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
#
# # === R² ===
# plot(MF2$R2, MF_pred_2$R2,
#      xlim = c(0, 1), ylim = c(0, 1),
#      xlab = "Explanatory Power (R²)",
#      ylab = "Predictive Power (R²)",
#      main = paste0("R²\nmean(expl) = ", round(mean(MF2$R2, na.rm = TRUE), 3),
#                    ", mean(pred) = ", round(mean(MF_pred_2$R2, na.rm = TRUE), 3)))
# abline(0, 1, col = "blue", lty = 2)
# abline(h = 0, v = 0, lty = 3)
# # • Points near diagonal = good generalization
# # • Points below = worse on test set → overfitting
# # • R² near 0 = poor model fit
#
# # === RMSE ===
# plot(MF2$RMSE, MF_pred_2$RMSE,
#      xlim = range(MF2$RMSE, MF_pred_2$RMSE, na.rm = TRUE),
#      ylim = range(MF2$RMSE, MF_pred_2$RMSE, na.rm = TRUE),
#      xlab = "Explanatory Power (RMSE)",
#      ylab = "Predictive Power (RMSE)",
#      main = paste0("RMSE\nmean(expl) = ", round(mean(MF2$RMSE, na.rm = TRUE), 3),
#                    ", mean(pred) = ", round(mean(MF_pred_2$RMSE, na.rm = TRUE), 3)))
# abline(0, 1, col = "blue", lty = 2)
# # • Lower = better
# # • Points above diagonal = worse generalization
# # • Points far from diagonal = fit mismatch
#
# # Optional: Add histogram of PSRF values for MCMC convergence
# hist(psrf.beta2, main="MCMC Convergence for Beta (Model 2)", xlab="PSRF (Beta)")
# abline(v=1.1, col="red", lty=2)
# # • PSRF < 1.1 suggests convergence
#
# # Optional: Plot residuals or distribution of prediction errors (if desired)
# hist(MF2$RMSE, main="RMSE Distribution (Model 2)", xlab="RMSE")

# Save plots and objects
# save.image(file = 'HMSC_MER/results/mcmc_convergence_plots_model2.RData')







# variance partitioning#####
library(ggplot2)
library(reshape2)

vp1 <- computeVariancePartitioning(models[[1]])
# See how much variance each fixed effect explains per species
head(vp1$vals)
# View mean variance explained by each covariate across species
print(vp1$R2T) # total R2 per species

# Prepare data
vp1_df <- as.data.frame(t(vp1$vals)) # species in rows, variables in columns
vp1_df$species <- rownames(vp1_df)
vp1_long <- melt(vp1_df, id.vars = "species")

# Plot
vp_PA <- ggplot(vp1_long, aes(x = species, y = value, fill = variable)) +
     geom_bar(stat = "identity", color = "black", linewidth = 0.05) +
     labs(y = "Proportion of explained variance", x = "OTUs", fill = "Covariate") +
     theme_classic() +
     theme(
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 16, face = "bold"),
          axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 16, face = "bold"),
          axis.line = element_line(linewidth = 2, colour = "black"),
          legend.title = element_text(size = 18, face = "bold")
     )
print(vp_PA)

ggsave("Plots/vp_PA.png", plot = vp_PA, width = 15, height = 7.5, bg = "white")

# ABUNDANCE
# vp2 <- computeVariancePartitioning(models[[2]])
#
# # Prepare data
# vp2_df <- as.data.frame(t(vp2$vals))  # species in rows, variables in columns
# vp2_df$species <- rownames(vp2_df)
# vp2_long <- melt(vp2_df, id.vars = "species")
#
# # Plot
# ggplot(vp2_long, aes(x = species, y = value, fill = variable)) +
#   geom_bar(stat = "identity") +
#   labs(y = "Proportion of explained variance", x = "Species", fill = "Covariate") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
#



## get posterior estimates

# beta
postBeta1 <- getPostEstimate(models[[1]], parName = "Beta", start = 1)
me <- as.data.frame(t(postBeta1$mean))
me <- cbind(models[[1]]$spNames, me)
colnames(me) <- c("Species", models[[1]]$covNames)
po <- as.data.frame(t(postBeta1$support))
po <- cbind(models[[1]]$spNames, po)
colnames(po) <- c("Species", models[[1]]$covNames)
ne <- as.data.frame(t(postBeta1$supportNeg))
ne <- cbind(models[[1]]$spNames, ne)
colnames(ne) <- c("Species", models[[1]]$covNames)
Beta1 <- list(mean = me, pos = po, neg = ne)

# postBeta2 <- getPostEstimate(models[[2]], parName='Beta', start=1)
# me = as.data.frame(t(postBeta2$mean))
# me = cbind(models[[2]]$spNames,me)
# colnames(me) = c("Species",models[[2]]$covNames)
# po = as.data.frame(t(postBeta2$support))
# po = cbind(models[[2]]$spNames,po)
# colnames(po) = c("Species",models[[2]]$covNames)
# ne = as.data.frame(t(postBeta2$supportNeg))
# ne = cbind(models[[2]]$spNames,ne)
# colnames(ne) = c("Species",models[[2]]$covNames)
# Beta2 <- list(mean = me, pos= po, neg= ne)

# omega
# names(models[[2]]$ranLevels)
# OmegaCor1 <- computeAssociations(models[[2]], start=200, thin=2)
# names(OmegaCor1) <- names(models[[2]]$ranLevels)#assign names to lists
# me = as.data.frame(t(OmegaCor1[["Site"]]$mean))
# support = as.data.frame(t(OmegaCor1[["Site"]]$support))
# OmegaCor1_Site<-list(mean=me,support=support)


# names(models[[2]]$ranLevels)
# OmegaCor2 <- computeAssociations(models[[2]], start=200, thin=2)
# names(OmegaCor2) <- names(models[[2]]$ranLevels)#assign names to lists
# me = as.data.frame(t(OmegaCor2[["Site"]]$mean))
# support = as.data.frame(t(OmegaCor2[["Site"]]$support))
# OmegaCor2_Site<-list(mean=me,support=support)

#
## save results
save(Beta1, file = "HMSC_MER/results/Beta1_host_freq.RData")
#
save(OmegaCor1_Site, file = "HMSC_MER/results/Omega1_host_freq.RData")
