library(Hmsc)

#set wd
localDir = "HMSC_MER"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")


load(file=file.path(dataDir,"unfitted_models.RData"))

# fit models -- one in each screen
message("Starting MCMC sampling...")
start_time <- Sys.time()

models2 <- sampleMcmc(models[[2]], samples = samples, thin = thin, transient = transient, nChains = nChains,
                      nParallel = nParallel, verbose = TRUE)

end_time <- Sys.time()
message("MCMC sampling completed in ", round(difftime(end_time, start_time, units = "hours"), 2), " hours")

filename = file.path(modelDir,paste("MER_model_thin_", as.character(thin),
                                    "_samples_", as.character(samples),
                                    "_chains_",as.character(nChains),
                                    "_2.Rdata",sep = ""))

save(models2,thin, samples, nChains, transient,file=filename)

models<-list()

models[[2]] <- models2

#check model convergence####

mpost2 <- convertToCodaObject(models[[2]])
psrf.beta2 <- gelman.diag(mpost2$Beta, multivariate=FALSE)$psrf
summary(psrf.beta2)
quantile(psrf.beta2[, 1], seq(0, 1, 0.05))

#traceplot(mpost2$Beta, ask = FALSE, col = 1:6, main = "Traceplot for Beta Parameters")

predY2 <- computePredictedValues(models[[2]], expected=FALSE, nParallel=2)
MF2 <- evaluateModelFit(hM=models[[2]], predY=predY2)
mean(MF2$R2, na.rm=T)  # 0.953

hist(psrf.beta2, main="COP", xlab="psrf (Beta)")

# # two-fold cross validation.
partition_2 <- createPartition(models[[2]], nfolds=2)
preds_new_2 <- computePredictedValues(models[[2]], partition=partition_2, nParallel=2)
MF_pred_2 <- evaluateModelFit(hM=models[[2]], predY=preds_new_2)

#plot in 2x2 layout
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # optional: adjust margins

# Plot explanatory vs predictive power using Tjur R²
plot(MF1$R2, MF_pred_2$R2,
     xlim = c(-1, 1), ylim = c(-1, 1),
     xlab = "Explanatory Power (Tjur R²)",
     ylab = "Predictive Power (Tjur R²)",
     main = paste0("Tjur R²\nmean(expl) = ", round(mean(MF1$TjurR2, na.rm = TRUE), 3),
                   ", mean(pred) = ", round(mean(MF_pred_1$TjurR2, na.rm = TRUE), 3)))
abline(0, 1, col = "blue", lty = 2)
abline(h = 0, v = 0, lty = 3)
# • Points near the diagonal = good generalization
# • Points above the line: model predicts better on test data (rare)
# • Points below the line: model performs worse in prediction → overfitting
# • Tjur R² values near or below 0 suggest poor discrimination between presence and absence

# Plot explanatory vs predictive power using AUC
plot(MF2$AUC, MF_pred_2$AUC,
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Explanatory Power (AUC)",
     ylab = "Predictive Power (AUC)",
     main = paste0("AUC\nmean(expl) = ", round(mean(MF1$AUC, na.rm = TRUE), 3),
                   ", mean(pred) = ", round(mean(MF_pred_1$AUC, na.rm = TRUE), 3)))
abline(0, 1, col = "blue", lty = 2)
abline(h = 0.5, v = 0.5, lty = 3)
# • AUC values > 0.7 = acceptable; > 0.8 = good; > 0.9 = excellent
# • Points near y = x = good match between explanatory and predictive performance
# • Points below diagonal = model performs worse on test data → possible overfitting
# • Points near 0.5 on y-axis = model not better than random in prediction


# Plot RMSE (lower = better)
plot(MF2$RMSE, MF_pred_2$RMSE,
     xlim = range(MF1$RMSE, MF_pred_1$RMSE, na.rm = TRUE),
     ylim = range(MF1$RMSE, MF_pred_1$RMSE, na.rm = TRUE),
     xlab = "Explanatory Power (RMSE)",
     ylab = "Predictive Power (RMSE)",
     main = paste0("RMSE\nmean(expl) = ", round(mean(MF1$RMSE, na.rm = TRUE), 3),
                   ", mean(pred) = ", round(mean(MF_pred_1$RMSE, na.rm = TRUE), 3)))
# Add reference line (perfect match between training and prediction error)
abline(0, 1, col = "blue", lty = 2)
# --- What to look for ---
# • Lower RMSE is better (closer to 0)
# • Points below the diagonal = lower RMSE in prediction → excellent generalization (rare)
# • Points above diagonal = model fits training data better than test → possible overfitting
# • Large RMSE values = poor fit (could indicate very rare or noisy OTUs)
