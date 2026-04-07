library(Hmsc)

#set wd
localDir = "HMSC_MER/severity_HMSC"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")

# model parameters
thin <- 50
samples <- 2300
nChains <- 2
transient <- 0.3*thin*samples
iterations<-thin*samples
nParallel<-nChains

load(file=file.path(dataDir,"unfitted_models_severity.RData"))

# fit models -- one in each screen
message("Starting MCMC sampling...")
start_time <- Sys.time()

models1 <- sampleMcmc(models[[1]], samples = samples, thin = thin, transient = transient, nChains = nChains,
                      nParallel = nParallel, verbose = TRUE)

end_time <- Sys.time()
message("MCMC sampling completed in ", round(difftime(end_time, start_time, units = "hours"), 2), " hours")

filename = file.path(modelDir,paste("MER_severity_model_thin_", as.character(thin),
                                    "_samples_", as.character(samples),
                                    "_chains_",as.character(nChains),
                                    "_1.Rdata",sep = ""))
save(models1,file=filename)

load(file=filename)

models<-list()

models[[1]] <- models1

## check model convergence####
mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate=FALSE)$psrf
summary(psrf.beta1)
quantile(psrf.beta1[, 1], seq(0, 1, 0.05))

traceplot(mpost1$Beta, ask = FALSE, col = 1:6, main = "Traceplot for Beta Parameters")


predY1 <- computePredictedValues(models[[1]], expected=FALSE, nParallel=2)
MF1 <- evaluateModelFit(hM=models[[1]], predY=predY1)
mean(MF1$TjurR2, na.rm=T)  # 0.1112463
mean(MF1$AUC, na.rm=T)  # 0.9213555

hist(psrf.beta1, main="Probit", xlab="psrf (Beta)")

# # two-fold cross validation.
partition_1 <- createPartition(models[[1]], nfolds=2)
preds_new_1 <- computePredictedValues(models[[1]], partition=partition_1, nParallel=2)
MF_pred_1 <- evaluateModelFit(hM=models[[1]], predY=preds_new_1)


#plot in 2x2 layout
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # optional: adjust margins

# Plot explanatory vs predictive power using Tjur R²
plot(MF1$TjurR2, MF_pred_1$TjurR2,
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
plot(MF1$AUC, MF_pred_1$AUC,
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
plot(MF1$RMSE, MF_pred_1$RMSE,
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

ggsave("Plots/vp_PA_severity.png", plot = vp_PA, width = 15, height = 7.5, bg = "white")



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


# omega
names(models[[1]]$ranLevels)
OmegaCor1 <- computeAssociations(models[[1]])
names(OmegaCor1) <- names(models[[1]]$ranLevels)#assign names to lists
me = as.data.frame(t(OmegaCor1[["Plot"]]$mean))
support = as.data.frame(t(OmegaCor1[["Plot"]]$support))
OmegaCor1_plot<-list(mean=me,support=support)

#
## save results
save(Beta1, file = "HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")
save(OmegaCor1_plot, file = "HMSC_MER/severity_HMSC/results/omega1_fire_severity.RData")

