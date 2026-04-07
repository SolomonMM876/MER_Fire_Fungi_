
#set local wd

localDir = "HMSC_MER"
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

#prepare models

# random effect structure -- spatial latent variable
rL0 <- HmscRandomLevel(sData=XY %>% 
                         group_by(Plot_ID) %>% 
                         slice(1) %>% 
                         column_to_rownames('Plot_ID') %>% 
                         select(X, Y), longlat=TRUE)


####Run models

## prepare models
models <- list()

# presence-absence model
models[[1]] <- Hmsc(Y = 1*(Y>0), XData = X, XFormula = XFormula, 
                    studyDesign = XY %>% select(Plot_ID, Subplot) %>% mutate(across(everything(), as_factor)), 
                    ranLevels = list(Plot_ID=rL1, Subplot=rL2, Plot_ID=rL0),
                    distr="probit",
                    YScale = TRUE)


# abundance model
Y2 <- Y
Y2[Y2==0] <- NA
models[[2]] <- Hmsc(Y = log(Y2), XData = X, XFormula = XFormula, 
                    studyDesign = XY %>% select(Plot_ID, Subplot) %>% mutate(across(everything(), as_factor)), 
                    ranLevels = list(Plot_ID=rL1, Subplot=rL2, Plot_ID=rL0),
                    distr="normal",
                    YScale = TRUE)

# model parameters for testing
# thin <- 1
# samples <- 2
# nChains <- 2
# transient <- 1

# model parameters
thin <- 20
samples <- 500
nChains <- 2
transient <- 0.3*thin*samples

# export workspace
save(models, thin, samples, nChains, transient, file='output/workspace/hmsc_reference/mcmc_input.RData')




library(Hmsc)
set.seed(1)

#load models
thin <- 20
samples <- 500
nChains <- 2
filename=file.path(modelDir, paste0("models_thin_",as.character(thin),"_samples_",as.character(samples),"_chains_",as.character(nChains),'_1.Rdata'))
load(filename)

thin <- 20
samples <- 500
nChains <- 2
filename=file.path(modelDir, paste0("models_thin_",as.character(thin),"_samples_",as.character(samples),"_chains_",as.character(nChains),'_2.Rdata'))
load(filename)

models<-list()

models[[1]] <- models1
models[[2]] <- models2
rm(models1, models2)

#save as one model
save(models, file='HMSC_ABS/models_updated/updated_model_thin_20_samples_500_Chains_2.RData')


filename_m=file.path(modelDir,'updated_model_thin_20_samples_500_Chains_2.RData')

file.exists(filename_m)  # Should return TRUE if path is correct


load(filename_m)
## check model convergence

mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate=FALSE)$psrf
summary(psrf.beta1)
quantile(psrf.beta1[, 1], seq(0, 1, 0.05))

mpost2 <- convertToCodaObject(models[[2]])
psrf.beta2 <- gelman.diag(mpost2$Beta, multivariate=FALSE)$psrf
summary(psrf.beta2)
quantile(psrf.beta2[, 1], seq(0, 1, 0.05))

predY1 <- computePredictedValues(models[[1]], expected=FALSE, nParallel=2)
MF1 <- evaluateModelFit(hM=models[[1]], predY=predY1)
mean(MF1$TjurR2, na.rm=T)  # 0.07
mean(MF1$AUC, na.rm=T)  # 0.96
#WAIC <- computeWAIC(models[[1]])

predY2 <- computePredictedValues(models[[2]], expected=FALSE, nParallel=2)
MF2 <- evaluateModelFit(hM=models[[2]], predY=predY2)
mean(MF2$R2, na.rm=T)  # 0.96


hist(psrf.beta1, main="Probit", xlab="psrf (Beta)")
hist(psrf.beta2, main="Probit", xlab="psrf (Beta)")


# # two-fold cross validation.
partition <- createPartition(models[[2]], nfolds=2)
preds_new <- computePredictedValues(models[[2]], partition=partition, nParallel=2)
MF_pred <- evaluateModelFit(hM=models[[2]], predY=preds_new)
# 
save.image(file='HMSC_ABS/results_updated/mcmc_convergence_pred.RData')


## get posterior estimates

# beta
postBeta1 <- getPostEstimate(models[[1]], parName='Beta', start=1)
me = as.data.frame(t(postBeta1$mean))
me = cbind(models[[1]]$spNames,me)
colnames(me) = c("Species",models[[1]]$covNames)
po = as.data.frame(t(postBeta1$support))
po = cbind(models[[1]]$spNames,po)
colnames(po) = c("Species",models[[1]]$covNames)
ne = as.data.frame(t(postBeta1$supportNeg))
ne = cbind(models[[1]]$spNames,ne)
colnames(ne) = c("Species",models[[1]]$covNames)
Beta1 <- list(mean = me, pos= po, neg= ne)

postBeta2 <- getPostEstimate(models[[2]], parName='Beta', start=1)
me = as.data.frame(t(postBeta2$mean))
me = cbind(models[[2]]$spNames,me)
colnames(me) = c("Species",models[[2]]$covNames)
po = as.data.frame(t(postBeta2$support))
po = cbind(models[[2]]$spNames,po)
colnames(po) = c("Species",models[[2]]$covNames)
ne = as.data.frame(t(postBeta2$supportNeg))
ne = cbind(models[[2]]$spNames,ne)
colnames(ne) = c("Species",models[[2]]$covNames)
Beta2 <- list(mean = me, pos= po, neg= ne)

# omega
OmegaCor1 <- computeAssociations(models[[1]], start=1)
me = as.data.frame(t(OmegaCor1[[1]]$mean))
me = cbind(models[[2]]$spNames,me) %>% rename(Species='models[[2]]$spNames')
po = as.data.frame(t(OmegaCor1[[1]]$support))
po = cbind(models[[2]]$spNames,po) %>% rename(Species='models[[2]]$spNames')


toPlot = ((OmegaCor1$support>support.level.omega) + (OmegaCor1$support<(1-support.level.omega))>0)*sign(OmegaCor1$mean)

OmegaCor2 <- computeAssociations(models[[2]], start=1)
me = as.data.frame(t(OmegaCor2[[1]]$mean))
me = cbind(models[[2]]$spNames,me) %>% rename(Species='models[[2]]$spNames')
po = as.data.frame(t(OmegaCor2[[1]]$support))
po = cbind(models[[2]]$spNames,po) %>% rename(Species='models[[2]]$spNames')

## save results
save( Beta1, Beta2, file='HMSC_ABS/results_updated/mcmc_output.RData')

