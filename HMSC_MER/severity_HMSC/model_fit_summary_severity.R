library(tidyverse)
library(Hmsc)

# --- Variance partitioning ---
# Only computing for model 1 (PA)
vp1 <- computeVariancePartitioning(models[[1]])

vp_PA <- rowMeans(vp1$vals, na.rm = TRUE)

# --- Model fit metrics ---
compare_models <- function(hM1) {
  p1 <- computePredictedValues(hM1)
  
  mf1 <- evaluateModelFit(hM1, predY = p1)
  
  r2_1 <- if (!is.null(mf1$R2)) mean(mf1$R2, na.rm = TRUE) else NA
  
  tibble(
    Model = c("PA"),
    Mean_R2 = c(r2_1),
    AUC = c(mean(mf1$AUC, na.rm = TRUE)),
    TjurR2 = c(mean(mf1$TjurR2, na.rm = TRUE)),
    RMSE = c(mean(mf1$RMSE, na.rm = TRUE))
  )
}

fit_tbl <- compare_models(models[[1]])

# --- Extract PSRF values (mean, min, max) ---
extract_psrf <- function(psrf_df) {
  df <- as.data.frame(psrf_df)
  df <- df %>% filter(Var2 == "  Point est.")
  
  # Helper to grab numbers out of strings like "Mean   :1.0014"
  grab_num <- function(pattern) {
    val <- df$Freq[grepl(pattern, df$Freq)]
    round(as.numeric(sub(".*:", "", val)), 4)
  }
  
  tibble(
    Mean = grab_num("Mean"),
    Min  = grab_num("Min"),
    Max  = grab_num("Max")
  )
}

## check model convergence####
mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate=FALSE)$psrf

psrf.b1 <- summary(psrf.beta1)

psrf1 <- extract_psrf(psrf.b1)

scale_red <- tibble(
  Model = c("PA"),
  Beta_mean=psrf1$Mean,
  Beta_scale_reduction = c(
    sprintf("(%.4f–%.4f)", psrf1$Min, psrf1$Max)
  )
)

# --- Build final table (PA ONLY) ---
final_tbl_cover_ITS <- tibble(
  Metric = c(
    "μ VP FireSeverity",
    "μ VP Frequency AM Hosts",
    "μ VP Frequency EcM Hosts",
    "μ VP Total Reads",
    "μ VP Random (Site)",
    "μ VP Random Spatial (Plot)",
    "AUC (discrimination ability)",
    "Tjur R2 (explanatory power)",
    "RMSE (accuracy)",
    "μ Scale reduction of Beta parameters (Range min–max)",
    "Range (min–max) of scale reduction of Beta parameters" 
  ),
  `Fire Severity` = c(
    round(vp_PA["Severity"], 3), 
    round(vp_PA["freq_AM"], 3),
    round(vp_PA["freq_ECM"], 3),
    round(vp_PA["total_reads"], 3), 
    round(vp_PA["Random: Spatial"], 3), 
    round(vp_PA["Random: Site"], 3), 
    round(fit_tbl$AUC, 3), 
    round(fit_tbl$TjurR2, 3), 
    round(fit_tbl$RMSE, 3), 
    scale_red$Beta_mean,
    scale_red$Beta_scale_reduction
  )
)

# View or export
print(final_tbl_cover_ITS)

# Ensure the directory exists before saving, or just save
saveRDS(final_tbl_cover_ITS, file= 'HMSC_MER/severity_HMSC/results/model_eval_PA_severity.rds')
