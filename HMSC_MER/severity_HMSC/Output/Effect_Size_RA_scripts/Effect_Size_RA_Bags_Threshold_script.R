# Load libraries
library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(performance)

# Read and process input data
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
myco_dat_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')


#load most up to date PA model
load('HMSC_MER/results/Beta1.RData')

Beta_estimates_PA <- Beta1$mean
prob_x_greater_0 <- Beta1$pos
prob_x_less_0 <- Beta1$neg

# Total sequence counts
hyph_count <- sum(myco_dat_hyph$count, na.rm = TRUE)
soil_count <- sum(myco_dat_soil$count, na.rm = TRUE)

# Relative abundance
RA_hyph <- myco_dat_hyph %>% group_by(OTU) %>% summarise(RA_hyph = sum(count, na.rm = TRUE)/hyph_count)
RA_soil <- myco_dat_soil %>% group_by(OTU) %>% summarise(RA_soil = sum(count, na.rm = TRUE)/soil_count)

# Merge and compute pseudocount
RA_joined <- full_join(RA_hyph, RA_soil, by = "OTU") %>%
  mutate(across(everything(), ~replace_na(.x, 0)))

min_ratio <- RA_joined %>%
  filter(RA_hyph > 0 & RA_soil > 0) %>%
  mutate(ratio = RA_hyph / RA_soil) %>%
  summarise(pseudo = min(ratio) / 2) %>%
  pull(pseudo)

# Output containers
thresholds <- c(0.9, 0.8, 0.7, 0.6,0.5)
results_all <- list()
plots_all <- list()
combined_plot_data <- list()

# Loop across thresholds
for (thresh in thresholds) {
  
  # Filter significant positive and negative effects
  beta <- "Severity"
  
  sig_pos <- Beta_estimates_PA %>%
    select(Species, !!beta) %>%
    left_join(prob_x_greater_0 %>% select(Species, !!beta), by = "Species", suffix = c("_mean", "_pos")) %>%
    filter(.data[[paste0(beta, "_pos")]] >= thresh) %>%
    rename(Beta_fire = ends_with("_mean")) %>%
    mutate(Sign = "Positive")
  
  sig_neg <- Beta_estimates_PA %>%
    select(Species, !!beta) %>%
    left_join(prob_x_less_0 %>% select(Species, !!beta), by = "Species", suffix = c("_mean", "_neg")) %>%
    filter(.data[[paste0(beta, "_neg")]] >= thresh) %>%
    rename(Beta_fire = ends_with("_mean")) %>%
    mutate(Sign = "Negative")
  
  Sig_OTUs <- bind_rows(sig_pos, sig_neg) %>%
    rename(OTU = Species)
  
  # Join with RA data
  RA_data <- RA_joined %>%
    mutate(ratio = (RA_hyph + min_ratio) / (RA_soil + min_ratio),
           log_ratio = log10(ratio))
  
  RA_ratio <- inner_join(RA_data, Sig_OTUs, by = "OTU") 
  
  # Run models and store stats
  model_types <- list(
    Absolute = RA_ratio %>% mutate(Response = abs(Beta_fire)),
    Positive = RA_ratio %>% filter(Beta_fire > 0) %>% mutate(Response = Beta_fire),
    Negative = RA_ratio %>% filter(Beta_fire < 0) %>% mutate(Response = Beta_fire)
  )
  
  model_results <- lapply(names(model_types), function(type) {
    dat <- model_types[[type]]
    if (nrow(dat) < 3) return(NULL)  # skip if too few points
    model <- lm(Response ~ log_ratio, data = dat)
    smry <- summary(model)
    anov <- Anova(model, Type='F')
    tibble(
      Threshold = thresh,
      Type = type,
      N = nrow(dat),
      Slope = coef(smry)[2, 1],
      SE = coef(smry)[2, 2],
      F = anov$`F value`[1],
      DF_num = anov$Df[1],
      DF_denom = anov$Df[2],
      P = coef(smry)[2, 4],
      R2 = smry$r.squared
    )
  })
  
  model_results <- lapply(names(model_types), function(type) {
    dat <- model_types[[type]]
    if (nrow(dat) < 3) return(NULL)  # skip if too few points
    model <- lm(Response ~ log_ratio, data = dat)
    smry <- summary(model)
    anov <- Anova(model, Type='F')
    tibble(
      Threshold = thresh,
      Type = type,
      Slope = coef(smry)[2, 1],
      SE = coef(smry)[2, 2],
      F = anov$`F value`[1],
      DF_num = anov$Df[1],
      DF_denom = anov$Df[2],
      P = coef(smry)[2, 4],
      R2 = smry$r.squared
    )
  })
  
  results_all[[as.character(thresh)]] <- bind_rows(model_results)
  
  
  # Plots
  plot_list <- lapply(names(model_types), function(type) {
    dat <- model_types[[type]]
    if (nrow(dat) < 3) return(NULL)
    
    smry <- results_all[[as.character(thresh)]] %>% filter(Type == type)
    
    plt <-ggplot(dat, aes(x = log_ratio, y = Response)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_smooth(method = "lm", color = "blue",se=FALSE, linewidth=2) +
      labs(
        #title = paste(type, "- Fire Effect vs log10(RA Hyphae / RA Soil)"),
        subtitle = paste("Threshold =", thresh),
        x = "log10(Hyphal / Soil RA)",
        y = switch(type,
                   "Absolute" = "|Fire effect size|",
                   "Positive" = "Fire effect size (positive)",
                   "Negative" = "Fire effect size (negative)")
      ) +
      annotate("text", x = max(dat$log_ratio, na.rm = TRUE)-1,
               y = max(abs(dat$Response), na.rm = TRUE),
               hjust = 0,
               label = paste0("Slope = ", round(smry$Slope, 3),
                              "\nSE = ", round(smry$SE, 2),
                              "\nF = ", round(smry$F, 2),
                              "\nP = ", signif(smry$P, 2),
                              "\nR² = ", round(smry$R2, 3))) +
      theme_minimal()
    
    
    # Store plot and metadata
    combined_plot_data[[length(combined_plot_data) + 1]] <<- tibble(
      Threshold = thresh,
      Type = type,
      Plot = list(plt)  # <- this wraps the ggplot in a list-column
    )
    
    return(plt)
  })
  
  names(plot_list) <- names(model_types)
  plots_all[[as.character(thresh)]] <- plot_list
}

# Combine and save results
final_results <- bind_rows(results_all) %>% 
  arrange(Type,Threshold)

final_results

saveRDS(final_results, 'HMSC_MER/Output/Processed_Data/threshold_Effect_Size_RA.RDS')



# Group plots by type across all thresholds
library(patchwork)

# Convert list to data frame
plot_df <- bind_rows(combined_plot_data)

library(patchwork)

view_combined_plots_by_type <- function(type_name, ncol = 2) {
  plot_df %>%
    filter(Type == type_name) %>%
    arrange(Threshold) %>%
    mutate(Plot = map2(Plot, Threshold, ~ .x + labs(subtitle = paste("Threshold:", .y)))) %>%
    pull(Plot) %>%
    wrap_plots(ncol = ncol) +
    plot_annotation(title = paste("Fire Effect Models by Threshold –", type_name),
                    theme = theme(plot.title = element_text(size = 16, face = "bold")))
}

# Example: view plots by type
view_combined_plots_by_type("Absolute")
view_combined_plots_by_type("Positive")
view_combined_plots_by_type("Negative")
