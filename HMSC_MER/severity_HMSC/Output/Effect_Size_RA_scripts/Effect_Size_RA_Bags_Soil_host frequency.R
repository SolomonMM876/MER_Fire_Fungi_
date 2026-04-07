library(tidyverse)
library(lme4)
library(lmerTest)
library(car)
library(performance)
library(broom.mixed)
library(patchwork)
library(cowplot)

#--------------------------
# Data import
#--------------------------
myco_tax_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
myco_dat_hyph <- readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
load('HMSC_MER/results/Beta1.RData')
myco_freq_final_sum_plot <- readRDS('Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% 
  dplyr::select(Site, Plot, Severity)

#--------------------------
# Step 1: Calculate log_ratio
#--------------------------
hyph_count <- sum(myco_dat_hyph$count, na.rm = TRUE)
soil_count <- sum(myco_dat_soil$count, na.rm = TRUE)

hyph_RA <- myco_dat_hyph %>%
  group_by(OTU) %>%
  summarise(RA_hyph = sum(count, na.rm = TRUE) / hyph_count, .groups = "drop")

soil_RA <- myco_dat_soil %>%
  group_by(OTU) %>%
  summarise(RA_soil = sum(count, na.rm = TRUE) / soil_count, .groups = "drop")

RA_joined <- full_join(hyph_RA, soil_RA, by = "OTU") %>%
  mutate(across(c(RA_hyph, RA_soil), ~ replace_na(., 0)))

min_ratio <- RA_joined %>%
  filter(RA_hyph > 0 & RA_soil > 0) %>%
  mutate(ratio = RA_hyph / RA_soil) %>%
  summarise(min_ratio = min(ratio, na.rm = TRUE)) %>%
  mutate(pseudocount = min_ratio/2) %>%
  pull(pseudocount)

RA_ratio <- RA_joined %>%
  mutate(ratio = RA_hyph / RA_soil,
         log_ratio = log10(ratio + min_ratio))

#--------------------------
# Step 2: Join metadata
#--------------------------
Site_OTU_RA <- myco_dat_soil %>%
  distinct(Site, Plot, OTU, .keep_all = TRUE) %>%
  left_join(RA_ratio %>% dplyr::select(OTU, log_ratio)) %>%
  left_join(myco_freq_final_sum_plot, by = c("Site", "Plot")) %>%
  left_join(All_Sites, by = c("Site", "Plot"))

#--------------------------
# Step 3: Model function
#--------------------------
run_logratio_host_model <- function(host_col, guild_filter = NA, subset_type = "all") {
  
  dat <- Site_OTU_RA
  if (!is.na(guild_filter)) {
    dat <- dat %>% filter(guild2 == guild_filter)
  }
  
  # Subset by log_ratio direction
  if (subset_type == "positive") dat <- dat %>% filter(log_ratio > 0)
  if (subset_type == "negative") dat <- dat %>% filter(log_ratio < 0)
  
  formula <- as.formula(paste0("log_ratio ~ ", host_col, " + (1|Site/Plot)"))
  model <- lmer(formula, data = dat)
  anov <- Anova(model, test = "F")
  tidy_model <- tidy(model, effects = "fixed") %>% filter(term == host_col)
  R2_vals <- performance::r2(model)
  
  # ANOVA table
  anova_tbl <- anov %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(
      Estimate = tidy_model$estimate,
      Std_Error = tidy_model$std.error,
      R2_marginal = R2_vals$R2_marginal,
      R2_conditional = R2_vals$R2_conditional,
      Subset = subset_type
    ) %>%
    dplyr::select(Subset, Factor, Estimate, Std_Error, Df, Df_res = Df.res, F, P = `Pr(>F)`,
           R2_marginal, R2_conditional)
  
  # Plot label
  label_text <- if (anova_tbl$P[1] > 0.10) {
    "ns"
  } else if (anova_tbl$P[1] < 0.001) {
    paste0("P < 0.001\nR²: ", round(R2_vals$R2_marginal, 2))
  } else {
    paste0("P: ", round(anova_tbl$P[1], 3), "\nR²: ", round(R2_vals$R2_marginal, 2))
  }
  
  plot <- ggplot(dat, aes(x = .data[[host_col]], y = log_ratio)) +
    geom_point(size = 6, alpha = 0.7, color = 'black') +
    geom_smooth(method = "lm", se = FALSE, linewidth = 5, color = 'purple4') +
    labs(
      x = case_when(
        host_col == "freq_ECM" ~ "Frequency of EcM hosts",
        host_col == "freq_AM" ~ "Frequency of AM hosts",
        host_col == "freq_total_myco" ~ "Frequency of all mycorrhizal hosts"
      ),
      y = "log10(RA Hyphae / RA Soil)"
    ) +
    annotate("text",
             x = min(dat[[host_col]], na.rm = TRUE) + 0.05,
             y = max(dat$log_ratio, na.rm = TRUE) * 0.88,
             label = label_text,
             hjust = 0, size = 16, color = "black") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 45, face = "bold"),
      axis.text.y = element_text(size = 45, face = "bold"),
      axis.title.x = element_text(size = 50, face = "bold"),
      axis.title.y = element_text(size = 50, face = "bold"),
      axis.line = element_line(linewidth = 3, colour = "black"),
      plot.title = element_text(size = 40, face = "bold")
    )
  
  return(list(anova = anova_tbl, plot = plot))
}

#--------------------------
# Step 4: Run all combinations
#--------------------------
hosts <- c("freq_ECM", "freq_AM", "freq_total_myco")
guilds <- c("EctoMycorrhizal", "Arbuscular Mycorrhizal", NA)
subsets <- c("all", "positive", "negative")

results <- list()

for (s in subsets) {
  temp <- map2(hosts, guilds, ~ run_logratio_host_model(.x, .y, s))
  results[[s]] <- temp
}

#--------------------------
# Step 5: Combine outputs
#--------------------------
anova_table_all <- bind_rows(
  map(results$all, "anova"),
  map(results$positive, "anova"),
  map(results$negative, "anova")
)

# Tag example plots (subset: all)
plot_a <- results$all[[1]]$plot + labs(tag = "a)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_b <- results$all[[2]]$plot + labs(tag = "b)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_c <- results$all[[3]]$plot + labs(tag = "c)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

final_plot_all <- (plot_a | plot_b) / plot_c



# Tag example plots (subset: all)
plot_a <- results$positive[[1]]$plot + labs(tag = "a)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_b <- results$positive[[2]]$plot + labs(tag = "b)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_c <- results$positive[[3]]$plot + labs(tag = "c)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

final_plot_positive <- (plot_a | plot_b) / plot_c

# Tag example plots (subset: all)
plot_a <- results$negative[[1]]$plot + labs(tag = "a)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_b <- results$negative[[2]]$plot + labs(tag = "b)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))
plot_c <- results$negative[[3]]$plot + labs(tag = "c)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

final_plot_negative <- (plot_a | plot_b) / plot_c

# Print
anova_table_all
final_plot_positive
final_plot_negative

#final_plot


# Save
save(anova_table_all, file = 'HMSC_MER/Output/Processed_Data/Logratio_Host_ANOVA_output_with_subsets.RDS')
ggsave("plots/Logratio_Host_Plot_all_positive_negative.png", final_plot, width = 30, height = 30, dpi = 300)



