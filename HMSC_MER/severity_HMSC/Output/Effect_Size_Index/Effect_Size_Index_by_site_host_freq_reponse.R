library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(car)
library(performance)
library(ggpubr)
library(emmeans)
library(broom.mixed)



#read in soil comm data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')

#read in hyph comm data
myco_dat_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_RA_Hyph.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')

#load most up to date PA model
load('HMSC_MER/results/Beta1.RData')
load('HMSC_MER/results/Beta1_host_freq.RData')


#unburnt host freq
host_freq_unburnt<-readRDS('Processed_data/metadata/veg/myco_freq_Site_unburnt_plots.Rdata')

#plot level data
myco_freq_final_sum_plot<-readRDS('Processed_data/metadata/veg/myco_freq_plot_level_df.Rdata')

All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% dplyr::select(Site,Plot,Severity)

#First lets just use PA HMSC
Beta_estimates_PA<-Beta1$mean

Beta_estimates_PA<-Beta_estimates_PA %>% 
  dplyr::select(OTU=Species,Beta=Severity)



# Define hosts and their corresponding mycorrhizal types (or NA if no filtering)
hosts <- c("freq_ECM", "freq_AM","freq_total_myco")
types  <- c("EctoMycorrhizal", "Arbuscular Mycorrhizal", NA)

# Fit models and collect ANOVA tables

model_outputs <- purrr::map2(hosts, types, function(host, type) {
  
  Site_OTU <- myco_dat_soil
  if (!is.na(type)) {
    relevant_otus <- myco_tax_soil %>%
      filter(guild2 == type) %>%
      pull(OTU)
    
    Site_OTU <- Site_OTU %>%
      filter(OTU %in% relevant_otus) 
  }
  
  Site_OTU <- Site_OTU %>%
    distinct(Site, Plot, OTU, .keep_all = TRUE) %>%
    left_join(Beta_estimates_PA, by = "OTU") %>%
    left_join(myco_freq_final_sum_plot)%>% 
    left_join(myco_tax_soil)
  
  # Correct formula
  formula <- as.formula(paste0("Beta ~ ", host, " + (1|Site/Plot)"))
  model <- lmer(formula, data = Site_OTU)
  anov <- Anova(model, test.statistic="F")
  tidy_model <- tidy(model, effects = "fixed")%>% 
    filter(term == host) 
  
  # Compute R²
  R2_marg <- performance::r2(model)$R2_marginal
  R2_cond <- performance::r2(model)$R2_conditional
  # Return ANOVA and plot
  anova_tbl <- anov %>%
    as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(
      Estimate = tidy_model$estimate,
      Std_Error = tidy_model$std.error,
      R2_marginal = R2_marg,
      R2_conditional= R2_cond
    ) %>%
    dplyr::select(Factor,Estimate,Std_Error, Df,Df_res=Df.res, F, P = `Pr(>F)`,R2_marginal,R2_conditional)
  
  # Main effect result
  p_val <- anov$`Pr(>F)`[1]
  # Label for plot
  label_text <- if (anova_tbl$P[1] > 0.10) {
    "ns"
  } else if (anova_tbl$P[1] < 0.001) {
    paste0("P < 0.001\nR²: ", round(R2_marg, 2))
  } else {
    paste0("P: ", round(anova_tbl$P[1], 3), "\nR²: ", round(R2_marg, 2))
  }
  
  # Define y-axis label separately using if/else
  y_label <- if (host == "freq_ECM") {
    expression(paste("Effect size of EcM taxa (", beta[fire], ")"))
  } else if (host == "freq_AM") {
    expression(paste("Effect size of AM taxa (", beta[fire], ")"))
  } else if (host == "freq_total_myco") {
    expression(paste("Effect size of all mycorrhizal taxa (", beta[fire], ")"))
  }
  
  # Main plot
  plot_main <- ggplot(Site_OTU, aes(x = .data[[host]], y = Beta)) +
    geom_point(size = 6,  alpha=0.7,
      #aes(color=family),
               color='black'
               ) +
    geom_smooth(method = 'lm', se = FALSE, linewidth = 5, color = 'purple4') +
    labs(
      y = y_label,
         x = case_when(
      host == "freq_ECM" ~ "Proportion of EcM hosts",
      host == "freq_AM" ~ "Proportion of AM hosts",
      host == "freq_total_myco" ~ "Proportion of all mycorrhizal hosts")) +
    annotate("text",
             x = 0.1,
             y = max(Site_OTU$Beta, na.rm = TRUE) * 0.88,
             label = label_text,
             hjust = 0, size = 16, color = "black") +
    theme_classic() +
    theme(
      #legend.position = "none",
      axis.text.x = element_text(hjust = 0.5, size = 45, face = "bold"),
      axis.text.y = element_text(size = 45, face = "bold"),
      axis.title.x = element_text(size = 50, face = "bold"),
      axis.title.y = element_text(size = 50, face = "bold"),
      axis.line = element_line(linewidth = 3, colour = "black")
    )
  
  return(list(
    anova = anova_tbl,
    plot_main = plot_main,
    reps=Site_OTU %>%
      mutate(host= paste0( host),
             n=n()) %>% 
      distinct(host,n)
  ))
})

# Combine all ANOVA tables
anova_table <- bind_rows(purrr::map(model_outputs, "anova"))
# Combine all ANOVA tables
rep_table <- bind_rows(purrr::map(model_outputs, "reps"))

# Extract plots
plots <- purrr::map(model_outputs, "plot_main")

library(patchwork)
library(cowplot)

# Tag and style each subplot individually
plot_a <- plots[[1]] +
  labs(tag = "a)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

plot_b <- plots[[2]] +
  labs(tag = "b)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

plot_c <- plots[[3]] +
  labs(tag = "c)") +
  theme(plot.tag = element_text(size = 65, face = "bold"))

# Arrange using patchwork: (a | b) / c
combined_patchwork <- (plot_a | plot_b) / plot_c +
  plot_layout(heights = c(1, 1), widths = c(1, 1))

final_plot<-combined_patchwork
# Add shared y-axis label using cowplot
# final_plot <- ggdraw() +
#   draw_label("Effect Size", x = 0.03, angle = 90, size = 50, fontface = "bold")+
#   draw_plot(combined_patchwork, x = 0.05, width = 0.95)

# Print final plot
final_plot

anova_table


save(anova_table, file='HMSC_MER/Output/Processed_Data/Myco_veg_host_ANOVA_output.RDS')

ggsave("plots/Host_freq_Effect_Size_Plot.png", final_plot, width = 30, height = 30, dpi = 300)
