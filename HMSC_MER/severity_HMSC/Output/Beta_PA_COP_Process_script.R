library(tidyverse)
library(ggplot2)
library(lme4)
library(car)
library(emmeans)
library(ggsignif)


select<-dplyr::select
# load Beta1 and 2

myco_tax <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds")
myco_dat <- readRDS("Processed_data/Seq_dat/Soil/myco_RA_soil.rds")

# load most up to date PA model
load("HMSC_MER/severity_HMSC/results/Beta1_fire_severity.RData")


############# PA####################

# Beta estimates
Beta_estimates <- Beta1$mean

prob_x_greater_0 <- Beta1$pos
prob_x_less_0 <- Beta1$neg

# Column 1 is "Species"; exclude Species and the Intercept by name
variable <- colnames(Beta_estimates) %>% setdiff(c("Species", "(Intercept)"))

# Get beta coefficient names (excluding Species and Intercept)
beta_coeffs <- variable

threshold <- 0.90

# Create lists to store results
significant_positive <- list()
significant_negative <- list()
# Initialize empty dataframe
significant_OTUs <- data.frame()

# Loop through each beta coefficient
for (beta in beta_coeffs) {
  # Filter for significant positive associations
  sig_pos <- Beta_estimates %>%
    select(Species, all_of(beta)) %>%
    left_join(prob_x_greater_0 %>% select(Species, all_of(beta)), by = "Species", suffix = c("_mean", "_pos")) %>%
    filter(.data[[paste0(beta, "_pos")]] >= threshold)

  # Filter for significant negative associations
  sig_neg <- Beta_estimates %>%
    select(Species, all_of(beta)) %>%
    left_join(prob_x_less_0 %>% select(Species, all_of(beta)), by = "Species", suffix = c("_mean", "_neg")) %>%
    filter(.data[[paste0(beta, "_neg")]] >= threshold)

  # Store in lists
  significant_positive[[beta]] <- sig_pos
  significant_negative[[beta]] <- sig_neg

  # Combine positive and negative associations
  significant_OTUs <- bind_rows(significant_OTUs, sig_pos, sig_neg)
}

# Check results
print(significant_positive) # List of dataframes for significantly positive OTUs
print(significant_negative) # List of dataframes for significantly negative OTUs

# ── Threshold sensitivity summary (Severity only) ───────────────────────
# Counts OTUs with significant positive, negative, or any response to
# Severity at each probability threshold.

thresholds <- c(0.95, 0.90, 0.80, 0.70, 0.60, 0.5)

threshold_summary <- map_dfr(thresholds, function(thr) {
  pos_species <- Beta_estimates %>%
    select(Species, Severity) %>%
    left_join(prob_x_greater_0 %>% select(Species, Severity),
      by = "Species", suffix = c("_mean", "_pos")
    ) %>%
    filter(Severity_pos >= thr) %>%
    pull(Species)

  neg_species <- Beta_estimates %>%
    select(Species, Severity) %>%
    left_join(prob_x_less_0 %>% select(Species, Severity),
      by = "Species", suffix = c("_mean", "_neg")
    ) %>%
    filter(Severity_neg >= thr) %>%
    pull(Species)

  tibble(
    Posterior_Prob = thr,
    Covariate = "Severity",
    N_Positive = length(pos_species),
    N_Negative = length(neg_species),
    N_Total = length(union(pos_species, neg_species))
  )
})

# Print the summary table
print(threshold_summary)


# Bind all significant positive and negative associations into one dataframe
Sig_OTUs_PA_ <- bind_rows(
  lapply(names(significant_positive), function(name) {
    if (nrow(significant_positive[[name]]) > 0) {
      df_subset <- significant_positive[[name]] %>%
        mutate(Covariate = name, Sign = "Positive") %>%
        rename(Mean = ends_with("_mean"), Prob = ends_with("_pos"))
      return(df_subset)
    }
  }),
  lapply(names(significant_negative), function(name) {
    if (nrow(significant_negative[[name]]) > 0) {
      df_subset <- significant_negative[[name]] %>%
        mutate(Covariate = name, Sign = "Negative") %>%
        rename(Mean = ends_with("_mean"), Prob = ends_with("_neg"))
      return(df_subset)
    }
  })
)
# Left join with myco_tax to get genus names
Sig_OTUs_PA <- Sig_OTUs_PA_ %>%
  rename(OTU = Species) %>%
  left_join(myco_tax, by = "OTU") %>%
  mutate(
    OTU_ID = case_when(
      !is.na(genus) ~ genus, # Keep genus if available
      is.na(genus) & !is.na(family) ~ family,
      is.na(genus) & is.na(family) & !is.na(order) ~ order,
      is.na(genus) & is.na(family) & is.na(order) & !is.na(class) ~ class,
      is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) ~ phylum
    ),
    OTU_phylo = paste(OTU, " (", OTU_ID, ")", sep = ""), # Create new OTU_phylo column  )%>%
    Color = case_when(
      Mean > 0 ~ "#2F3D7D", # Blue (positive)
      Mean < 0 ~ "#7B2525" # Red (negative)
    )
  ) %>%
  filter(!str_detect(Covariate, "total_reads")) %>% # Correctly filters out these values
  left_join(myco_dat %>%
    group_by(OTU) %>%
    summarise(total_count = sum(count))) %>%
  mutate(
    OTU_phylo = str_remove(OTU_phylo, "ITSall_"),
    OTU_phylo = fct_reorder(OTU_phylo, total_count, max),
    SignSymbol = case_when(
      Sign == "Positive" ~ "+",
      Sign == "Negative" ~ "–",
      TRUE ~ ""
    )
  )

# --- ADD THIS BLOCK after the existing Sig_OTUs_PA pipeline ---
Sig_OTUs_PA_gen <- Sig_OTUs_PA %>%
  group_by(genus, Covariate) %>%
  mutate(
    SignSymbol = str_c(unique(SignSymbol[SignSymbol != ""]), collapse = "/"),
    Color = case_when(
      SignSymbol == "+"  ~ "#2F3D7D",   # Blue (positive)
      SignSymbol == "–"  ~ "#7B2525",   # Red (negative)
      TRUE               ~ "#A9A9A9"    # Grey (both)
    )
  ) %>%
  ungroup()
# --------------------------------------------------------------



covariate_labels <- c(
  "Severity" = "Fire Severity"
)

# Assuming Sig_OTUs already contains the necessary columns
Beta_Cor_color_PA <- Sig_OTUs_PA %>%
  filter(Covariate == "Severity") %>%
  arrange(guild, OTU_ID) %>%
  mutate(OTU_phylo = factor(OTU_phylo, levels = unique(OTU_phylo))) %>%
  ggplot(aes(x = Covariate, y = OTU_phylo, fill = Color)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = SignSymbol), size = 5, color = "black") +
  scale_fill_identity() +
  scale_x_discrete(labels = covariate_labels) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
    axis.text.y = element_text(face = "italic", size = 16, color = "black"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  labs(y = NULL, x = "Fixed Effects")


print(Beta_Cor_color_PA)

plotly::ggplotly(Beta_Cor_color_PA)

ggsave(filename = "Plots/Beta_PA.png", plot = Beta_Cor_color_PA, dpi = 300, device = "png", width = 8, height = 12)


# ####################################################### ABUNDANCE##########################################
# NOTE: Beta2 (abundance model) is currently not fitted/available.
# The abundance section below is commented out until models[[2]] is refitted.
# When ready, uncomment and run from here.

# Beta_estimates_COP <- Beta2$mean
# prob_x_greater_0_COP <- Beta2$pos
# prob_x_less_0_COP <- Beta2$neg
# ... (loop and Sig_OTUs_Abun_ construction mirrors the PA section above)

# Placeholder so downstream JDSM code still has a Sig_OTUs_Abun_ object
# Sig_OTUs_Abun_ <- Sig_OTUs_PA_ %>% mutate(Covariate = NA_character_)[0, ] # empty tibble

# ggsave(filename = "Plots/Beta_Abundance.png", plot = Beta_Abund, dpi=300, device = "eps", width = 70, height = 50, units = "cm")

################## PA ONLY (Abundance model not yet available) ######################################


# PA <- Sig_OTUs_PA_ %>%
#   rename(OTU = Species) %>%
#   left_join(myco_tax, by = "OTU") %>%
#   mutate(
#     OTU_ID = case_when(
#       !is.na(genus) ~ genus, # Keep genus if available
#       is.na(genus) & !is.na(family) ~ family,
#       is.na(genus) & is.na(family) & !is.na(order) ~ order,
#       is.na(genus) & is.na(family) & is.na(order) & !is.na(class) ~ class,
#       is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) ~ phylum
#     ),
#     SignSymbol = case_when(
#       Sign == "Positive" ~ "+",
#       Sign == "Negative" ~ "–",
#       TRUE ~ ""
#     )
#   ) %>%
#   group_by(OTU_ID, Covariate) %>%
#   summarise(
#     SignSymbol = str_c(unique(SignSymbol[SignSymbol != ""]), collapse = "/"),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Color = case_when(
#       SignSymbol == "+" ~ "#2F3D7D", # Blue (positive)
#       SignSymbol == "–" ~ "#7B2525", # Red (negative)
#       SignSymbol == "+/–" ~ "#A9A9A9" # grey
#     )
#   ) %>%
#   filter(!str_detect(Covariate, "total_reads")) %>% # Correctly filters out these values
#   mutate(Model = "Presence-Absence")



# Abun <- Sig_OTUs_Abun_ %>%
#   rename(OTU = Species) %>%
#   left_join(myco_tax, by = "OTU") %>%
#   mutate(
#     OTU_ID = case_when(
#       !is.na(genus) ~ genus, # Keep genus if available
#       is.na(genus) & !is.na(family) ~ family,
#       is.na(genus) & is.na(family) & !is.na(order) ~ order,
#       is.na(genus) & is.na(family) & is.na(order) & !is.na(class) ~ class,
#       is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) ~ phylum
#     ),
#     SignSymbol = case_when(
#       Sign == "Positive" ~ "+",
#       Sign == "Negative" ~ "–",
#       TRUE ~ ""
#     )
#   ) %>%
#   group_by(OTU_ID, Covariate) %>%
#   summarise(
#     SignSymbol = str_c(unique(SignSymbol[SignSymbol != ""]), collapse = "/"),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Color = case_when(
#       SignSymbol == "+" ~ "#2F3D7D", # Blue (positive)
#       SignSymbol == "–" ~ "#7B2525", # Red (negative)
#       SignSymbol == "+/–" ~ "#A9A9A9" # grey
#     )
#   ) %>%
#   filter(!str_detect(Covariate, "total_reads")) %>% # Correctly filters out these values
#   mutate(Model = "Abundance")

# covariate_order <- names(covariate_labels) # To reorder x-axis


# JDSM <- PA %>%
#   bind_rows(Abun) %>%
#   left_join(
#     myco_dat %>%
#       mutate(
#         OTU_ID = case_when(
#           !is.na(genus) ~ genus, # Keep genus if available
#           is.na(genus) & !is.na(family) ~ family,
#           is.na(genus) & is.na(family) & !is.na(order) ~ order,
#           is.na(genus) & is.na(family) & is.na(order) & !is.na(class) ~ class,
#           is.na(genus) & is.na(family) & is.na(order) & is.na(class) & !is.na(phylum) ~ phylum
#         )
#       ) %>%
#       group_by(OTU_ID) %>%
#       summarise(total_count = sum(count))
#   ) %>%
#   mutate(
#     OTU_ID = fct_reorder(OTU_ID, total_count, max),
#     Covariate = factor(Covariate, levels = covariate_order) # reorder for plotting
#   )


# # Preserve original facet titles
# # JDSM$Model <- factor(JDSM$Model, levels = c("Presence-Absence", "Abundance"), labels = c("Presence-Absence Model", "Abundance Model"))

# # Add facet tags to the Model variable
# # JDSM <- JDSM %>%
# #   mutate(Model = factor(Model,
# #     levels = unique(Model),
# #     labels = c(
# #       "(a)            Presence-Absence Model",
# #       "(b)            Abundance Model"
# #     )
# #   ))

# # Create the base plot
# JDSM_plot <- JDSM %>%
#   ggplot(aes(x = Covariate, y = OTU_ID, fill = Color)) +
#   geom_tile(color = "black", linewidth = 0.5) +
#   geom_text(aes(label = SignSymbol), size = 12, color = "black") +
#   scale_fill_identity() +
#   scale_x_discrete(labels = covariate_labels) +
#   facet_wrap(~Model, nrow = 1, strip.position = "top") +
#   theme_minimal() +
#   theme(
#     strip.text = element_text(size = 36),
#     strip.placement = "outside",
#     panel.spacing = unit(3, "lines"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = .25, size = 32, color = "black"),
#     axis.text.y = element_text(face = "italic", size = 34, color = "black"),
#     axis.title.x = element_blank(),
#     axis.ticks.x = element_blank()
#   ) +
#   labs(y = NULL, x = "Fixed Effects")
# JDSM_plot

# ggsave(filename = "plots/JDSM_ABUN_PA.png", plot = JDSM_plot, dpi = 300, device = "png", width = 70, height = 50, units = "cm")
