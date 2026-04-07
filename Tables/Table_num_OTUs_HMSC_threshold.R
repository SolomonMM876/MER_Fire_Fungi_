library(tidyverse)
library(flextable)
library(officer)


# Load your data
myco_tax <- readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_dat <- readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
load('HMSC_MER/results/Beta1.RData')

# Extract components
Beta_estimates <- Beta1$mean
prob_x_greater_0 <- Beta1$pos
prob_x_less_0 <- Beta1$neg

# Get beta coefficient names (excluding Species and Intercept)
beta_coeffs <- colnames(Beta_estimates)[-c(1,2,4)]

# Define thresholds
thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95)

# Initialize results table
summary_table <- data.frame(
  Threshold = numeric(),
  Positive_OTUs = integer(),
  Negative_OTUs = integer(),
  Total_OTUs = integer()
)

# Loop over thresholds
for (threshold in thresholds) {
  sig_pos_all <- data.frame()
  sig_neg_all <- data.frame()
  
  for (beta in beta_coeffs) {
    # Significant positive
    sig_pos <- Beta_estimates %>%
      select(Species, !!beta) %>%
      left_join(prob_x_greater_0 %>% select(Species, !!beta), by = "Species", suffix = c("_mean", "_pos")) %>%
      filter(.data[[paste0(beta, "_pos")]] >= threshold)
    
    # Significant negative
    sig_neg <- Beta_estimates %>%
      select(Species, !!beta) %>%
      left_join(prob_x_less_0 %>% select(Species, !!beta), by = "Species", suffix = c("_mean", "_neg")) %>%
      filter(.data[[paste0(beta, "_neg")]] >= threshold)
    
    sig_pos_all <- bind_rows(sig_pos_all, sig_pos)
    sig_neg_all <- bind_rows(sig_neg_all, sig_neg)
  }
  
  # Unique species counts
  n_pos <- n_distinct(sig_pos_all$Species)
  n_neg <- n_distinct(sig_neg_all$Species)
  n_total <- n_distinct(bind_rows(sig_pos_all, sig_neg_all)$Species)
  
  summary_table <- summary_table %>%
    add_row(
      Threshold = threshold,
      Positive_OTUs = n_pos,
      Negative_OTUs = n_neg,
      Total_OTUs = n_total
    ) 
}
summary_table<-summary_table%>% 
  rename(`Positive OTUs`=Positive_OTUs,`Negative OTUs`= Negative_OTUs, `Total OTUs`= Total_OTUs)

# View result
print(summary_table)



# Build flextable
ft_summary <- summary_table %>%
  flextable() %>%
  theme_booktabs() %>%
  add_header_lines("Table SX.Number of OTUs (operational taxonomic units) at each posterior probability") %>%
  align(align = "center", part = "all") %>%
  set_table_properties(layout = "autofit") 

# View in RStudio Viewer
ft_summary

# Create and save Word document
doc <- read_docx()
doc <- doc %>%
  body_add_flextable(ft_summary)

print(doc, target = "Tables/Output/Table_OTUs_per_threshold.docx")
