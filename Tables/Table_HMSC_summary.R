library(tidyverse)
library(flextable)
library(officer)

# -------------------------------
# 1. Read and clean input tables
# -------------------------------

MER_fire_host_HMSC <- readRDS("HMSC_MER/results/model_eval_PA_host.rds")

# -------------------------------
# 5. Build flextable with multi-headers
# -------------------------------
ft_model_eval <- flextable(MER_fire_host_HMSC) %>%
  theme_booktabs() %>%
  add_header_lines("Supp Table X. Evaluation metrics for HMSC model fit to soil mycorrhizal fungal communities. Shown are variance partitioning results (VP), discrimination ability (AUC), explanatory power (Tjur’s R²), accuracy (RMSE), and scale reduction factors for Beta parameters across models.") %>%
  theme_booktabs() %>%
  merge_h(part = "header") %>%
  align(align = "center", part = "header") %>%
  autofit() %>%
  align(j = 1, align = "left", part = "all") # keep Metric left-aligned

ft_model_eval

# -------------------------------
# 6. Export to Word
# -------------------------------
doc <- read_docx() %>%
  body_add_flextable(ft_model_eval)

print(doc, target = "Tables/Output/HMSC_model_eval_Table.docx")
