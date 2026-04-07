# =============================================================================
# EcM Exploration Type Summary Table
# Input:  Processed_data/Seq_dat/Soil/myco_tax_soil.rds
# Output: Tables/ecm_exploration_type_table.docx
# =============================================================================
library(tidyverse)
library(flextable)
library(officer)

# --- Load data ----------------------------------------------------------------
myco_tax_soil <- readRDS("Processed_data/Seq_dat/Soil/myco_tax_soil.rds") 

# --- Define exploration type order --------------------------------------------
exploration_levels <- c(
  "contact",
  "short-distance_coarse",
  "short-distance_delicate",
  "mat",
  "medium-distance_fringe",
  "medium-distance_smooth",
  "long-distance"
)
unique(myco_tax_soil$exploration_type)
# --- Summarise exploration type per EcM genus ---------------------------------
ecm_exploration <- myco_tax_soil %>%
  filter(str_detect(guild, regex("ectomycorrhizal", ignore_case = TRUE))) %>%
  filter(!is.na(exploration_type)) %>%
  group_by(genus, exploration_type, Ecm_lineage) %>%
  summarise(
    n_OTUs = n_distinct(OTU),
    .groups = "drop"
  ) %>%
  mutate(
    exploration_type = factor(exploration_type, levels = exploration_levels),
    Ecm_lineage = str_remove(Ecm_lineage, "/")
  ) %>%
  arrange(exploration_type, Ecm_lineage, genus) %>%
  select(exploration_type, Ecm_lineage, genus, n_OTUs)

# --- Build flextable ----------------------------------------------------------
ft <- flextable(ecm_exploration) %>%
  
  # Column labels
  set_header_labels(
    exploration_type = "Exploration Type",
  #  Ecm_lineage      = "EcM Lineage",
    genus            = "Genus",
    n_OTUs           = "No. OTUs"
  ) %>%
  
  # Merge repeated exploration type cells for readability
  merge_v(j = "exploration_type") %>%
  valign(j = "exploration_type", valign = "top") %>%
  
  # Style header
  bold(part = "header") %>%
  align(part = "header", align = "center") %>%
  
  # Style body
  italic(j = "genus") %>%
  align(j = "n_OTUs", align = "right", part = "all") %>%
  
  # Borders
  border_inner_h(part = "all", border = fp_border(color = "#CCCCCC", width = 0.5)) %>%
  
  # Font and width
  font(fontname = "Aptos", part = "all") %>%
  fontsize(size = 11, part = "all") %>%
  autofit() %>%
  set_table_properties(layout = "autofit")

ft  # preview in RStudio Viewer

# --- Export to Word -----------------------------------------------------------
dir.create("Tables", showWarnings = FALSE)

doc <- read_docx() %>%
  body_add_par("EcM Exploration Types by Genus", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_add_par("")

doc <- read_docx() %>%
  body_add_par("EcM Exploration Types by Genus", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_add_par("")   # trailing blank line

print(doc, target = "Tables/output/ecm_exploration_type_table.docx")

message("Table saved to Tables/ecm_exploration_type_table.docx")