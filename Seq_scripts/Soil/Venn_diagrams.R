library(tidyverse)
library(vegan)
library(ggrepel)

# Load data
wide_myco_soil <- readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)


# Load data
wide_myco_hyph<-readRDS('Processed_data/Seq_dat/Hyph/wide_myco.rds')
myco_tax_Hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)

site_list <- unique(wide_myco_hyph$Site)

# Color palette for plotting
burn_colors <- c("B" = "darkred", "U" = "orange")





# Load libraries
library(tidyverse)
library(ggvenn)

# Load data
wide_myco_soil <- readRDS('Processed_data/Seq_dat/Soil/wide_myco.rds')
wide_myco_hyph <- readRDS('Processed_data/Seq_dat/Hyph/wide_myco.rds')
All_Sites <- readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>%
  select(Site, Plot, Fire_Treatment)

# Merge metadata
wide_myco_soil <- All_Sites %>%
  left_join(wide_myco_soil) %>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

wide_myco_hyph <- All_Sites %>%
  left_join(wide_myco_hyph) %>%
  filter(!if_all(starts_with("ITSall"), ~ . == 0))

### Function to get OTU sets by treatment
get_otu_sets <- function(df, treatment_col = "Fire_Treatment") {
  otu_cols <- df %>% select(starts_with("ITSall")) %>% colnames()
  list(
    Unburnt = df %>% filter(!!sym(treatment_col) == "U") %>%
      select(all_of(otu_cols)) %>%
      summarise(across(everything(), ~ any(. > 0))) %>%
      pivot_longer(everything(), names_to = "OTU", values_to = "present") %>%
      filter(present) %>% pull(OTU),
    
    Burnt   = df %>% filter(!!sym(treatment_col) == "B") %>%
      select(all_of(otu_cols)) %>%
      summarise(across(everything(), ~ any(. > 0))) %>%
      pivot_longer(everything(), names_to = "OTU", values_to = "present") %>%
      filter(present) %>% pull(OTU)
  )
}

# Get OTU sets
otu_sets_soil <- get_otu_sets(wide_myco_soil) 

otu_sets_hyph <- get_otu_sets(wide_myco_hyph)

library(patchwork)

# Soil Venn diagram
soil_venn <-  ggvenn(
  otu_sets_soil,
  fill_color = c( "#cccccc", "#2F4F4F"),
  show_percentage = TRUE,
  digits=0,
  set_name_size = 6,
  stroke_size = 2,
  #@auto_scale = TRUE,
  text_size = 7
) +
  labs(
    tag= "b)")+
  theme(    plot.tag = element_text(size = 16, face = "bold"),
)

soil_venn
  

# Hyph Venn diagram
hyph_venn <- ggvenn(
  otu_sets_hyph,
  fill_color = c( "#5497B6","#9F2121"),
  show_percentage = TRUE,
  set_name_size = 8,
  stroke_size = 2,
  text_size = 7
) +
  labs(
    tag= "b)")+
  theme(plot.title = element_text(size = 16, face = "bold"))

# Combine with patchwork
combined_plot <- soil_venn + hyph_venn

# Display the plot
print(combined_plot)

# Save the plot
ggsave("plots/venn_diagrams.png", combined_plot, width = 15, height = 11, dpi = 300)

# Extract OTU column names
otu_cols_soil <- wide_myco_soil %>% select(starts_with("ITSall")) %>% colnames()
otu_cols_hyph <- wide_myco_hyph %>% select(starts_with("ITSall")) %>% colnames()

# Helper function to get present OTUs
present_otus <- function(df, treatment, otu_cols) {
  df %>%
    filter(Fire_Treatment == treatment) %>%
    select(all_of(otu_cols)) %>%
    summarise(across(everything(), ~ any(. > 0))) %>%
    pivot_longer(everything(), names_to = "OTU", values_to = "present") %>%
    filter(present) %>%
    pull(OTU)
}

# Define each group
venn_sets <- list(
  Soil_Burnt   = present_otus(wide_myco_soil, "B", otu_cols_soil),
  Hyph_Burnt   = present_otus(wide_myco_hyph, "B", otu_cols_hyph),
  Hyph_Unburnt = present_otus(wide_myco_hyph, "U", otu_cols_hyph),
  Soil_Unburnt = present_otus(wide_myco_soil, "U", otu_cols_soil)
)

# Plot 4-set Venn diagram
ggvenn(venn_sets,label = "count") +
  ggtitle("OTU Overlap: Soil vs Hyph & Burnt vs Unburnt") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
