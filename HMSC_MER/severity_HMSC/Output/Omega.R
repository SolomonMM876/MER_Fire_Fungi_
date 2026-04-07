library(tidyverse)
library(circlize)

#load data
load(file='HMSC_MER/results/Omega.RData')

myco_dat<-readRDS('Processed_data/Seq_dat/Soil/myco_RA_soil.rds')
myco_tax<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')



# 1. Threshold for significance (adjustable)
threshold <- 0.95
# 2. Extract matrices
omega_mean <- OmegaCor1_Site$mean
omega_support <- OmegaCor1_Site$support
# === STEP 1: Create a mask for significant correlations ===
sig_mat <- (omega_support > threshold | omega_support < (1 - threshold))













# === STEP 2: Extract OTU correlations into long-format data ===
# This gets all OTU pairs that are significantly correlated
sig_df <- which(sig_mat, arr.ind = TRUE) %>%
  as.data.frame() %>%
  mutate(
    OTU1 = rownames(omega_mean)[row],
    OTU2 = colnames(omega_mean)[col],
    correlation = mapply(function(i, j) omega_mean[i, j], row, col),
    direction = ifelse(correlation > 0, "positive", "negative")
  ) %>%
  select(OTU1, OTU2, correlation, direction)

# === STEP 3: Merge in Taxonomy for Labeling and Color Grouping ===
# You use genus > family > order > class > phylum hierarchy
grouping_var <- "guild2"

otu_colors <- myco_tax %>%
  mutate(
    OTU_ID = case_when(
      !is.na(genus) ~ genus,
      is.na(genus) & !is.na(family) ~ family,
      is.na(genus) & is.na(family) & !is.na(order) ~ order,
      is.na(genus) & is.na(family) & is.na(order) & !is.na(class) ~ class,
      TRUE ~ phylum
    ),
    OTU_phylo = paste(OTU, " (", OTU_ID, ")", sep = "")
  ) %>%
  select(OTU, OTU_phylo, !!sym(grouping_var)) %>%
  distinct()

# Generate color palette
phylo_groups <- unique(otu_colors[[grouping_var]])
color_palette <- setNames(rainbow(length(phylo_groups)), phylo_groups)
otu_colors$color <- color_palette[otu_colors[[grouping_var]]]

# Vector of OTU_phylo → color
col_vector <- setNames(otu_colors$color, otu_colors$OTU_phylo)

# === STEP 4: Filter and Annotate Significant Correlations ===
sig_df <- sig_df %>%
  filter(OTU1 %in% otu_colors$OTU & OTU2 %in% otu_colors$OTU) %>%
  left_join(otu_colors, by = c("OTU1" = "OTU")) %>%
  rename(OTU1_phylo = OTU_phylo) %>%
  left_join(otu_colors, by = c("OTU2" = "OTU")) %>%
  rename(OTU2_phylo = OTU_phylo)

# Use phylogenetic names for plotting
sig_plot_df <- sig_df %>%
  filter(!OTU1_phylo==OTU2_phylo) %>% #remove auto correlations
  select(OTU1,OTU2,OTU1_phylo, OTU2_phylo, correlation, direction)


# Join taxonomy to both OTUs
sig_plot_df <- sig_plot_df %>%
  left_join(myco_tax%>% select(OTU,genus,family,phylum), by = c("OTU1" = "OTU")) %>%
  rename(Genus1 = genus, Family1=family,Phylum1 = phylum) %>%
  left_join(myco_tax%>% select(OTU,genus,family,phylum), by = c("OTU2" = "OTU")) %>%
  rename(Genus2 = genus, Family2=family,Phylum2 = phylum) 
  
# OPTIONAL: Limit to top 100 for legibility
# sig_plot_df_filt <- sig_plot_df %>%
#   arrange(desc(abs(correlation))) %>%
#   slice(1:100)

# Replace Genus with Phylum in group_by() below if needed
grouped_genus <- sig_plot_df %>%
  filter(!is.na(Genus1), !is.na(Genus2), Genus1 != Genus2) %>%
  group_by(Genus1, Genus2) %>%
  summarise(
  total_strength = sum(abs(correlation)),
mean_correlation = mean(correlation),
.groups = "drop")

# Create matrix of interaction strengths
link_mat_genus <- grouped_genus %>%select(-mean_correlation) %>% 
  pivot_wider(names_from = Genus2, values_from = total_strength, values_fill = 0) %>%
  column_to_rownames("Genus1") %>%
  as.matrix()

all_genera <- union(rownames(link_mat_genus), colnames(link_mat_genus))
set.seed(123)
col_vector <- setNames(rainbow(length(all_genera)), all_genera)
#plot
circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(
  link_mat_genus,
  grid.col = col_vector,
  transparency = 0.4,
  directional = 1,
  annotationTrack = c("name", "grid"),  # Adds group labels
  preAllocateTracks = list(track.height = 0.05)
)

#Now by phylum
grouped_phylum <- sig_plot_df %>%
  filter(!is.na(Phylum1), !is.na(Phylum2), Phylum1 != Phylum2) %>%
  group_by(Phylum1, Phylum2) %>%
  summarise(
    total_strength = sum(abs(correlation)),
    mean_correlation = mean(correlation),
    .groups = "drop")

# Create matrix of interaction strengths
link_mat_phylum <- grouped_phylum %>%select(-mean_correlation) %>% 
  pivot_wider(names_from = Phylum2, values_from = total_strength, values_fill = 0) %>%
  column_to_rownames("Phylum1") %>%
  as.matrix()

all_phylum <- union(rownames(link_mat_phylum), colnames(link_mat_phylum))
set.seed(123)
col_vector <- setNames(rainbow(length(all_phylum)), all_phylum)

circos.clear()
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(
  link_mat_phylum,
  grid.col = col_vector,
  transparency = 0.4,
  directional = 2,
  annotationTrack = c("name", "grid"),  # Adds group labels
  preAllocateTracks = list(track.height = 0.05)
)


# PDF
pdf("Plots/ChordDiagram_Genus_or_Phylum.pdf", width = 12, height = 12)
chordDiagram(
  link_mat,
  grid.col = col_vector,
  transparency = 0.4,
  directional = 1,
  annotationTrack = c("name", "grid"),
  preAllocateTracks = list(track.height = 0.05)
)
dev.off()

# PNG
png("Plots/ChordDiagram_Genus.png", width = 1000, height = 1000, res = 150)
chordDiagram(
  link_mat,
  grid.col = col_vector,
  transparency = 0.4,
  directional = 0,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.05)
)
dev.off()

#################################################
#correlation by species

# Step 1: Define threshold and create significance mask
threshold <- 0.95
sig_mat <- (omega_support > threshold | omega_support < (1 - threshold))

# Step 2: Keep only significant values in omega_mean
omega_sig <- omega_mean
omega_sig[!sig_mat] <- NA  # Set non-significant values to NA
# Step 3: Remove self-correlations
diag(omega_sig) <- NA
# Remove rows with all NA
omega_sig <- omega_sig[rowSums(!is.na(omega_sig)) > 0, ]
# Remove columns with all NA
omega_sig <- omega_sig[, colSums(!is.na(omega_sig)) > 0]

otus_in_mat<-colnames(omega_sig)
################

tax_df <- myco_tax %>%
  filter(OTU %in% otus_in_mat) %>%
  mutate(
  genus_ID = make.unique(paste(phylum, genus, sep = "_"))
  ) %>%
  arrange(factor(phylum))

otu_order <- tax_df$OTU
omega_sig_ordered <- omega_sig[otu_order, otu_order]

# === 3. Create color vector for OTUs by phylum ===
phylum_colors <- c(
  "Glomeromycota" = "blue",
  "Basidiomycota" = "green",
  "Ascomycota" = "red"
)

col_vector <- tax_df %>%
  mutate(color = phylum_colors[phylum]) %>%
  pull(color)  # names are OTUs

# === 4. Set matrix row/col names to phylum only (for labeling) ===
phylum_labels <- tax_df$phylum
names(phylum_labels) <- tax_df$OTU

names(col_vector) <- phylum_labels[names(col_vector)]  # color vector must match new names


# === 5. Compute gap vector (breaks only between phyla) ===
gap_vec <- tax_df %>%
  mutate(phylum_shift = phylum != lag(phylum, default = first(phylum))) %>%
  mutate(gap = if_else(phylum_shift, 5, 0)) %>%
  pull(gap)
gap_vec[length(gap_vec)] <- 5


# === 6. Plot the chord diagram ===
circos.clear()
circos.par(gap.after = gap_vec)

chordDiagram(
  omega_sig_ordered,
  grid.col = col_vector,
  transparency = 0.4,
  directional = 2,
  annotationTrack = "name",  # shows sector labels
  preAllocateTracks = list(track.height = 0.05)
)

