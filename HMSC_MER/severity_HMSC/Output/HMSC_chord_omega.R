################################################################################
# HMSC Associations Chord Diagram
# Single model, pre-computed OmegaCor1_plot (mean + support as data frames)
# - Static PNG/PDF: one combined plot (positive = blue, negative = red edges)
# - Interactive HTML: one combined chord diagram
# - OTU labels: genus > family > order > class > phylum
################################################################################

library(circlize)
library(tidyverse)
library(htmlwidgets)
library(chorddiag)
library(base64enc)

################################################################################
# 1. SETTINGS — edit these
################################################################################

support_threshold <- 0.95
plot_title        <- "Fire Severity — Fungal Associations"
output_dir        <- "HMSC_MER/severity_HMSC/results/plots"
omega_path        <- "HMSC_MER/severity_HMSC/results/omega1_fire_severity.RData"

# Path to myco_tax RDS file
taxonomy_path <- "Processed_data/Seq_dat/Soil/myco_tax_soil.rds"

################################################################################
# 2. LOAD PRE-COMPUTED OMEGA
################################################################################

cat("Loading OmegaCor1_plot from:", omega_path, "\n")
load(omega_path)   # loads OmegaCor1_plot — list with $mean and $support

mean_mat    <- as.matrix(OmegaCor1_plot$mean)
support_mat <- as.matrix(OmegaCor1_plot$support)

stopifnot(nrow(mean_mat) == ncol(mean_mat))
stopifnot(all(dim(mean_mat) == dim(support_mat)))

################################################################################
# 3. LOAD TAXONOMY AND BUILD OTU LABELS
################################################################################

cat("Loading taxonomy from:", taxonomy_path, "\n")
myco_tax <- readRDS(taxonomy_path)

# Standardise column names to lower-case for safety
colnames(myco_tax) <- tolower(colnames(myco_tax))

# Ensure all required rank columns exist (fill missing ranks with NA)
for (.col in c("genus", "family", "order", "class", "phylum")) {
  if (!.col %in% colnames(myco_tax)) myco_tax[[.col]] <- NA_character_
}

# Build display labels using the finest available taxonomic rank
myco_tax <- myco_tax %>%
  mutate(
    OTU_ID = case_when(
      !is.na(genus)                                                        ~ genus,
      is.na(genus) & !is.na(family)                                        ~ family,
      is.na(genus) & is.na(family) & !is.na(order)                        ~ order,
      is.na(genus) & is.na(family) & is.na(order) & !is.na(class)        ~ class,
      is.na(genus) & is.na(family) & is.na(order) & is.na(class) &
        !is.na(phylum)                                                     ~ phylum,
      TRUE                                                                  ~ otu   # fallback to raw OTU ID
    ),
    OTU_phylo = paste(otu, " (", OTU_ID, ")", sep = "")
  )

# Named vector: raw OTU ID -> display label (e.g. "OTU1 (Russula)")
otu_label_map <- setNames(myco_tax$OTU_phylo, myco_tax$otu)

# Warn about any OTUs in the omega matrix that are absent from the taxonomy
otu_names <- rownames(mean_mat)
missing_tax <- setdiff(otu_names, names(otu_label_map))
if (length(missing_tax) > 0) {
  cat("  Warning:", length(missing_tax),
      "OTU(s) not found in myco_tax — raw names will be used:\n")
  cat("  ", paste(missing_tax, collapse = ", "), "\n")
  otu_label_map[missing_tax] <- missing_tax   # identity fallback
}

cat("  Taxonomy labels built for", length(otu_label_map), "OTUs.\n")

# Helper: apply label map to both dimensions of a matrix
relabel_mat <- function(m, label_map) {
  rn <- label_map[rownames(m)]
  cn <- label_map[colnames(m)]
  rn[is.na(rn)] <- rownames(m)[is.na(rn)]   # safety fallback
  cn[is.na(cn)] <- colnames(m)[is.na(cn)]
  rownames(m) <- rn
  colnames(m) <- cn
  m
}
################################################################################
# 4a FILTER TO AM TAXA ONLY
################################################################################

# Identify AM OTUs from myco_tax using the guild2 column
if (!"guild2" %in% colnames(myco_tax)) {
  stop("Column 'guild2' not found in myco_tax. Check column names with colnames(myco_tax).")
}

am_otus <- myco_tax %>%
  filter(guild2 == "Arbuscular Mycorrhizal") %>%
  pull(otu)

cat("AM OTUs in taxonomy:", length(am_otus), "\n")

# Restrict to AM OTUs that are actually present in the omega matrix
am_in_omega <- intersect(am_otus, rownames(mean_mat))
cat("AM OTUs present in omega matrix:", length(am_in_omega), "\n")

if (length(am_in_omega) < 2) {
  stop("Fewer than 2 AM OTUs found in the omega matrix — cannot draw associations.")
}

# Subset both matrices to AM-only rows and columns
mean_mat    <- mean_mat[am_in_omega, am_in_omega]
support_mat <- support_mat[am_in_omega, am_in_omega]
################################################################################
# 4. FILTER TO ECTOMYCORRHIZAL TAXA ONLY
################################################################################

# Identify EcM OTUs from myco_tax using the guild2 column
if (!"guild2" %in% colnames(myco_tax)) {
  stop("Column 'guild2' not found in myco_tax. Check column names with colnames(myco_tax).")
}

ecm_otus <- myco_tax %>%
  filter(guild2 == "Ectomycorrhizal") %>%
  pull(otu)

cat("EcM OTUs in taxonomy:", length(ecm_otus), "\n")

# Restrict to EcM OTUs that are actually present in the omega matrix
ecm_in_omega <- intersect(ecm_otus, rownames(mean_mat))
cat("EcM OTUs present in omega matrix:", length(ecm_in_omega), "\n")

if (length(ecm_in_omega) < 2) {
  stop("Fewer than 2 EcM OTUs found in the omega matrix — cannot draw associations.")
}

# Subset both matrices to EcM-only rows and columns
mean_mat    <- mean_mat[ecm_in_omega, ecm_in_omega]
support_mat <- support_mat[ecm_in_omega, ecm_in_omega]

################################################################################
# 5. COMPUTE SIGNIFICANT ASSOCIATION MATRIX
################################################################################

toPlot <- ((support_mat > support_threshold) +
             (support_mat < (1 - support_threshold)) > 0) * mean_mat
diag(toPlot) <- 0

n_sig <- sum(toPlot != 0) / 2
cat("Significant EcM-EcM associations at threshold", support_threshold, ":", n_sig, "\n")

if (n_sig == 0) stop("No significant EcM-EcM associations found at threshold ", support_threshold)

# Drop taxa with no significant associations
has_link <- rowSums(toPlot != 0) > 0
toPlot   <- toPlot[has_link, has_link]
cat("EcM taxa retained after filtering:", nrow(toPlot), "\n")

# Apply taxonomic labels
toPlot <- relabel_mat(toPlot, otu_label_map)

################################################################################
# 6. COLOR PALETTE
################################################################################

n_taxa     <- nrow(toPlot)
base_cols  <- c("#4e79a7","#f28e2b","#59a14f","#e15759","#76b7b2",
                "#edc948","#b07aa1","#ff9da7","#9c755f","#bab0ac",
                "#499894","#86bcb6")
sector_cols <- setNames(
  colorRampPalette(base_cols)(n_taxa),
  rownames(toPlot)
)

################################################################################
# 7. DRAW COMBINED CHORD DIAGRAM (positive = blue, negative = red edges)
################################################################################

draw_chord_combined <- function(mat, title_text) {
  
  if (!any(mat != 0)) {
    plot.new()
    mtext(title_text, side = 3, line = 0.5, cex = 0.9, font = 2)
    text(0.5, 0.5, "(No significant associations)", cex = 1)
    return(invisible(NULL))
  }
  
  circos.par(canvas.xlim = c(-1.12, 1.12), canvas.ylim = c(-1.12, 1.12))
  
  chordDiagram(
    mat,
    col              = function(x) ifelse(x < 0, "#cc3333", "#3366cc"),
    grid.col         = sector_cols[rownames(mat)],
    annotationTrack  = "grid",
    preAllocateTracks = list(
      track.height = max(strwidth(rownames(mat))) * 0.55
    )
  )
  
  mtext(title_text, side = 3, line = 1.5, cex = 0.95, font = 2)
  
  circos.track(
    track.index = 1,
    track.margin = c(0.06, 0),
    panel.fun = function(x, y) {
      circos.text(
        CELL_META$xcenter, CELL_META$ylim[1],
        CELL_META$sector.index,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.2), cex = 0.65
      )
    },
    bg.border = NA
  )
  
  circos.clear()
}

################################################################################
# 8. SAVE STATIC FIGURE (PNG + PDF) — single combined panel
################################################################################

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

png_path <- file.path(output_dir, "hmsc_chord_omega.png")
pdf_path <- file.path(output_dir, "hmsc_chord_omega.pdf")

# PNG
png(png_path, width = 7, height = 7, units = "in", res = 300)
par(mar = c(1, 1, 3, 1))
draw_chord_combined(toPlot, plot_title)
dev.off()
cat("Saved PNG:", png_path, "\n")

# PDF
pdf(pdf_path, width = 7, height = 7)
par(mar = c(1, 1, 3, 1))
draw_chord_combined(toPlot, plot_title)
dev.off()
cat("Saved PDF:", pdf_path, "\n")

################################################################################
# 9. INTERACTIVE HTML — single combined chord diagram
################################################################################

# Build a square symmetric matrix using absolute values for chorddiag display;
# chord colour is determined per-edge by sign at render time via JS tooltip,
# but visually we encode sign via a two-colour scheme below.
#
# Strategy: split into pos/neg matrices, give each a distinct chord colour,
# then embed both as separate overlapping widgets inside one HTML page with
# a shared legend — OR use a single widget with a blended colour.
#
# Here we use a single widget and encode:
#   positive edges -> blue   (#3366cc)
#   negative edges -> red    (#cc3333)
# chorddiag doesn't support per-chord colours natively, so we render TWO
# overlapping semi-transparent widgets inside a flex container.

make_chorddiag_widget <- function(mat, sign_type) {
  # sign_type: "positive" or "negative"
  sub <- mat
  if (sign_type == "positive") sub[sub < 0] <- 0
  if (sign_type == "negative") sub[sub > 0] <- 0
  
  # Drop taxa with no links for this sign
  keep <- rowSums(sub != 0) > 0
  sub  <- sub[keep, keep, drop = FALSE]
  
  if (nrow(sub) == 0 || !any(sub != 0)) return(NULL)
  
  # Symmetric absolute-value matrix for display
  sq <- abs(sub)
  sq <- pmax(sq, t(sq))
  
  gcols <- sector_cols[rownames(sq)]
  gcols[is.na(gcols)] <- "#aaaaaa"
  
  chord_col <- if (sign_type == "positive") "#3366cc" else "#cc3333"
  
  chorddiag::chorddiag(
    sq,
    groupColors           = unname(gcols),
    chordedgeColor        = chord_col,
    groupnamePadding      = 30,
    groupnameFontsize     = 11,
    showTicks             = FALSE,
    margin                = 120,
    tooltipGroupConnector = " \u2194 ",
    tooltipUnit           = "",
    precision             = 3,
    fadeLevel             = 0.3
  )
}

widget_to_iframe <- function(w, height_px = 600) {
  if (is.null(w)) {
    return('<p style="color:#999;font-style:italic;text-align:center;">(None)</p>')
  }
  tmp <- tempfile(fileext = ".html")
  on.exit(unlink(tmp))
  htmlwidgets::saveWidget(w, file = tmp, selfcontained = TRUE)
  b64 <- base64enc::base64encode(
    charToRaw(paste(readLines(tmp, warn = FALSE), collapse = "\n"))
  )
  paste0(
    '<iframe src="data:text/html;base64,', b64,
    '" style="width:100%;height:', height_px, 'px;border:none;"></iframe>'
  )
}

# Build positive and negative widgets
w_pos <- make_chorddiag_widget(toPlot, "positive")
w_neg <- make_chorddiag_widget(toPlot, "negative")

html_path <- file.path(output_dir, "hmsc_chord_omega_interactive.html")

html <- paste(c(
  '<!DOCTYPE html><html><head><meta charset="utf-8">',
  paste0('<title>', plot_title, '</title>'),
  '<style>',
  'body { max-width: 1400px; margin: 0 auto; padding: 20px;',
  '       font-family: Arial, sans-serif; background: #fff; }',
  'h1 { font-size: 1.4em; }',
  '.legend { display: flex; gap: 24px; margin: 10px 0 20px; align-items: center; }',
  '.legend-item { display: flex; align-items: center; gap: 8px; font-size: 0.95em; }',
  '.swatch { width: 28px; height: 10px; border-radius: 3px; display: inline-block; }',
  '.row { display: flex; gap: 10px; flex-wrap: wrap; justify-content: center; }',
  '.col { flex: 1; min-width: 480px; }',
  '.col h2 { text-align: center; font-size: 1.05em; margin-bottom: 4px; }',
  '</style></head><body>',
  paste0('<h1>', plot_title, '</h1>'),
  paste0('<p>Posterior support threshold: <strong>', support_threshold, '</strong></p>'),
  '<div class="legend">',
  '  <div class="legend-item">',
  '    <span class="swatch" style="background:#3366cc;"></span>',
  '    <span>Positive associations</span>',
  '  </div>',
  '  <div class="legend-item">',
  '    <span class="swatch" style="background:#cc3333;"></span>',
  '    <span>Negative associations</span>',
  '  </div>',
  '</div>',
  '<div class="row">',
  # Positive panel
  '<div class="col">',
  '<h2 style="color:#1a3a6b;">&#x2b; Positive Associations</h2>',
  widget_to_iframe(w_pos, height_px = 620),
  '</div>',
  # Negative panel
  '<div class="col">',
  '<h2 style="color:#8b1a1a;">&#x2212; Negative Associations</h2>',
  widget_to_iframe(w_neg, height_px = 620),
  '</div>',
  '</div>',
  '</body></html>'
), collapse = "\n")

writeLines(html, con = html_path)
cat("Saved interactive HTML:", html_path, "\n")

cat("\nDone. Outputs in:", output_dir, "\n")
cat("Legend: Blue edges = positive, Red edges = negative\n")
cat("OTU labels: genus > family > order > class > phylum\n")

