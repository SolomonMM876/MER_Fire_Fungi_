library(tidyverse)
library(patchwork)

# Source all files that start with 'Effect_size'
files <- list.files(
  path = "HMSC_MER/severity_HMSC/Output/Effect_Size_Most_up_to_date",
  pattern = "^Effect_Size", 
  full.names = TRUE
)

# Source each file quietly
lapply(files, source)

# --- Define a consistent axis theme ---
axis_theme <- theme(
  plot.tag = element_text(size=20, face='bold'),
  axis.title.x = element_blank(),
  axis.text.x = element_text(size=16, face='bold')
)

axis_theme_b <- theme(
  plot.tag = element_text(size=20, face='bold'),
  axis.title.x = element_blank(),
  axis.text.x = element_text(size=14, face='bold')
)

# --- Define a consistent axis theme ---
axis_theme_c <- theme(
  plot.tag = element_text(size=20, face='bold'),
  axis.text.x = element_text(size=16, face='bold')
)

axis_theme_d <- theme(
  plot.tag = element_text(size=20, face='bold')
)

# --- Apply to all plots ---
phy_mod          <- phy + axis_theme
explo_mod        <- explo + axis_theme_b
AM_rhiz_guild_mod <- AM_rhiz_guild + axis_theme
repo_mod         <- repo + axis_theme
Fire_interval_mod <- Fire_interval + axis_theme_c
site_mod <- site+ axis_theme_d
aridity_mod<- aridity + axis_theme_c
# --- Combine with aligned axes ---
top<-(phy_mod +repo_mod )  +plot_layout(
  widths = c(1, 1),     # Equal column widths
  axis_titles = "collect_y"   ) # Align shared legends/axes

middle<-(AM_rhiz_guild_mod + explo_mod )+plot_layout(
  widths = c(1, 1),     # Equal column widths
  axes = "collect_y"   ) # Align shared legends/axes

bottom<-( site_mod/Fire_interval_mod/aridity_mod)

p_fungal<-(top/middle) +
  plot_layout(
    widths = c(1, 1),     # Equal column widths
    heights = c(1, 1))
p_fungal  

p_site<-(bottom) +
  plot_layout(
    widths = c(1),     # Equal column widths
    heights = c(1, 1, 1))
p_site  

# --- Save ---
ggsave("plots/combined_fungal_size_Plot.png", p_fungal, width = 15, height = 12, dpi = 300)

ggsave("plots/combined_site_size_Plot.png", p_site, width = 13, height = 17, dpi = 300)

