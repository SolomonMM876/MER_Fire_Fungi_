library(tidyverse)
library(ggplot2)

#soil tax data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds') 

AM_tax<-myco_tax_soil %>% 
  filter(phylum=='Glomeromycota') %>% 
  rename(Family=family,Order=order)

AM_genera<-unique(AM_tax$genus)
# lookup table of genus → family → order
taxonomy_lookup <- tribble(
  ~genus,             ~Family,               ~Order,
  "Acaulospora",      "Acaulosporaceae",     "Diversisporales",
  "Ambispora",        "Ambisporaceae",       "Diversisporales",
  "Archaeospora",     "Archaeosporaceae",    "Archaeosporales",
  "Bulbospora",       "Glomeraceae",         "Glomerales",
  "Cetraspora",       "Gigasporaceae",       "Diversisporales",
  "Claroideoglomus",  "Claroideoglomeraceae","Glomerales",
  "Corymbiglomus",    "Glomeraceae",         "Glomerales",
  "Dentiscutata",     "Gigasporaceae",       "Diversisporales",
  "Diversispora",     "Diversisporaceae",    "Diversisporales",
  "Dominikia",        "Glomeraceae",         "Glomerales",
  "Endogone",         "NA",                  "NA",
  "Entrophospora",    "Entrophosporaceae",   "Diversisporales",
  "Funneliformis",    "Glomeraceae",         "Glomerales",
  "Fuscutata",        "Gigasporaceae",       "Diversisporales",
  "Geosiphon",        "Geosiphonaceae",      "Geosiphonales",   # check, sometimes treated separately
  "Gigaspora",        "Gigasporaceae",       "Diversisporales",
  "Glomus",           "Glomeraceae",         "Glomerales",
  "Septoglomus",      "Glomeraceae",         "Glomerales",
  "Intraornatospora", "Gigasporaceae",       "Diversisporales",
  "Kamienskia",       "Glomeraceae",         "Glomerales",
  "Otospora",         "Glomeraceae",         "Glomerales",
  "Pacispora",        "Pacisporaceae",       "Diversisporales",
  "Palaeospora",      "Glomeraceae",         "Glomerales",
  "Paradentiscutata", "Gigasporaceae",       "Diversisporales",
  "Paraglomus",       "Paraglomeraceae",     "Paraglomerales",
  "Paurocotylis_fulva var. zealandica", "none (invalidly published)", "none (invalidly published)",
  "Racocetra",        "Gigasporaceae",       "Diversisporales",
  "Redeckera",        "Glomeraceae",         "Glomerales",
  "Rhizoglomus",      "Glomeraceae",         "Glomerales",
  "Rhizophagus",      "Glomeraceae",         "Glomerales",
  "Sacculospora",     "Glomeraceae",         "Glomerales",
  "Sclerocystis",     "Glomeraceae",         "Glomerales",
  "Scutellospora",    "Gigasporaceae",       "Diversisporales"
)


#read in Carlos AMF data
#read in hyph comm data
AMF_repo_type<-read.csv('Processed_data/Seq_dat/datasets_external/AMF_Spore_Database_Volume.csv') %>% 
  separate(good.names, into = c("genus", "species"), sep = "_", remove = FALSE)

AMF_repo_type_fixed <- AMF_repo_type %>%
  left_join(taxonomy_lookup, by = "genus") %>%
  mutate(
    Family = if_else(Family.x == "" | is.na(Family.x), Family.y, Family.x),
    Order  = if_else(Order.x == ""  | is.na(Order.x),  Order.y,  Order.x)
  ) %>%
  select(-Family.y, -Order.y,-Family.x,-Order.x)


summarise_trait_distribution <- function(df, tax_level = "xxx") {
  
  # check tax_level column exists
  if (!tax_level %in% colnames(df)) {
    stop(paste("Column", tax_level, "not found in dataframe"))
  }
  
  traits <- c("Spore_formation", "Spore_placement", "Sporocarp", "Sporocarp_placement")
  
  # reshape into long format
  trait_long <- df %>%
    mutate(across(all_of(traits), ~as.factor(.x))) %>% 
    pivot_longer(
      cols = all_of(traits),
      names_to = "Trait",
      values_to = "Value"
    )
  
  AM_phylo<-AM_tax %>% distinct(.data[[tax_level]]) %>% pull()
  
  # summarise counts and percentages per tax group × trait × value
  trait_summary <- trait_long %>%
    filter(.data[[tax_level]] %in% AM_phylo) %>% 
    group_by(.data[[tax_level]], Trait, Value) %>%
    summarise(
      n_obs = n(),
      .groups = "drop_last"
    ) %>%
    group_by(.data[[tax_level]], Trait) %>%
    mutate(
      total_obs = sum(n_obs),
      perc = round((n_obs / total_obs) * 100, 0)
    ) %>%
    ungroup() %>%
    rename(phylo = !!tax_level)
  
  # summary table for annotating total_obs
  total_labels <- trait_summary %>%
    distinct(phylo, total_obs) %>% 
    mutate(Trait='Spore_formation')
  
  # plot: stacked bar chart with value + % labels, no legend
  p <- ggplot(trait_summary, aes(x = Trait, y = perc, fill = Value)) +
    geom_col(position = "stack") +
    geom_text(aes(label = paste0(Value, " (", perc, "%)")),
              position = position_stack(vjust = 0.5),
              size = 3, color = "black") +
    geom_text(data = total_labels,
              aes(x = Trait, y = 110, label = paste0("n=", total_obs)),
              inherit.aes = FALSE,
              size = 3.5) +
    facet_wrap(~ phylo, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      #axis.text.x = element_text(angle = 45, hjust = 0.5,vjust=.7,size=10),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    labs(title = paste("Trait distributions across", tax_level),
         y = "% of observations",
         x = "Trait")
  
  return(list(summary = trait_summary, plot = p))
}




# Genus level
out_genus <- summarise_trait_distribution(AMF_repo_type_fixed, "genus")
genus<-out_genus$summary   # table
genus
out_genus$plot      # ggplot

# Family level
out_family <- summarise_trait_distribution(AMF_repo_type_fixed, "Family")
family <-out_family$summary   # table
family 
out_family$plot   

# Order level
out_order <- summarise_trait_distribution(AMF_repo_type_fixed, "Order")
order <-out_order$summary   # table
order
out_order$plot 

