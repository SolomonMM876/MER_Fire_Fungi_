library(tidyverse)
library(vegan)
library(lme4)
library(car)
library(purrr)
library(APCalign)
library(broom.mixed)

#read in Auscribe year 1 data already extract from veg_import_clean
auscribe_veg<-readRDS('Processed_data/metadata/veg/auscribe_veg_yr1.Rdata')
#read in all site locations
All_Sites<-readRDS("raw_data/MER_Site_Data/All_Sites.RDS") %>% select(Site,Plot,Fire_Treatment)


#load for austraits
tax_resources <- load_taxonomic_resources()

# Extract distinct species names as a character vector
spp_names <- auscribe_veg %>%
  distinct(herbarium_determination) %>%
  filter(!is.na(herbarium_determination)) %>% 
  pull(herbarium_determination) %>%  # extract as a character vector
  iconv(from = "", to = "UTF-8", sub = "byte")  # Convert to valid UTF-8


# Align the species names using APC resources
aligned_spp <- align_taxa(spp_names, resources = tax_resources)

upd_spp_list<-update_taxonomy(aligned_spp, taxonomic_splits = "most_likely_species",resources = tax_resources)%>%
  select(original_name,suggested_name,aligned_name,genus:taxon_rank)%>%
  filter(!is.na(suggested_name))#remove indeterminate IDs

#read table of myco IDs
brundrett_2017 <- read_csv("Processed_data/metadata/veg/datasets_published/brundrett_2017_myco.csv")%>%
  #correct family names so joins work
  mutate(Family = str_replace(Family, "^(Euphorbaceae|Fabaceae)\\b.*", "\\1"),
         Family = str_replace(Family, "Euphorbaceae", "Euphorbiaceae"))

#read in North american classification dataset
Brundrett_Tedersoo_2020.csv<- read_csv("Processed_data/metadata/veg/datasets_published/Brundrett_Tedersoo_2020.csv")

library(readxl)
#read in global fungalroot db for myco association Nadejda A. Soudzilovskaia et al 2020
Soudzilovskaia2020<-read_excel('Processed_data/metadata/veg/datasets_published/Soudzilovskaia_etal_2020.xlsx', sheet = 'Table S2',skip = 1)


# First join: Match trait data based on Genus, some taxa are dual, but we just care if they are or are not
genus_myco <- upd_spp_list %>%
  left_join(
    brundrett_2017 %>%
      rename(Myco_type = Primary.Obs) %>%
      group_by(Genus) %>%
      summarise(Myco_type = paste(unique(Myco_type), collapse = " | "), .groups = "drop"),  
    by = c("genus" = "Genus")
  ) %>%
  mutate(
    Myco_type = case_when(
      genus %in% c("Thysanotus", "Cassytha", "Olax",'Anarthria', 'Crassula') ~ "NM",
      genus %in% c('Doryanthes',"Brunoniella",'Caesia') ~ "Unkown",
      TRUE ~ Myco_type),
    Notes = case_when(
      genus == "Thysanotus" ~ "Unique sub-epidermal association - McGee 1988, Brundrett & Abbott 1991.",
      genus == "Cassytha" ~ "Can have members of family that are mycorrhizal, but all species observed are hemi-parasitic growing on stems and are NM as such.",
      genus == "Brunoniella" ~ "Type not determined due to lack of records, removed from analysis.",
      #taxon_name == "Doryanthes excelsa" ~ "Type not determined due to lack of records, classified as NM for analysis.",
      #taxon_name == "Caesia sp." ~ "Type not determined due to lack of records, classified as NM for analysis.",
      genus == "Olax" ~ "Non-mycorrhizal Bellgard et al. 1994",
      genus == "Anarthria" ~ "Non-myorrhizal Brundrett 2009 (https://mycorrhizas.info/nmplants.html)",
      genus == "Crassula" ~ "Non-myorrhizal Brundrett 2009 Most samples have NM roots, some few examples of VAM *More sampling required.(https://mycorrhizas.info/nmplants.html)",
      is.na(Myco_type) ~ "Classification based on other genera in family",
      TRUE ~ ""))

#join remaining based on family
family_myco <- genus_myco %>%
  # Join to Brundrett data again, this time by Family
  left_join(
    brundrett_2017 %>%
      rename(Myco_type_fam = Primary.Obs) %>%
      group_by(Family) %>%
      summarise(
        Myco_type_fam = paste(unique(Myco_type_fam), collapse = " | "),
        .groups = "drop"
      ),
    by = c("family" = "Family")  # Join on family names
  ) %>%
  # If genus-level Myco_type is still NA, use the family-level classification
  mutate(Myco_type = coalesce(Myco_type, Myco_type_fam)) %>%
  # Drop the intermediate family trait column
  select(-Myco_type_fam) 

#combine classifications
myco_tax<-family_myco%>%
  left_join(genus_myco) %>% 
  rename(herbarium_determination=original_name)

#ID unknown plant taxa
unk_myco_plants<-myco_tax%>%
  filter(is.na(Myco_type))%>%
  distinct()

#ID unknown myco types by Soudzilovskaia et al 2020
Soudzilovskaia_genus_myco<-unk_myco_plants %>% 
  left_join(Soudzilovskaia2020, by=c('genus'='Genus')) %>% 
  filter(!is.na(`Mycorrhizal type`)) %>% 
  mutate(Myco_type = coalesce(Myco_type, `Mycorrhizal type`),
         Myco_type=if_else(Myco_type=='uncertain','Unkown',Myco_type),
         Myco_type=if_else(Myco_type=='AM','VAM',Myco_type),
         Notes= 'Classification based Soudzilovskaia et al 2020 FungalRoot databse') %>% 
  select(-`Mycorrhizal type`)


# Remove original unknown rows and add updated identifications
myco_tax_final <- myco_tax %>%
  filter(!(is.na(Myco_type)) ) %>%
  bind_rows(Soudzilovskaia_genus_myco)

#quick check
myco_tax_final%>%filter(is.na(Myco_type))%>%distinct()
#######################count Frequency myco hosts##############################################################

# Step 1: Join vegetation with mycorrhizal classifications, remove 'Unknown' types
auscribe_myco <- auscribe_veg %>% 
  left_join(myco_tax_final) %>% 
  filter(!str_detect(Myco_type, "Unkown"))%>% 
  select(-plot, Plot = plot_join) 

# Define filters
all_myco <- auscribe_myco %>% filter(str_detect(Myco_type, "VAM|Ericoid|ECM"))
ecm_only <- auscribe_myco %>% filter(str_detect(Myco_type, "ECM"))
vam_only <- auscribe_myco %>% filter(str_detect(Myco_type, "VAM"))
non_myco <- auscribe_myco %>% filter(!str_detect(Myco_type, "VAM|Ericoid|ECM"))

#create function for calculating alpha
calculate_alpha_div <- function(df, all_sites_df, group_label) {
  # 1. Community matrix
  veg_comm <- df %>%
    filter(!is.na(suggested_name)) %>% 
    count(Plot, herbarium_determination) %>%
    pivot_wider(names_from = herbarium_determination, values_from = n, values_fill = 0)
  
  # 2. Calculate alpha diversity
  alpha_div <- veg_comm %>%
    rowwise() %>%
    mutate(
      Shannon = diversity(c_across(-Plot), index = "shannon"),
      Simpson = diversity(c_across(-Plot), index = "simpson"),
      Chao1   = fossil::chao1(c_across(-Plot)),
      Observed = sum(c_across(-Plot) > 0),
      Pielou  = ifelse(Observed > 1, Shannon / log(Observed), NA)
    ) %>%
    ungroup()
  
  # 3. Join with site metadata
  alpha_div_full <- all_sites_df %>%
    select(Site, Plot, Fire_Treatment) %>%
    left_join(alpha_div %>% 
                select(Plot,Shannon,Simpson,Chao1,Observed,Pielou), by = "Plot") %>%
    mutate(Myco_group = group_label)
  
  return(alpha_div_full)
}


alpha_all_myco <- calculate_alpha_div(all_myco, All_Sites, "All_Mycorrhizal")
alpha_ecm      <- calculate_alpha_div(ecm_only, All_Sites, "ECM_only")
alpha_vam      <- calculate_alpha_div(vam_only, All_Sites, "VAM_only")
alpha_non      <- calculate_alpha_div(non_myco, All_Sites, "Non_Mycorrhizal")




alpha_all_groups <- bind_rows(alpha_all_myco, alpha_ecm, alpha_vam, alpha_non)

# Define diversity metrics
metrics <- c("Shannon", "Simpson", "Chao1", "Pielou")

# Define mycorrhizal groups
myco_groups <- c("All_Mycorrhizal", "ECM_only", "VAM_only", "Non_Mycorrhizal")

# Filter to remove rows with missing values in metrics
alpha_all_groups_filtered <- alpha_all_groups %>%
  filter(if_all(all_of(metrics), ~ !is.na(.)))

#set factor order
alpha_all_groups_filtered <- alpha_all_groups_filtered %>%
  mutate(Fire_Treatment = factor(Fire_Treatment, levels = c("U", "B")))


# Fit LMMs for each combination of Myco_group and Metric
# Loop over Metric and Myco_group combinations
anova_table_all_myco <- cross_df(list(Metric = metrics, Myco_group = myco_groups)) %>%
  mutate(
    results = map2(Metric, Myco_group, function(metric, group) {
      # Filter for group
      df_sub <- alpha_all_groups_filtered %>%
        filter(Myco_group == group)
      
      # Skip if underpowered
      if (n_distinct(df_sub$Fire_Treatment) < 2 | length(unique(df_sub$Site)) < 2) {
        return(NULL)
      }
      
      # Model formula
      formula <- as.formula(paste0(metric, " ~ Fire_Treatment + (1|Site)"))
      model <- tryCatch(lmer(formula, data = df_sub), error = function(e) return(NULL))
      if (is.null(model)) return(NULL)
      
      # ANOVA
      anov <- tryCatch(Anova(model, test = "F"), error = function(e) return(NULL))
      if (is.null(anov)) return(NULL)
      
      # Extract Fire_Treatment effect
      tidy_model <- tryCatch(
        tidy(model, effects = "fixed") %>%
          filter(str_detect(term, "Fire_Treatment")),
        error = function(e) tibble(estimate = NA, std.error = NA)
      )
      
      # Assemble summary table
      tibble(
        Group = group,
        Metric = metric,
        Factor = rownames(anov)[1],
        Estimate = tidy_model$estimate[1],
        Std_Error = tidy_model$std.error[1],
        DF_num = anov$Df[1],
        DF_denom = anov$Df.res[1],
        F = anov$F[1],
        P = anov$`Pr(>F)`[1]
        )
    })
  ) %>%
  select(results) %>%
  unnest(results)

# Final table
print(anova_table_all_myco)

save(anova_table_all_myco, file='HMSC_MER/Output/Processed_Data/Vegitation_myco_alpha_metrics.RDS')


alpha_long_all <- alpha_all_groups %>%
  pivot_longer(cols = c("Shannon", "Simpson", "Chao1", "Pielou"), 
               names_to = "Metric", 
               values_to = "Value") %>%
  select(Site, Plot, Fire_Treatment, Group=Myco_group, Metric, Value) %>% 
  filter(!is.na(Value))




library(ggpubr)

# 1. Create bracket_df including both Metric and Myco_group
bracket_df <- anova_table_all_myco %>%
  mutate(
    p_label = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      P < 0.1   ~ "+"
    ),
    xmin = 1,
    xmax = 2
  ) %>%
  left_join(
    alpha_long_all %>%
      group_by(Metric, Group) %>%
      summarise(y.position = max(Value, na.rm = TRUE)*.9, .groups = "drop"),
    by = c("Metric", "Group")
  ) %>%
  select(Metric, Group, xmin, xmax, y.position, p_label)

#plot

ggplot(alpha_long_all, aes(x = Fire_Treatment, y = Value)) +
  geom_boxplot(aes(fill = Fire_Treatment),outliers=FALSE,alpha = 0.6, linewidth=1) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  geom_text(
    data = bracket_df,
    aes(x = .5, y = y.position, label = p_label),  # position label near top-left
    hjust = 0, vjust = 0,
    size = 9
  ) +
  facet_grid(Metric ~ Group, scales = "free_y") +
  labs(
    title = "Alpha Diversity by Mycorrhizal Host Group and Fire Treatment",
    x = "Fire Treatment",
    y = "Alpha Diversity Value"
  ) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( hjust = 0.5, size = 25, face = "bold"),
        axis.text.y = element_text(size = 23, face = "bold"),
        axis.title.x = element_text(size = 30, face = "bold"),
        axis.title.y = element_text(size = 27, face = "bold"),
        strip.text = element_text(size = 27, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 3),
        axis.line = element_line(linewidth = 2, colour = "black")
  )



