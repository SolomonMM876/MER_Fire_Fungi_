
library(tidyverse)

#DAte accessed: 16/7/25

fungal_repro_strategy <- tribble(
  ~Clade_No, ~Family_or_Group, ~Epigeous, ~hypogeous_semi_hypogeous, ~Resupinate,
  "1 (34)", "Lyophyllum", "Lyophyllum*", NA, NA,
  "2 (37-38)", "Tricholomataceae", "Leucopaxillus, Tricholoma*", NA, NA,
  "3 (43/44)", "Entolomataceae", "Entoloma* (not all species)", "Rhodogaster, Richoniella", NA,
  "4 (55)", "Amanitaceae", "Amanita, Limacella", "Ammarrendia, Torrendia", NA,
  "5 (65)", "Hygrophoraceae", "Gliophorus, Humidicutis, Hygrophorus*", NA, NA,
  "6 (72)", "Hydnangiaceae (Laccaria clade)", "Laccaria", "Hydnangium*, Podohydnangium, Gigaspera", NA,
  "7 (73)", "Cortinariaceae", "Cortinarius*, Dermocybe, Naucoria*, Rozites*", "Cortinarius, Destuntzia, Protoglossum, Quadrispora, Setchelliogaster*, Stephanopus*, Thaxterogaster*", NA,
  "8 (74)", "Phaeocollybia clade in Hymenogastraceae", "Phaeocollybia*", NA, NA,
  "9 (95)", "Hebeloma clade in Hymenogastraceae", "Hebeloma*", "Hymenogaster (in part)", NA,
  "10 (78/79)", "Descolea clade ", "Descolea*", "Descomyces*, Setchelliogaster", NA,
  "11 (104)", "Inocybaceae", "Inocybe*, Auritella", "Auritella", NA,
  "12A (Outgroup 1)", "Boletes (Boletales) Includes many families", 
  "Aureoboletus, Austroboletus, Austropaxillus*, Boletellus*, Boletochaete, Boletus*, Chroogomphus*, Fuscoboletinus, Heimielia, Gomphidius*, Gyroporus*, Leccinum*, Paragyrodon, Paxillus*, Phlebopus, Phylloporus, Poryphyrellus*, Psiloboletinus, Rubinoboletus, Strobilomyces, Suillus*, Tylopilus*, Xanthoconium, Xerocomus*",
  "Alpova*, Austrogaster, Austrogautieria, Chamonixia*, Gastroboletus, Gastrotylopilus, Gymnogaster, Gymnopaxillus, Horakiella, Hysterogaster, Melanogaster*, Mycoamaranthus, Octaviania, Rhizopogon*, Royoungia, Sclerogaster, Truncocolumella*, Timgrovea, Wakefieldia",
  "Serpula (are saprophytes)",
  "12B (Outgroup 1)", "Boletales - Sclerodermatales", "Pisolithus*, Scleroderma*, Calostoma*", "Astaeus*, Scleroderma, Velligaster", NA,
  "13 (Outgroup 2)", "Russulales", "Lactarius*, Russula*", "Arcangeliella*, Cystangium, Gymnomyces, Leucogaster, Macowanites, Octaviania*, Stephanospora, Zelleromyces", "Byssoporia, Albatrellus*, Polyporoletus* (Stereum etc. are saprophytes)",
  "14 (Outgroup 3)", "Gomphales & Hysterangiales", "Bankera*, Boletopsis*, Clavaridelphus*, Gomphus*, Hydnum*, Hydnellum*, Phellodon*, Sarcodon*, Ramaria* (some NM?)",
  "Hysterangium clades: Aroramyces, Chondrogaster, Hysterangium* , Mesophelliaceae: Andebbia, Castoreum, Gummiglobus, Malajczukia, Mesophellia, Nothocastoreum, Gallacaceaceae: Austrogautiera, Gallacea, Phallogaster clade: Protrubera, Phallogaster, Trappea, Ramaria clade: Gautieria*",
  NA,
  "15 (Outgroup 4)", "Cantharellales", "Cantharellula, Cantharellus*, Craterellus*", NA, "Sistotrema?",
  "16AB (Outgroup(s) 5)", "Rhizoctonia alliance clades (Sebacinaceae, Ceratobasidiaceae)", NA, NA, "Thanatephorus*, Sebacina*",
  "17 (Outgroup 6)", "Thelephorales", "Bankera, Boletopsis, Thelephora*, Phellodon*, Sarcodon", NA, "Tomentella*, Pseudotomentella*, Tomentellopsis*",
  "18 (Outgroup 7)", "Clavariaceae", "Aphelaria, Clavaria, Clavariadelphus, Clavicorona, Clavulina, Clavulinopsis, Ramariopsis", NA, NA,
  "19 (Outgroup 8)", "Atheliales", NA, NA, "Amphinema*, Byssocorticium*, Byssosporia*, Piloderma*, Tylospora*",
  "20 (Ascomycete-1)", "Pezizales (inc. Tuberaceae, Helvellaceae)", "Genea*, Geopora, Humaria*, Hydnotrya, Helvella, Leucangium*, Pachyphloeus, Pulvinula, Sarcosphaera, Sphaerosporella*, Sphaerozone*, Tirmania, Tricharina*, Wilcoxina*", "Choiromyces, Dingleya, Eremiomyces, Kalaharituber, Labyrinthomyces, Pachyphloeus, Reddellomyces, Tuber*, Terfezia, Turmania", "Phialophora* (anamorph)",
  "21 (Ascomycete-2)", "Eurotiales (Elaphomycetaceae)", NA, "Elaphomyces* Pseudotulostoma*", NA,
  "22 (Ascomycete-3)", "Melanommatales", NA, NA, "Cenococcum (sterile fungus)",
  "23 (Zygomycete-1)", "Endogonaceae", NA, "Densospora, Endogone, Peridiospora, Youngiomyces", NA
)

# 1. Reshape to long format
strategy_long <- fungal_repro_strategy %>%
  pivot_longer(cols = c(Epigeous, hypogeous_semi_hypogeous, Resupinate),
               names_to = "repo_strategy", values_to = "taxa") %>%
  filter(!is.na(taxa)) %>%
  separate_rows(taxa, sep = ",\\s*") %>%
  mutate(
    # Keep only the part after the colon if present
    taxa = str_replace(taxa, "^.*?:\\s*", ""),
    
    # Remove asterisks
    taxa = str_remove(taxa, "\\*"),
    
    # Extract family names: any word ending in "aceae"
    family = str_extract(Family_or_Group, "\\b\\w*aceae\\b"),
    
    # Extract group/order names: any word ending in "ales"
    group = str_extract(Family_or_Group, "\\b\\w*ales\\b"),
    
    notes = str_extract(taxa, "\\(([^()]*)\\)"),     # Extract text in parentheses
    taxa = str_remove_all(taxa, "\\s*\\([^)]*\\)"),  # Remove the parenthetical notes from genus
    # Clean up taxa (e.g., remove * and trim spaces)
    taxa = str_trim(taxa)
  ) %>%
  select(genus = taxa, repo_strategy, family, group,notes) %>% 
  #remove spelling errors or incorrect classification
  filter(
    !(genus == "Austrogautiera" & group == "Gomphales"),
    !(genus %in% c("Auritella", "Cortinarius", "Scleroderma") & repo_strategy == "hypogeous_semi_hypogeous"),
    !(genus %in% c("Bankera", "Boletopsis", "Phellodon", "Sarcodon") & group == "Gomphales"),
    !(genus == "Cantharella"),
    !(genus == "Octaviania" & group == "Russulales"),
    !(genus == "Setchelliogaster" & is.na(family))
  )

#load myco tax data
myco_tax_soil<-readRDS('Processed_data/Seq_dat/Soil/myco_tax_soil.rds')
myco_tax_hyph<-readRDS('Processed_data/Seq_dat/Hyph/myco_tax_Hyph.rds')


myco_tax<-myco_tax_soil %>% bind_rows(myco_dat_hyph) %>% 
  distinct() %>% 
  filter(!phylum=='Glomeromycota')


# 2. Create fallback family-level strategy lookup
family_strategy <- fungal_repro_strategy %>%
  pivot_longer(cols = c(Epigeous, hypogeous_semi_hypogeous, Resupinate),
               names_to = "repo_strategy", values_to = "taxa") %>%
  group_by(Family_or_Group, repo_strategy) %>%
  filter(!is.na(taxa)) %>%
  summarise(.groups = "drop") %>%
  rename(family = Family_or_Group)


# 3. Join to `myco_tax` by genus
myco_tax_labeled_genus <- myco_tax %>%
  distinct(genus,family) %>% 
  left_join(strategy_long, by = c("genus")) %>%
  rename(repo_strategy_genus = repo_strategy,family=family.x) %>% 
  select(-family.y)

#combine the two
final_repo_strat<-myco_tax_labeled_genus%>%
  mutate(
    # Manual fixes 
    #all of these were duplicated or needed to be adjusted to line up with webpage
    repo_strategy_genus = case_when(
      genus == "Ruhlandiella" ~ "hypogeous_semi_hypogeous",
      genus == "Coltricia" ~ "hypogeous_semi_hypogeous",
      genus == 'Bonomyces'~ "hypogeous_semi_hypogeous",
      family == "Tricholomataceae" & is.na(genus)~ 'Epigeous',
      family == "Thelephoraceae" & is.na(genus)~ 'Epigeous',
      family =='Inocybaceae' & genus =='Tubariomyces' ~ 'Epigeous',
      family == 'Pluteaceae'& is.na(genus)~ 'Epigeous',
      genus %in% c('Serendipita', 'Pezoloma','Myxotrichum') ~ NA, 
      family %in% c('Boletaceae', 'Hysterangiaceae') & is.na(genus) ~ NA,
      TRUE ~ repo_strategy_genus
    ),
    
    notes = case_when(
      genus == "Ruhlandiella" ~ "Epigeous mushrooms (Warcup 1991)",
      genus == 'Coltrcia' ~ 'Epigeous mushrooms (Bian et al. 2022)',
      genus == 'Bonomyces'~ "Epigeous mushrooms (Mao et al. 2023)",
      family == "Tricholomataceae" & is.na(genus)~ 'Epigeous based on other genera (Brundrett, 2009)',
      family == "Thelephoraceae" & is.na(genus) ~ 'coral-like (form higly variable) (Brundrett, 2009)',
      family =='Inocybaceae' & genus =='Tubariomyces' ~ 'Epigeous mushrooms (Alvardo and Manjon 2010)',
      family == 'Pluteaceae'~ 'Epigeous mushrooms (Malysheva et al. 2023)',
      genus %in% c('Serendipita', 'Pezoloma','Myxotrichum') ~ 'Unable to determine reproductive structures from available literature.',
      family %in% c('Boletaceae', 'Hysterangiaceae') & is.na(genus)  ~ 'Taxonomic identification limited, conflicting fruiting structures within family.',
      TRUE ~ notes
    )
  ) %>%
  rename(repo_strategy=repo_strategy_genus) %>% 
  select(-group,family,genus,repo_strategy,notes) %>% 
  mutate( repo_strategy=str_replace(repo_strategy,'hypogeous_semi_hypogeous', 'Hypogeous or Semi-hypogeous') )





saveRDS(final_repo_strat,'Processed_data/Seq_dat/datasets_external/repro_strat_df.Rds')





