#Read in Suzzanne Prober  AusEcoModel spreadsheet
AusEcoModel<-read_excel("Processed_data/metadata/veg/MER eco_type fire details 3Jul25_SP.xlsx") %>% 
  dplyr::select(Site,site_nickname=Manuscript_names, Veg_type=Relevant_archetype_model,expected_fire_frequency=fire_interval)

AusEcoModel %>% 
  mutate(Veg_type=str_remove(Veg_type,'Related to Coastal floodplain eucalypt forests (see Casuarina expression)'))
