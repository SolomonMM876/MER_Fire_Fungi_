
# Load libraries
library(officer)
library(flextable)
library(dplyr)

#SOIL PERM OUTPUT
# Round numeric columns to 2 decimals
model_results_rounded <- model_results_soil %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Create flextable with rounded values
ft <- model_results_rounded %>%
  select(variable, estimate, std_error, p_value, adj_r2, f_statistic,
         df1, df2, r2, AIC) %>%
  flextable() %>%
  autofit() %>%
  hline(i=11) %>% 
  set_caption("Univariate Linear Model Results") 

ft

# Create Word doc and add table
doc <- read_docx() %>%
  body_add_par("Model Summary Table", style = "heading 1") %>%
  body_add_flextable(ft)

# Save Word file
print(doc, target = "Tables/Output/permanova_soil_metadata.docx")



################Soil veg


#SOIL PERM OUTPUT
# Round numeric columns to 2 decimals
model_results_rounded <- model_results_soil_veg %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

# Create flextable with rounded values
ft <- model_results_rounded %>%
  select(variable, estimate, std_error, p_value, adj_r2, f_statistic,
         df1, df2, r2, AIC) %>%
  flextable() %>%
  autofit() %>%
  set_caption("Univariate Linear Model Results") 

ft

# Create Word doc and add table
doc <- read_docx() %>%
  body_add_par("Model Summary Table", style = "heading 1") %>%
  body_add_flextable(ft)

# Save Word file
print(doc, target = "Tables/Output/permanova_soil_metadata.docx")

