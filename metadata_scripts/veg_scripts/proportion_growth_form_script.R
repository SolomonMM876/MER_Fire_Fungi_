library(tidyverse)

auscribe_veg_both<-readRDS('Processed_data/metadata/veg/auscribe_veg_both.Rdata')



library(tidyverse)

# Count frequency of growth_form per site
growth_form_counts <- auscribe_veg_both %>%
  filter(!is.na(herbarium_determination)) %>% 
  group_by(Site, growth_form) %>%
  summarise(count = n(), .groups = "drop")

# Check top rows
head(growth_form_counts)

# Stacked barplot
ggplot(growth_form_counts, aes(x = Site, y = count, fill = growth_form)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Growth Form Frequency per Site",
    x = "Site",
    y = "Frequency",
    fill = "Growth Form"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )


ggplot(growth_form_counts, aes(x = Site, y = count, fill = growth_form)) +
  geom_bar(stat = "identity", position = "fill") +  # proportion instead of absolute
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Proportional Growth Form Frequency per Site",
    x = "Site",
    y = "Proportion",
    fill = "Growth Form"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )






library(tidyverse)

# Calculate growth form proportions per site
growth_form_site <- auscribe_veg_both %>%
  filter(!is.na(herbarium_determination)) %>% 
  group_by(Site, growth_form) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Site) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Check results
head(growth_form_site)

ggplot(growth_form_site, aes(x = Site, y = prop, fill = growth_form)) +
  geom_bar(stat = "identity", position = "stack",color='black') +
  labs(
    title = "Proportion of Growth Forms per Site",
    x = "Site",
    y = "Proportion",
    fill = "Growth Form"
  ) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Calculate growth form proportions per plot within each site
growth_form_plot <- auscribe_veg_both %>%
  filter(!is.na(herbarium_determination)) %>% 
  group_by(Site, plot, growth_form) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Site, plot) %>%
  mutate(prop = count / sum(count)) %>%
  ungroup()

# Order plots within each site by total count for better visualization
growth_form_plot <- growth_form_plot %>%
  group_by(Site, plot) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(plot = fct_reorder(plot, total_count))

# Facetted stacked barplot
ggplot(growth_form_plot, aes(x = plot, y = prop, fill = growth_form)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(labels = scales::percent) +
  facet_wrap(~ Site, scales = "free_x") +   # switched to facet_wrap for better control
  labs(
    title = "Proportion of Growth Forms per Plot within Each Site",
    x = "Plot",
    y = "Proportion",
    fill = "Growth Form"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

growth_form_plot %>% distinct(growth_form)
