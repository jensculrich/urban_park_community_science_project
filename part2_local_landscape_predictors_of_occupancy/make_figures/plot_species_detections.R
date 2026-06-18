library(tidyverse)

city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Denver",
  "Des_moines",
  "Detroit",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "SD",
  "SF",
  "St_louis",
  "Tampa"
)

city_names_labels <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "Washington D.C.",
  "Denton",
  "Denver",
  "Des Moines",
  "Detroit",
  "Houston",
  "Los Angeles",
  "Minneapolis",
  "New York City",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "San Diego",
  "San Francisco",
  "St. Louis",
  "Tampa"
)

city_names_labels <- as.data.frame(cbind(city_names_labels, seq(1:length(city_names_labels)))) %>%
  rename("city_number" = "V2") %>%
  mutate(city_number = as.integer(city_number))

# source the prep function
source("./part2_local_landscape_predictors_of_occupancy/run_model/get_species_ranges.R")

ranges_raw <- get_species_ranges(city_names)

temp <- select(ranges_raw, species, city) %>%
  group_by(city) %>%
  mutate(city_number = cur_group_id()) %>%
  ungroup() %>% 
  left_join(., city_names_labels) %>%
  select(-city) %>%
  pivot_wider(names_from = city_number, values_from = city_names_labels) %>%
  unite(city_combined, 2:23, na.rm = TRUE, remove = TRUE) %>%
  mutate(city_combined = gsub("_", ", ", city_combined))

df <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))$species_info

df <- left_join(df, temp)
df <- select(df, species, family, n_binary_detections, city_combined)

#write.csv(df, "./part2_local_landscape_predictors_of_occupancy/plot_results/species_detections.csv", row.names = FALSE)

df <- df %>%
  arrange(family) %>%
  rename("Species" = "species",
         "Family" = "family",
         "n" = "n_binary_detections",
         "Cities" = "city_combined") 

df <- select(df, Family, Species, Cities, n)

library(gt)

DT <- df
DT1 <- df[1:46,]
DT2 <- df[47:92,]

tab <- DT %>%
  gt() %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = Species) # Target specific column
  ) %>%
  data_color(
    columns = n,
    method = "numeric",
    palette = "viridis",
    domain = c(1, 3938)
  )  %>%
  cols_width(
    n ~ px(50),
    Cities ~ px(500),
    everything() ~ px(200)
  ) %>%
  tab_header(
    title = md("**Butterfly Species Included in Occupancy Models**. \\
  Our analyses included 183 different butterfly species from five different butterfly families.
                           Species are listed by alphabetical order within family. We modeled butterfly species occupancy within
                           cities that were considered to be within the geographical range. We considered a city to be within range
                           if the species was detected by iNaturalist users at least once in a park in the city, or if the species
                           was detected by iNaturalist users at least once anywhere within a 20 km buffer of the city administrative
                           boundary. Cities determined to be within range of each species are listed under \'Cities\' column. 
                           *n* indicates the number of unique (park by month by year) detections of each species in our dataset.
                           Duplicate detections of the same species from the same time and place were not counted. The *n* detection
                           counts are coloured from lowest (dark colours) to highest (light colours).")
  ) 

tab

#tab %>% gtsave("./part2_local_landscape_predictors_of_occupancy/figures/species_detections_tab_1.tex")
tab %>% gtsave("./part2_local_landscape_predictors_of_occupancy/figures/species_detections_tab_1.docx")
