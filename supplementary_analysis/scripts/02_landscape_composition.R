#This script is for calculating the landscape composition in each city using the urbanwatch data
library(tidyverse)
library(data.table)
library(terra)
library(sf)


#list all the files in the repo
data_path <- "E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data"
file_list<-list.files(path = data_path, full.names = TRUE, pattern = "\\.tif$") 
data_list <- lapply(file_list, rast)

list_names <- sapply(file_list, function(x) {
  str_match(x, "/([^/]+)_aggregated_land_cover\\.tif")[,2]
})

data_list <-set_names(data_list, list_names)


# Get frequency table for each land cover type
freq_table <- lapply(names(data_list), function(city_name) {
  x <- data_list[[city_name]]
  cell_size <- res(x)[1] * res(x)[2]
  freq_df <- freq(x)
  freq_df$area_m2 <- freq_df$count * cell_size
  
  # Filter out unclassified and calculate percentage
  freq_df <- freq_df %>%
    filter(value != 0) %>%  # Exclude unclassified
    mutate(percentage = (count / sum(count)) * 100,
           land_cover_type = case_when(
             value == 1 ~ "road",
             value == 2 ~ "grass_shrub",
             value == 3 ~ "tree",
             value == 4 ~ "building",
             value == 5 ~ "parking_lot",
             value == 6 ~ "water",
             value == 7 ~ "barren",
             value == 8 ~ "agriculture",
             value == 9 ~ "other"),
           city = city_name) %>%
    select(city, land_cover_type, percentage) %>%
    pivot_wider(names_from = land_cover_type, values_from = percentage)
  
  return(freq_df)
})


land_cover_percentage<-bind_rows(freq_table)

write.csv(land_cover_percentage, "E:/phd_study/urban_park_community_science_project/supplementary_analysis/output/land_cover_percentage.csv", row.names = FALSE)

read.csv("E:/phd_study/urban_park_community_science_project/supplementary_analysis/output/land_cover_percentage.csv")


plot(x)



