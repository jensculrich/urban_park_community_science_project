library(tidyverse)

# list of city names
n_cites <- (length(city_names <- c(
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
)))


## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

calculate_shannon <- function(proportions) {
  # Ensure input is a vector of proportions
  if(round(sum(proportions), 5) != 1) {
    stop("Input data must be proportions that sum to 1.")
  }
  
  # Remove zero proportions to avoid issues with log(0)
  proportions <- proportions[proportions > 0]
  
  # Calculate H' = -sum(pi * ln(pi))
  H <- -sum(proportions * log(proportions))
  
  return(H)
}

#-----------------------------------------------------------------------------
# get park site land cover and isolation data
# which is stored in separate files for each city

site_data <- data.frame(matrix(ncol = 17, nrow = 0))

for(i in 1:length(city_names)){
  
  city <- city_names[i]
  
  # first read the data 
  temp <- cbind(city, read.csv(paste0(
    "./data/detections_by_city/", city, "/04_0m_", city,
    "_isolation_non_water_only.csv"
  ))) %>%
    
    # standardize infitsimely small greenspace areas
    mutate(total_green_space_area = ifelse(total_green_space_area < 0.1, 0.1, total_green_space_area),
           connectivity = isolation * (-1)) %>%
    #   get scaled covariate values
    select(city, new_id, total_area_sqm, total_green_space_area, connectivity, tree_percent_cover, grass_shrub__percent_cover)

  # add park flower data
  # first read the data 
  flower_data <- cbind(city, read.csv(paste0(
    "./data/detections_by_city/", city, "/03_0m_", city,
    "_flowers_classified_park.csv"
  ))) %>%
    
    # get number of flowering plant genera per site
    group_by(new_id, genus) %>%
    # only need one row per genus in each site to count number of genera at each site
    slice(1) %>%
    ungroup() %>%
    group_by(city, new_id) %>%
    summarise(n_plant_genera = n()) %>%
    ungroup()
  
  # and join the variables we want with the site data by site id
  temp <- temp %>%
    left_join(., flower_data, by=c("city", "new_id")) %>%
    mutate(n_plant_genera = replace_na(n_plant_genera, 0),
           log_n_plant_genera = log(n_plant_genera+1))
  
  site_data <- rbind(site_data, temp)
  
}

#-----------------------------------------------------------------------------
# get landscape data
landscape_data <- data.frame(matrix(ncol = 12, nrow = 0))

for(i in 1:length(city_names)){
  
  city <- city_names[i]
  
  # first read the data 
  temp <- cbind(city, read.csv(paste0(
    "./data/buffer_around_parks/2km_park_buffer/01_", city, 
    "_2km_buffer_area_around_park.csv"
  ))) %>%
    
    # get landscape vegetation metrics  
    mutate(
      # proportion of landscape that is a natural land class
      proportion_landscape_vegetation = 
        rowSums(.[,c(13, 14, 15, 16, 17, 18, 19, 21, 22)]) / total_sur_area_sqm,
      # proportion of landscape that is open developed area
      proportion_landscape_open_developed = 
        .[,9] / total_sur_area_sqm,
      # proportion of landscape that is med-highly developed area
      proportion_landscape_medhigh_developed = 
        rowSums(.[,c(11, 12)]) / total_sur_area_sqm,
      # proportion of landscape that is woody veg
      proportion_landscape_woody = 
        rowSums(.[,c(14, 15, 16, 17, 21)]) / total_sur_area_sqm,
      # proportion of landscape that is grassland/herbaceous
      proportion_landscape_grassherb = 
        rowSums(.[,c(18, 19, 22)]) / total_sur_area_sqm,
    ) %>%
    # get shannon diversity of viable landscape types
    mutate(sum_viable_lands = rowSums(.[,c(9, 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)]),
           p_9 = .[,9] / sum_viable_lands,
           p_10 = .[,10] / sum_viable_lands,
           p_13 = .[,13] / sum_viable_lands,
           p_14 = .[,14] / sum_viable_lands,
           p_15 = .[,15] / sum_viable_lands,
           p_16 = .[,16] / sum_viable_lands,
           p_17 = .[,17] / sum_viable_lands,
           p_18 = .[,18] / sum_viable_lands,
           p_19 = .[,19] / sum_viable_lands,
           p_20 = .[,20] / sum_viable_lands,
           p_21 = .[,21] / sum_viable_lands,
           p_22 = .[,22] / sum_viable_lands,
    ) 
  shannon_diversity <- apply(temp[,seq(30,41)], 1, calculate_shannon)
  temp <- cbind(temp, shannon_diversity) %>%
    
    #   get scaled covariate values
    select(city, new_id, landcover_type_diversity, shannon_diversity,
           proportion_landscape_vegetation, proportion_landscape_open_developed, 
           proportion_landscape_medhigh_developed, proportion_landscape_woody, proportion_landscape_grassherb)
  
  
  landscape_data <- rbind(landscape_data, temp)
  
}

site_data <- left_join(site_data, landscape_data)

write.csv(site_data, "./data/derived_park_site_covariate_data/derived_park_site_covariate_data.csv", row.names = FALSE)
