library(tidyverse)

## --------------------------------------------------
## Operation Functions
## predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

get_site_data <- function(city_names) {
  
  site_data <- data.frame(matrix(ncol = 15, nrow = 0))
  
  for(i in 1:length(city_names)){
    
    city <- city_names[i]
    
    # first read the data 
    temp <- cbind(city, read.csv(paste0(
      "./data/detections_by_city/", city, "/04_50m_", city,
      "_isolation.csv"
    ))) %>%
    
    # get scaled (log) green space area  
    filter(total_green_space_area > 0) %>%
    select(city, new_id, total_green_space_area, isolation, tree_percent_cover, grass_shrub__percent_cover) %>%
    mutate(log_total_green_space_area = log(total_green_space_area),
           log_total_green_space_area_scaled = center_scale(log_total_green_space_area),
           log_isolation_scaled = center_scale(log(isolation)),
           tree_cover_scaled = center_scale(tree_percent_cover),
           grass_shrub_cover_scaled = center_scale(grass_shrub__percent_cover))
    
    # add park flower data
      # first read the data 
      flower_data <- cbind(city, read.csv(paste0(
        "./data/detections_by_city/", city, "/03_50m_", city,
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
             log_n_plant_genera = log(n_plant_genera+1),
             plant_genera_density = log_n_plant_genera / log_total_green_space_area) %>%
      mutate(plant_genera_density_scaled = center_scale(plant_genera_density))
    
    cor <- cor(temp$log_total_green_space_area_scaled, temp$log_isolation_scaled)
    print(paste0(city, " - correlation between park size and isolation = ", cor))
    
    site_data <- rbind(site_data, temp)
    
  }
  
  ## --------------------------------------------------
  # Return stuff
  return(site_data)
  
}