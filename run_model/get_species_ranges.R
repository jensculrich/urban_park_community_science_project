library(tidyverse)

get_species_ranges <- function(city_names) {
  
  ranges <- as.data.frame(matrix(nrow = 0, ncol = 71))
  
  for(i in 1:length(city_names)){
    
    city <- city_names[i]
    
    # first read the data 
    temp <- cbind(city, read.csv(paste0(
      "./data/detections_by_city/", city, "/02_0m_", city,
      "_regional_species_pool.csv"))
    )
    
    if(colnames(temp[2]) == "X.1"){
      temp <- select(temp, -X.1)
    }
    
    ranges <- rbind(ranges, temp)
    
  }
  
  butterfly_families <- c("Hesperiidae", "Lycaenidae", "Nymphalidae", 
                          "Papilionidae", "Pieridae", "Riodinidae")
  
  ranges <- ranges %>%
    filter(family %in% butterfly_families) %>%
    filter(species != "") %>%
    group_by(city, species) %>%
    slice(1) %>%
    ungroup() %>%
    select(city, family, species) %>%
    mutate(in_range = 1)
  
  ## --------------------------------------------------
  # Return stuff
  return(ranges)
  
}