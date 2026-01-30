##------------------------------------------------------------------------------
# class/landscape level connectivity vector metrics

library(sf)
library(lconnect)
library(tidyverse)

##------------------------------------------------------------------------------
# prepare the shapefile

# list of city names
n_cities <- length(city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Raleigh",
  "SD",
  "SF"
))


for(i in 1:n_cities){
  # now choose a city (enter the number of the city)
  city <- city_names[i]
  
  # lconnect requires us to have a habitat column in the shapefile
  # here we assign any greenspaces as "habitat" (habitat = 1)
  sf <- readRDS(paste0(
    "./data/city_shapefiles/park_classification_", city, 
    ".Rdata")) 
  
  sf <- sf[[1]]
  
  #plot(sf$geometry)
  
  sf <- sf %>%
    # filter(type == "classified") %>%
    mutate(habitat = 1) %>%
    select(habitat, ParkID)
  
  # save the modified file for lconnect 
  sf::write_sf(sf, paste0(
    "./data/city_shapefiles/", city, 
    "_classified.shp"))
}



##------------------------------------------------------------------------------
# calculate landscape connectivity metrics shapefile

IIC <- vector(length=n_cities)

for(i in 1:n_cities){
  # now choose a city (enter the number of the city)
  city <- city_names[i]
  
  # Load the landscape data
  land <- upload_land(paste0(
    "./data/city_shapefiles/", city, 
    "_classified.shp"), 
    habitat = 1, max_dist = 2000)
  
  # Confirm the class
  class(land)
  # Plot the landscape aggregate by clusters defined by the “max_dist” argument
  #plot(land, main = "Landscape clusters")
  
  # Compute the connectivity metrics
  IIC[i] <- con_metric(land, metric = c("IIC"))
}

df <- as.data.frame(cbind(city_names, IIC)) %>%
  mutate(IIC = as.numeric(IIC))

ggplot(df) +
  geom_point(aes(x = as.factor(city_names), y = log(IIC)))

# Save these outputs as shapefiles, using the sf package
write.csv(df, paste0(
  "./data/city_wide_data/landscape_metrics.csv"), row.names=FALSE)