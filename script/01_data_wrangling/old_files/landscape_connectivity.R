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
  "NY",     
  "Philadelphia",
  "Raleigh",
  "SD",
  "SF",
  "tampa",
  "denver"
))

summary<-data.frame()

for(i in 1:n_cities){
  # now choose a city (enter the number of the city)
  city <- city_names[i]
  
  print(paste("Working on", city, ", reading in shapefile"))
  # lconnect requires us to have a habitat column in the shapefile
  # here we assign any greenspaces as "habitat" (habitat = 1)
  
sf <- st_read(paste0(
    "/Volumes/sea_angel/iNat_urbanwatch/data/shapefile_0m_buffered_park_2km_regional_pool/",city, "_0_buffered_park_2km_regional_pool.shp")) 
  
  
  sf <- sf %>%
    filter(type == "classified") %>%
    mutate(habitat = 1) %>%
    select(habitat, new_id)
  
  # Create a temporary shapefile path
  temp_shp <- tempfile(fileext = ".shp")
  
  sf::write_sf(sf, temp_shp)
  


##------------------------------------------------------------------------------
# calculate landscape connectivity metrics shapefile
  print(paste("Working on", city, ", calculating landscape connectivity metrics"))
  
  # Load the landscape data
  land <- upload_land(temp_shp, 
                      habitat = 1, 
                      max_dist = 2000)

  # requires us to set a max dist between which parks can be connected
  # I chose 2000 metres to be consistent with our park site isolation metric 
  
  # Confirm the class
  class(land)
  # Plot the landscape aggregate by clusters defined by the “max_dist” argument
  #plot(land, main = "Landscape clusters")
  
  # Compute the connectivity metrics
  # IIC is fast to calculate
  # other metrics of interest might include AWF
  # https://www.r-bloggers.com/2019/03/lconnect-connectivity-metrics/
  result <- con_metric(land, metric = c("IIC"))

  summary<-rbind(summary, data.frame(city = city, IIC = result))
  
}


ggplot(summary) +
  geom_point(aes(x = as.factor(city, y = log(IIC)))

# Save these outputs as a csv
write.csv(summary, paste0(
  "/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/add_on_parameters/02_urbanwatch_city_wide_connectivity_metrics_classified_parks_only.csv"), row.names=FALSE)
