##------------------------------------------------------------------------------
# class/landscape level connectivity vector metrics

library(lconnect)
library(tidyverse)

##------------------------------------------------------------------------------
# prepare the shapefile

# list of city names
city_names <- c(
  # list in alphabetical order
  "Atlanta",
  "Boston",
  "Charlotte",
  "Dallas",
  "DC",
  "Denton",
  "Houston",
  "NYC",
  "Philadelphia",
  "Raleigh"
)

# now choose a city (enter the number of the city)
city <- city_names[6]

# lconnect requires us to have a habitat column in the shapefile
# here we assign any greenspaces as "habitat" (habitat = 1)
land_sf <- sf::read_sf(paste0(
  "./data/city_shapefiles/", city, 
  "/", city, 
  "_50_buffered_park_2km_regional_pool.shp")) %>%
  # filter(type == "classified") %>%
  mutate(habitat = 1) %>%
  select(habitat, ParkID, ParkCnt, type, new_id)

# save the modified file for lconnect 
sf::write_sf(land_sf, paste0(
  "./data/city_shapefiles/", city, 
  "/", city, 
  "_50_buffered_park_2km_regional_pool_with_habitat_column.shp"))


##------------------------------------------------------------------------------
# calculate landscape connectivity metrics shapefile

# Load the landscape data
land <- upload_land(paste0(
  "./data/city_shapefiles/", city, 
  "/", city, 
  "_50_buffered_park_2km_regional_pool_with_habitat_column.shp"), 
  habitat = 1, max_dist = 1000)

# Confirm the class
class(land)
# Plot the landscape aggregate by clusters defined by the “max_dist” argument
plot(land, main = "Landscape clusters")

# Compute the connectivity metrics
metrics <- con_metric(land, metric = c("IIC"))

# Visualize the metrics
print(as.data.frame(metrics))
df <- as.data.frame(metrics) %>%
  mutate(city = city) %>%
  rename("IIC" = "metrics")
rownames(df) <- NULL

# Save these outputs as shapefiles, using the sf package
write.csv(df, paste0(
  "./data/city_shapefiles/", city, 
  "/", city, 
  "_landscape_metrics.csv"), row.names=FALSE)
#sf::st_write(land$landscape, "./data/city_shapefiles/NYC/land.shp")
#sf::st_write(importance$landscape, "./data/city_shapefiles/NYC/importance.shp")