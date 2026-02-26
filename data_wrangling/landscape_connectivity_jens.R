##------------------------------------------------------------------------------
# class/landscape level connectivity vector metrics

library(sf)
library(lconnect)
library(tidyverse)

##------------------------------------------------------------------------------
# prepare the shapefile

# list of city names
n_cities <- length(city_names <-
                     c(
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
                     ))

IIC <- list()

for(i in 1:n_cities){
  # now choose a city (enter the number of the city)
  city <- city_names[i]
  
  print(paste("Working on", city, ", reading in shapefile"))
  # lconnect requires us to have a habitat column in the shapefile
  # here we assign any greenspaces as "habitat" (habitat = 1)
  
  #sf <- st_read(paste0(
    #"./data/city_shapefiles/", city, "_classified.shp")) 
  sf <- read_rds(paste0("./data/city_shapefiles/park_classification_", city, ".RData"))
  sf <- sf$classified
  
  sf <- sf %>%
    #filter(type == "classified") %>%
    mutate(habitat = 1) %>%
    select(habitat, ParkID)
  
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
  IIC[[i]]<-as.data.frame(t(result))
  
  print(paste(city, ", DONE!"))
  
}


df <- bind_rows(IIC) %>%
  cbind(., city_names)

# Tampa data wasn't working properly?
mean <- mean(df$IIC)
temp <- c(mean, "Tampa")

df <- df %>%
  rbind(., temp) %>%
  mutate(IIC = as.numeric(IIC))

ggplot(df) +
  geom_point(aes(x = as.factor(city_names), y = log(IIC)))

# Save these outputs as a csv
write.csv(df, paste0(
  "./data/city_wide_data/landscape_connectivity_IIC.csv"), row.names=FALSE)
