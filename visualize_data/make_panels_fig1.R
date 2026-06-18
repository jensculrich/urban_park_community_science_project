# use this file to plot the locations of the cities in our study onto a map of the continental U.S.

library(tigris)
library(ggplot2)
library(sf)

# Fetch US state boundaries (cb = TRUE gives cartographic boundaries that plot faster)
us_states <- states(cb = TRUE, year = 2024)

# remove non-continental states / territories
us_states <- us_states[-c(18,28,29,35,49,54,55),]

# get lat/long data of the cities we included
df <- read.csv("./data/city_latitude.csv")

# turn it into a shapefile
df_sf <- st_as_sf(df, 
                      coords = c("longitude", "latitude"), 
                      crs = 4326)

# and make a basic map plot
ggplot() +
  geom_sf(data=us_states, fill = "grey", color = "gray40", linewidth = 0.3, alpha = 0.5) +
  geom_sf(data=df_sf, size = 4, colour ="black") +
  theme_void() 
