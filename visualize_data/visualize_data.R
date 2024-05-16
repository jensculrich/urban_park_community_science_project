library(tidyverse)
library(sf) # spatial data processing

#-----------------------------------------------------
# summarize the detection data

# first read the data 
# downloaded via https://www.citynaturechallenge.org/participating-cities
# first read the data 
df <- rbind(
  cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2020-los-angeles-county.csv"), year = 1),
  cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2021-los-angeles-county.csv"), year = 2),
  cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2022-los-angeles-county.csv"), year = 3),
  cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2023-los-angeles-county.csv"), year = 4),
  cbind(read.csv("./data/city_nature_challenge_2020_2024_los_angeles_county/city-nature-challenge-2024-los-angeles-county.csv"), year = 5)
)

# and perform some initial filters
df <- df %>%
  
  # for now, let's only look at detections with a species ID
  filter(str_count(scientific_name, "\\w+") == 2) %>%
  
  # for now, let's filter out geo obscured detections
  filter(geoprivacy != "obscured") %>%
  
  # and additionally, let's get rid of any points with huge coordinate uncertainty
  # say > 1000m for now
  filter(positional_accuracy <= 1000) %>%
  
  # for now, to speed up the inference, let's also filter to the more common species
  group_by(scientific_name) %>%
  add_tally() %>%
  filter(n >= min_species_detections) %>%
  dplyr::select(-n) %>%
  ungroup()

# how many species were detected?
nrow(species_names <- df %>%
  # group by species ID
  group_by(scientific_name) %>%
  add_tally() %>%
  # and take one record
  slice(1) %>%
  select(scientific_name, common_name, n))

#-----------------------------------------------------
# map the data on a spatial file

county_shp <- read_sf("./data/los_angeles_county_boundary_shapefile/County_Boundary.SHP")

# USA_Contiguous_Albers_Equal_Area_Conic
crs <- 5070
grid_size <- 2000 # grid size in metres (1000 = 1km)

county_shp <- county_shp  %>% 
  st_transform(., crs) %>% # USA_Contiguous_Albers_Equal_Area_Conic
  filter(!OBJECTID %in% c(1,2,3,5,6,7))

# view the state shapefile
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none") 

# create "grid_size" km grid over the area
grid <- st_make_grid(county_shp, cellsize = c(grid_size, grid_size)) %>% 
  st_sf(grid_id = 1:length(.))

clipped_grid <- st_intersection(grid, county_shp)

# view the grid on the polygons
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = clipped_grid, fill = 'transparent', lwd = 0.3) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none") + 
  ggtitle("2km x 2km grid cells")

# include the detection data on the map

# make the detection data a spatial file
(df_sf <- st_as_sf(df,
                   coords = c("longitude", "latitude"), 
                   crs = 4326))

# and then transform it to the crs
df_sf <- st_transform(df_sf, crs = crs) %>%
  st_join(clipped_grid, join = st_intersects) %>% as.data.frame %>%
  # filter out records from outside of the urban grid
  filter(!is.na(grid_id)) %>%
  # now rejoin the lat/long data for each point
  left_join(., dplyr::select(
    df, id, latitude, longitude), by="id") 
  
# reproject
(df_sf <- st_as_sf(df_sf,
                   coords = c("longitude", "latitude"), 
                   crs = 4326))

# transform crs for plotting
df_sf <- st_transform(df_sf, crs = crs)

# plot by species
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = clipped_grid, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = df_sf, aes(colour = scientific_name)) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="") + 
  ggtitle("Lepidoptera detections from the LA city-nature-challenge\n(2020-2024) (coloured by species)")

# plot by year
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = clipped_grid, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = df_sf, aes(colour = (as.factor(year+2019)))) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  #theme(legend.position="") + 
  ggtitle("Lepidoptera detections from the LA city-nature-challenge\n(2020-2024) (coloured by year)")


#-----------------------------------------------------
# for now, let's also just include only grid cells with 
# one or more detections. This may be important to correct 
# later because it's quite likely people looked for  
# butterflies in other places but just didn't see any.
# could add in all the sites by rejoining the clipped grid?

nrow(df_sf %>%
       group_by(grid_id) %>%
       slice(1))
