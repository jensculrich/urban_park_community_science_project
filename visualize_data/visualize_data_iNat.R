library(tidyverse)
library(sf) # spatial data processing

min_species_detections = 20
grid_size = 2000 # meters
min_park_size_acres = 10 # acres
max_park_size_acres = 500 # acres
buffer_distance <- 250 # meters

#-----------------------------------------------------
# summarize the detection data

# first read the data 
# downloaded via https://www.citynaturechallenge.org/participating-cities
# first read the data 
df <- read.csv(
  "./data/all_inat_lep_records_2020-2023_los_angeles_county/all_inat_lep_records_2020-2023_los_angeles_county.csv")
# 29,276 records

# and perform some initial filters
df <- df %>%
  
  # for now, let's only look at detections with a species ID
  filter(str_count(species, "\\w+") == 2) %>%
  
  # for now, let's filter out geo obscured detections
  #filter(geoprivacy != "obscured") %>%
  
  # and additionally, let's get rid of any points with huge coordinate uncertainty
  # say > 1000m for now
  filter(coordinateUncertaintyInMeters <= 1000) %>%
  filter(!is.na(coordinateUncertaintyInMeters)) %>%
  
  # for now, to speed up the inference, let's also filter to the more common species
  group_by(species) %>%
  add_tally() %>%
  filter(n >= min_species_detections) %>%
  dplyr::select(-n) %>%
  ungroup()
# 21,417 records

# how many species were detected?
nrow(species_names <- df %>%
  # group by species ID
  group_by(species) %>%
  add_tally() %>%
  # and take one record
  slice(1) %>%
  select(species, family, n))

#-----------------------------------------------------
# map the data on a spatial file

county_shp <- read_sf("./data/los_angeles_county_boundary_shapefile/County_Boundary.SHP")

# USA_Contiguous_Albers_Equal_Area_Conic
crs <- 5070
grid_size <- 2000 # grid size in metres (1000 = 1km)

county_shp <- county_shp  %>% 
  st_transform(., crs) %>% # USA_Contiguous_Albers_Equal_Area_Conic
  filter(!OBJECTID %in% c(1,2,3,5,6,7)) # filter out catalina islands and some other random disparate chunks

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
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# and then transform it to the crs
df_sf <- st_transform(df_sf, crs = crs) %>%
  st_join(clipped_grid, join = st_intersects) %>% as.data.frame %>%
  # filter out records from outside of the urban grid
  filter(!is.na(grid_id)) %>%
  # now rejoin the lat/long data for each point
  left_join(., dplyr::select(
    df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 
  
# reproject
(df_sf <- st_as_sf(df_sf,
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# transform crs for plotting
df_sf <- st_transform(df_sf, crs = crs)

# plot by species
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = clipped_grid, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = df_sf, aes(colour = species)) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="") + 
  ggtitle("Lepidoptera detections from iNat\n(Los Angeles, 2020-2023) (coloured by species)")

# plot by year
ggplot() +
  geom_sf(data = county_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = clipped_grid, fill = 'transparent', lwd = 0.3) +
  geom_sf(data = df_sf, aes(colour = (as.factor(year)))) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  #theme(legend.position="") + 
  ggtitle("Lepidoptera detections from iNat\n(Los Angeles, 2020-2023) (coloured by year)")



#-----------------------------------------------------
# map the data on a spatial file of urban parks

parks_shp <- read_sf("./data/los_angeles_county_parks_shapefile/Regional_Site_Inventory.shp")

# USA_Contiguous_Albers_Equal_Area_Conic
crs <- 5070

parks_shp <- parks_shp  %>% 
  st_transform(., crs) # USA_Contiguous_Albers_Equal_Area_Conic

# view the state shapefile
ggplot() +
  geom_sf(data = parks_shp, fill = 'white', lwd = 0.05) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="none") 

# for now, I filtered out some things that were clearly not urban parks 
parks_shp <- parks_shp  %>% 
  filter(!str_detect(PARK_NAME, "Angeles National Forest")) %>% 
  filter(!str_detect(PARK_NAME, "Los Padres National Forest")) %>% 
  filter(!str_detect(PARK_NAME, "Edwards AFB")) %>%
  filter(!str_detect(PARK_NAME, "Hungry Valley")) %>%
  filter(!str_detect(PARK_NAME, "Air Force")) %>%
  filter(!str_detect(PARK_NAME, "State Recreation Area")) %>%
  filter(!str_detect(PARK_NAME, "National Recreation")) %>%
  filter(!str_detect(PARK_NAME, "State Park")) %>%
  filter(!str_detect(PARK_NAME, "State Beach")) %>%
  filter(!str_detect(PARK_NAME, "County Beach")) %>%
  filter(!str_detect(PARK_NAME, "Santa Catalina Island")) %>%
  
  # for now I also filtered out small parks just to speed up the estimation times (fewer sites)
  filter(RRE_ACRES > min_park_size_acres) %>%
  # and also the really really big ones 
  filter(RRE_ACRES < max_park_size_acres)

# let's add a buffer around each park and then merge parks that are touching or overlapping
parks_shp <- st_buffer(parks_shp, buffer_distance)

# include the detection data on the map

# make the detection data a spatial file
(df_sf <- st_as_sf(df,
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# and then transform it to the crs
df_sf <- st_transform(df_sf, crs = crs) %>%
  st_join(parks_shp, join = st_intersects) %>% as.data.frame %>%
  # filter out records from outside of the urban grid
  filter(!is.na(PARK_NAME)) %>%
  # now rejoin the lat/long data for each point
  left_join(., dplyr::select(
    df, gbifID, decimalLatitude, decimalLongitude), by="gbifID") 

# reproject
(df_sf <- st_as_sf(df_sf,
                   coords = c("decimalLongitude", "decimalLatitude"), 
                   crs = 4326))

# transform crs for plotting
df_sf <- st_transform(df_sf, crs = crs)

# plot by species
ggplot() +
  geom_sf(data = parks_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = df_sf, aes(colour = species)) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  theme(legend.position="") + 
  ggtitle("Lepidoptera detections from iNat\n(Los Angeles, 2020-2023) (coloured by species)")

# plot by year
ggplot() +
  geom_sf(data = parks_shp, fill = 'white', lwd = 0.05) +
  geom_sf(data = df_sf, aes(colour = (as.factor(year)))) +
  coord_sf(datum = NA)  +
  labs(x = "") +
  labs(y = "") +
  #theme(legend.position="") + 
  ggtitle("Lepidoptera detections from iNat\n(Los Angeles, 2020-2023) (coloured by year)")


#-----------------------------------------------------
# for now, let's also just include only grid cells with 
# one or more detections. This may be important to correct 
# later because it's quite likely people looked for  
# butterflies in other places but just didn't see any.
# could add in all the sites by rejoining the clipped grid?

nrow(df_sf %>%
       group_by(PARK_NAME) %>%
       slice(1))
