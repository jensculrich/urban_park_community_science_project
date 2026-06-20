# Load required libraries
library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(prettymapr)
library(cowplot)


# Load the LA park shape file (LA, NY, Seattle, dalla, houston, sf, riverside, sd)
city <-"st_louis"
state <- "MO"
buffer_size<-0
parks_data <- readRDS(paste0("E:/phd_study/urban_park_community_science_project/data/parks/", buffer_size , "m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))

#need to check if parks_data is in meters, denton is in ft
st_crs(parks_data)

#if it is in ft, then use this function to fix it
# Simple function: transform to metric UTM based on location
# to_meters <- function(sf_object) {
#   center <- st_transform(st_centroid(st_union(sf_object)), 4326)
#   coords <- st_coordinates(center)
#   utm_zone <- floor((coords[1] + 180) / 6) + 1
#   epsg <- ifelse(coords[2] >= 0, 32600 + utm_zone, 32700 + utm_zone)
#   st_transform(sf_object, epsg)
# }

#parks_data <- to_meters(parks_data)

st_area(parks_data[1,])

# Load the biodiversity data
leps_data <- read.csv(paste0("data/inat_data/02_filtered_data/Leps/leps_", state, "_data_coordUncertainty_cleaned.csv"))
#For Detroit we also need data from Canada
#leps_data_essex <-read.csv("data/inat_data/02_filtered_data/Leps/leps_Essex_data_coordUncertainty_cleaned.csv")

#leps_data<- rbind(leps_data, leps_data_essex)

table(leps_data$stateProvince)
# Convert the biodiversity data to spatial points 
observation_coords <-st_as_sf(leps_data, 
                              coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
  st_transform(crs = st_crs(parks_data))


####Creating the 2km buffer surrounding all the classified parks to form a small regional pool
#Create a 2km buffer around classified parks
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km
# Combine all buffers into a single multipolygon
buffer_combined <- st_union(buffer_2km)
# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer for further connectivity calculation
# st_intersection will keep only the parts that overlap within the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)
#ensure the unclassified parks do not overlap with the classified ones
unclassified_parks_clean <- st_difference(unclassified_parks_clipped, st_buffer( st_union(classified_parks), 0.001))

#Check to see if they are still overlapping
st_intersection(unclassified_parks_clean , classified_parks)

# combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clean)

# Visualize the result 
ggplot() +
  # Add the buffer outline
  geom_sf(data = buffer_combined, fill = NA, color = "red", size = 1) +
  # Add the classified parks
  geom_sf(data = classified_parks, fill = "darkgreen", color = NA) +
  # Add the clipped unclassified parks
  geom_sf(data = unclassified_parks_clipped, fill = "lightgreen", color = NA) +
  # Add a title and theme
  labs(title = "Parks Analysis",
       subtitle = "Classified parks (dark green) with 2km buffer (red outline)\nand clipped unclassified parks (light green)") +
  theme_minimal() 

#Looks the combined data
ggplot() +
  geom_sf(data = result_parks, aes(fill=type))+
  theme_minimal() 

st_write(result_parks, paste0("E:/phd_study/urban_park_community_science_project/data/shapefile_", buffer_size, "m_buffered_park_2km_regional_pool/", city, "_", buffer_size, "_buffered_park_2km_regional_pool.shp"))


#join the clipped parks with the iNat Observation
#observations_with_parks <- st_join(observation_coords, result_parks, join = st_intersects) #this is a left join, which preserve the observations in the space between parks

# Create 50m detection range for each park
result_parks_buffered <- st_buffer(result_parks, dist = 50)


# Join
observations_with_parks <- st_join(observation_coords, result_parks_buffered, join = st_intersects)



#filter out observations by the 2km regional buffer
observations_with_parks_2km_clipped <- st_intersection(observations_with_parks, buffer_combined)%>%
  #group_by(gbifID) %>%
  #slice(1) %>%
  ungroup()


# See which values are duplicated
#observations_with_parks_2km_clipped$gbifID[duplicated(observations_with_parks_2km_clipped$gbifID)]

sum(duplicated(observations_with_parks_2km_clipped$gbifID))

ggplot(result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf(data= observations_with_parks_2km_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = paste0(city,", Parks by Type")) +
  theme(legend.position = "right")

nrow(observations_with_parks_2km_clipped)

as.data.frame(st_drop_geometry(observations_with_parks_2km_clipped))%>%filter(type=="classified")%>%dplyr::select(new_id)%>%unique()%>%pull%>%length()

View(as.data.frame(st_drop_geometry(observations_with_parks_2km_clipped)))

write.csv(as.data.frame(st_drop_geometry(observations_with_parks_2km_clipped%>%dplyr::select(gbifID:type))), paste0("E:/phd_study/urban_park_community_science_project/data/final_merged_data/01_", buffer_size, "m_", city, "_observations_parkID_2km_clipped.csv"), row.names = FALSE)


#### create a large regional species pool with a 20-km buffer around the classified park
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_20km <- st_buffer(classified_parks, dist = 20000)  # 20000 meters = 20km
# Combine all buffers into a single multipolygon
buffer_20km_combined <- st_union(buffer_20km)
# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]

unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_20km_combined)

unclassified_parks_clean <- st_difference(unclassified_parks_clipped, st_buffer( st_union(classified_parks), 0.001))
#Check to see if they are still overlapping
st_intersection(unclassified_parks_clean , classified_parks)

# Visualize the result

plot(st_geometry(buffer_20km_combined), border = "red", lwd = 2, main = "20-km regional species pool")

plot(st_geometry(classified_parks), col = "darkgreen", main = "20-km regional species pool", add = TRUE)

plot(st_geometry(unclassified_parks_clean), col = "orange", main = "20-km regional species pool", add = TRUE)

result_parks <- rbind(
  classified_parks,
  unclassified_parks_clean
)

# Create 50m detection range for each park
result_parks_buffered <- st_buffer(result_parks, dist = 50)


# Join
observations_with_parks <- st_join(observation_coords, result_parks_buffered, join = st_intersects)


observations_with_parks_20km_clipped <- st_intersection(observations_with_parks, buffer_20km_combined)%>%
  #group_by(gbifID) %>%
  #slice(1) %>%
  ungroup()

# See which values are duplicated
#observations_with_parks_20km_clipped$gbifID[duplicated(observations_with_parks_20km_clipped$gbifID)]
sum(duplicated(observations_with_parks_20km_clipped$gbifID))

#plot the result
ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data = result_parks%>%filter(type== "classified"), aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= observations_with_parks_20km_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "20km Regional Species Pool") +
  theme(legend.position = "right")+
  theme_cowplot()

id_lookup<-observations_with_parks_20km_clipped%>%
  st_drop_geometry()%>%
  filter(type== "classified")%>%
  group_by(new_id)%>%
  summarise(nrow=n())%>%
  slice_max(order_by = nrow, n=1, with_ties = FALSE)%>%
  dplyr::select(new_id)%>%
  pull()


ggplot()+
  geom_sf(data = result_parks, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data = result_parks%>%filter(new_id == id_lookup), inherit.aes = FALSE, fill = "green") +
  geom_sf(data= observations_with_parks_20km_clipped%>%filter(new_id == id_lookup), color ="red", size=1)

nrow(observations_with_parks_20km_clipped)


write.csv(as.data.frame(st_drop_geometry(observations_with_parks_20km_clipped)), paste0("data/final_merged_data/02_",buffer_size,"m_", city, "_regional_species_pool.csv"), row.names = FALSE)


#### Load the flowering plants data
flowers_df<-read.csv(paste0("data/inat_data/02_filtered_data/Plants/flowers_", state, "_data_coordUncertainty_cleaned.csv"))

table(flowers_df$stateProvince)
#Convert the flower observations into spatial points

flowers_coords <-st_as_sf(flowers_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs = st_crs(parks_data))

#join the clipped parks with the flowers Observation
flowers_with_parks_clipped <- st_join(flowers_coords, parks_data%>%filter(type == "classified"), join = st_intersects)

flowers_with_parks_clipped<-flowers_with_parks_clipped%>%
  filter(!is.na(new_id))

ggplot(result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf(data= flowers_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type with iNat flowering plant observations") +
  theme(legend.position = "right")

as.data.frame(st_drop_geometry(flowers_with_parks_clipped))%>%filter(type=="classified")%>%dplyr::select(new_id)%>%unique()%>%pull%>%length()


write.csv(as.data.frame(st_drop_geometry(flowers_with_parks_clipped)), paste0("E:/phd_study/urban_park_community_science_project/data/final_merged_data/03_", buffer_size , "m_flowers_", city, "_classified_park.csv"), row.names = FALSE)  


