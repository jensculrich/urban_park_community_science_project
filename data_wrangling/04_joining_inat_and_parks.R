# Load required libraries
library(sf)
library(dplyr)
library(raster)
library(ggplot2)
library(prettymapr)


# Load the LA park shape file
la_parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_LA.rds")

ny_parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_NY.rds")

seattle_parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_Seattle.rds")

# Load the biodiversity data
leps_CA_data <- read.csv("data/inat_data/02_filtered_data/Leps/leps_CA_data_coordUncertainty_cleaned.csv")

leps_NY_data <- read.csv("data/inat_data/02_filtered_data/Leps/leps_NY_data_coordUncertainty_cleaned.csv")

leps_WA_data <- read.csv("data/inat_data/02_filtered_data/Leps/leps_WA_data_coordUncertainty_cleaned.csv")

# Convert the biodiversity data to spatial points 
LA_observation_coords <-st_as_sf(leps_CA_data, 
                              coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
                              st_transform(crs = st_crs(la_parks_data))

NY_observation_coords <-st_as_sf(leps_NY_data, 
                                 coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
                                st_transform(crs = st_crs(ny_parks_data))

WA_observation_coords <-st_as_sf(leps_WA_data, 
                                 coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
                                st_transform(crs = st_crs(seattle_parks_data))


####Creating the 2km buffer surrounding all the classified parks to form a small regional pool
##LA####
#Create a 2km buffer around classified parks
la_classified_parks <- la_parks_data[la_parks_data$type == "classified", ]
la_buffer_2km <- st_buffer(la_classified_parks, dist = 2000)  # 2000 meters = 2km
# Combine all buffers into a single multipolygon
la_buffer_combined <- st_union(la_buffer_2km)
# Select unclassified parks
la_unclassified_parks <- la_parks_data[la_parks_data$type == "unclassified", ]
# Clip unclassified parks to only include areas within the buffer for further connectivity calculation 
# st_intersection will keep only the parts that overlap with the buffer
la_unclassified_parks_clipped <- st_intersection(la_unclassified_parks, la_buffer_combined)
la_unclassified_parks_clean <- st_difference(la_unclassified_parks_clipped, st_union(la_classified_parks))
#Check to see if they are still overlapping
st_intersection(la_unclassified_parks_clean , la_classified_parks)
# combine both datasets back together
la_result_parks <- rbind(
  la_classified_parks,
  la_unclassified_parks_clean
)

# Visualize the result 
ggplot() +
  # Add the buffer outline
  geom_sf(data = la_buffer_combined, fill = NA, color = "red", size = 1) +
  # Add the classified parks
  geom_sf(data = la_classified_parks, fill = "darkgreen", color = NA) +
  # Add the clipped unclassified parks
  geom_sf(data = la_unclassified_parks_clipped, fill = "lightgreen", color = NA) +
  # Add a title and theme
  labs(title = "Parks Analysis",
       subtitle = "Classified parks (dark green) with 2km buffer (red outline)\nand clipped unclassified parks (light green)") +
  theme_minimal() 

ggplot() +
  geom_sf(data = la_result_parks, aes(fill=type))+
  theme_minimal() 


#join the clipped parks with the iNat Observation
#la_observations_with_parks_clipped <- st_join(LA_observation_coords, la_buffer_combined, join = st_intersects)
la_observations_with_parks <- st_join(LA_observation_coords, la_result_parks, join = st_intersects)

la_observations_with_parks_clipped <- st_intersection(la_observations_with_parks, la_buffer_combined)

la_observations_park_clipped_sf<-la_observations_with_parks_clipped

la_observations_park_clipped<-la_observations_with_parks_clipped%>%
  as.data.frame()

# See which values are duplicated
la_observations_park_clipped$gbifID[duplicated(la_observations_park_clipped$gbifID)]

ggplot(la_result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf(data= la_observations_park_clipped_sf, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type") +
  theme(legend.position = "right")


write.csv(la_observations_park_clipped,"E:/phd_study/urban_park_community_science_project/data/final_merged_data/01_50m_LA_observations_parkID_2km_clipped.csv")


#### create a large regional species pool with a 20-km buffer around the classified park
la_classified_parks <- la_parks_data[la_parks_data$type == "classified", ]
la_buffer_20km <- st_buffer(la_classified_parks, dist = 20000)  # 20000 meters = 20km
# Combine all buffers into a single multipolygon
la_buffer_20km_combined <- st_union(la_buffer_20km)
# Select unclassified parks
la_unclassified_parks <- la_parks_data[la_parks_data$type == "unclassified", ]

la_unclassified_parks_clipped <- st_intersection(la_unclassified_parks, la_buffer_20km_combined)

la_unclassified_parks_clean <- st_difference(la_unclassified_parks_clipped, st_union(la_classified_parks))
#Check to see if they are still overlapping
st_intersection(la_unclassified_parks_clean , la_classified_parks)

# Visualize the result

plot(st_geometry(la_buffer_20km_combined), border = "red", lwd = 2, main = "20-km regional species pool")

plot(st_geometry(la_classified_parks), col = "darkgreen", main = "20-km regional species pool", add = TRUE)

plot(st_geometry(la_unclassified_parks_clean), col = "orange", main = "20-km regional species pool", add = TRUE)

la_result_parks <- rbind(
  la_classified_parks,
  la_unclassified_parks_clean
)

la_observations_with_parks <- st_join(LA_observation_coords, la_result_parks, join = st_intersects)

la_observations_with_parks_clipped <- st_intersection(la_observations_with_parks, la_buffer_20km_combined)


ggplot() +
  #annotation_map_tile(type = "osm")+
geom_sf(data = la_parks_data, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= la_observations_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Regional Species Pool") +
  theme(legend.position = "right")

write.csv(as.data.frame(la_observations_with_parks_clipped), "data/final_merged_data/02_50m_LA_regional_species_pool.csv")


#New York----
# Create a 2km buffer around classified parks
ny_classified_parks <- ny_parks_data[ny_parks_data$type == "classified", ]
ny_buffer_2km <- st_buffer(ny_classified_parks, dist = 2000)  # 2000 meters = 2km

# Combine all buffers into a single multipolygon
ny_buffer_combined <- st_union(ny_buffer_2km)

# Select unclassified parks
ny_unclassified_parks <- ny_parks_data[ny_parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
ny_unclassified_parks_clipped <- st_intersection(ny_unclassified_parks, ny_buffer_combined)
ny_unclassified_parks_clean <- st_difference(ny_unclassified_parks_clipped, st_union(ny_classified_parks))


st_intersection(ny_unclassified_parks_clean , ny_classified_parks)
# If you want to combine both datasets back together

ny_result_parks <- rbind(
  ny_classified_parks,
  ny_unclassified_parks_clean
)

# Visualize the result (optional)
ggplot() +
  # Add the buffer outline
  geom_sf(data = ny_buffer_combined, fill = NA, color = "red", size = 1) +
  # Add the classified parks
  geom_sf(data = ny_classified_parks, fill = "darkgreen", color = NA) +
  # Add the clipped unclassified parks
  geom_sf(data = ny_unclassified_parks_clipped, fill = "orange", color = NA) +
  # Add a title and theme
  labs(title = "NY Parks Analysis",
       subtitle = "Classified parks (dark green) with 2km buffer (red outline)\nand clipped unclassified parks (orange)") +
  theme_minimal() 

ggplot() +
  # Add the classified parks
  geom_sf(data = ny_result_parks, aes(fill=type))+
  theme_minimal() 


#join the clipped parks with the iNat Observation
ny_observations_with_parks<- st_join(NY_observation_coords, ny_result_parks, join = st_intersects)
ny_observations_with_parks_clipped <- st_intersection(ny_observations_with_parks, ny_buffer_combined)


# See which values are duplicated
ny_observations_with_parks_clipped$gbifID[duplicated(ny_observations_with_parks_clipped$gbifID)]


ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data=ny_result_parks, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= ny_observations_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "NYC Parks by Type (observations filtered out by the 2km-buffer around classified parks") +
  theme(legend.position = "right")

write.csv(as.data.frame(ny_observations_with_parks_clipped),"E:/phd_study/urban_park_community_science_project/data/final_merged_data/01_50m_NYC_observations_parkID_2km_clipped.csv")

#### create a regional species pool with a 20-km buffer around the NYC classified park
ny_classified_parks <- ny_parks_data[ny_parks_data$type == "classified", ]
ny_buffer_20km <- st_buffer(ny_classified_parks, dist = 20000)  # 20000 meters = 20km
# Combine all buffers into a single multipolygon
ny_buffer_20km_combined <- st_union(ny_buffer_20km)
# Select unclassified parks
ny_unclassified_parks <- ny_parks_data[ny_parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
ny_unclassified_parks_clipped <- st_intersection(ny_unclassified_parks, ny_buffer_20km_combined)
ny_unclassified_parks_clean <- st_difference(ny_unclassified_parks_clipped, st_union(ny_classified_parks))

st_intersection(ny_unclassified_parks_clean , ny_classified_parks)
#combine both datasets back together

ny_result_parks <- rbind(
  ny_classified_parks,
  ny_unclassified_parks_clean
)


# Visualize the result
ggplot(ny_result_parks) +
  geom_sf(aes(fill = type)) +
  #geom_sf(data= observations_park_sf, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type",
  fill = "New ID") +
  theme(legend.position = "right")


ny_observations_within_20km_buffer <- st_join(NY_observation_coords, ny_result_parks, join = st_intersects)

ny_observations_within_20km_buffer_clipped <- st_intersection(ny_observations_within_20km_buffer, ny_buffer_20km_combined)

ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data = ny_result_parks, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= ny_observations_within_20km_buffer_clipped, color ="red", size=0.01) +
    scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "NYC Regional Species Pool") +
  theme(legend.position = "right")

# See which values are duplicated
ny_observations_with_parks_clipped$gbifID[duplicated(ny_observations_with_parks_clipped$gbifID)]


write.csv(as.data.frame(ny_observations_within_20km_buffer_clipped), "data/final_merged_data/02_NYC_regional_species_pool.csv")


#Seattle----
# Create a 2km buffer around classified parks
seattle_classified_parks <- seattle_parks_data[seattle_parks_data$type == "classified", ]
seattle_buffer_2km <- st_buffer(seattle_classified_parks, dist = 2000)  # 2000 meters = 2km

# Combine all buffers into a single multipolygon
seattle_buffer_combined <- st_union(seattle_buffer_2km)

# Select unclassified parks
seattle_unclassified_parks <- seattle_parks_data[seattle_parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
seattle_unclassified_parks_clipped <- st_intersection(seattle_unclassified_parks, seattle_buffer_combined)
seattle_unclassified_parks_clean <- st_difference(seattle_unclassified_parks_clipped, st_union(seattle_classified_parks))


st_intersection(seattle_unclassified_parks_clean , seattle_classified_parks)
# If you want to combine both datasets back together

seattle_result_parks <- rbind(
  seattle_classified_parks,
  seattle_unclassified_parks_clean
)

# Visualize the result (optional)
ggplot() +
  # Add the buffer outline
  geom_sf(data = seattle_buffer_combined, fill = NA, color = "red", size = 1) +
  # Add the classified parks
  geom_sf(data = seattle_classified_parks, fill = "darkgreen", color = NA) +
  # Add the clipped unclassified parks
  geom_sf(data = seattle_unclassified_parks_clipped, fill = "orange", color = NA) +
  # Add a title and theme
  labs(title = "Seattle Parks Analysis",
       subtitle = "Classified parks (dark green) with 2km buffer (red outline)\nand clipped unclassified parks (orange)") +
  theme_minimal() 

ggplot() +
  # Add the classified parks
  geom_sf(data = seattle_result_parks, aes(fill=type))+
  theme_minimal() 


#join the clipped parks with the iNat Observation
seattle_observations_with_parks<- st_join(WA_observation_coords, seattle_result_parks, join = st_intersects)
seattle_observations_with_parks_clipped <- st_intersection(seattle_observations_with_parks, seattle_buffer_combined)


# See which values are duplicated
seattle_observations_with_parks_clipped$gbifID[duplicated(seattle_observations_with_parks_clipped$gbifID)]


ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data=seattle_result_parks, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= seattle_observations_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "NYC Parks by Type (observations filtered out by the 2km-buffer around classified parks") +
  theme(legend.position = "right")

write.csv(as.data.frame(seattle_observations_with_parks_clipped),"E:/phd_study/urban_park_community_science_project/data/final_merged_data/01_50m_Seattle_observations_parkID_2km_clipped.csv")

#### create a regional species pool with a 20-km buffer around the NYC classified park
seattle_classified_parks <- seattle_parks_data[ny_parks_data$type == "classified", ]
seattle_buffer_20km <- st_buffer(seattle_classified_parks, dist = 20000)  # 20000 meters = 20km
# Combine all buffers into a single multipolygon
seattle_buffer_20km_combined <- st_union(seattle_buffer_20km)
# Select unclassified parks
seattle_unclassified_parks <- seattle_parks_data[seattle_parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
seattle_unclassified_parks_clipped <- st_intersection(seattle_unclassified_parks, seattle_buffer_20km_combined)
seattle_unclassified_parks_clean <- st_difference(seattle_unclassified_parks_clipped, st_union(seattle_classified_parks))

st_intersection(seattle_unclassified_parks_clean , seattle_classified_parks)
#combine both datasets back together

seattle_result_parks <- rbind(
  seattle_classified_parks,
  seattle_unclassified_parks_clean
)


# Visualize the result
ggplot(seattle_result_parks) +
  geom_sf(aes(fill = type)) +
  #geom_sf(data= observations_park_sf, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type",
       fill = "New ID") +
  theme(legend.position = "right")


seattle_observations_within_20km_buffer <- st_join(WA_observation_coords, seattle_result_parks, join = st_intersects)

seattle_observations_within_20km_buffer_clipped <- st_intersection(seattle_observations_within_20km_buffer, seattle_buffer_20km_combined)

ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data = seattle_result_parks, aes(fill = type), inherit.aes = FALSE) +
  geom_sf(data= seattle_observations_within_20km_buffer_clipped, color ="red", size=0.01) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "NYC Regional Species Pool") +
  theme(legend.position = "right")

write.csv(as.data.frame(seattle_observations_within_20km_buffer_clipped), "data/final_merged_data/02_50m_Seattle_regional_species_pool.csv")


#### Load the flowering plants data
#Los Angeles
flowers_CA_df<-read.csv("data/inat_data/02_filtered_data/Plants/inat_flowers_CA.csv")

ggplot(la_parks_data) +
  geom_sf(aes(fill = type)) +
  #geom_sf(data= observations_park_clipped_sf, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type") +
  theme(legend.position = "right")


#Convert the flower observations into spatial points

flowers_CA_coords <-st_as_sf(flowers_CA_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs = st_crs(la_parks_data))

#join the clipped parks with the flowers Observation
flowers_with_parks_clipped <- st_join(flowers_CA_coords, la_parks_data%>%filter(type == "classified"), join = st_intersects)

flowers_with_parks_clipped<-flowers_with_parks_clipped%>%
  filter(!is.na(new_id))
  
  ggplot(la_result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf(data= flowers_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type with iNat flowering plant observations") +
  theme(legend.position = "right")
  
  
  total_no_genus_LA_flowers<-flowers_with_parks_clipped%>%
    group_by(new_id)%>%
    summarize(no_of_genus=length(unique(genus)))%>%
    ungroup() 
  
  ggplot(la_result_parks) +
    geom_sf(aes(fill = type)) +
    geom_sf_text(data = total_no_genus_LA_flowers, aes(label = no_of_genus), color="purple", size = 3) +scale_fill_viridis_d() +  # Use viridis color palette
    theme_minimal() +
    labs(title = "Parks by Type with Flowering Plant Genus Counts (Whole studying periods)") +
    theme(legend.position = "right")

write.csv(flowers_with_parks_clipped, "E:/phd_study/urban_park_community_science_project/data/final_merged_data/03_50m_flowers_LA_classified_park.csv")  

#New York
flowers_NY_df<-read.csv("data/inat_data/02_filtered_data/Plants/inat_flowers_NY.csv")
  
#Convert the flower observations into spatial points
  
  flowers_NY_coords <-st_as_sf(flowers_NY_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs = st_crs(ny_parks_data))
  
  #join the clipped parks with the flowers Observation
  flowers_with_parks_clipped <- st_join(flowers_NY_coords, ny_result_parks%>%filter(type == "classified"), join = st_intersects)
  
  flowers_with_parks_clipped<-flowers_with_parks_clipped%>%
    filter(!is.na(new_id))
  
  ggplot(ny_result_parks) +
    geom_sf(aes(fill = type)) +
    geom_sf(data= flowers_with_parks_clipped, color ="red", size=0.1) +
    scale_fill_viridis_d() +  # Use viridis color palette
    theme_minimal() +
    labs(title = "Parks by Type with iNat flowering plant observations") +
    theme(legend.position = "right")
  
  total_no_genus_NY_flowers<-flowers_with_parks_clipped%>%
    group_by(new_id)%>%
    summarize(no_of_genus=length(unique(genus)))%>%
    ungroup() 
  
  ggplot(ny_result_parks) +
    geom_sf(aes(fill = type)) +
    geom_sf_text(data = total_no_genus_NY_flowers, aes(label = no_of_genus), color="red", size = 3) +scale_fill_viridis_d() +  # Use viridis color palette
    theme_minimal() +
    labs(title = "NYC Parks by Type with Flowering Plant Genus Counts (Whole studying periods)") +
    theme(legend.position = "right")
  
write.csv(flowers_with_parks_clipped, "E:/phd_study/urban_park_community_science_project/data/final_merged_data/03_50m_flowers_NY_classified_park.csv")  
  
#Washington  
flowers_WA_df<-read.csv("data/inat_data/02_filtered_data/Plants/inat_flowers_WA.csv")

#Convert the flower observations into spatial points

flowers_WA_coords <-st_as_sf(flowers_WA_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(crs = st_crs(seattle_parks_data))

#join the clipped parks with the flowers Observation
flowers_with_parks_clipped <- st_join(flowers_WA_coords, seattle_result_parks%>%filter(type == "classified"), join = st_intersects)

flowers_with_parks_clipped<-flowers_with_parks_clipped%>%
  filter(!is.na(new_id))

ggplot(seattle_result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf(data= flowers_with_parks_clipped, color ="red", size=0.1) +
  scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "Parks by Type with iNat flowering plant observations") +
  theme(legend.position = "right")

total_no_genus_seattle_flowers<-flowers_with_parks_clipped%>%
  group_by(new_id)%>%
  summarize(no_of_genus=length(unique(genus)))%>%
  ungroup() 

ggplot(seattle_result_parks) +
  geom_sf(aes(fill = type)) +
  geom_sf_text(data = total_no_genus_seattle_flowers, aes(label = no_of_genus), color="red", size = 3) +scale_fill_viridis_d() +  # Use viridis color palette
  theme_minimal() +
  labs(title = "NYC Parks by Type with Flowering Plant Genus Counts (Whole studying periods)") +
  theme(legend.position = "right")

write.csv(flowers_with_parks_clipped, "E:/phd_study/urban_park_community_science_project/data/final_merged_data/03_50m_flowers_Seattle_classified_park.csv")  
  
  
  
  
  
  
  
  
  
  
  
  
  # species_not_match_inat<-observations_with_parks%>%
  #   as.data.frame()%>%
  #   filter(!is.na(new_id))%>%
  #   filter(species %in% as.matrix(observations_per_species_final$species))%>%
  #   mutate(species_match=(species==verbatimScientificName))%>%
  #   filter(species_match==FALSE)%>%
  #   select(species, verbatimScientificName)%>%
  #   unique()%>%
  #   rename("gbif_name"=species, "inat_name"=verbatimScientificName)%>%
  #   arrange(gbif_name)
  
  # write.csv(species_not_match_inat,"E:/phd_study/urban_park_community_science_project/data/final_merged_data/flagged_species_taxonomy.csv")
  