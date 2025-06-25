library(sf)
library(terra)
library(tidyverse)
library(lconnect)
library(landscapemetrics)

# Load the LA park shape file
parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_LA.rds")

#handling the 2km buffer surrounding all the classified parks

# Create a 2km buffer around classified parks
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km

# Combine all buffers into a single multipolygon
buffer_combined <- st_union(buffer_2km)

# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)

#check overkapping between classified and unclassified parks
#id_to_be_removed<-st_intersection(unclassified_parks_clipped, classified_parks)%>%
#  dplyr::select(new_id)%>%
#  pull(new_id)

#remove unclassified parks that overlaps with the classified ones, we dont want to include this in the connectivity calculation
#unclassified_parks_clean <- unclassified_parks_clipped %>%
#  filter( !new_id %in% id_to_be_removed)

#combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clipped)

ggplot(result_parks) +
  geom_sf(aes(fill = type))

# Calculate distance matrix between all parks
dist_matrix <- st_distance(result_parks)
dist_matrix_numeric <- units::drop_units(dist_matrix)

# Calculate nearest neighbor distance for each park
diag(dist_matrix_numeric) <- Inf  # Exclude self-distances
result_parks$nearest_dist <- apply(dist_matrix_numeric, 1, min)

# Calculate average distance to all other parks
result_parks$avg_dist_all <- rowMeans(replace(dist_matrix_numeric, 
                                              cbind(1:nrow(dist_matrix_numeric), 
                                                    1:nrow(dist_matrix_numeric)), 
                                              NA), na.rm = TRUE)

# Function to calculate average distance within specific buffers
calc_avg_dist_in_buffer <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000)) {
  for (buffer_dist in buffer_dists) {
    col_name <- paste0("avg_dist_", buffer_dist, "m")
    parks_sf[[col_name]] <- NA_real_ #NA_real_ explicitly creates a numeric/double type NA value
    
    for (i in 1:nrow(parks_sf)) {
      # Create buffer around current park
      buffer <- st_buffer(parks_sf[i,], buffer_dist)
      
      # Find parks within buffer (excluding current park)
      in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
      in_buffer[i] <- FALSE  # Exclude the current park
      
      # If parks found in buffer, calculate average distance
      if (sum(in_buffer) > 0) {
        parks_sf[[col_name]][i] <- mean(dist_matrix_numeric[i, in_buffer])
      }
    }
  }
  return(parks_sf)
}

# Apply the function to calculate all buffer distances
parks_with_metrics <- calc_avg_dist_in_buffer(result_parks)

final_connectivity_df<-parks_with_metrics%>%
  st_drop_geometry()%>%
  filter(type=="classified")%>%
  dplyr::select(new_id, nearest_dist, avg_dist_all, avg_dist_100m, avg_dist_500m, avg_dist_1000m, avg_dist_2000m)

write.csv(final_connectivity_df, "data/final_merged_data/04_50m_LA_connectivity.csv")



#New York
# Load the NY park shape file
parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_NY.rds")
#handling the 2km buffer surrounding all the classified parks

# Create a 2km buffer around classified parks
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km

# Combine all buffers into a single multipolygon
buffer_combined <- st_union(buffer_2km)

# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)

#check overkapping between classified and unclassified parks
st_intersection(unclassified_parks_clipped, classified_parks)

#combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clipped)

#Check the visualization
ggplot(result_parks) +
  geom_sf(aes(fill = type))


# Calculate distance matrix between all parks
dist_matrix <- st_distance(result_parks)
dist_matrix_numeric <- units::drop_units(dist_matrix)

# Calculate nearest neighbor distance for each park
diag(dist_matrix_numeric) <- Inf  # Exclude self-distances
result_parks$nearest_dist <- apply(dist_matrix_numeric, 1, min)

# Calculate average distance to all other parks
result_parks$avg_dist_all <- rowMeans(replace(dist_matrix_numeric, 
                                              cbind(1:nrow(dist_matrix_numeric), 
                                                    1:nrow(dist_matrix_numeric)), 
                                              NA), na.rm = TRUE)

# Function to calculate average distance within specific buffers
calc_avg_dist_in_buffer <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000)) {
  for (buffer_dist in buffer_dists) {
    col_name <- paste0("avg_dist_", buffer_dist, "m")
    parks_sf[[col_name]] <- NA_real_ #NA_real_ explicitly creates a numeric/double type NA value
    
    for (i in 1:nrow(parks_sf)) {
      # Create buffer around current park
      buffer <- st_buffer(parks_sf[i,], buffer_dist)
      
      # Find parks within buffer (excluding current park)
      in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
      in_buffer[i] <- FALSE  # Exclude the current park
      
      # If parks found in buffer, calculate average distance
      if (sum(in_buffer) > 0) {
        parks_sf[[col_name]][i] <- mean(dist_matrix_numeric[i, in_buffer])
      }
    }
  }
  return(parks_sf)
}

# Apply the function to calculate all buffer distances
parks_with_metrics <- calc_avg_dist_in_buffer(result_parks)

final_connectivity_df<-parks_with_metrics%>%
  st_drop_geometry()%>%
  filter(type=="classified")%>%
  dplyr::select(new_id, nearest_dist, avg_dist_all, avg_dist_100m, avg_dist_500m, avg_dist_1000m, avg_dist_2000m)

write.csv(final_connectivity_df, "data/final_merged_data/04_50m_NYC_connectivity.csv")


#Seattle
# Load the Seattle park shape file
parks_data <- readRDS("E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_seattle.rds")
#handling the 2km buffer surrounding all the classified parks

# Create a 2km buffer around classified parks
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km

# Combine all buffers into a single multipolygon
buffer_combined <- st_union(buffer_2km)

# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]

# Clip unclassified parks to only include areas within the buffer
# st_intersection will keep only the parts that overlap with the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)

#check overkapping between classified and unclassified parks
st_intersection(unclassified_parks_clipped, classified_parks)

#combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clipped)

#Check the visualization
ggplot(result_parks) +
  geom_sf(aes(fill = type))


# Calculate distance matrix between all parks
dist_matrix <- st_distance(result_parks)
dist_matrix_numeric <- units::drop_units(dist_matrix)

# Calculate nearest neighbor distance for each park
diag(dist_matrix_numeric) <- Inf  # Exclude self-distances
result_parks$nearest_dist <- apply(dist_matrix_numeric, 1, min)

# Calculate average distance to all other parks
result_parks$avg_dist_all <- rowMeans(replace(dist_matrix_numeric, 
                                              cbind(1:nrow(dist_matrix_numeric), 
                                                    1:nrow(dist_matrix_numeric)), 
                                              NA), na.rm = TRUE)

# Function to calculate average distance within specific buffers
calc_avg_dist_in_buffer <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000)) {
  for (buffer_dist in buffer_dists) {
    col_name <- paste0("avg_dist_", buffer_dist, "m")
    parks_sf[[col_name]] <- NA_real_ #NA_real_ explicitly creates a numeric/double type NA value
    
    for (i in 1:nrow(parks_sf)) {
      # Create buffer around current park
      buffer <- st_buffer(parks_sf[i,], buffer_dist)
      
      # Find parks within buffer (excluding current park)
      in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
      in_buffer[i] <- FALSE  # Exclude the current park
      
      # If parks found in buffer, calculate average distance
      if (sum(in_buffer) > 0) {
        parks_sf[[col_name]][i] <- mean(dist_matrix_numeric[i, in_buffer])
      }
    }
  }
  return(parks_sf)
}

# Apply the function to calculate all buffer distances
parks_with_metrics <- calc_avg_dist_in_buffer(result_parks)

final_connectivity_df<-parks_with_metrics%>%
  st_drop_geometry()%>%
  filter(type=="classified")%>%
  dplyr::select(new_id, nearest_dist, avg_dist_all, avg_dist_100m, avg_dist_500m, avg_dist_1000m, avg_dist_2000m)

write.csv(final_connectivity_df, "data/final_merged_data/04_50m_Seattle_connectivity.csv")



