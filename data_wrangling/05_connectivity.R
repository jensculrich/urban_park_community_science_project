library(sf)
library(terra)
library(tidyverse)
library(lconnect)
library(landscapemetrics)

city<-"des_moines"
buffer_size<-50
# Load the park shape file
parks_data <- readRDS(paste0("E:/phd_study/urban_park_community_science_project/data/parks/", buffer_size, "m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))
# Function required to calculate the connectivity metrics
# FUNCTION 1: DISTANCE-BASED metrics within buffers
calc_avg_dist_in_buffer <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000)) {
  for (buffer_dist in buffer_dists) {
    col_name <- paste0("avg_dist_", buffer_dist, "m")
    parks_sf[[col_name]] <- NA_real_
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
# FUNCTION 2: IFM connectivity within BUFFERS
calc_ifm_connectivity_in_buffers <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000),
                                             alpha = 0.005, b = 0.3, c = 0.3) {
  # Get areas for IFM calculation - UPDATED TO USE CORRECT COLUMN
  areas <- parks_sf$total_area_sqm
  for (buffer_dist in buffer_dists) {
    # Column name for IFM connectivity within this buffer
    ifm_conn_col <- paste0("ifm_connectivity_", buffer_dist, "m")
    # Initialize column
    parks_sf[[ifm_conn_col]] <- NA_real_
    for (i in 1:nrow(parks_sf)) {
      # Create buffer around current park
      buffer <- st_buffer(parks_sf[i,], buffer_dist)
      # Find parks within buffer (excluding current park)
      in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
      in_buffer[i] <- FALSE  # Exclude the current park
      # If parks found in buffer, calculate IFM connectivity
      if (sum(in_buffer) > 0) {
        # Focal patch area term: A_i^c
        focal_area_term <- areas[i]^c
        # Calculate connectivity from all source patches within buffer
        source_connectivity_sum <- 0
        neighbor_indices <- which(in_buffer)
        for (j in neighbor_indices) {
          # Distance decay: exp(-alpha * d_ij)
          distance_decay <- exp(-alpha * dist_matrix_numeric[i, j])
          # Source area scaling: A_j^b
          source_area_term <- areas[j]^b
          # Add to sum: exp(-alpha * d_ij) * A_j^b
          source_connectivity_sum <- source_connectivity_sum +
            (distance_decay * source_area_term)
        }
        # Complete IFM formula: I_i = A_i^c * sum
        parks_sf[[ifm_conn_col]][i] <- focal_area_term * source_connectivity_sum
      }
    }
  }
  return(parks_sf)
}
# FUNCTION 3: GLOBAL IFM connectivity (no distance limit)
calc_global_ifm_connectivity <- function(parks_sf, alpha = 0.005, b = 0.3, c = 0.3) {
  # UPDATED TO USE CORRECT COLUMN
  areas <- parks_sf$total_area_sqm
  n_patches <- nrow(parks_sf)
  connectivity <- numeric(n_patches)
  for (i in 1:n_patches) {
    # Focal patch area term: A_i^c
    focal_area_term <- areas[i]^c
    # Sum connectivity from ALL other patches
    source_connectivity_sum <- 0
    for (j in 1:n_patches) {
      if (i != j) {  # exclude self
        # Distance decay: exp(-alpha * d_ij)
        distance_decay <- exp(-alpha * dist_matrix_numeric[i, j])
        # Source area scaling: A_j^b
        source_area_term <- areas[j]^b
        # Add to sum
        source_connectivity_sum <- source_connectivity_sum +
          (distance_decay * source_area_term)
      }
    }
    # Complete IFM formula: I_i = A_i^c * sum
    connectivity[i] <- focal_area_term * source_connectivity_sum
  }
  parks_sf$ifm_connectivity_global <- connectivity
  return(parks_sf)
}
# FUNCTION 4: Simple area-weighted average distance within buffers
calc_area_weighted_avg_dist <- function(parks_sf, buffer_dists = c(100, 500, 1000, 2000)) {
  areas <- parks_sf$total_area_sqm
  for (buffer_dist in buffer_dists) {
    col_name <- paste0("area_weighted_avg_dist_", buffer_dist, "m")
    parks_sf[[col_name]] <- NA_real_
    for (i in 1:nrow(parks_sf)) {
      buffer <- st_buffer(parks_sf[i,], buffer_dist)
      in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
      in_buffer[i] <- FALSE
      if (sum(in_buffer) > 0) {
        distances <- dist_matrix_numeric[i, in_buffer]
        neighbor_areas <- areas[in_buffer]
        # Area-weighted average: Σ(distance_i × area_i) / Σ(area_i)
        parks_sf[[col_name]][i] <- sum(distances * neighbor_areas) / sum(neighbor_areas)
      }
    }
  }
  return(parks_sf)
}
# FUNCTION 5: Simple area-weighted average distance (no distance limit)
calc_global_area_weighted_avg_dist <- function(parks_sf) {
  areas <- parks_sf$total_area_sqm
  parks_sf$area_weighted_avg_dist_global <- NA_real_
  for (i in 1:nrow(parks_sf)) {
    # All other parks (no buffer limit)
    other_parks <- setdiff(1:nrow(parks_sf), i)
    if (length(other_parks) > 0) {
      distances <- dist_matrix_numeric[i, other_parks]
      neighbor_areas <- areas[other_parks]
      # Global area-weighted average
      parks_sf$area_weighted_avg_dist_global[i] <-
        sum(distances * neighbor_areas) / sum(neighbor_areas)
    }
  }
  return(parks_sf)
}
# FUNCTION 6: Simple area to distance ratio (no distance limit)
connectivity_fun <- function(parks_sf) {
  areas <- parks_sf$total_area_sqm
  parks_sf$connectivity <- NA_real_
  for (i in 1:nrow(parks_sf)) {
    # All other parks (no buffer limit)
    other_parks <- setdiff(1:nrow(parks_sf), i)
    if (length(other_parks) > 0) {
      distances <- dist_matrix_numeric[i, other_parks]
      neighbor_areas <- areas[other_parks]
      # Global area-weighted average
      parks_sf$connectivity[i] <-
        sum((log(neighbor_areas + 1) / (distances+1)))
    }
  }
  return(parks_sf)
}
#handling the 2km buffer surrounding all the classified parks
# Create a 2km buffer around classified parks
classified_parks <- parks_data[parks_data$type == "classified", ]
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km
# Combine all buffers into a single multipolygon
buffer_combined <- st_union(buffer_2km)
# Select unclassified parks
unclassified_parks <- parks_data[parks_data$type == "unclassified", ]
# Clip unclassified parks to only include areas within the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)
#combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clipped)
ggplot(result_parks) +
  geom_sf(aes(fill = type))
# Calculate distance matrix between all parks
#dist_matrix <- st_distance(result_parks)
dist_matrix <- st_distance(st_centroid(result_parks)) #distance is calculated based on centroid
dist_matrix_numeric <- units::drop_units(dist_matrix)
#distance-based metrics
parks_with_connectivity <- connectivity_fun(result_parks)
# Create final connectivity dataframe
final_connectivity_df <- parks_with_connectivity%>%
  st_drop_geometry() %>%
  filter(type == "classified") %>%
  dplyr::select(
    new_id, total_area_sqm, connectivity)
# View results
head(final_connectivity_df)
summary(final_connectivity_df)
View(final_connectivity_df)
write.csv(final_connectivity_df, paste0("data/final_merged_data/04_", buffer_size , "m_", city, "_connectivity.csv"), row.names = FALSE)
