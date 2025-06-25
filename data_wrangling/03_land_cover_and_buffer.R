## Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(daymetr)
library(exactextractr)
library(ggplot2)
library(units)
library(igraph)
## Load the park shapefile
parks <- st_read("E:/phd_study/urban_park_community_science_project/data/Parkserve_Shapefiles_05212024/ParkServe_Parks.shp")

# Filter for LA, Ny and WA County parks
parks_la <- parks %>% filter(Park_Count == "Los Angeles County")

parks_ny <- parks %>% filter(Park_State == "New York" | Park_State == "New Jersey") #New Jersey is added because it is very close to Manhattan

parks_wa <- parks %>% filter(Park_State == "Washington")

plot(st_geometry(parks_ny))
## Load the land cover raster using terra
# Land cover data with classified values 0-9 (0 = Unclassified, 2 = Grass/Shrub, 3 = Tree Canopy; )
land_cover_la <- rast("data/urbanwatch_data/03_classified_land_cover_data/la_classified_land_cover.tif")

land_cover_ny <- rast("data/urbanwatch_data/03_classified_land_cover_data/ny_classified_land_cover.tif")

land_cover_seattle <-rast("data/urbanwatch_data/03_classified_land_cover_data/seattle_classified_land_cover.tif")


# Reproject parks to match the CRS of the land cover raster if needed
parks_la_reproj <- st_transform(parks_la, crs = crs(land_cover_la))

parks_ny_reproj <- st_transform(parks_ny, crs = crs(land_cover_ny))

parks_wa_reproj <- st_transform(parks_wa, crs = crs(land_cover_seattle))

# Aggregate the land cover raster
calculate_mode <- function(x) {
  unique_vals <- unique(na.omit(x))  # Remove NA values and get unique values
  if (length(unique_vals) == 0) {
    return(NA)  # Return NA if no valid values
  }
  return(unique_vals[which.max(tabulate(match(x, unique_vals)))])  # Find the mode
}

# We use a factor of 10 to merge 10x10 cells and `modal` function to keep the most frequent category in each block
aggregated_land_cover <- aggregate(land_cover_la, fact = 10, fun = calculate_mode)
aggregated_land_cover_ny <- aggregate(land_cover_ny, fact = 10, fun = calculate_mode)
aggregated_land_cover_seattle <- aggregate(land_cover_seattle, fact = 10, fun = calculate_mode)

# Plot the aggregated raster to check the result
png(filename = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/plot/aggregated_land_cover_LA.png", width = 500, height = 800)
dev.off()

plot(land_cover_la, main = "Land Cover Raster ") # Grass/Shrub=3; Tree Canopy=4

plot(aggregated_land_cover, main = "Aggregated Land Cover Raster (10-meter resolution)")


# Save the aggregated raster as a GeoTIFF file
writeRaster(aggregated_land_cover, "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/processing_data/10_meters_aggregated_urbanwatch_data/aggregated_land_cover.tif", overwrite = TRUE)

writeRaster(aggregated_land_cover_ny, "data/processing_data/10_meters_aggregated_urbanwatch_data/ny_aggregated_land_cover.tif", overwrite = TRUE)

writeRaster(aggregated_land_cover_seattle, "data/processing_data/10_meters_aggregated_urbanwatch_data/seattle_aggregated_land_cover.tif", overwrite = TRUE)

# To load it back into R
aggregated_land_cover_la <- rast("E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data/LA_aggregated_land_cover.tif")

aggregated_land_cover_ny <- rast("E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data/ny_aggregated_land_cover.tif")
  
aggregated_land_cover_seattle <- rast("E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data/seattle_aggregated_land_cover.tif")

#Define a function to filter parks based on the mode of land cover values (classified vs. unclassified)
filter_parks_by_mode <- function(parks, land_cover_raster) {
  classified_parks <- list()  # List to store parks on classified land (categories 1-9)
  unclassified_parks <- list()  # List to store parks on unclassified land (category 0)
  
  # Loop through each park individually
  for (i in 1:nrow(parks)) {
    park <- parks[i, ]  # Extract one park at a time
    
    # Extract land cover values for this park using exact_extract
    land_cover_values <- exact_extract(land_cover_raster, park, include_cell = FALSE)[[1]]$value
    
    # Calculate the mode of land cover values for this park
    mode_val <- calculate_mode(land_cover_values)
    
    # If the mode is greater than 0, add to classified_parks (on classified land), else add to unclassified_parks
    if (!is.na(mode_val) && mode_val > 0) {
      classified_parks[[length(classified_parks) + 1]] <- park  # Add to classified parks
    } else {
      unclassified_parks[[length(unclassified_parks) + 1]] <- park  # Add to unclassified parks
    }
  }
  
  # Combine classified and unclassified parks into separate data frames
  classified_parks <- do.call(rbind, classified_parks)
  unclassified_parks <- do.call(rbind, unclassified_parks)
  
  return(list(classified = classified_parks, unclassified = unclassified_parks))
}

#Apply the function to classify parks (classified vs. unclassified) based on land cover value
la_park_classification <- filter_parks_by_mode(parks_la_reproj, aggregated_land_cover_la)

ny_park_classification <- filter_parks_by_mode(parks_ny_reproj, aggregated_land_cover_ny)

seattle_park_classification <- filter_parks_by_mode(parks_wa_reproj, aggregated_land_cover_seattle)

# Save the entire park_classification list as an RData file
#saveRDS(la_park_classification, "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/processing_data/classified_parks_urbanwatch/park_classification_LA.RData")

saveRDS(ny_park_classification, "data/processing_data/classified_parks_urbanwatch/park_classification_NY.RData")

saveRDS(seattle_park_classification, "data/processing_data/classified_parks_urbanwatch/park_classification_seattle.RData")

la_park_classification <-readRDS("E:/phd_study/urban_park_community_science_project/data/processing_data/classified_parks_urbanwatch/park_classification_LA.RData")


# Plot the land cover raster
plot(aggregated_land_cover_ny, main = "Merged Parks in NYC", col = c("black", "grey", "lightgreen", "darkgreen", "red", "magenta", "blue", "brown", "yellow", "white"))

# Plot the park classification
plot(st_geometry(la_park_classification$classified), col = "blue", pch = 20, cex = 1.5)
plot(st_geometry(la_park_classification$unclassified), add = TRUE, col = "red", pch = 20, cex = 1.5)


plot(st_geometry(ny_park_classification$classified), col = "blue", pch = 20, cex = 1.5)
plot(st_geometry(ny_park_classification$unclassified), add = TRUE, col = "red", pch = 20, cex = 1.5)

plot(st_geometry(seattle_park_classification$classified), col = "blue", pch = 20, cex = 1.5)
plot(st_geometry(seattle_park_classification$unclassified), add = TRUE, col = "red", pch = 20, cex = 1.5)


#Create 250-meter buffer around all the classified parks
la_parks_buffered <- st_buffer(la_park_classification$classified, dist = 250)
ny_parks_buffered <- st_buffer(ny_park_classification$classified, dist = 250)
seattle_parks_buffered <- st_buffer(seattle_park_classification$classified, dist = 250)

#Create 100-meter buffer around all the classified parks
la_parks_buffered_100 <- st_buffer(la_park_classification$classified, dist = 100)
ny_parks_buffered_100 <- st_buffer(ny_park_classification$classified, dist = 100)
seattle_parks_buffered_100 <- st_buffer(seattle_park_classification$classified, dist = 100)

#Create 50-meter buffer around all the classified parks
la_parks_buffered_50 <- st_buffer(la_park_classification$classified, dist = 50)
ny_parks_buffered_50 <- st_buffer(ny_park_classification$classified, dist = 50)
seattle_parks_buffered_50 <- st_buffer(seattle_park_classification$classified, dist = 50)


##Merge the buffer between classified park
# Identify overlapping buffers using st_intersects
#100m
la_buffer_intersections_100 <- st_intersects(la_parks_buffered_100, sparse = FALSE)
ny_buffer_intersections_100 <- st_intersects(ny_parks_buffered_100, sparse = FALSE)
seattle_buffer_intersections_100 <- st_intersects(seattle_parks_buffered_100, sparse = FALSE)

# Convert intersection matrix to an igraph object
la_graphs_100 <- graph_from_adjacency_matrix(la_buffer_intersections_100, mode = "undirected")
ny_graphs_100 <- graph_from_adjacency_matrix(ny_buffer_intersections_100, mode = "undirected")
seattle_graphs_100 <- graph_from_adjacency_matrix(seattle_buffer_intersections_100, mode = "undirected")

# Find connected components (groups of overlapping buffers)
la_components_100 <- igraph::components(la_graphs_100)
ny_components_100 <- igraph::components(ny_graphs_100)
seattle_components_100 <- igraph::components(seattle_graphs_100)

#50m
la_buffer_intersections_50 <- st_intersects(la_parks_buffered_50, sparse = FALSE)
ny_buffer_intersections_50 <- st_intersects(ny_parks_buffered_50, sparse = FALSE)
seattle_buffer_intersections_50 <- st_intersects(seattle_parks_buffered_50, sparse = FALSE)

# Convert intersection matrix to an igraph object
la_graphs_50 <- graph_from_adjacency_matrix(la_buffer_intersections_50, mode = "undirected")
ny_graphs_50 <- graph_from_adjacency_matrix(ny_buffer_intersections_50, mode = "undirected")
seattle_graphs_50 <- graph_from_adjacency_matrix(seattle_buffer_intersections_50, mode = "undirected")

# Find connected components (groups of overlapping buffers)
la_components_50 <- igraph::components(la_graphs_50)
ny_components_50 <- igraph::components(ny_graphs_50)
seattle_components_50 <- igraph::components(seattle_graphs_50)

# Merge buffers by component and track ParkID
buffer_merge_fun<- function(components, parks_buffered, park_classification){

merged_buffers <- list()
merged_park_ids <- list()

for (i in unique(components$membership)) {
  # Get the indices of buffers in the current component
  group_indices <- which(components$membership == i)
  
  # Merge the buffers in this group using st_union()
  merged_buffer <- st_union(parks_buffered[group_indices, ])
  
  # Convert to single polygon if it's a multipolygon
  # # Use convex hull to get a single polygon enclosing all parts
  #if (inherits(merged_buffer, "sfc_MULTIPOLYGON")) {
  #   merged_buffer <- st_convex_hull(merged_buffer)
  # }
  
  # Store the merged buffer
  merged_buffers[[i]] <- merged_buffer
  
  # Collect the ParkIDs of the merged parks
  merged_park_ids[[i]] <- park_classification$classified$ParkID[group_indices]
}

# Ensure all merged geometries are valid sf objects
merged_buffers <- lapply(merged_buffers, function(geom) {
  if (inherits(geom, "sfc")) {
    return(geom)  # Already in correct format
  } else {
    return(st_geometry(geom))  # Convert to geometry if needed
  }
})

# Flatten the list of geometries
merged_buffers <- do.call(c, merged_buffers)

# Create the sfc object with the correct CRS
merged_buffers_sfc <- st_transform(merged_buffers, st_crs(parks_buffered))

# Create the sf object with ParkIDs and merged geometries
merged_buffers_sf <- st_sf(
  geometry = merged_buffers_sfc,
  ParkIDs = sapply(merged_park_ids, function(ids) paste(ids, collapse = ", ")),
  ParkCount = sapply(merged_park_ids, length)
)

return(merged_buffers_sf)

}

la_merged_buffers_100 <- buffer_merge_fun(la_components_100, la_parks_buffered_100, la_park_classification)
ny_merged_buffers_100 <- buffer_merge_fun(ny_components_100, ny_parks_buffered_100, ny_park_classification)
seattle_merged_buffers_100 <- buffer_merge_fun(seattle_components_100, seattle_parks_buffered_100, seattle_park_classification)


la_merged_buffers_50 <- buffer_merge_fun(la_components_50, la_parks_buffered_50, la_park_classification)
ny_merged_buffers_50 <- buffer_merge_fun(ny_components_50, ny_parks_buffered_50, ny_park_classification)
seattle_merged_buffers_50 <- buffer_merge_fun(seattle_components_50, seattle_parks_buffered_50, seattle_park_classification)

#Look at the merged buffers
#100m
png("data/processing_data/la_parks_map_100m_buffered.png", width = 800, height = 800)
plot(st_geometry(la_park_classification$classified), col = "blue", main="100m-buffer, Los Angeles, buffered park (purple, n=425), unclassified park (orange)")
plot(st_geometry(la_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(la_merged_buffers_100$geometry,add = TRUE, col = "purple")
dev.off()

png("data/processing_data/ny_parks_map_100m_buffered.png", width = 800, height = 800)
plot(st_geometry(ny_park_classification$classified), col = "blue", main="100m-buffer, New york City, buffered park (purple, n=779), unclassified park (orange)")
plot(st_geometry(ny_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(ny_merged_buffers_100$geometry,add = TRUE, col = "purple")
dev.off()


png("data/processing_data/seattle_parks_map_100m_buffered.png", width = 800, height = 800)
plot(st_geometry(seattle_park_classification$classified), col = "blue", main="100m-buffer, Seattle, buffered park (purple, n=211), unclassified park (orange)")
plot(st_geometry(seattle_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(seattle_merged_buffers_100$geometry,add = TRUE, col = "purple")
dev.off()


#50m
png("data/processing_data/la_parks_map_50m_buffered.png", width = 800, height = 800)
plot(st_geometry(la_park_classification$classified), col = "blue", main="50m-buffer, Los Angeles, buffered park (purple, n=475), unclassified park (orange)")
plot(st_geometry(la_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(la_merged_buffers_50$geometry,add = TRUE, col = "purple")
dev.off()

png("data/processing_data/ny_parks_map_50m_buffered.png", width = 800, height = 800)
plot(st_geometry(ny_park_classification$classified), col = "blue", main="50m-buffer, New york City, buffered park (purple, n=1353), unclassified park (orange)")
plot(st_geometry(ny_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(ny_merged_buffers_50$geometry,add = TRUE, col = "purple")
dev.off()

png("data/processing_data/seattle_parks_map_50m_buffered.png", width = 800, height = 800)
plot(st_geometry(seattle_park_classification$classified), col = "blue", main="50m-buffer, Seattle, buffered park (purple, n=390), unclassified park (orange)")
plot(st_geometry(seattle_park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(seattle_merged_buffers_50$geometry,add = TRUE, col = "purple")
dev.off()

saveRDS(la_merged_buffers_50, "E:/phd_study/urban_park_community_science_project/data/parks/parks_merged_50_LA.rds")

saveRDS(ny_merged_buffers_50, "E:/phd_study/urban_park_community_science_project/data/parks/parks_merged_50_NY.rds")

saveRDS(seattle_merged_buffers_50, "E:/phd_study/urban_park_community_science_project/data/parks/parks_merged_50_seattle.rds")


####Join the data of the newly merged classified park to the unclassified parks data
###50m
##LA
la_merged_parks<-la_merged_buffers_50%>%
  mutate(type="classified")%>%
rename(ParkID=ParkIDs)

la_unclassified_parks<-la_park_classification$unclassified%>%
  dplyr::select(ParkID, geometry)%>%
  dplyr::mutate(ParkCount=1, type="unclassified")%>%
  dplyr::select(ParkID, ParkCount, type, geometry)

la_classified_unclassified_parks<-rbind(la_merged_parks, la_unclassified_parks)
la_classified_unclassified_parks <- st_make_valid(la_classified_unclassified_parks)
# Ensure CRS alignment
st_crs(la_classified_unclassified_parks) == crs(aggregated_land_cover_la)
la_classified_unclassified_parks <- st_transform(la_classified_unclassified_parks, crs = crs(aggregated_land_cover_la))

plot(st_geometry(la_classified_unclassified_parks$geometry))


##NYC
ny_merged_parks<-ny_merged_buffers_50%>%
  mutate(type="classified")%>%
  rename(ParkID=ParkIDs)

ny_unclassified_parks<-ny_park_classification$unclassified%>%
  dplyr::select(ParkID, geometry)%>%
  dplyr::mutate(ParkCount=1, type="unclassified")%>%
  dplyr::select(ParkID, ParkCount, type, geometry)

ny_classified_unclassified_parks<-rbind(ny_merged_parks, ny_unclassified_parks)
ny_classified_unclassified_parks <- st_make_valid(ny_classified_unclassified_parks)
# Ensure CRS alignment
st_crs(ny_classified_unclassified_parks) == crs(aggregated_land_cover_ny)
ny_classified_unclassified_parks <- st_transform(ny_classified_unclassified_parks, crs = crs(aggregated_land_cover_ny))

plot(st_geometry(ny_classified_unclassified_parks$geometry))

##Seattle
seattle_merged_parks<-seattle_merged_buffers_50%>%
  mutate(type="classified")%>%
  rename(ParkID=ParkIDs)

seattle_unclassified_parks<-seattle_park_classification$unclassified%>%
  dplyr::select(ParkID, geometry)%>%
  dplyr::mutate(ParkCount=1, type="unclassified")%>%
  dplyr::select(ParkID, ParkCount, type, geometry)

seattle_classified_unclassified_parks<-rbind(seattle_merged_parks, seattle_unclassified_parks)
seattle_classified_unclassified_parks <- st_make_valid(seattle_classified_unclassified_parks)
# Ensure CRS alignment
st_crs(seattle_classified_unclassified_parks) == crs(aggregated_land_cover_seattle)
seattle_classified_unclassified_parks <- st_transform(seattle_classified_unclassified_parks, crs = crs(aggregated_land_cover_seattle))

plot(st_geometry(seattle_classified_unclassified_parks$geometry))

#### Calculate areas for every label using exact extraction
###Calculated the total area within a merged park using the raster data
##LA
la_exact_areas <- exact_extract(
  aggregated_land_cover_la, 
  la_classified_unclassified_parks, 
  include_cols = NULL, 
  progress = FALSE
)

la_classified_unclassified_parks$total_area_sqm <- sapply(la_exact_areas, function(df) {
  sum(df$coverage_fraction * 100, na.rm = TRUE)})
la_classified_unclassified_parks$unclassified_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 0, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$road_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 1, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$grass_shrub_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 2, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$tree_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 3, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$building_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 4, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$parking_lot_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 5, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$water_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 6, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$barren_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 7, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$agriculture_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 8, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
la_classified_unclassified_parks$other_area <- sapply(la_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 9, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
  

#Calculating the percentage covers of some land cover types
la_classified_unclassified_parks<-la_classified_unclassified_parks%>%
  dplyr::rowwise() %>% 
  dplyr::mutate(total_classified_area=sum(c_across(road_area:other_area)))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(total_green_space_area=tree_area+grass_shrub_area,
  tree_percent_cover=100*tree_area/total_classified_area, 
       grass_shrub__percent_cover=100*grass_shrub_area/total_classified_area, 
       impervious_surface_percent_cover=100*(road_area+building_area+parking_lot_area+other_area)/total_classified_area)%>%
  dplyr::mutate(new_id=row_number())%>%
  dplyr::select(ParkID, new_id, everything())


saveRDS(la_classified_unclassified_parks, "E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_LA.rds")

##NYC
ny_exact_areas <- exact_extract(
  aggregated_land_cover_ny, 
  ny_classified_unclassified_parks, 
  include_cols = NULL, 
  progress = FALSE
)

ny_classified_unclassified_parks$total_area_sqm <- sapply(ny_exact_areas, function(df) {
  sum(df$coverage_fraction * 100, na.rm = TRUE)})
ny_classified_unclassified_parks$unclassified_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 0, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$road_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 1, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$grass_shrub_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 2, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$tree_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 3, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$building_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 4, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$parking_lot_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 5, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$water_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 6, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$barren_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 7, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$agriculture_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 8, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
ny_classified_unclassified_parks$other_area <- sapply(ny_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 9, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})


#Calculating the percentage covers of some land cover types
ny_classified_unclassified_parks<-ny_classified_unclassified_parks%>%
  dplyr::rowwise() %>% 
  dplyr::mutate(total_classified_area=sum(c_across(road_area:other_area)))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(total_green_space_area=tree_area+grass_shrub_area,
                tree_percent_cover=100*tree_area/total_classified_area, 
                grass_shrub__percent_cover=100*grass_shrub_area/total_classified_area, 
                impervious_surface_percent_cover=100*(road_area+building_area+parking_lot_area+other_area)/total_classified_area)%>%
  dplyr::mutate(new_id=row_number())%>%
  dplyr::select(ParkID, new_id, everything())


saveRDS(ny_classified_unclassified_parks, "E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_NY.rds")

##Seattle
seattle_exact_areas <- exact_extract(
  aggregated_land_cover_seattle, 
  seattle_classified_unclassified_parks, 
  include_cols = NULL, 
  progress = FALSE
)

seattle_classified_unclassified_parks$total_area_sqm <- sapply(seattle_exact_areas, function(df) {
  sum(df$coverage_fraction * 100, na.rm = TRUE)})
seattle_classified_unclassified_parks$unclassified_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 0, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$road_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 1, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$grass_shrub_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 2, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$tree_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 3, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$building_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 4, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$parking_lot_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 5, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$water_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 6, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$barren_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 7, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$agriculture_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 8, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
seattle_classified_unclassified_parks$other_area <- sapply(seattle_exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 9, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})


#Calculating the percentage covers of some land cover types
seattle_classified_unclassified_parks<-seattle_classified_unclassified_parks%>%
  dplyr::rowwise() %>% 
  dplyr::mutate(total_classified_area=sum(c_across(road_area:other_area)))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(total_green_space_area=tree_area+grass_shrub_area,
                tree_percent_cover=100*tree_area/total_classified_area, 
                grass_shrub__percent_cover=100*grass_shrub_area/total_classified_area, 
                impervious_surface_percent_cover=100*(road_area+building_area+parking_lot_area+other_area)/total_classified_area)%>%
  dplyr::mutate(new_id=row_number())%>%
  dplyr::select(ParkID, new_id, everything())


saveRDS(seattle_classified_unclassified_parks, "E:/phd_study/urban_park_community_science_project/data/parks/50m_merged_classified_parks_with_unclassified_parks_sqm_area_Seattle.rds")