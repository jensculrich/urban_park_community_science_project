## Load necessary libraries
library(terra)
library(sf)
library(dplyr)
library(daymetr)
library(exactextractr)
library(ggplot2)
library(units)
library(igraph)
library(purrr)
## Load the park shapefile
parks <- st_read("E:/phd_study/urban_park_community_science_project/data/Parkserve_Shapefiles_05212024/ParkServe_Parks.shp")

#park shapefile of the Essex County in Ontario, Canada (needed for Detroit)
#essex_parks <- st_read("E:/phd_study/urban_park_community_science_project/data/essex_county_Recreational_Areas/Recreational_Areas.shp")
#plot(st_geometry(essex_parks))

##Define the city (dallas, houston, sf, riverside, sd)
city <- "denver"
state <- "Colorado"

#filter the park shapefile by state, check to see if the city is close to the state line. If so, may need to include the state next to it.
#parks_filtered <- parks %>% filter(Park_State == state | Park_State == "Maryland"  | Park_State == "Virginia" )
parks_filtered <- parks %>% filter(Park_State == state)

table(parks_filtered$Park_State)

# Select only ID and geometry columns with differentiation
#parks_filtered <- parks_filtered %>%
  #select(ParkID, geometry, state = Park_State)

#essex_parks <- essex_parks %>%
  #select(ParkID = OBJECTID, geometry) %>%
  #mutate(state = "Essex_Ontario")

#align CRS
#essex_parks <- st_transform(essex_parks, st_crs(parks_filtered))

#parks_filtered <- rbind(parks_filtered, essex_parks)


#Make sure the shapefile matches with the city/state
plot(st_geometry(parks_filtered))

## Load the land cover raster using terra
# Land cover data with classified values 0-9 (0 = Unclassified, 2 = Grass/Shrub, 3 = Tree Canopy; )
land_cover <- rast(paste0("data/urbanwatch_data/03_classified_land_cover_data/", city, "_classified_land_cover.tif"))

# Reproject parks to match the CRS of the land cover raster of the city, except dallas and denton which some multipolygons cannot be converted
parks_reproj <- st_transform(parks_filtered, crs = crs(land_cover))%>%
    st_cast(to = "POLYGON")

#for dallas
# safe_cast_row <- function(sf_row) {
#   original_area <- st_area(sf_row)
#   
#   # Try to cast
#   tryCatch({
#     cast_result <- st_cast(sf_row, "POLYGON", do_split = TRUE)
#     cast_area <- sum(st_area(cast_result))
#     
#     # Check if area is preserved (within 0.1% tolerance)
#     if (abs(cast_area - original_area) / original_area < 0.001) {
#       return(cast_result)  # Return cast version
#     } else {
#       return(sf_row)  # Keep original if area doesn't match
#     }
#   }, error = function(e) {
#     return(sf_row)  # Keep original if error
#   })
# }
# 
#parks_reproj <- st_transform(parks_filtered, crs = crs(land_cover)) %>%
#st_make_valid()

#parks_reproj <- do.call(rbind, lapply(1:nrow(parks_reproj), function(i) {
#  safe_cast_row(parks_reproj[i, ])
#}))

nrow(parks_filtered)
nrow(parks_reproj)

table(st_geometry_type(parks_filtered))
table(st_geometry_type(parks_reproj))

sum(st_area(parks_filtered))
sum(st_area(parks_reproj))

plot(st_geometry(parks_filtered))
plot(st_geometry(parks_reproj))


# Aggregate the land cover raster
calculate_mode <- function(x) {
  unique_vals <- unique(na.omit(x))  # Remove NA values and get unique values
  if (length(unique_vals) == 0) {
    return(NA)  # Return NA if no valid values
  }
  return(unique_vals[which.max(tabulate(match(x, unique_vals)))])  # Find the mode
}

#I use a factor of 10 to merge 10x10 cells and `modal` function to keep the most frequent category in each block
#aggregated_land_cover <- aggregate(land_cover, fact = 10, fun = calculate_mode)


# Save the aggregated raster as a GeoTIFF file
#writeRaster(aggregated_land_cover, paste0("data/processing_data/10_meters_aggregated_urbanwatch_data/", city, "_aggregated_land_cover.tif"), overwrite = TRUE)

aggregated_land_cover <- rast(paste0("E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data/", city, "_aggregated_land_cover.tif"))

plot(aggregated_land_cover, main=paste0(city, ", Land Cover Raster"), col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white")) # Grass/Shrub=2; Tree Canopy=3


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
park_classification <- filter_parks_by_mode(st_make_valid(parks_reproj), aggregated_land_cover)

##For DC, we have to convert one classified park to unclassified because it is in Maryland
#park_classification$unclassified <- rbind(park_classification$unclassified, park_classification$classified%>%filter(ParkID=="2485100-0002"))

#park_classification$classified<-park_classification$classified%>%filter(ParkID!="2485100-0002")
#table(park_classification$classified$Park_State)

#Check to see if that park is merged correctly
#park_classification$unclassified%>%filter(ParkID=="2485100-0002")


##For Denton County, we have to convert one classified park to unclassified because it is in Maryland
#park_classification$unclassified <- rbind(park_classification$unclassified, park_classification$classified%>%filter(ParkID %in% c("69400-0002", "4872530-0037")))

#park_classification$classified<-park_classification$classified%>%filter(!ParkID %in% c("69400-0002", "4872530-0037"))
#table(park_classification$classified$Park_State)

#Check to see if that park is merged correctly
#park_classification$unclassified%>%filter(ParkID=="2485100-0002")

# Plot the park classification
plot(st_geometry(park_classification$classified), col = "blue", pch = 20, cex = 1.5, main=paste0(city,", raw, classified park (blue), unclassified park (red)"))
plot(st_geometry(park_classification$unclassified), add = TRUE, col = "red", pch = 20, cex = 1.5)

saveRDS(park_classification, paste0("data/processing_data/classified_parks_urbanwatch/park_classification_", city, ".RData"))

park_classification <-readRDS(paste0("data/processing_data/classified_parks_urbanwatch/park_classification_", city, ".RData"))

#Create buffer around all the classified parks
buffer_size<-0
parks_buffered <- st_buffer(park_classification$classified, dist = buffer_size)

##Merge the buffer between classified park
# Identify overlapping buffers using st_intersects
buffer_intersections <- st_intersects(parks_buffered, sparse = FALSE)

# Convert intersection matrix to an igraph object
graphs <- graph_from_adjacency_matrix(buffer_intersections, mode = "undirected")

# Find connected components (groups of overlapping buffers)
components <- igraph::components(graphs)

# Merge buffers by component and track ParkID (only classified parks are merged)
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

merged_buffers <- buffer_merge_fun(components, parks_buffered, park_classification)

plot(st_geometry(park_classification$classified), col = "blue", main=paste0(city,", ", buffer_size, "m-buffer, buffered park (purple), unclassified park (orange)"))
plot(st_geometry(park_classification$unclassified), add = TRUE, col = "orange")
#plot(st_geometry(seattle_parks_buffered_100), add = TRUE, col = NA, border = "orange", lwd = 1)
plot(merged_buffers$geometry, add = TRUE, col = "purple")


saveRDS(merged_buffers, paste0("E:/phd_study/urban_park_community_science_project/data/parks/parks_merged_", buffer_size , "_", city, ".rds"))

merged_buffers<-readRDS(paste0("E:/phd_study/urban_park_community_science_project/data/parks/parks_merged_", buffer_size , "_", city, ".rds"))

####Join the data of the newly merged classified park to the unclassified parks data
merged_parks<-merged_buffers%>%
  mutate(type="classified")%>%
  rename(ParkID=ParkIDs)

unclassified_parks<-park_classification$unclassified%>%
  dplyr::select(ParkID, geometry)%>%
  dplyr::mutate(ParkCount=1, type="unclassified")%>%
  dplyr::select(ParkID, ParkCount, type, geometry)

classified_unclassified_parks<-rbind(merged_parks, unclassified_parks)
classified_unclassified_parks <- st_make_valid(classified_unclassified_parks)
# Ensure CRS alignment
st_crs(classified_unclassified_parks) == crs(aggregated_land_cover)
classified_unclassified_parks <- st_transform(classified_unclassified_parks, crs = crs(aggregated_land_cover))

plot(st_geometry(classified_unclassified_parks$geometry))
#### Calculate areas for every label using exact extraction
###Calculated the total area within a merged park using the raster data
exact_areas <- exact_extract(
  aggregated_land_cover, 
  classified_unclassified_parks, 
  include_cols = NULL, 
  progress = FALSE
)

classified_unclassified_parks$total_area_sqm <- sapply(exact_areas, function(df) {
  sum(df$coverage_fraction * 100, na.rm = TRUE)})
classified_unclassified_parks$unclassified_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 0, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$road_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 1, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$grass_shrub_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 2, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$tree_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 3, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$building_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 4, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})



classified_unclassified_parks$parking_lot_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 5, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$water_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 6, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})
classified_unclassified_parks$barren_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 7, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})


classified_unclassified_parks$agriculture_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 8, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})

classified_unclassified_parks$other_area <- sapply(exact_areas, function(df) {
  # Filter for only the cells with value 0 (unclassified)
  unclassified_df <- df[df$value == 9, ]
  # Sum up the covered areas
  sum(unclassified_df$coverage_fraction * 100, na.rm = TRUE)
})


#Calculating the percentage covers of some land cover types
classified_unclassified_parks<-classified_unclassified_parks%>%
  dplyr::rowwise() %>% 
  dplyr::mutate(total_classified_area=sum(c_across(road_area:other_area)))%>%
  dplyr::ungroup()%>%
  dplyr::mutate(total_green_space_area=tree_area+grass_shrub_area,
                tree_percent_cover=100*tree_area/total_classified_area, 
                grass_shrub__percent_cover=100*grass_shrub_area/total_classified_area, 
                impervious_surface_percent_cover=100*(road_area+building_area+parking_lot_area+other_area)/total_classified_area)%>%
  dplyr::mutate(new_id=row_number())%>%
  dplyr::select(ParkID, new_id, everything())

plot(st_geometry(classified_unclassified_parks[classified_unclassified_parks$type=="classified",]), 
     col = "darkgreen", 
     main = "Classified Parks",
     border = "white")

ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data = classified_unclassified_parks, aes(fill = type), inherit.aes = FALSE)

length(unique(classified_unclassified_parks[classified_unclassified_parks$type=="classified",]$ParkID))
nrow(classified_unclassified_parks[classified_unclassified_parks$type=="classified",])

saveRDS(classified_unclassified_parks, paste0("E:/phd_study/urban_park_community_science_project/data/parks/", buffer_size, "m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))

ggplot() +
  #annotation_map_tile(type = "osm")+
  geom_sf(data = classified_unclassified_parks%>%filter(type=="classified")%>%slice_max(order_by = total_area_sqm, n =1), inherit.aes = FALSE)

