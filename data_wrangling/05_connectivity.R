#Here we calculate the isolation metrics for each classified parks within its 2km buffer
library(sf)
library(terra)
library(tidyverse)
library(lconnect)

city<-"denver"
buffer_size<-0
# Load the park shape file
parks_data <- readRDS(paste0("data/parks/", buffer_size, "m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))

#need to check if parks_data is in meters, denton is in ft
st_area(parks_data[1,])

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

#st_area(parks_data[1,])

#load the land cover data
aggregated_land_cover <- rast(paste0("E:/phd_study/urban_park_community_science_project/data/processing_data/10_meters_aggregated_urbanwatch_data/", 
                                     city, "_aggregated_land_cover.tif"))

#align projection
aggregated_land_cover <- project(aggregated_land_cover, crs(parks_data))

#create water raster
water_raster <- aggregated_land_cover == 6


# Get indices of classified parks only
classified_idx <- which(parks_data$type == "classified")

# Store original areas for comparison
parks_data$original_area <- as.numeric(st_area(parks_data))

# Initialize tracking
successful_subtractions <- 0
failed_subtractions <- 0
no_water_found <- 0

for (i in classified_idx) {
  
  tryCatch({
    # Convert park to SpatVector
    park_vect <- vect(parks_data[i, ])
    
    # Crop and mask water raster to this park
    water_crop <- crop(water_raster, park_vect)
    water_mask <- mask(water_crop, park_vect)
    
    # Check if there's any water (value = 1) in this park
    water_values <- values(water_mask)
    if (sum(water_values == 1, na.rm = TRUE) == 0) {
      # No water found in this park
      no_water_found <- no_water_found + 1
      next
    }
    
    # Convert water pixels to polygons
    water_poly <- as.polygons(water_mask)
    
    # Convert to sf
    water_sf <- st_as_sf(water_poly)
    
    #Keep only water pixels (value = 1)
    water_sf <- water_sf[water_sf$label == 1, ]

    # Ensure same CRS
    water_sf <- st_transform(water_sf, st_crs(parks_data))

    # Subtract water from park geometry
    permeable_geom <- st_difference(parks_data[i, ], st_union(water_sf))
    
    # Replace geometry
    st_geometry(parks_data)[i] <- st_geometry(permeable_geom)
    
    successful_subtractions <- successful_subtractions + 1
    
  }, error = function(e) {
    # If error occurs, keep original geometry
    warning(paste("Could not subtract water from park", parks_data$ParkID[i], ":", e$message))
    failed_subtractions <- failed_subtractions + 1
  })
}


# Function required to calculate the isolation metrics
isolation_fun <- function(parks_sf) {
  
  parks_sf$isolation <- NA_real_
  
  # Get indices of classified parks only
  classified_idx <- which(parks_sf$type == "classified")
  
  for (i in classified_idx) {
    
    # Create buffer around current park
    buffer <- st_buffer(parks_sf[i,], 2000)
    
    # Find parks within buffer (excluding current park)
    in_buffer <- st_intersects(buffer, parks_sf, sparse = FALSE)[1,]
    in_buffer[i] <- FALSE  # Exclude the current park

      #Index the other parks
      other_parks <-which(in_buffer)
      
      # clip parks within buffer (excluding current park)
      other_parks_sf <- st_intersection(parks_sf[other_parks,], buffer)
      
      #calculate the clipped area of each park
      other_parks_area <- as.numeric(st_area(other_parks_sf))
      
      #calculate the distance between i and each other park
      other_parks_distance <- dist_matrix_numeric[i, other_parks]
      
      if (length(other_parks)== 0 ){
        other_parks_area = 0.0
        other_parks_distance = 0.0
      }
      
      parks_sf$isolation[i] <- 1/sum(1+log(other_parks_area + 1)/(other_parks_distance + 1))
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
# Clip unclassified parks to only include unclassified within the buffer
unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_combined)

#combine both datasets back together
result_parks <- rbind(
  classified_parks,
  unclassified_parks_clipped)


#A visullization to make sure we have th correct data
ggplot(result_parks) + geom_sf(aes(fill = type))

# Calculate distance matrix between all parks
dist_matrix <- st_distance(result_parks) #distance is calculated based on edge to edge
#dist_matrix <- st_distance(st_centroid(parks_data)) 
dist_matrix_numeric <- units::drop_units(dist_matrix)

#distance-based metrics
parks_with_isolation <- isolation_fun(result_parks)


ggplot() +
  geom_sf(data = parks_with_isolation%>%
            filter(type=="classified"))+
    geom_sf(data=parks_with_isolation%>%
              filter(type=="classified")%>%
              slice_max(order_by=isolation, n=1)%>%st_buffer(50), color = "darkgreen", linewidth = 1)+
    geom_sf(data=parks_with_isolation%>%
              filter(type=="classified")%>%
              slice_min(order_by=isolation, n=1)%>%st_buffer(50), color = "blue", linewidth = 1)+

    geom_sf(data=parks_with_isolation%>%
              filter(type=="classified")%>%
              slice_max(order_by=total_area_sqm, n=5), aes(fill = rank(total_area_sqm)))+
    scale_fill_gradientn(colours = colorspace::heat_hcl(5))+
    geom_sf_text(data=parks_with_isolation%>%
                   filter(type=="classified")%>%
                   slice_max(order_by=isolation, n=1), aes(label = round(isolation, digits = 4)), cex=4,  vjust = 2, color = "darkgreen")+
    geom_sf_text(data=parks_with_isolation%>%
                   filter(type=="classified")%>%
                   slice_min(order_by=isolation, n=1), aes(label = round(isolation, digits = 4)), cex=4, hjust = -0.1, vjust = 2, color = "blue")+
    labs(title = paste0(city, " (", buffer_size, "m-buffer), Ranked by park sizes (top 5), labeled with isolation values (max outlined in green, min outlined in blue)"), size=1)+
    theme(plot.title = element_text(size =8))
    
  

# Create final connectivity dataframe
final_connectivity_df <- parks_with_isolation%>%
  st_drop_geometry() %>%
  filter(type == "classified") %>%
  dplyr::select(
    new_id, total_area_sqm:isolation)



# View results
View(final_connectivity_df)
sum(is.na(final_connectivity_df$isolation))
sum(final_connectivity_df$isolation==Inf)
write.csv(final_connectivity_df, paste0("data/final_merged_data/04_", buffer_size , "m_", city, "_isolation_non_water_only.csv"), row.names = FALSE)

final_connectivity_df<-read.csv(paste0("data/final_merged_data/04_", buffer_size , "m_", city, "_isolation_non_water_only.csv"))

