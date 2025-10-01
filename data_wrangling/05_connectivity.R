#Here we calculate the isolation metrics for each classified parks within its 2km buffer
library(sf)
library(terra)
library(tidyverse)
library(lconnect)

city<-"philadelphia"
buffer_size<-50
# Load the park shape file
parks_data <- readRDS(paste0("data/data_for_calculating_connectivity/", buffer_size, "m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))
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
    if (sum(in_buffer) > 0) {
      
      #Index the other parks
      other_parks <-which(in_buffer)
      
      # clip parks within buffer (excluding current park)
      other_parks_sf <- st_intersection(parks_sf[other_parks,], buffer)
      
      #calculate the clpped area of each park
      other_parks_area <- as.numeric(st_area(other_parks_sf))
      
      #calculate the distance between i and each other park
      other_parks_distance <- dist_matrix_numeric[i, other_parks]
      
      parks_sf$isolation[i] <- 1/sum(1+log(other_parks_area + 1)/(other_parks_distance + 1))
    }
  }
  return(parks_sf)
}

#A visullization to make sure we have th correct data
ggplot(parks_data) +
  geom_sf(aes(fill = type))

# Calculate distance matrix between all parks
#dist_matrix <- st_distance(result_parks)
dist_matrix <- st_distance(st_centroid(parks_data)) #distance is calculated based on centroid
dist_matrix_numeric <- units::drop_units(dist_matrix)

#distance-based metrics
parks_with_isolation <- isolation_fun(parks_data)

parks_with_isolation%>%
  filter(type=="classified")%>%
  ggplot() +
  geom_sf(fill ="lightblue")+
  geom_sf(data=parks_with_isolation%>%
            filter(new_id==340), fill="red" )+ #highest isolation
  geom_sf(data=parks_with_isolation%>%
            filter(new_id==16), fill="blue" ) #lowest isolation

# Create final connectivity dataframe
final_connectivity_df <- parks_with_isolation%>%
  st_drop_geometry() %>%
  filter(type == "classified") %>%
  dplyr::select(
    new_id, total_area_sqm:isolation)



# View results
head(final_connectivity_df)
summary(final_connectivity_df)
View(final_connectivity_df)
sum(is.na(final_connectivity_df$isolation))
write.csv(final_connectivity_df, paste0("data/final_merged_data/04_", buffer_size , "m_", city, "_isolation.csv"), row.names = FALSE)



