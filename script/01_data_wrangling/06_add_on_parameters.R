#this script is for calculating the add-on parameters
library(dplyr)
library(sf)
library(terra)
library(purrr)
library(exactextractr)
library(vectormetrics) 

#Setting up the data repo path
setwd("/Volumes/sea_angel/iNat_urbanwatch/data")
##Intermediate-Area of land cover types around the park (2km-buffer)
#First we need the sf object of the classified park and unclassified parks
cities <-c("detroit", "chicago", "minneapolis", "des_moines", "st_louis", "ny", "boston", "philadelphia", "dc", "dallas", "houston", "raleigh", "charlotte", "denton", "atlanta", "la", "sf", "riverside", "sd", "phoenix")

for (city in cities){
  print(city)
#read in the sf object
parks <- readRDS(paste0("parks/0m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))

#loop to extract both the attributes and geometries
for (i in 1:nrow(parks)) {
  current_row <- parks[i, ]
  geom <- st_geometry(current_row)
}

#Have a look at the geometries to make sure it is the correct city
#plot(st_geometry(parks))

#read in the land cover data from NLCD
land_cover<- rast("Annual_NLCD_LndCov_2020_CU_C1V1/Annual_NLCD_LndCov_2020_CU_C1V1.tif")

#look at the land cover types
levels(land_cover)

#Create a 2km buffer around classified parks
classified_parks <- parks[parks$type == "classified", ]
classified_parks <- st_transform(classified_parks, crs = crs(land_cover))
buffer_2km <- st_buffer(classified_parks, dist = 2000)  # 2000 meters = 2km

#substracting the parks from the buffers
buffer_donut <- buffer_2km%>%
  dplyr::select(ParkID, new_id, ParkCount, type, geometry)
buffer_donut$geometry <- map2(
  buffer_2km$geometry, 
  classified_parks$geometry,
  ~st_difference(.x, .y)
) %>% st_sfc(crs = st_crs(buffer_2km))


#check to see if the substraction is correct
plot(st_geometry(classified_parks%>%filter(new_id == 1)))
plot(st_geometry(buffer_2km%>%filter(new_id == 1)))
plot(st_geometry(buffer_donut%>%filter(new_id == 1)))


exact_areas <- exact_extract(
  land_cover, 
  buffer_donut, 
  include_cols = NULL, 
  progress = TRUE
)

# Get cell area
cell_area <- prod(res(land_cover))

buffer_donut$total_sur_area_sqm <- sapply(exact_areas, function(df) {
  sum(df$coverage_fraction * cell_area, na.rm = TRUE)})

buffer_donut$open_water_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 11, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$perennial_ice_snow_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 12, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$developed_open_space_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 21, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$developed_low_intensity_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 22, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$developed_med_intensity_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 23, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$developed_high_intensity_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 24, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$barren_land_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 31, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$deciduous_forest_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 41, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$evergreen_forest_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 42, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$mixed_forest_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 43, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$shrub_scrub_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 52, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$grassland_herbaceous_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 71, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$pasture_hay_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 81, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$cultivated_crops_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 82, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$woody_wetlands_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 90, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

buffer_donut$emergent_herbaceous_wetlands_sqm <- sapply(exact_areas, function(df) {
  result_df <- df[df$value == 95, ]
  # Sum up the covered areas
  sum(result_df$coverage_fraction * cell_area, na.rm = TRUE)
})

#double check if the sum of all land types equals to the total area
buffer_donut %>%
  st_drop_geometry() %>%
  mutate(check_sum = rowSums(across(open_water_sqm:emergent_herbaceous_wetlands_sqm))
)%>%
  dplyr::select(total_sur_area_sqm, check_sum)


final_donut<-buffer_donut %>%
  st_drop_geometry()%>%
  mutate(landcover_type_diversity = rowSums(across(open_water_sqm:emergent_herbaceous_wetlands_sqm) > 0, na.rm = TRUE))

write.csv(final_donut, paste0("final_merged_data/add_on_parameters/01_", city, "_2km_buffer_area_around_park.csv"), row.names = FALSE)

}

#now calculate the the average park size across the city  (classified parks only) 
result_df <- data.frame()

for (city in cities){
  
  print(city)
  
  #read in the sf object
  parks <- readRDS(paste0("parks/0m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))
  
  #average size
  average_park_size<-parks%>%
    st_drop_geometry()%>%
    filter(type == "classified")%>%
    summarise(mean_park_size_sqm = mean(total_area_sqm),
              median_park_size_sqm = median(total_area_sqm),
              mean_log_park_size = mean(log(total_area_sqm)),
              median_log_park_size = median(log(total_area_sqm)))%>%
    mutate(city = city)%>%
    select(city, everything())
  
  result_df<-rbind(result_df, average_park_size)
} 

write.csv(result_df, paste0("final_merged_data/add_on_parameters/all_cities_average_park_size_classified_parks_only.csv"), row.names = FALSE)

#Region 1, city wide park connectivity
connectivity_list<-list()

for(city in cities){
  
  print(city)
  
  #read in the sf object
  parks <- readRDS(paste0("parks/0m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))
  
  parks<- parks %>%
    filter(type == "classified")
  
  # Skip if no classified parks
  if(nrow(parks) == 0) {
    warning(paste("No classified parks found for", city))
    next
  }
  
  aggregated_land_cover <- rast(paste0("processing_data/10_meters_aggregated_urbanwatch_data/", 
                                       city, "_aggregated_land_cover.tif"))
  
  #align projection
  aggregated_land_cover <- project(aggregated_land_cover, crs(parks))
  
  
  #create water raster
  water_raster <- aggregated_land_cover == 6
  
  # Store original areas for comparison
  parks$original_area <- as.numeric(st_area(parks))
  
  print("Subtracting water from parks...")
  
  for (i in seq_len(nrow(parks))) {
    
    tryCatch({
      park_vect <- vect(parks[i, ])
      
      # Crop and mask - combined operation is faster
      water_mask <- mask(crop(water_raster, park_vect), park_vect)
      
      # Early exit if no water
      if (!any(values(water_mask) == 1, na.rm = TRUE)) next
      
      print(i)
      # Convert to polygons and sf
      water_sf <- st_as_sf(as.polygons(water_mask))
      water_sf <- water_sf[water_sf$label == 1, ]
      #water_sf <- st_transform(water_sf, st_crs(parks))
      
      # Do the difference
      park_minus_water <- st_difference(parks[i, ], st_union(water_sf))
      
      #check if the whole park is water
      if (!st_is_empty(park_minus_water)) {
        st_geometry(parks)[i] <- st_geometry(park_minus_water)
      } else {
        # Mark for removal
        parks_to_keep[i] <- FALSE
        print(paste("Park", parks$new_id[i], "completely covered by water - will be removed"))
      }
      
    }, error = function(e) {
      stop(paste("Could not subtract water from park", parks$new_id[i], ":", e$message))
    })
  }
  print("Cleaning geometries...")
  
  
  #Remove empty geometries
  parks <- parks[!st_is_empty(parks), ]
  
  # Extract only polygon/multipolygon parts from geometry collections
  parks <- st_collection_extract(parks, "POLYGON")
  
  #Cast all to MULTIPOLYGON 
  parks <- st_cast(parks, "MULTIPOLYGON")
  
  #Remove any that are still invalid
  parks <- parks[st_is_valid(parks), ]
  
  
  # Calculate metrics
  print("Calculating connectivity metrics...")
  
  connectivity_list[[city]] <- data.frame(
    city = city,
    average_proximity_index = vm_l_proxim_mn(parks, n = 100)$value,  #Proximity_index
    average_perimeter_area_ratio_index = vm_l_perarea_mn(parks)$value, #The mean value of all patches Perimeter-Area ratio index at landscape level
    average_elongation = vm_l_elong_mn(parks)$value #The mean elongation
  )
  
  # Cleanup
  rm(parks, aggregated_land_cover, water_raster)
  gc()
}

# Combine results
connectivity_df <- do.call(rbind, connectivity_list)


write.csv(connectivity_df, paste0("final_merged_data/add_on_parameters/02_urbanwatch_city_wide_connectivity_metrics_classified_parks_only.csv"), row.names = FALSE)

##regional 1 - Area of land cover types in the whole city and the land cover type diversity
#read in the urbanwatch land cover data
rg1_lc_area_div_list<-list()
 for (city in cities){
   print(city)
   urbanwatch_rast<-rast(paste0("/Volumes/sea_angel/iNat_urbanwatch/data/processing_data/10_meters_aggregated_urbanwatch_data/",city, "_aggregated_land_cover.tif"))
   
   land_cover<- rast("Annual_NLCD_LndCov_2020_CU_C1V1/Annual_NLCD_LndCov_2020_CU_C1V1.tif")
   
   # Crop land_cover to urbanwatch extent
   urbanwatch_ext_aea <- project(ext(urbanwatch_rast), 
                                 from = crs(urbanwatch_rast),
                                 to = crs(land_cover))
   
   land_cover_cropped <- crop(land_cover, urbanwatch_ext_aea)
   

   # Get frequency counts for each category
   label_counts <- freq(land_cover_cropped)
   res_land <- res(land_cover_cropped)
   pixel_area_m2 <- res_land[1] * res_land[2]
   label_counts$area_m2 <- label_counts$count * pixel_area_m2
   total_area_sqm = sum(label_counts$area_m2)
   
   rg1_lc_area_div_list[[city]]<-label_counts%>%
     mutate(city = city,
            value =paste0(str_to_lower(str_replace_all(str_remove_all(value, ",|/"), " ", "_")), "_sqm"))%>%
     select(city, everything(), -layer, -count)%>%
  rename(land_cover_type = value)%>%
     pivot_wider(names_from = land_cover_type, values_from = area_m2)%>%
    mutate(landcover_type_diversity = rowSums(across(open_water_sqm:emergent_herbaceous_wetlands_sqm) > 0, na.rm = TRUE),
           total_area_sqm = total_area_sqm)
}

rg1_lc_area_combined <- bind_rows(rg1_lc_area_div_list) %>%
  mutate(across(everything(), ~replace_na(., 0)))%>%
  select(city, total_area_sqm, everything())

write.csv(rg1_lc_area_combined, "/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/add_on_parameters/02_urbanwatch_city_wide_land_cover_area_diversity.csv")

#Regional 2, Area of land cover types in the whole city and Land cover type diversity of the whole city
rg2_lc_area_div_list<-list()

for (city in cities){
  print(city)
parks <- readRDS(paste0("parks/0m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))

#### create a large regional species pool with a 20-km buffer around the classified park
classified_parks <- parks[parks$type == "classified", ]
buffer_20km <- st_buffer(classified_parks, dist = 20000)  # 20000 meters = 20km
buffer_20km_combined <- st_union(buffer_20km)
buffer_20km_vect <- vect(buffer_20km_combined)

land_cover<- rast("Annual_NLCD_LndCov_2020_CU_C1V1/Annual_NLCD_LndCov_2020_CU_C1V1.tif")

# Crop land_cover to urbanwatch extent
buffer_20km_ext_area <- project(ext(buffer_20km_vect), 
                              from = crs(buffer_20km_vect),
                              to = crs(land_cover))

land_cover_cropped <- crop(land_cover, buffer_20km_ext_area)

# Get frequency counts for each category
label_counts <- freq(land_cover_cropped)
res_land <- res(land_cover_cropped)
pixel_area_m2 <- res_land[1] * res_land[2]
label_counts$area_m2 <- label_counts$count * pixel_area_m2
total_area_sqm = sum(label_counts$area_m2)

rg2_lc_area_div_list[[city]]<-label_counts%>%
  mutate(city = city,
         value =paste0(str_to_lower(str_replace_all(str_remove_all(value, ",|/"), " ", "_")), "_sqm"))%>%
  select(city, everything(), -layer, -count)%>%
  rename(land_cover_type = value)%>%
  pivot_wider(names_from = land_cover_type, values_from = area_m2)%>%
  mutate(landcover_type_diversity = rowSums(across(open_water_sqm:emergent_herbaceous_wetlands_sqm) > 0, na.rm = TRUE),
         total_area_sqm = total_area_sqm)
}

rg2_lc_area_combined <- bind_rows(rg2_lc_area_div_list) %>%
  mutate(across(everything(), ~replace_na(., 0)))%>%
  select(city, total_area_sqm, everything())

write.csv(rg2_lc_area_combined, "/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/add_on_parameters/03_20km_buffer_city_wide_land_cover_area_diversity.csv")



#Regional 2, City wide park connectivity
land_cover<- rast("Annual_NLCD_LndCov_2020_CU_C1V1/Annual_NLCD_LndCov_2020_CU_C1V1.tif")

rg2_connectivity_list<-list()

for(city in cities){
  
  print(city)
  
  #read in the sf object
  parks <- readRDS(paste0("parks/0m_merged_classified_parks_with_unclassified_parks_sqm_area_", city, ".rds"))
  
  classified_parks <- parks[parks$type == "classified", ]
  buffer_20km <- st_buffer(classified_parks, dist = 20000)  # 20000 meters = 20km
  # Combine all buffers into a single multipolygon
  buffer_20km_combined <- st_union(buffer_20km)
  # Select unclassified parks
  unclassified_parks <- parks[parks$type == "unclassified", ]
  
  unclassified_parks_clipped <- st_intersection(unclassified_parks, buffer_20km_combined)
  
  unclassified_parks_clean <- st_difference(unclassified_parks_clipped, st_buffer(st_union(classified_parks), 0.001))

  result_parks <- rbind(
    classified_parks,
    unclassified_parks_clean
  )
  
  #result_parks <- result_parks[!st_is_empty(result_parks), ]
  
  #align projection
  result_parks <- st_transform(result_parks, crs(land_cover))
  
  #create water raster
  parks_bbox <- st_bbox(st_buffer(result_parks, 5000))  # 5km buffer
  land_cover_cropped <- crop(land_cover, parks_bbox)
  water_raster <- land_cover_cropped == 11
  
  # Store original areas for comparison
  result_parks$original_area <- as.numeric(st_area(result_parks))
  
  print("Subtracting water from parks...")
  
  # Track which parks to keep
  parks_to_keep <- rep(TRUE, nrow(result_parks))
  
  for (i in seq_len(nrow(result_parks))) {
    
    tryCatch({
      park_vect <- vect(result_parks[i, ])
      
      # Crop and mask - combined operation is faster
      water_mask <- mask(crop(water_raster, park_vect), park_vect)
      
      # Early exit if no water
      if (!any(values(water_mask) == 1, na.rm = TRUE)) next
      
      print(i)
      # Convert to polygons and sf
      water_sf <- st_as_sf(as.polygons(water_mask))
      water_sf <- water_sf[water_sf$`NLCD Land Cover Class` == 1, ]

      # Do the difference
      park_minus_water <- st_difference(result_parks[i, ], st_union(water_sf))
      
      #check if the whole park is water
      if (nrow(park_minus_water) > 0) {
        st_geometry(result_parks)[i] <- st_geometry(park_minus_water)
      } else {
        # Mark for removal
        parks_to_keep[i] <- FALSE
        print(paste("Park", result_parks$new_id[i], "completely covered by water - will be removed"))
      }
      
    }, error = function(e) {
      stop(paste("Could not subtract water from park", result_parks$new_id[i], ":", e$message))
    })
  }
  
  print("Cleaning geometries...")
  
  # Remove water-covered parks
  result_parks <- result_parks[parks_to_keep, ]
  
  #Remove empty geometries
  result_parks <- result_parks[!st_is_empty(result_parks), ]
  
  # Extract only polygon/multipolygon parts from geometry collections
  result_parks <- st_collection_extract(result_parks, "POLYGON")
  
  #Cast all to MULTIPOLYGON 
  result_parks <- st_cast(result_parks, "MULTIPOLYGON")
  
  #Remove any that are still invalid
  result_parks <- result_parks[st_is_valid(result_parks), ]
  
  
  # Calculate metrics
  print("Calculating connectivity metrics...")
  
  rg2_connectivity_list[[city]] <- data.frame(
    city = city,
    average_proximity_index = vm_l_proxim_mn(result_parks, n = 100)$value,  #Proximity_index
    average_perimeter_area_ratio_index = vm_l_perarea_mn(result_parks)$value, #The mean value of all patches Perimeter-Area ratio index at landscape level
    average_elongation = vm_l_elong_mn(result_parks)$value #The mean elongation
  )
  
  # Cleanup
  rm(parks, result_parks, water_raster)
  gc()
}

# Combine results
rg2_connectivity_df <- do.call(rbind, rg2_connectivity_list)


write.csv(rg2_connectivity_df, paste0("final_merged_data/add_on_parameters/03_20km_buffer_city_wide_connectivity_metrics_classified_parks_only.csv"), row.names = FALSE)

dark_mask <- (park[[1]] == 0) & (park[[2]] == 0) & (park[[3]] == 0)
park_mask <- mask(park, dark_mask, maskvalue = 1)
plotRGB(park_mask, na.color = "white")
