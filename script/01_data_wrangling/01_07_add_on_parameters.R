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
cities <-c("detroit", "chicago", "minneapolis", "des_moines", "st_louis", "nyc", "boston", "philadelphia", "dc", "dallas", "houston", "raleigh", "charlotte", "denton", "atlanta", "la", "sf", "riverside", "sd", "phoenix", "denver", "tampa")

#calculate the the average park size across the city  (classified parks only) 
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

#city wide park connectivity
summary<-data.frame()

for(i in 1:n_cities){
  # now choose a city (enter the number of the city)
  city <- city_names[i]
  
  print(paste("Working on", city, ", reading in shapefile"))
  # lconnect requires us to have a habitat column in the shapefile
  # here we assign any greenspaces as "habitat" (habitat = 1)
  
  sf <- st_read(paste0(
    "/Volumes/sea_angel/iNat_urbanwatch/data/shapefile_0m_buffered_park_2km_regional_pool/",city, "_0_buffered_park_2km_regional_pool.shp")) 
  
  
  sf <- sf %>%
    filter(type == "classified") %>%
    mutate(habitat = 1) %>%
    select(habitat, new_id)
  
  # Create a temporary shapefile path
  temp_shp <- tempfile(fileext = ".shp")
  
  sf::write_sf(sf, temp_shp)
  
  
  
  ##------------------------------------------------------------------------------
  # calculate landscape connectivity metrics shapefile
  print(paste("Working on", city, ", calculating landscape connectivity metrics"))
  
  # Load the landscape data
  land <- upload_land(temp_shp, 
                      habitat = 1, 
                      max_dist = 2000)
  
  # requires us to set a max dist between which parks can be connected
  # I chose 2000 metres to be consistent with our park site isolation metric 
  
  # Confirm the class
  class(land)
  # Plot the landscape aggregate by clusters defined by the “max_dist” argument
  #plot(land, main = "Landscape clusters")
  
  # Compute the connectivity metrics
  # IIC is fast to calculate
  # other metrics of interest might include AWF
  # https://www.r-bloggers.com/2019/03/lconnect-connectivity-metrics/
  result <- con_metric(land, metric = c("IIC"))
  
  summary<-rbind(summary, data.frame(city = city, IIC = result))
  
}

ggplot(summary) +
  geom_point(aes(x = as.factor(city, y = log(IIC)))
             
# Save these outputs as a csv
write.csv(summary, paste0(
 "/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/add_on_parameters/02_urbanwatch_city_wide_connectivity_metrics_classified_parks_only.csv"), row.names=FALSE)
 
#city wide mean isolation
isolation_summary <- data.frame()

for (city in cities){
  isolation <- read.csv(paste0("/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/04_0m_", city, "_isolation_non_water_only.csv"))
  mean_isolation <- mean(isolation$isolation)
  
  isolation_summary<- rbind(isolation_summary, data.frame(city = city, mean_isolation = mean_isolation))
}

write.csv(isolation_summary, paste0(
  "/Volumes/sea_angel/iNat_urbanwatch/data/final_merged_data/add_on_parameters/04_city_wide_isolation_metrics.csv"), row.names=FALSE)


##Area of land cover types in the whole city and the land cover type diversity
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



