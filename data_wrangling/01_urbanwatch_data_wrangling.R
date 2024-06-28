#Load the reuqired oackages ----
library(terra)
library(dplyr)
library(sf)
#library(raster)
library(ggplot2)

# The UrbanWatch data ----
##Seattle ----
### List all GeoTIFF files in the directory----
seattle_file_paths <- list.files(path = "~/Documents/phd_study/iNat_project/data/UrbanWatch/seattle", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
seattle_rasters <- lapply(seattle_file_paths, rast)
# Merging rasters
seattle_merged_raster <- do.call(merge, seattle_rasters)
# Visual inspection to confirm the merge has been successful
plot(seattle_merged_raster, main="Merged Seattle Raster")
#### Save the Merged Raster ----
writeRaster(seattle_merged_raster, filename = "~/Documents/phd_study/iNat_project/data/UrbanWatch/02_merged_raster_data/seattle_merged_raster.tif", overwrite = TRUE)

##San Francisco
### List all GeoTIFF files in the directory----
sf_file_paths <- list.files(path = "~/Documents/phd_study/iNat_project/data/UrbanWatch/01_raw_data/san_francisco", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
sf_rasters <- lapply(sf_file_paths, rast)
# Merging rasters
sf_merged_raster <- do.call(merge, sf_rasters)
# Visual inspection to confirm the merge has been successful
plot(sf_merged_raster, main="Merged San Francisco Raster")
#### Save the Merged Raster ----
writeRaster(sf_merged_raster, filename = "~/Documents/phd_study/iNat_project/data/UrbanWatch/02_merged_raster_data/san_franc_merged_raster.tif", overwrite = TRUE)

#Joining with the iNat data ----
## Washington ----
leps_data_washington<-read.csv("~/Documents/phd_study/iNat_project/data/leps_data_washington.csv")
### Assign labels of land cover to the raster data ----
#Extract RGB values directly from the raster
rgb_values <- values(raster_1)  # This extracts a matrix where each column corresponds to a band
# Convert the matrix to a data frame
rgb_data <- data.frame(
  red = rgb_values[, 1],
  green = rgb_values[, 2],
  blue = rgb_values[, 3]
)
# Now find unique combinations that actually appear in your data
unique_rgb_combinations <- unique(rgb_data)
# View the resulting unique combinations
print(unique_rgb_combinations)

#A function to assign labels of land cover to the raster data
classify_rgb_to_code <- function(r, g, b) {
  if (r == 0 && g == 0 && b == 0) return(0)  # "Unclassified"
  if (r == 133 && g == 133 && b == 133) return(1)  # "Road"
  if (r == 128 && g == 236 && b == 104) return(2)  # "Grass/Shrub"
  if (r == 34 && g == 139 && b == 34) return(3)  # "Tree Canopy"
  if (r == 255 && g == 0 && b == 0) return(4)  # "Building"
  if (r == 255 && g == 0 && b == 192) return(5)  # "Parking lot" (Pink)
  if (r == 0 && g == 0 && b == 255) return(6)  # "Water"
  if (r == 128 && g == 0 && b == 0) return(7)  # "Barren"
  if (r == 255 && g == 193 && b == 37) return(8)  # "Agriculture"
  if (r == 255 && g == 255 && b == 255) return(9)  # "Others" (White)
  return(0)  # Default to "Unclassified"
}

# Apply the updated classification function to each pixel in the raster
seattle_classified_raster <- app(seattle_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))

# Create a data frame for category definitions
levels(seattle_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))

### Save the classified raster ----
writeRaster(seattle_classified_raster, "~/Documents/phd_study/iNat_project/data/UrbanWatch/03_classified_land_cover_data/seattle_classified_land_cover.tif", overwrite=TRUE)

#Compare plots between the raw and the classified ones
plot(seattle_merged_raster)
plot(seattle_classified_raster)

### Convert the biodiversity data to spatial points ----
washington_coords <- st_as_sf(leps_data_washington, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
#Convert the coordinate system of the biodiversity data from WGS84 to NAD83 / UTM zone 10N (EPSG:26910)
washington_coordinates_utm <- st_transform(washington_coords, 26910)

# Extract raster values at the data points of the biodiversity data
washington_coordinates_extracted_values <- extract(seattle_classified_raster, washington_coordinates_utm)

# Add these values back to your CSV data
washington_coordinates_utm$raster_value <- washington_coordinates_extracted_values

# Convert the geometry to longitude and latitude
washington_complete_geo <- st_transform(washington_coordinates_utm, crs = 4326)
# Extract longitude and latitude
washington_geo_coords <- st_coordinates(washington_complete_geo)

# Convert back to a regular data frame for CSV output or other non-spatial uses
seattle_data_complete <- as.data.frame(washington_coordinates_utm)%>%
  mutate(raster_label = raster_value$label,
         longitude = washington_geo_coords[, "X"],
         latitude = washington_geo_coords[, "Y"]) %>%
  dplyr::select(-raster_value) %>%
  filter(!is.na(raster_label))

### Save the updated CSV ----
write.csv(seattle_data_complete, "~/Documents/phd_study/iNat_project/data/biodiversity_data_with_land_classification/leps_seattle_data_with_land_classification.csv")

### Final Plots ----
png("~/Documents/phd_study/iNat_project/data/biodiversity_data_with_land_classification/Seattle_Leps_with_Land_Cover_Map.png", width = 2000, height = 2000, res = 300)

plot(seattle_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))
points(washington_coordinates_utm, col = 'darkorange', pch = 20, cex = 0.6)

dev.off()
## California ----
leps_data_cali<-read.csv("~/Documents/phd_study/iNat_project/data/inat_data/united_states/filtered_data/leps_data_cali.csv")
### Assign labels of land cover to the raster data ----
#Extract RGB values directly from the raster
rgb_values <- values(raster_1)  # This extracts a matrix where each column corresponds to a band
# Convert the matrix to a data frame
rgb_data <- data.frame(
  red = rgb_values[, 1],
  green = rgb_values[, 2],
  blue = rgb_values[, 3]
)
# Now find unique combinations that actually appear in your data
unique_rgb_combinations <- unique(rgb_data)
# View the resulting unique combinations
print(unique_rgb_combinations)

#A function to assign labels of land cover to the raster data
classify_rgb_to_code <- function(r, g, b) {
  if (r == 0 && g == 0 && b == 0) return(0)  # "Unclassified"
  if (r == 133 && g == 133 && b == 133) return(1)  # "Road"
  if (r == 128 && g == 236 && b == 104) return(2)  # "Grass/Shrub"
  if (r == 34 && g == 139 && b == 34) return(3)  # "Tree Canopy"
  if (r == 255 && g == 0 && b == 0) return(4)  # "Building"
  if (r == 255 && g == 0 && b == 192) return(5)  # "Parking lot" (Pink)
  if (r == 0 && g == 0 && b == 255) return(6)  # "Water"
  if (r == 128 && g == 0 && b == 0) return(7)  # "Barren"
  if (r == 255 && g == 193 && b == 37) return(8)  # "Agriculture"
  if (r == 255 && g == 255 && b == 255) return(9)  # "Others" (White)
  return(0)  # Default to "Unclassified"
}

# Apply the updated classification function to each pixel in the raster
sf_classified_raster <- app(sf_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))

# Create a data frame for category definitions
levels(sf_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))

### Save the classified raster ----
writeRaster(sf_classified_raster, "~/Documents/phd_study/iNat_project/data/UrbanWatch/03_classified_land_cover_data/san_fran_classified_land_cover.tif", overwrite=TRUE)

#Compare plots between the raw and the classified ones
plot(sf_merged_raster)
plot(sf_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))

### Convert the biodiversity data to spatial points ----


cali_coords <- st_as_sf(leps_data_cali, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
#Convert the coordinate system of the biodiversity data from WGS84 to NAD83 / UTM zone 10N (EPSG:26910)
washington_coordinates_utm <- st_transform(washington_coords, 26910)

# Extract raster values at the data points of the biodiversity data
washington_coordinates_extracted_values <- extract(seattle_classified_raster, washington_coordinates_utm)

# Add these values back to your CSV data
washington_coordinates_utm$raster_value <- washington_coordinates_extracted_values

# Convert the geometry to longitude and latitude
washington_complete_geo <- st_transform(washington_coordinates_utm, crs = 4326)
# Extract longitude and latitude
washington_geo_coords <- st_coordinates(washington_complete_geo)

# Convert back to a regular data frame for CSV output or other non-spatial uses
seattle_data_complete <- as.data.frame(washington_coordinates_utm)%>%
  mutate(raster_label = raster_value$label,
         longitude = washington_geo_coords[, "X"],
         latitude = washington_geo_coords[, "Y"]) %>%
  dplyr::select(-raster_value) %>%
  filter(!is.na(raster_label))

### Save the updated CSV ----
write.csv(seattle_data_complete, "~/Documents/phd_study/iNat_project/data/biodiversity_data_with_land_classification/leps_seattle_data_with_land_classification.csv")

### Final Plots ----
png("~/Documents/phd_study/iNat_project/data/biodiversity_data_with_land_classification/Seattle_Leps_with_Land_Cover_Map.png", width = 2000, height = 2000, res = 300)

plot(seattle_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))
points(washington_coordinates_utm, col = 'darkorange', pch = 20, cex = 0.6)

dev.off()





