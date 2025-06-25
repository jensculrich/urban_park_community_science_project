#Load the reuqired oackages ----
library(terra)
library(dplyr)
library(sf)
#library(raster)
library(ggplot2)

# The UrbanWatch data ----
##Seattle ----
### List all GeoTIFF files in the directory
seattle_file_paths <- list.files(path = "~/Documents/phd_study/iNat_project/data/UrbanWatch/seattle", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
seattle_rasters <- lapply(seattle_file_paths, rast)
# Merging rasters
seattle_merged_raster <- do.call(merge, seattle_rasters)
# Visual inspection to confirm the merge has been successful
plot(seattle_merged_raster, main="Merged Seattle Raster")
#### Save the Merged Raster
writeRaster(seattle_merged_raster, filename = "~/Documents/phd_study/iNat_project/data/UrbanWatch/02_merged_raster_data/seattle_merged_raster.tif", overwrite = TRUE)

##San Francisco----
### List all GeoTIFF files in the directory
sf_file_paths <- list.files(path = "~/Documents/phd_study/iNat_project/data/UrbanWatch/01_raw_data/san_francisco", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
sf_rasters <- lapply(sf_file_paths, rast)
# Merging rasters
sf_merged_raster <- do.call(merge, sf_rasters)
# Visual inspection to confirm the merge has been successful
plot(sf_merged_raster, main="Merged San Francisco Raster")
#### Save the Merged Raster 
writeRaster(sf_merged_raster, filename = "~/Documents/phd_study/iNat_project/data/UrbanWatch/02_merged_raster_data/san_franc_merged_raster.tif", overwrite = TRUE)
sf_merged_raster<-rast("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/san_franc_merged_raster.tif")
##Los Angeles----
# List all GeoTIFF files in the directory
la_file_paths <- list.files(path = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/01_raw_data/los_angeles/los_angeles", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
la_rasters <- lapply(la_file_paths, rast)
# Merging rasters
la_merged_raster <- do.call(merge, la_rasters)

# Visual inspection to confirm the merge has been successful
plot(la_merged_raster, main="Merged San Francisco Raster")
##Riverisde----
# List all GeoTIFF files in the directory
riverside_file_paths <- list.files(path = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/01_raw_data/riverside", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
riverside_rasters <- lapply(riverside_file_paths, rast)
# Merging rasters
riverside_merged_raster <- do.call(merge, riverside_rasters)

# Visual inspection to confirm the merge has been successful
plot(riverside_merged_raster, main="Merged San Francisco Raster")
#### Save the Merged Raster
writeRaster(riverside_merged_raster, filename = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/riverside_merged_raster.tif", overwrite = TRUE)
riverside_merged_raster<-rast("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/riverside_merged_raster.tif")
##San Diego----
# List all GeoTIFF files in the directory
sd_file_paths <- list.files(path = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/01_raw_data/san_diego", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
sd_rasters <- lapply(sd_file_paths, rast)
# Merging rasters
sd_merged_raster <- do.call(merge, sd_rasters)

# Visual inspection to confirm the merge has been successful
plot(sd_merged_raster)

# Save the Merged Raster
writeRaster(sd_merged_raster, filename = "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/san_diego_merged_raster.tif", overwrite = TRUE)

sd_merged_raster<-rast("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/san_diego_merged_raster.tif")


##New York City ----
### List all GeoTIFF files in the directory
ny_file_paths <- list.files(path = "E:/phd_study/urban_park_community_science_project/data/urbanwatch_data/01_raw_data/new_york_city/", pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
ny_rasters <- lapply(ny_file_paths, rast)
# Merging rasters
ny_merged_raster <- do.call(merge, ny_rasters)
# Visual inspection to confirm the merge has been successful
plot(ny_merged_raster, main="Merged New York City Raster")
#### Save the Merged Raster
writeRaster(ny_merged_raster, filename = "E:/phd_study/urban_park_community_science_project/data/urbanwatch_data/02_merged_raster_data/ny_merged_raster.tif", overwrite = TRUE)





#Extract the land cover type data----
##Extract RGB values directly from the raster----
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
  if (is.na(r) || is.na(g) || is.na(b)) return(0)  # "Unclassified"
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


## Assign labels of land cover to the raster data
##Seattle----
# Apply the updated classification function to each pixel in the raster
seattle_classified_raster <- app(seattle_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))

# Create a data frame for category definitions
levels(seattle_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))

### Save the classified raster 
writeRaster(seattle_classified_raster, "~/Documents/phd_study/iNat_project/data/UrbanWatch/03_classified_land_cover_data/seattle_classified_land_cover.tif", overwrite=TRUE)

#Compare plots between the raw and the classified ones
plot(seattle_merged_raster)
plot(seattle_classified_raster)

## California ----
### Assign labels of land cover to the raster data

# Apply the updated classification function to each pixel in the raster
sf_classified_raster <- app(sf_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))
la_classified_raster <- app(la_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))
riverside_classified_raster <- app(riverside_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))
# Create a data frame for category definitions
levels(sf_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))
levels(la_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))
levels(riverside_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))
#Compare plots between the raw and the classified ones
plot(sf_merged_raster)
plot(la_merged_raster)

plot(sf_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))

plot(la_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))
plot(riverside_classified_raster, main="Seattle Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))
### Save the classified raster 
writeRaster(sf_classified_raster, "~/Documents/phd_study/iNat_project/data/UrbanWatch/03_classified_land_cover_data/san_fran_classified_land_cover.tif", overwrite=TRUE)
writeRaster(la_classified_raster, "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/03_classified_land_cover_data/la_classified_land_cover.tif", overwrite=TRUE)
writeRaster(riverside_classified_raster, "C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/03_classified_land_cover_data/riverside_classified_land_cover.tif", overwrite=TRUE)


sf_classified_raster<-rast("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/03_classified_land_cover_data/san_fran_classified_land_cover.tif")
la_classified_raster<-rast("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/03_classified_land_cover_data/la_classified_land_cover.tif")




##New York City----
# Apply the updated classification function to each pixel in the raster
ny_classified_raster <- app(ny_merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))

# Create a data frame for category definitions
levels(ny_classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))

### Save the classified raster 
writeRaster(ny_classified_raster, "E:/phd_study/urban_park_community_science_project/data/urbanwatch_data/03_classified_land_cover_data/ny_classified_land_cover.tif", overwrite=TRUE)

#Compare plots between the raw and the classified ones
plot(ny_merged_raster)
plot(ny_classified_raster, main="NYC Classified Land Cover Map", col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))





