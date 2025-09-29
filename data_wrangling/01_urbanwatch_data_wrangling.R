#Load the reuqired oackages ----
library(terra)
library(dplyr)
library(sf)
#library(raster)
library(ggplot2)

###Merge the data from seperated tif files from Urbanwatch
city<-"denton" #The city I want to work on, all in lowercase letters

### List all GeoTIFF files in the directory
file_paths <- list.files(path = paste0("data/urbanwatch_data/01_raw_data/",city), pattern = "\\.tif$", full.names = TRUE)
# Load the GeoTIFF files as RasterLayers using terra
rasters <- lapply(file_paths, rast)
# Merging rasters
merged_raster <- do.call(merge, rasters)

# Visual inspection to confirm the merge has been successful
plot(merged_raster, main=paste0(city, ", Merged Raster"))
#### Save the Merged Raster
writeRaster(merged_raster, filename = paste0("data/urbanwatch_data/02_merged_raster_data/", city,"_merged_raster.tif"), overwrite = TRUE)


###Extract the land cover type data----
##Extract RGB values directly from the raster----

#rgb_values <- values(rasters[[1]])  # This extracts a matrix where each column corresponds to a band

# Convert the matrix to a data frame

#rgb_data <- data.frame(
#  red = rgb_values[, 1],
#  green = rgb_values[, 2],
#  blue = rgb_values[, 3]
#)

# Now find unique combinations that actually appear in your data
#unique_rgb_combinations <- unique(rgb_data)
#saveRDS(unique_rgb_combinations,"data/rgb_data.rds")
#unique_rgb_combinations <- readRDS("data/rgb_data.rds")
# View the resulting unique combinations

#print(unique_rgb_combinations)

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
# Apply the updated classification function to each pixel in the raster

merged_raster<-rast(paste0("data/urbanwatch_data/02_merged_raster_data/",city, "_merged_raster.tif"))

classified_raster <- app(merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))

# Create a data frame for category definitions
levels(classified_raster) <- data.frame(
  code = 0:9,
  label = c("Unclassified",  "Road", "Grass/Shrub","Tree Canopy", "Building","Parking Lot", 
            "Water",  "Barren", "Agriculture", "Others"))

#Compare plots between the raw and the classified ones
plot(merged_raster)
plot(classified_raster, col=c("black", "grey",  "lightgreen", "darkgreen", "red", "magenta",  "blue",  "brown", "yellow", "white"))

### Save the classified raster 
writeRaster(classified_raster, paste0("data/urbanwatch_data/03_classified_land_cover_data/", city, "_classified_land_cover.tif"), overwrite=TRUE)










