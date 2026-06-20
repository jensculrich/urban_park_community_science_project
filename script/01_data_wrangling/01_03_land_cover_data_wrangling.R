library(terra)
library(sf)
library(dplyr)

setwd("/Volumes/sea_angel/iNat_urbanwatch/data")


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

# A function to process the raster filtes and spatial join with the iNat data
process_geotiff_files <- function(city_name) {
  # Create the file paths based on the city name
  input_path <- paste0("urbanwatch_data/01_raw_data/", city_name)
  
  output_path_merged <- paste0("~/Documents/research/urban_park_community_science_project/data/Temp/urbanwatch_data/02_merged_raster_data/", city_name, "_merged_raster.tif")
  
  output_path_classified <- paste0("~/Documents/research/urban_park_community_science_project/data/Temp/urbanwatch_data/03_classified_land_cover_data/", city_name, "_classified_land_cover.tif")
  
  # List all GeoTIFF files in the directory
  file_paths <- list.files(path = input_path, pattern = "\\.tif$", full.names = TRUE)
  
  # Load the GeoTIFF files as RasterLayers using terra
  rasters <- lapply(file_paths, rast)
  
  # Merging rasters
  merged_raster <- do.call(merge, rasters)
  
  # Visual inspection to confirm the merge has been successful
  plot(merged_raster, main = paste("Merged", city_name, "Raster"))
  
  # Save the Merged Raster
  writeRaster(merged_raster, filename = output_path_merged, overwrite = TRUE)
  
  ## Assign labels of land cover to the raster data
  # Apply the classification function to each pixel in the raster
  classified_raster <- app(merged_raster, fun = function(x) classify_rgb_to_code(x[1], x[2], x[3]))
  
  # Create a data frame for category definitions
  levels(classified_raster) <- data.frame(
   code = 0:9,
    label = c("Unclassified", "Road", "Grass/Shrub", "Tree Canopy", "Building", "Parking Lot", 
              "Water", "Barren", "Agriculture", "Others")
  )
  
  ### Save the classified raster 
  writeRaster(classified_raster, output_path_classified, overwrite = TRUE)
  
  # Compare plots between the raw and the classified ones
  plot(merged_raster)
  plot(classified_raster, main = paste(city_name, "Classified Land Cover Map"), 
       col = c("black", "grey", "lightgreen", "darkgreen", "red", "magenta", "blue", "brown", "yellow", "white"))
  
}

# using the function to process and output the data
process_geotiff_files("sd")

