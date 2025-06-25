library(terra)
library(sf)
library(dplyr)
##Extract RGB values directly from the raster
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

# A function to process the raster filtes and spatial join with the iNat data
process_geotiff_files <- function(city_name, state_name, map_title_city_name) {
  # Create the file paths based on the city name
  input_path <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/01_raw_data/", city_name)
  output_path_merged <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/02_merged_raster_data/", city_name, "_merged_raster.tif")
  output_path_classified <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/urbanwatch_data/03_classified_land_cover_data/", city_name, "_classified_land_cover.tif")
  output_csv_path <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/biodiversity_data_with_land_cover_classification/leps_", city_name, "_data_with_land_classification.csv")
  output_plot_path <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/biodiversity_data_with_land_cover_classification/", city_name, "_Leps_with_Land_Cover_Map.png")
  biodiversity_data_path <- paste0("C:/Users/yyjen/Documents/urban_park_community_science_project/data_wrangling/data/inat_data/02_filtered_data/leps_data_", state_name, ".csv")
  
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
  
  ## Load the biodiversity data
  leps_data <- read.csv(biodiversity_data_path)
  
  ## Assign labels of land cover to the raster data
  
  # Apply the updated classification function to each pixel in the raster
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
  
  ### Convert the biodiversity data to spatial points 
  coords <- st_as_sf(leps_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, agr = "constant")
  
  # Transform the coordinate system of the biodiversity data to match the classified raster
  coords_utm <- st_transform(coords, crs(classified_raster))
  
  # Extract raster values at the data points of the biodiversity data
  extracted_values <- extract(classified_raster, coords_utm)
  
  # Add these values back to your CSV data
  coords_utm$raster_value <- extracted_values
  
  # Convert the geometry to longitude and latitude
  complete_geo <- st_transform(coords_utm, crs = 4326)
  # Extract longitude and latitude
  geo_coords <- st_coordinates(complete_geo)
  
  # Convert back to a regular data frame for CSV output or other non-spatial uses
  complete_data <- as.data.frame(coords_utm) %>%
    mutate(raster_label = raster_value$label,
           longitude = geo_coords[, "X"],
           latitude = geo_coords[, "Y"]) %>%
    dplyr::select(-raster_value) %>%
    filter(!is.na(raster_label))
  
  ### Save the updated CSV 
  write.csv(complete_data, output_csv_path)
  
  ### Final Plots 
  png(output_plot_path, width = 2500, height = 2000, res = 300)
  
  plot(classified_raster, main = paste(map_title_city_name, "Classified Land Cover Map"), 
       col = c("black", "grey", "lightgreen", "darkgreen", "red", "magenta", "blue", "brown", "yellow", "white"))
  points(coords_utm, col = 'darkorange', pch = 20, cex = 0.6)
  
  dev.off()
}

# using the function to merge our data
process_geotiff_files("san_diego", "cali")
process_geotiff_files("dallas", "texas")
process_geotiff_files("minneapolis", "minnesota")
process_geotiff_files("denver", "colorado", "Denver")
process_geotiff_files("des_moines", "colorado", "Des Moines")
process_geotiff_files("seattle", "washington", "Seattle")

