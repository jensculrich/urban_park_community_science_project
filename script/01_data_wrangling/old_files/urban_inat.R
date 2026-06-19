### Urban Greenspace Insect iNaturalist Project
### Chris Cosma
### 10/16/2024
### Connectivity analysis 

#### Prep ####

## WD
setwd("/Users/chriscosma/Desktop/Urban iNat Project")

## Packages
library(terra)
library(sf)
library(tidyverse)
library(zip)

#### ParkServe data ####

## Load data
parks = vect("/Users/chriscosma/Desktop/Urban iNat Project/Parkserve_Shapefiles_05212024/ParkServe_Parks.shp")

## Change CRS
parks = project(parks, "EPSG:4326")

# ## Examine
# print(parks)
# 
# ## Test plot
# plot(parks, xlim = c(-118.67, -118.15), ylim = c(33.70, 34.34), main = "Los Angeles")

#### UrbanWatch Data ####

## Since the UrbanWatch data downloads as multiple zipped folders with the content spread randomly, we need to download and consolidate everything into one folder for each city

## Run entire chunk below just once to create the consolidated city folders

# Define the path to your folder
project_path <- "/Users/chriscosma/Desktop/Urban iNat Project/UrbanWatch"

# Get all the zipped files in the folder
zipped_folders <- list.files(project_path, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)

# Destination for consolidated city folders
destination_folder <- file.path(project_path, "Consolidated_Cities")

# Ensure the destination folder exists
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

# Loop through all zipped files
# sometimes it fails to unzip the some of the folders, and you have to go in a specify which folders in zipped_folders to try again
for (zip_file in zipped_folders) {

  ##test
  #zip_file = zipped_folders[1]

  # Create a temporary folder for unzipping
  unzip_folder <- tempfile()

  # Unzip the current zip file
  system(paste("unzip", shQuote(zip_file), "-d", shQuote(unzip_folder)))

  # Navigate to the subfolder that contains the city folders
  subfolder <- list.dirs(unzip_folder, full.names = TRUE, recursive = FALSE)

  # If there's only one subfolder, go into that one to find the city folders
  if (length(subfolder) == 1) {
    city_folders <- list.dirs(subfolder, full.names = TRUE, recursive = FALSE)
  } else {
    city_folders <- subfolder  # Adjust if structure is different
  }

  # Iterate through each city folder
  for (city_folder in city_folders) {

    ##test
    #city_folder = city_folders[1]

    city_name <- basename(city_folder)  # Extract the city name
    target_city_folder <- file.path(destination_folder, city_name)  # Define the target city folder

    # Create the target city folder if it doesn't exist
    if (!dir.exists(target_city_folder)) {
      dir.create(target_city_folder, recursive = TRUE)
    }

    # Copy files from the current city folder to the target city folder
    city_files <- list.files(city_folder, full.names = TRUE)
    file.copy(city_files, target_city_folder, overwrite = TRUE)
  }

  # Clean up the temporary folder
  unlink(unzip_folder, recursive = TRUE)
}

#### UrbanWatch Processing and Joining with ParkServe####

## Join all of the UrbanWatch sub-.tifs for each city into one and merge with the ParkServe Data

###NOTE: The below was taking very long and there could be errors. It may be better to convert everything to spatial polygon data frame format first. Jenny's new spatial files may negate the need to do any of this. 

# Set the path to the directory with city subfolders
city_folder_path <- "/Users/chriscosma/Desktop/Urban iNat Project/UrbanWatch/Consolidated_Cities"

# List all city folders
city_folders <- list.dirs(city_folder_path, full.names = TRUE, recursive = FALSE)

# Create an empty list to store merged objects
city_shapes <- list()

for (city_folder in city_folders) {
  
  # test 
  city_folder = city_folders[8]
  
  # Get the city name from the folder path
  city_name <- basename(city_folder)
  
  print(city_name)
  
  # List all .tif files in the city folder
  tif_files <- list.files(city_folder, pattern = "\\.tif$", full.names = TRUE)

  # Open as rasters
  rasters <- lapply(tif_files, rast)
  
  # Merge the raster layers
  merged_raster <- do.call(merge, rasters)
  
  #average the values of the 3 layers
  averaged_raster <- mean(merged_raster, na.rm = TRUE)
  
  #convert to polygon outline
  outline <- as.polygons(averaged_raster > 0)
  
  #change CRS
  outline <- project(outline, "EPSG:4326")
  
  #crop parks data to the city extent
  parks_ext <- crop(parks, ext(outline))
  
  #crop parks data to the city outline
  parks_crop <- intersect(parks_ext, outline)
  
  parks_filtered <- parks_ext[, ]
  
  
  plot(parks_crop)
  
  # Store the city outlines in the list
  city_shapes[[city_name]] <- outline
    
  
}

# Combine all city shapes into one sf object
combined_city_shapes <- do.call(rbind, lapply(city_shapes, st_as_sf))

# Save
st_write(combined_city_shapes, "/Users/chriscosma/Desktop/Urban iNat Project/UrbanWatch/combined_city_shapes.shp")


##### Landscape Metrics ######

#once you have a shapefile or raster containing all the parks of interest with their associated park ID to link back to the rest of the data, use the below packages to calculate the patch-level and landscape-level metrics. Make sure to preserve park ID throughout!

#There may be other good R packages to use, but these seem to be the most common

#packages
#https://cran.r-project.org/web/packages/lconnect/index.html
#https://cran.r-project.org/web/packages/landscapemetrics/index.html 

install.packages(c("lconnect", "landscapemetrics"))

library(lconnect)
library(landscapemetrics)

#Below I'm just showing the functions that calculate the comprehensive suite of metrics. Look into package details for other functions

#Assuming that the raster/SF object is called "greenspaces"

#Landscape metrics 

landscape_metrics <- calculate_lsm(greenspaces, level = c("patch", "class", "landscape"))

#lconnect

#I think for this one you need to specify a max dispersal distance, so not sure if we can use this
vec_path <- "path/to/your/greensapces.shp"
greenspaces <- upload_land("path/to/your/shapefile.shp", habitat = 1, max_dist = 500)

# Calculate a suite of connectivity metrics
connectivity_metrics <- con_metric(greenspaces, metric = c("NC", "LNK", "SLC", "MSC", "CCP", "LCP", "CPL", "ECS", "AWF", "IIC"))

#Will have to compare landscapemetrics and lconnect to see which one is preferable. They likely do many of the same things. 

#To choose which metrics to focus on, this paper may be useful: 41.	https://link.springer.com/article/10.1007/s10531-024-02938-2 

#“Out of 68 connectivity indicators we found through a literature review, we identified a key-set of six indicators that align with the Essential Biodiversity Variables framework and are suitable to guide rapid action for connectivity and conservation targets in the KM-GBF. We provide an R-tool to support our general approach, which enables a comprehensive evaluation of connectivity for regional spatial planning for biodiversity in regions with moderate to high human disturbance.” 

