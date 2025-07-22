install.packages("lconnect")
library(lconnect)
library(tidyverse)

land_sf <- sf::read_sf("./visualize_data/visual_spatial_data/merged_classified_parks_with_unclassified_parks_sqm_area_LA.shp") %>%
  filter(type == "classified") %>%
  mutate(habitat = 1) %>%
  select(habitat, ParkID, ParkCnt, type, new_id)

sf::write_sf(land_sf, "./visualize_data/visual_spatial_data/merged_classified_parks_with_unclassified_parks_sqm_area_LA2.shp")

land <- upload_land("./visualize_data/visual_spatial_data/merged_classified_parks_with_unclassified_parks_sqm_area_LA2.shp", habitat = 1, max_dist = 500)

# Confirm the class
class(land)

# Plot the landscape aggregate by clusters defined by the “max_dist” argument
plot(land, main = "Landscape clusters")

#In this code, the function system.file retrieves the path to a file or directory which is part of an installed package. Here it is used to retrieve the sample shapefile.
#As an example, we will compute the landscape connectivity metrics for all the available connectivity metrics described in Table 3, and present them as a data frame:

# Compute the connectivity metrics
metrics <- con_metric(land, metric = c("NC", "LNK", "SLC", "MSC", "CCP", "LCP", "CPL", "ECS", "AWF", "IIC"))

metrics <- con_metric(land, metric = "AWF")

# Visualize the metrics
print(as.data.frame(metrics))

# Computing patch importance based on AWF
importance <- patch_imp(land, metric="AWF")
# Confirm the class
class(importance)
# Plot the landscape with patch importance for global connectivity
plot(importance, main="Patch Importance - AWF")