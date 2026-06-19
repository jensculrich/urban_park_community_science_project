#remotes::install_github("r-spatialecology/vectormetrics")
library(vectormetrics)
library(sf)
library(tidyverse)

# landscape and patch metrics are described here:
# https://r-spatialecology.github.io/vectormetrics/reference/index.html

##------------------------------------------------------------------------------
# load shapefile vector

# list of city names
city_names <- c(
  "Boston", # 1
  "Dallas", # 2
  "Houston", # 3
  "LA", # 4
  "NYC", # 5
  "Riverside", # 6
  "SF" # 7
)

# now choose a city (enter the number of the city)
city_name <- city_names[5]

sf <- sf::read_sf(paste0(
  "./data/city_shapefiles/",
  city_name, "/", city_name, 
  "_50_buffered_park_2km_regional_pool.shp")) 

# filter to classified parks 
sf <- sf %>%
  filter(type == "classified") %>%
  mutate(class = "1") %>%
  mutate(id = row_number())

ggplot(sf) +
  geom_sf(aes(fill=type), color = "white") +
  theme_void()

##------------------------------------------------------------------------------
# patch level vector metrics

## Perimeter to Area Ratio indexes

# perim_idx is the one we will use for the analysis
perarea_idx <- vm_p_perim_idx(sf, class_col = "class") %>%
  mutate(id = as.integer(id)) %>%
  rename("perarea_idx" = "value")

sf_new1 <- left_join(sf, perarea_idx, by = join_by(id == id)) 

q <- ggplot(sf_new1) +
  geom_sf(aes(fill=perarea_idx), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Perimeter Index")


## perimter:area
perarea <- vm_p_perarea(sf, class_col = "class") %>%
  mutate(id = as.integer(id)) %>%
  rename("perarea" = "value")

sf_new1 <- left_join(sf, perarea, by = join_by(id == id)) 

r <- ggplot(sf_new1) +
  geom_sf(aes(fill=perarea), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Perimeter:Area Ratio")

## Proximity Index
proximity <- vm_p_proxim(sf, class_col = "class") %>%
  mutate(id = as.integer(id)) %>%
  rename("proximity" = "value")

sf_new2 <- left_join(sf, proximity, by = join_by(id == id))

s <- ggplot(sf_new2) +
  geom_sf(aes(fill=proximity), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Proximity Index")

##------------------------------------------------------------------------------
# save outputs

perarea_idx <- perarea_idx %>%
  dplyr::select(id, perarea_idx)

df <- as.data.frame(cbind(perarea_idx, perarea$perarea, proximity$proximity)) %>%
  rename("perarea" = "perarea$perarea",
         "proximity" = "proximity$proximity")

write.csv(df, 
  paste0(
  "./data/detections_by_city/",
  city_name, "/",  
  "05_50m_", city_name,
  "_patch_shape.csv"))


