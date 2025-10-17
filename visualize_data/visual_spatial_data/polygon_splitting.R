library(sf)
library(tidyverse)

land_sf <- sf::read_sf("./data/city_shapefiles/NYC/NYC_50_buffered_park_2km_regional_pool.shp") %>%
  filter(type == "classified") %>%
  mutate(habitat = 1) %>%
  select(habitat, ParkID, ParkCnt, type, new_id) 

land_sf$area <- st_area(land_sf) #Take care of units

land_sf <- filter(land_sf, new_id == "16")

# Map the spatial data
ggplot() +
  geom_sf(data = land_sf, aes(fill = as.factor(new_id))) 

sdist = 1200
land_sf_splitnarrow <- splitnarrow(land_sf, sdist, 1e-3)

# Map the spatial data
ggplot() +
  geom_sf(data = bparts, 
          aes(fill = as.factor(new_id))
          ) 

splitnarrow <- function(pol, sdist, eps){
  ###
  ### split a polygon at its narrowest point.
  ###
  
  ### sdist is the smallest value for internal buffering that splits the
  ### polygon into a MULTIPOLYGON and needs computing before running this.
  
  ### eps is another tolerance that is needed to get the points at which the
  ### narrowest point is to be cut.
  
  ## split the polygon into two separate polygons
  bparts = st_buffer(pol, sdist)
  features = st_cast(st_as_sfc(bparts), "POLYGON")
  
  ## find where the two separate polygons are closest, this is where
  ## the internal buffering pinched off into two polygons.
  
  pinch = st_nearest_points(features[1],features[2])
  
  ## buffering the pinch point by a slightly larger buffer length should intersect with
  ## the polygon at the narrow point. 
  inter = st_intersection(
    st_cast(pol,"MULTILINESTRING"),
    st_buffer(pinch,-(sdist-(eps))
    )
  )
  join = st_cast(st_as_sfc(inter), "LINESTRING")
  
  ## join is now two small line segments of the polygon across the "waist".
  ## find the line of closest approach of them:
  splitline = st_nearest_points(join[1], join[2])
  
  ## that's our cut line. Now put that with the polygon and make new polygons:
  mm = st_union(splitline, st_cast(pol, "LINESTRING"))
  parts = st_collection_extract(st_polygonize(mm))
  parts
}
