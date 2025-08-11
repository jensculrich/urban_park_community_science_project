#remotes::install_github("r-spatialecology/vectormetrics")
library(vectormetrics)
library(sf)
library(tidyverse)

# landscape and patch metrics are described here:
# https://r-spatialecology.github.io/vectormetrics/reference/index.html

##------------------------------------------------------------------------------
# load shapefile vector

sf <- sf::read_sf(
  "./visualize_data/visual_spatial_data/New York/New York/NY_100m_buffered_park_2km_regional_pool.shp") 

ggplot(sf) +
  geom_sf(aes(fill=type), color = "white") +
  theme_void()

# filter to classified parks for simplicity
sf <- sf %>%
  mutate(class = "1") %>%
  mutate(id = row_number())

##------------------------------------------------------------------------------
# patch level vector metrics

## Area
area <- vm_p_area(sf, class_col = "class") %>%
  mutate(id = as.integer(id))

sf_new <- left_join(sf, area, by = join_by(id == id)) %>%
  rename("area" = "value")

p <- ggplot(sf_new) +
  geom_sf(aes(fill=log(area)), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Patch Area")


## Perimeter to Area Ratio index
perarea <- vm_p_perarea(sf, class_col = "class") %>%
  mutate(id = as.integer(id))

sf_new2 <- left_join(sf, perarea, by = join_by(id == id)) %>%
  rename("perarea" = "value")

q <- ggplot(sf_new2) +
  geom_sf(aes(fill=log(perarea)), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Perimeter:Area Ratio")


## Proximity Index
proximity <- vm_p_proxim(sf, class_col = "class") %>%
  mutate(id = as.integer(id))

sf_new <- left_join(sf, proximity, by = join_by(id == id)) %>%
  rename("proximity" = "value")

r <- ggplot(sf_new) +
  geom_sf(aes(fill=proximity), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Proximity Index")


## Elongation index
elong <- vm_p_elong(sf, class_col = "class") %>%
  mutate(id = as.integer(id))

sf_new <- left_join(sf, elong, by = join_by(id == id)) %>%
  rename("elong" = "value")

s <- ggplot(sf_new) +
  geom_sf(aes(fill=elong), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Elongation Index")

## plot grid
cowplot::plot_grid(p, q, r, s, ncol = 2)

##------------------------------------------------------------------------------
# class/landscape level vector metrics
# I've designated all park units as the same class


## Mean Euclidean Distance to Nearest Neighbor
sf_samples <- sample_n(sf, size = 10)
ggplot(sf_samples) +
  geom_sf(aes(fill=type), color = "white") +
  theme_void()

enn_mn <- vm_c_enn_mn(sf_samples, class_col = "class") %>%
  mutate(id = as.integer(id))

# get euclids for all patches then take the average of only the classified ones
p_enn <- vm_p_enn(sf, class_col = "class") %>%
  mutate(id = as.integer(id))

sf_new <- left_join(sf, enn_mn, by = join_by(id == id)) %>%
  rename("enn_mn" = "value")

s <- ggplot(sf_new) +
  geom_sf(aes(fill=enn_mn), color = "white") +
  scale_fill_viridis_c() +
  theme_bw() +
  ggtitle("Mean Euclidean Distance to Nearest Neighbor")


sf_1 <- filter(sf_new, area > 10)
sf_2 <- filter(sf_new, area < 10)

ggplot(sf_1) +
  geom_sf(aes(fill=type), color = "white") +
  theme_void()
ggplot(sf_2) +
  geom_sf(aes(fill=type), color = "white") +
  theme_void()

# landscape shape index
vm_l_lsi(sf_1)
vm_l_lsi(sf_2)
# core area index
vm_l_cai_mn(sf_1, edge_depth=100)
vm_l_cai_mn(sf_2, edge_depth=100)
# area index
vm_l_area_mn(sf_1)
vm_l_area_mn(sf_2)
# number of distinct cores per patch
vm_l_dcore_mn(sf_1, edge_depth=100)
vm_l_dcore_mn(sf_2, edge_depth=100)
# total core area 
vm_l_tca(sf_1, edge_depth=100)
vm_l_tca(sf_2, edge_depth=100)
# core area 
vm_l_perarea_mn(sf_1)
vm_l_perarea_mn(sf_2)

##------------------------------------------------------------------------------
# class/landscape level vector metrics
# I've designated all park units as the same class

library(lconnect)
library(tidyverse)

land_sf <- sf::read_sf("./data/city_shapefiles/NYC/NYC_50_buffered_park_2km_regional_pool.shp") %>%
  mutate(habitat = 1) %>%
  select(habitat, ParkID, ParkCnt, type, new_id) #%>%
 # sample_n(., 120)

sf::write_sf(land_sf, "./data/city_shapefiles/NYC/NYC_50_buffered_park_2km_regional_pool_lconnect.shp")

# Load the landscape data

land <- upload_land(
  "./data/city_shapefiles/NYC/NYC_50_buffered_park_2km_regional_pool_lconnect.shp", 
  habitat = 1, max_dist = 1000)

# Confirm the class
class(land)
# Plot the landscape aggregate by clusters defined by the “max_dist” argument
plot(land, main = "Landscape clusters")

# calculate landscape wide metrics
metrics <- con_metric(land, metric = "IIC")

# Computing patch importance based on IIC
importance <- patch_imp(land, metric="IIC")

# Confirm the class
class(importance)
# Plot the landscape with patch importance for global connectivity
plot(importance, main="Patch Importance - IIC")

# Save these outputs as shapefiles, using the sf package
sf::st_write(land$landscape, "./data/city_shapefiles/NYC/land.shp")
sf::st_write(importance$landscape, "./data/city_shapefiles/NYC/importance.shp")

##------------------------------------------------------------------------------
# class/landscape level vector metrics
# I've designated all park units as the same class

#library(devtools)
#install_github("oehrij/Reconnect",build_vignettes = TRUE)
library(Reconnect)
library(sf)

help(package = 'Reconnect', help_type = 'html')

set.seed(15)
### create a name of the simulation
## initialize parameters
nc1     = 100                                       ## define number of cells of landscape # 500 cells for dimx and dimy is reasonably good!!
hfr1    = c(0.15)                                   ## define fraction of habitat in landscape, in combination with npatches, this defines patch area
npatches= c(3) ## define number of patches (maximum nr should not be more than 10 times the number of cells..) in combination with hfrs, this defines patch area
sdpa1   = 0                                         ## define sd of patch area
cf1     = 1                                         ## clumping factor: 1 = non-overlapping patches, >1: clumped & overlapping patches
remedge1= TRUE                                      ## remove the edge of the landscape?

### prepare plot
par(mfrow=c(1,1))
### make a simple circle simulation
for(npatch1 in npatches){
  res = simpcirc(dimx=nc1,dimy=nc1,hfr=hfr1,npatch=npatch1,sdpa=sdpa1,cf=cf1,form="circle",return="all",remedge=remedge1)
  ### unify new shape
  sps1 = sf::st_union(res$sps,by_feature=FALSE) %>% st_cast("POLYGON") %>% st_sf # st_cast(sps1,"POLYGON")
  ## total new area
  #tna   = sum(res$rast[res$rast==1])
  tna    = sum(as.numeric(sf::st_area(sps1))) # new area
  ## total nr of new patches
  tnpatch = length(sps1$geometry)
  ##plot result
  plot(res$rast,main=sprintf("total area (ha) =%0.0f habitat fraction = %0.02f \n nr of patches = %0.0f clumping factor = %0.02f",tna,hfr1,npatch1,cf1),cex.main=0.8)
  plot(sps1$geometry,add=TRUE,border="darkgreen")
}
#> [1] "creating patch areas with mean pa=500.0 and sdpa=0.0"
#> [1] "sample landscape for final number of patches = 3"
#> [1] "minimum distance = 25.2 units"
#> [1] 0
#> [1] "3 samples finished after 1 rounds"
#> [1] "finished strata 0 going to next..."

sf <- sf::read_sf(
  "./visualize_data/visual_spatial_data/New York/New York/NY_100m_buffered_park_2km_regional_pool.shp") 

plot(st_geometry(sf) ,main="habitat shapefile")
