# select a region
regions <- c(
  "midwest",
  "northeast",
  "southeast",
  "southwest"
)

region <- regions[2]

# get data for northeast
my_data2 <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))

V2 <- my_data2$V_detections # detections (1==detected)
dim(V2)


species_info2 <- my_data2$species_info
n_species2 <- my_data2$n_species # number of species
species2 <- as.integer(as.factor(species_info2$species))
site_data2 <- my_data2$site_data
n_sites2 <- my_data2$n_sites # number of sites
city2 <- as.integer(as.factor(site_data2$city))
n_cities2 <- length(unique(city2))

# sum detections of species at sites across years/months
detections_per_site2 <- apply(V2, c(2, 1), sum)
detections_per_site2 <- as.data.frame(detections_per_site2)
rownames(detections_per_site2)

# get site metadata
metadata2 <- as.data.frame(cbind(1:n_sites2, city2))

# get data for southeast
region <- regions[3]

my_data3 <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))

V3 <- my_data3$V_detections # detections (1==detected)
dim(V3)


species_info3 <- my_data3$species_info
n_species3 <- my_data3$n_species # number of species
species3 <- as.integer(as.factor(species_info3$species))
site_data3 <- my_data3$site_data
n_sites3 <- my_data3$n_sites # number of sites
city3 <- as.integer(as.factor(site_data3$city))
n_cities3 <- length(unique(city3))

# sum detections of species at sites across years/months
detections_per_site3 <- apply(V3, c(2, 1), sum)
detections_per_site3 <- as.data.frame(detections_per_site3)
rownames(detections_per_site3)

# get site metadata
metadata3 <- as.data.frame(cbind(1:n_sites3, city3))

# now join data from all regions
region2 <- rep("northeast", times = n_sites2) 
region3 <- rep("southeast", times = n_sites3) 

metadata2 <- as.data.frame(cbind(metadata2, region2))
metadata3 <- as.data.frame(cbind(metadata3, region3))

detections_per_site <- rbind(detections_per_site2, detections_per_site3)

# calculate bray curtis differences using vegan
library(vegan) # For distance and ordination methods
library(ggplot2) # For plotting
library(ggforce) # For enhanced aesthetics

dist_mat <- vegdist(detections_per_site3, method = "bray")
pcoa_result <- cmdscale(dist_mat, eig = TRUE, k = 2)

points <- as.data.frame(pcoa_result$points)
colnames(points) <- c("PCoA1", "PCoA2")
points$city <- as.factor(metadata3$city)
variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

ggplot(points, aes(x = PCoA1, y = PCoA2, color = as.factor(city))) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3) +
  labs(
    title = "PCoA Analysis (Bray-Curtis)",
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal()


# reshape into long form with presence/absence

df <- tibble::rownames_to_column(detections_per_site2, "site") %>%
  mutate(site = as.numeric(site))

df <- pivot_longer(df, 
                   cols=-site,
                   names_to='species',
                   values_to='detections') 

metadata2 <- rename(metadata2, "site" = "V1")
test <- left_join(df, metadata2) %>%
  mutate(city2 = as.factor(city2))

test <- test %>%
  group_by(city2, species) %>%
  mutate(sum_detections = sum(detections)) %>%
  slice(1) %>%
  mutate(binary = ifelse(sum_detections > 0, 1, 0))

labs <- c("Boston", "DC", "NYC", "Philly")
ggplot(test, aes(city2, species, fill = factor(binary))) +
  geom_tile(width = .9, height = .9) +
  labs(fill = "Detection") + 
  scale_x_discrete(labels= labs)

  facet_wrap(~Period)