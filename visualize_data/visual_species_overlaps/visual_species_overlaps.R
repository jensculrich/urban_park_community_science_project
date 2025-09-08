# select a region
regions <- c(
  "midwest",
  "northeast",
  "southeast",
  "southwest"
)

region <- regions[2]

# OR do all at once

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data_multicity.R")

# list of city names

city_names <- c(
  # list in alphabetical order
  "Boston",
  "Charlotte",
  "Dallas",
  "DC",
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",
  "Philadelphia",
  "Raleigh",
  "SD",
  "SF"
)

min_species_detections <- 2 # binary park/year/species detections
min_species_for_community_sampling_event = 1 
family_sampling = TRUE # Should enter either TRUE or FALSE 
# family_sampling:
# if false infer sampling event for all butterflies if any butterflies detected
# if true only infer sampling event for butterflies in same family as any butterflies detected

my_data <- prep_data(city_names,
                     min_species_detections,
                     min_species_for_community_sampling_event,
                     family_sampling
)

saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_all_cities.rds"))
my_data <- readRDS( paste0("./run_model/prepped_data/prepped_data_all_cities.rds"))

V <- my_data$V_detections # detections (1==detected)
dim(V)

species_info <- my_data$species_info
n_species <- my_data$n_species # number of species
species <- as.integer(as.factor(species_info$species))
site_data <- my_data$site_data
n_sites <- my_data$n_sites # number of sites
city <- as.integer(as.factor(site_data$city))
n_cities <- length(unique(city))

# sum detections of species at sites across years/months
detections_per_site <- apply(V, c(2, 1), sum)
detections_per_site <- as.data.frame(detections_per_site)
rownames(detections_per_site)

# get site metadata
metadata <- as.data.frame(cbind(1:n_sites, city))

region_df <- as.data.frame(cbind(city = seq(1:n_cities), 
                   region = c(
                     "northeast",
                     "southeast",
                     "southeast",
                     "northeast",
                     "southeast",
                     "southwest",
                     "midwest",
                     "northeast",
                     "northeast",
                     "southeast",
                     "southwest",
                     "southwest")
)) %>%
  mutate(city = as.integer(city))

metadata <- left_join(metadata, region_df)

# reshape into long form with presence/absence

df <- tibble::rownames_to_column(detections_per_site, "site") %>%
  mutate(site = as.numeric(site))

df <- pivot_longer(df, 
                   cols=-site,
                   names_to='species',
                   values_to='detections') 

metadata <- rename(metadata, "site" = "V1")
df <- left_join(df, metadata) %>%
  mutate(species = substr(species, 2, 4))

df <- df %>%
  mutate(species = as.integer(species)) %>%
  group_by(city, species) %>%
  mutate(sum_detections = sum(detections)) %>%
  slice(1) %>%
  mutate(binary = ifelse(sum_detections > 0, 1, 0)) %>%
  ungroup()

city_df <- as.data.frame(cbind(city = seq(1:n_cities),
                 city_character = city_names)) %>%
  mutate(city = as.integer(city))

df <- left_join(df, city_df)

ggplot(df, aes(city_character, as.factor(species), fill = factor(binary))) +
  geom_tile(width = .9, height = .9) +
  labs(fill = "Detection") + 
  scale_y_discrete(labels=species_info$species) +
  theme_bw() +
  scale_fill_manual(breaks = c("0", "1"),
                       values = c("white", "black")) +
  facet_grid(~ region, scale="free", space="free_x")


# calculate bray curtis differences using vegan
library(vegan) # For distance and ordination methods
library(ggplot2) # For plotting
library(ggforce) # For enhanced aesthetics

dist_mat <- vegdist(detections_per_site, method = "bray")
pcoa_result <- cmdscale(dist_mat, eig = TRUE, k = 2)

points <- as.data.frame(pcoa_result$points)
colnames(points) <- c("PCoA1", "PCoA2")
points$region <- as.factor(metadata$region)
variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

ggplot(points, aes(x = PCoA1, y = PCoA2, color = as.factor(region))) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.3) +
  labs(
    title = "PCoA Analysis (Bray-Curtis)",
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  theme_minimal()


# The metaMDS function automatically transforms data and checks solution
# robustness
library(picante)
comm <- detections_per_site[-(c(20, 138, 762, 486, 492)),]
metadata2 <- metadata[-(c(20, 138, 762, 486, 492)),]
comm.bc.mds <- metaMDS(comm, dist = "bray")
# Assess goodness of ordination fit (stress plot)
stressplot(comm.bc.mds)

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
# plot site scores as text
ordiplot(comm.bc.mds, display = "sites", type = "text")
mds.fig <- ordiplot(comm.bc.mds, type = "none")
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = metadata2$region == 
         "northeast")
points(mds.fig, "sites", pch = 19, col = "red", select = metadata2$region == 
         "southwest")
points(mds.fig, "sites", pch = 19, col = "blue", select = metadata2$region == 
         "southeast")
points(mds.fig, "sites", pch = 19, col = "cyan", select = metadata2$region == 
         "midwest")
# add confidence ellipses around habitat types
ordiellipse(comm.bc.mds, metadata2$region, conf = 0.95, label = TRUE)
s# overlay the cluster results we calculated earlier
ordicluster(comm.bc.mds, comm.bc.clust, col = "gray")

