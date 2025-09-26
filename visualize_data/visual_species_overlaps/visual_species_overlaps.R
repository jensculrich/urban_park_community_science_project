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
  "Atlanta",
  "Boston",
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Denver",
  "Des_Moines",
  "Detroit", 
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "SD",
  "SF",
  "St_Louis"
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
                     "southeast",
                     "northeast",
                     "southeast",
                     "midwest",
                     "southeast",
                     "northeast",
                     "southeast",
                     "midwest",
                     "midwest",
                     "midwest", 
                     "southeast",
                     "southwest",
                     "midwest",
                     "northeast",
                     "northeast",
                     "southwest",
                     "southeast",
                     "southwest",
                     "southwest",
                     "southwest",
                     "midwest")
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
  ungroup() %>%
  select(-site)

city_df <- as.data.frame(cbind(city = seq(1:n_cities),
                 city_character = city_names)) %>%
  mutate(city = as.integer(city))

df <- left_join(df, city_df)

ggplot(df, aes(city_character, as.factor(species), fill = factor(binary))) +
  geom_tile(width = .9, height = .9) +
  labs(fill = "Detection") + 
  scale_y_discrete(labels=species_info$species) +
  theme_bw() +
  xlab("") +
  ylab("") +
  scale_fill_manual(breaks = c("0", "1"),
                       values = c("white", "black")) +
  facet_grid(~ region, scale="free", space="free_x")


# calculate bray curtis differences using vegan
library(vegan) # For distance and ordination methods
library(ggplot2) # For plotting
library(ggforce) # For enhanced aesthetics

detections_by_city <- df %>%
  select(city_character, species, sum_detections) %>%
  pivot_wider(names_from = species, values_from = sum_detections)

detections_by_city <- as.matrix(detections_by_city[,2:ncol(detections_by_city)])

dist_mat <- vegdist(detections_by_city, method = "bray")
pcoa_result <- cmdscale(dist_mat, eig = TRUE, k = 2)

points <- as.data.frame(pcoa_result$points)
colnames(points) <- c("PCoA1", "PCoA2")

region = c(
  "southeast",
  "northeast",
  "southeast",
  "midwest",
  "southeast",
  "northeast",
  "southeast",
  "midwest",
  "midwest",
  "midwest", 
  "southeast",
  "southwest",
  "midwest",
  "northeast",
  "northeast",
  "southwest",
  "southeast",
  "southwest",
  "southwest",
  "southwest",
  "midwest")

points$region <- as.factor(region)
variance <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

# Then 'relabel' the cities of interest
label <- c(
  # list in alphabetical order
  "Atlanta",
  "Boston",
  "Charlotte",
  "Chicago",
  "Dallas",
  "Washington DC",
  "Denton",
  "Denver",
  "Des Moines",
  "Detroit", 
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "SD",
  "SF",
  "St. Louis"
)

library(ggrepel)

ggplot(points, aes(x = PCoA1, y = PCoA2, color = as.factor(region))) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, alpha = 0.9) +
  labs(
    title = "PCoA Analysis (Bray-Curtis)",
    x = paste0("PCoA1 (", variance[1], "%)"),
    y = paste0("PCoA2 (", variance[2], "%)")
  ) +
  labs(color='Region') +
  geom_text_repel(aes(label = label), size = 5,
                  seed = 42, box.padding = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold"),
        legend.title=element_text(size=18),
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
