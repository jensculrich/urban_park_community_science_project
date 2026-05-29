library(tidyverse)
library(cowplot)
library(corrplot)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

n_cities <- length(city_names <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "DC",
  "Denton",
  "Denver",
  "Des_moines",
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
  "St_louis",
  "Tampa"
))

city_names_labels <- c(
  "Atlanta",
  "Boston", 
  "Charlotte",
  "Chicago",
  "Dallas",
  "Washington D.C.",
  "Denton",
  "Denver",
  "Des Moines",
  "Detroit",
  "Houston",
  "Los Angeles",
  "Minneapolis",
  "New York City",     
  "Philadelphia",
  "Phoenix",
  "Raleigh",
  "Riverside",
  "San Diego",
  "San Francisco",
  "St. Louis",
  "Tampa"
)

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

#-------------------------------------------------------------------------------
# get the city covariate data

park_size_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv") 
connectivity_data <- read.csv("./data/city_wide_data/04_city_wide_isolation_metrics.csv") 
IIC_connectivity_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_connectivity_metrics_classified_parks_only.csv") 
landcover_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_land_cover_area_diversity.csv") %>%
  rowwise() %>%
  mutate(semi_natural = sum(grass_shrub, tree)) %>%
  mutate(percent_grass_shrub = grass_shrub / total_area_sqm,
         percent_tree = tree / total_area_sqm,
         percent_semi_natural = semi_natural / total_area_sqm, 
         percent_agriculture = agriculture / total_area_sqm) %>%
  ungroup() %>%
  select(city, percent_grass_shrub, percent_tree, percent_semi_natural, total_area_sqm)
regional_landcover_data <-  read.csv("./data/city_wide_data/03_20km_buffer_city_wide_land_cover_area_diversity.csv")  %>%
  rowwise() %>%
  mutate(semi_natural = sum(deciduous_forest_sqm, evergreen_forest_sqm, mixed_forest_sqm,
                            grasslandherbaceous_sqm, woody_wetlands_sqm, emergent_herbaceous_wetlands_sqm,
                            shrubscrub_sqm)) %>%
  mutate(percent_semi_natural_20km = semi_natural / total_area_sqm) %>%
  ungroup() %>%
  select(city, percent_semi_natural_20km)
latitude <- read.csv("./data/city_latitude.csv")

# join the data
city_data <- park_size_data %>%
  left_join(., connectivity_data, by ="city") %>%
  left_join(., IIC_connectivity_data, by ="city") %>%
  left_join(., landcover_data, by = "city") %>%
  left_join(., latitude, by = "city")
#left_join(., regional_landcover_data, by = "city") 

city_data <- city_data[order(city_data$city), ]

city_data <- city_data %>%
  cbind(., city_factor = seq(1:n_cities)) %>% 
  mutate(log_IIC = log(IIC), 
         log_park_size_scaled = center_scale(median_log_park_size),
         log_IIC_scaled = center_scale(log_IIC),
         isolation_scaled = center_scale(mean_isolation),
         percent_tree_scaled = center_scale(percent_tree),
         percent_grassshrub_scaled = center_scale(percent_grass_shrub),
         log_total_area = log(total_area_sqm),
         log_total_area_scaled = center_scale(log_total_area),
         latitude_scaled = center_scale(latitude),
         longitude_scaled = center_scale(longitude)
         #,
         #percent_semi_natural_20km_scaled = center_scale(percent_semi_natural_20km)
  ) 



#-------------------------------------------------------------------------------
# get city-wide park area and add this to the corr plots

log_total_park_area <- vector(length=n_cities)

for(i in 1:length(city_names)){
  
  city <- city_names[i]
  
  # first read the data 
  temp <- cbind(city, read.csv(paste0(
    "./data/detections_by_city/", city, "/04_0m_", city,
    "_isolation_non_water_only.csv"
  ))) 
  
  log_total_park_area[i] <- log(sum(temp$total_green_space_area))
  
}

city_data <- cbind(city_data, log_total_park_area)

#-------------------------------------------------------------------------------
# corrplot

city_data <- city_data %>% 
  mutate(prop_park_area = log_total_park_area / log_total_area)

city_data_plot <- city_data %>%
  rename(#"mean log(park size)" = "mean_log_park_size",
         "median log(park size)" = "median_log_park_size",
         "% herbaceous" = "percent_grass_shrub",
         "% tree" = "percent_tree",
         "% park" = "prop_park_area",
         "mean log(isolation)" = "log_isolation",
         "log(city spatial area)" = "log_total_area"
         )
M = cor(select(city_data_plot, -city, -isolation_scaled, -log_park_size_scaled,
               -median_park_size_sqm, -mean_park_size_sqm, -percent_semi_natural,
               -percent_semi_natural_scaled, -percent_grassshrub_scaled, -percent_tree_scaled,
               -log_IIC_scaled, -log_IIC, -city_factor, -mean_isolation, -mean_log_park_size,
               -total_area_sqm, -log_total_park_area,
               -latitude_scaled))
corrplot::corrplot(M, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)


#-------------------------------------------------------------------------------
# View the citywide data as scatterplots

city_data <- cbind(city_data, city_names_labels)

p <- ggplot(city_data) +
  geom_point(aes(
    x=median_log_park_size, y=log_IIC, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "Median log(Park Size)") +
  xlab("Median log(Park Size)") +
  ylab("City-wide Park Connectivity (IIC)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

q <- ggplot(city_data) +
  geom_point(aes(
    x=percent_tree, y=percent_grass_shrub, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  xlab("City-wide Tree Cover") +
  ylab("City-wide Herbaceous Cover") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

r <- ggplot(city_data) +
  geom_point(aes(
    x=latitude, y=log_total_area, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  xlab("Latitude") +
  ylab("log(City Area (m^2))") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

legend_plot <- ggplot(city_data) +
  geom_point(aes(
    x=latitude, y=log_total_area, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  xlab("Latitude") +
  ylab("log(City Area (m^2))") +
  theme_classic() +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

legend <- get_legend(legend_plot)


cowplot::plot_grid(p, q, r, legend, ncol = 2)


#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity_may25.RDS")

median_richness <- simmed_diversity[[2]] # 2 is median not the mean

#  calculate Means and CI's for the diversity metrics for each city
median_richness_quantiles <- apply(median_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
median_richness_quantiles_df <- as.data.frame(t(median_richness_quantiles))
colnames(median_richness_quantiles_df) <- c("lower90", "lower50",
                                            "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, median_richness_quantiles_df)

#-------------------------------------------------------------------------------
# plot median

p <- ggplot(city_data) +
  geom_point(aes(
    x=city_names_labels, y=mean, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  geom_errorbar(aes(x = city_names_labels, ymin=lower90, ymax=upper90,
                    colour = city_names_labels), width=.2,
                position=position_dodge(.9)) +
  xlab("") +
  ylab("Median Park\nSpecies Richness") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

city_data <- city_data[,-(27:31)]

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity_may25.RDS")

mean_disturbance_or_edge_avoidant <- simmed_diversity[[14]] 

#  calculate Means and CI's for the diversity metrics for each city
mean_disturbance_or_edge_avoidant_quantiles <- apply(mean_disturbance_or_edge_avoidant, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_disturbance_or_edge_avoidant_quantiles_df <- as.data.frame(t(mean_disturbance_or_edge_avoidant_quantiles))
colnames(mean_disturbance_or_edge_avoidant_quantiles_df) <- c("lower90", "lower50",
                                                              "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, mean_disturbance_or_edge_avoidant_quantiles_df)

#-------------------------------------------------------------------------------
# plot median

q <- ggplot(city_data) +
  geom_point(aes(
    x=city_names_labels, y=mean, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  geom_errorbar(aes(x = city_names_labels, ymin=lower90, ymax=upper90,
                    colour = city_names_labels), width=.2,
                position=position_dodge(.9)) +
  xlab("") +
  ylab("Avg. Disturbance or Edge\nAvoidant Species Richness") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

city_data <- city_data[,-(27:31)]

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity_may25.RDS")

# elements of simmed_diversity list are as follows:
beta_diversity <- simmed_diversity[[7]] # 7 is beta diversity

#  calculate Means and CI's for the diversity metrics for each city
beta_diversity_quantiles <- apply(beta_diversity, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
beta_diversity_quantiles_df <- as.data.frame(t(beta_diversity_quantiles))
colnames(beta_diversity_quantiles_df) <- c("lower90", "lower50",
                                           "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, beta_diversity_quantiles_df)

#-------------------------------------------------------------------------------
# plot median

r <- ggplot(city_data) +
  geom_point(aes(
    x=city_names_labels, y=mean, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  geom_errorbar(aes(x = city_names_labels, ymin=lower90, ymax=upper90,
                    colour = city_names_labels), width=.2,
                position=position_dodge(.9)) +
  xlab("") +
  ylab("Jaccard Index of\nSpecies Dissimilarity") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

city_data <- city_data[,-(27:31)]

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity_may25.RDS")

# elements of simmed_diversity list are as follows:
total_richness <- simmed_diversity[[10]] 

#  calculate mean and CI's for the diversity metrics for each city
total_richness_quantiles <- apply(total_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
total_richness_quantiles_df <- as.data.frame(t(total_richness_quantiles))
colnames(total_richness_quantiles_df) <- c("lower90", "lower50",
                                           "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, total_richness_quantiles_df) 

#-------------------------------------------------------------------------------
# plot median

s <- ggplot(city_data) +
  geom_point(aes(
    x=city_names_labels, y=mean, colour=city_names_labels), size=6) +
  scale_colour_manual(values=my_palette, name = "City") +
  geom_errorbar(aes(x = city_names_labels, ymin=lower90, ymax=upper90,
                    colour = city_names_labels), width=.2,
                position=position_dodge(.9)) +
  xlab("") +
  ylab("Total City-Wide\nSpecies Richness") +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size=18),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16))

city_data <- city_data[,-(27:31)]

cowplot::plot_grid(p, q, r, s, ncol = 1, labels = c("a)", "b)", "c)", "d)"),
                   label_size = 18)
