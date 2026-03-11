library(tidyverse)

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
         percent_semi_natural = semi_natural / total_area_sqm) %>%
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

# join the data
city_data <- park_size_data %>%
  left_join(., connectivity_data, by ="city") %>%
  left_join(., IIC_connectivity_data, by ="city") %>%
  left_join(., landcover_data, by = "city") 
  #left_join(., regional_landcover_data, by = "city") 

city_data <- city_data[order(city_data$city), ]

city_data <- city_data %>%
  cbind(., city_factor = seq(1:n_cities)) %>% 
  mutate(log_IIC = log(IIC), 
         log_park_size_scaled = center_scale(median_log_park_size),
         log_IIC_scaled = center_scale(log_IIC),
         isolation_scaled = center_scale(mean_isolation),
         percent_semi_natural_scaled = center_scale(percent_semi_natural)
         #,
         #percent_semi_natural_20km_scaled = center_scale(percent_semi_natural_20km)
  ) 

# View the citywide data
ggplot(city_data) +
  geom_point(aes(
    x=log_IIC_scaled, y=percent_semi_natural_scaled, 
    size=log_park_size_scaled, colour=city)) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "Median log(Park Size)") +
  ylab("City-wide Semi-Natural Vegetation Cover") +
  xlab("City-wide Park Connectivity (IIC Metric)") +
  theme_classic() +
  theme(axis.title = element_text(size=16))

ggplot(city_data) +
  geom_point(aes(
    x=log_IIC_scaled, y=isolation_scaled, 
    size=log_park_size_scaled, colour=city)) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "Median log(Park Size)") +
  ylab("City_wide Mean Isolation") +
  xlab("City-wide Park Connectivity (IIC Metric)") +
  theme_classic() +
  theme(axis.title = element_text(size=16))

ggplot(city_data) +
  geom_point(aes(
    x=log_park_size_scaled, y=percent_semi_natural_scaled, 
    size=log_IIC_scaled, colour=city)) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "City-wide Park Connectivity (IIC Metric)") +
  ylab("City-wide Semi-Natural Vegetation Cover") +
  xlab("Median log(Park Size)") +
  theme_classic() +
  theme(axis.title = element_text(size=16))

ggplot(city_data) +
  geom_point(aes(
    x=percent_semi_natural_scaled, y=percent_semi_natural_20km_scaled, 
    size=log_park_size_scaled, colour=city)) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "Median log(Park Size)") +
  ylab("Regional Semi-Natural Vegetation Cover") +
  xlab("City-wide Semi-Natural Vegetation Cover") +
  theme_classic() +
  theme(axis.title = element_text(size=16))

city_data <- select(city_data, city, city_factor, 
                    log_park_size_scaled, 
                    #log_IIC_scaled,
                    percent_semi_natural_scaled, percent_semi_natural_20km_scaled)

# corrplot
M = cor(select(city_data, -city, -isolation_scaled, -log_park_size_scaled, -percent_semi_natural_scaled, -log_IIC_scaled, -log_IIC, -city_factor))
corrplot::corrplot(M, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_richness <- simmed_diversity[[5]] # 2 is median not the mean

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, mean_richness_quantiles_df)

#-------------------------------------------------------------------------------
# get a basic sense of the relationships by plotting

a1 <- ggplot(city_data, aes(log(total_area_sqm) , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Median log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(#legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

a2 <- ggplot(city_data, aes(isolation_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Connectivity (IIC)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

a3 <- ggplot(city_data, aes(percent_semi_natural_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Vegetation Cover") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(a1, a2, a3, ncol = 3)
