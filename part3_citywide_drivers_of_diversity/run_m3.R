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
  "Houston",
  "LA",
  "Minneapolis",
  "NYC",     
  "Philadelphia",
  "Raleigh",
  "SD",
  "SF"
))


simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_richness <- simmed_diversity[[1]]
mean_prop_disturbance_avoidant <- simmed_diversity[[2]]
beta_diversity <- simmed_diversity[[3]]
gamma_diversity <- simmed_diversity[[4]]

#-------------------------------------------------------------------------------
# get the city covariate data
park_size_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")
connectivity_data <- read.csv("./data/city_wide_data/landscape_metrics.csv") %>%
  rename("city" = "city_names")
landcover_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_land_cover_area_diversity.csv")

city_data <- park_size_data %>%
  left_join(., connectivity_data)
%>%
  left_join(., landcover_data) 
%>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:n_cities))



city_data <- city_data %>%
  mutate(log_avg_park_size = log(average_park_size_sqm),
         log_park_size_scaled = center_scale(log_avg_park_size))


