library(tidyverse)

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

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

#-------------------------------------------------------------------------------
# get the city covariate data

park_size_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")
connectivity_data <- read.csv("./data/city_wide_data/landscape_connectivity_metrics.csv") %>%
  rename("city" = "city_names")
landcover_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_land_cover_area_diversity.csv")

city_data <- park_size_data %>%
  left_join(., connectivity_data) %>%
  left_join(., landcover_data, by = "city") %>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:n_cities)) %>%
  rowwise() %>%
  mutate(semi_natural = sum(deciduous_forest_sqm, evergreen_forest_sqm, mixed_forest_sqm,
                            grasslandherbaceous_sqm, woody_wetlands_sqm, emergent_herbaceous_wetlands_sqm,
                            shrubscrub_sqm)) %>%
  mutate(percent_semi_natural = semi_natural / total_area_sqm) %>%
  ungroup()

city_data <- city_data %>%
  mutate(log_avg_park_size = log(average_park_size_sqm),
         log_IIC = log(IIC), 
         log_park_size_scaled = center_scale(log_avg_park_size),
         log_IIC_scaled = center_scale(log_IIC),
         percent_semi_natural_scaled = center_scale(percent_semi_natural)
         )

# View the citywide data
ggplot(city_data) +
  geom_point(aes(
    x=log_IIC_scaled, y=percent_semi_natural_scaled, 
    size=log_park_size_scaled, colour=city)) +
  scale_colour_manual(values=my_palette, name = "City") +
  scale_size(name = "Mean Park Size") +
  ylab("City-wide Semi-Natural Vegetation Cover") +
  xlab("City-wide Park Connectivity (IIC Metric)") +
  theme_classic() +
  theme(axis.title = element_text(size=16))

city_data <- select(city_data, city, city_factor, 
                    log_park_size_scaled, log_IIC_scaled,
                    percent_semi_natural_scaled)

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_richness <- simmed_diversity[[1]]
mean_prop_disturbance_avoidant <- simmed_diversity[[2]]
beta_diversity <- simmed_diversity[[3]]
gamma_diversity <- simmed_diversity[[4]]

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data_richness <- cbind(city_data, mean_richness_quantiles_df)

#-------------------------------------------------------------------------------
# get a basic sense of the relationships by plotting

a1 <- ggplot(city_data_richness, aes(log_park_size_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Mean log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

a2 <- ggplot(city_data_richness, aes(log_IIC_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Connectivity (IIC)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

a3 <- ggplot(city_data_richness, aes(percent_semi_natural_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Vegetation Cover") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

cowplot::plot_grid(a1, a2, a3, ncol = 3)

#-------------------------------------------------------------------------------
# model relationship between city-wide predictors and mean richness

#-------------------------------------------------------------------------------
# first just model relationship between city-wide predictors and MEAN ESTIMATE of mean richness

## IMPORTANTLY this does not account for the uncertainty in our simulated richness predictions

library(rstanarm)

#---------------------------------------------------------
m3.1 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + log_IIC_scaled + percent_semi_natural_scaled, 
                           data = city_data_richness)

plot(m3.1)
pp_check(m3.1) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1, pars = "log_park_size_scaled", prob = 0.9)
plot(m3.1, "areas", pars = "log_park_size_scaled", prob = 0.9)

#---------------------------------------------------------
m3.1.1 <- rstanarm::stan_glm(mean ~ log_park_size_scaled, 
                             data = city_data_richness)

plot(m3.1.1)
pp_check(m3.1.1) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1.1, pars = "log_park_size_scaled", prob = 0.9)
plot(m3.1.1, "areas", pars = "log_park_size_scaled", prob = 0.9)

#---------------------------------------------------------
m3.1.2 <- rstanarm::stan_glm(mean ~ log_IIC_scaled, 
                             data = city_data_richness)

plot(m3.1.2)
pp_check(m3.1.2) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1.2, pars = "log_IIC_scaled", prob = 0.9)
plot(m3.1.2, "areas", pars = "log_IIC_scaled", prob = 0.9)

#---------------------------------------------------------
m3.1.3 <- rstanarm::stan_glm(mean ~ percent_semi_natural_scaled, 
                             data = city_data_richness)

plot(m3.1.3)
pp_check(m3.1.3) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1.3, pars = "percent_semi_natural_scaled", prob = 0.9)
plot(m3.1.3, "areas", pars = "percent_semi_natural_scaled", prob = 0.9)
