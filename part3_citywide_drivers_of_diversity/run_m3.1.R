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

my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

#-------------------------------------------------------------------------------
# get the city covariate data

park_size_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv") 
#connectivity_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_connectivity_metrics_classified_parks_only.csv") %>%
connectivity_data <- read.csv("./data/city_wide_data/landscape_connectivity_metrics.csv") %>%
  rename("city" = "city_names")
landcover_data <- read.csv("./data/city_wide_data/02_urbanwatch_city_wide_land_cover_area_diversity.csv")  %>%
  rowwise() %>%
  mutate(semi_natural = sum(deciduous_forest_sqm, evergreen_forest_sqm, mixed_forest_sqm,
                            grasslandherbaceous_sqm, woody_wetlands_sqm, emergent_herbaceous_wetlands_sqm,
                            shrubscrub_sqm)) %>%
  mutate(percent_semi_natural = semi_natural / total_area_sqm) %>%
  ungroup() %>%
  select(city, percent_semi_natural)
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
  left_join(., landcover_data, by = "city") %>%
  left_join(., regional_landcover_data, by = "city") %>%
  filter(city %in% city_names)

city_data <- city_data[order(city_data$city), ]

city_data <- city_data %>%
  cbind(., city_factor = seq(1:n_cities)) %>% 
  mutate(log_IIC = log(IIC), 
         log_park_size_scaled = center_scale(median_log_park_size),
         log_IIC_scaled = center_scale(log_IIC),
         percent_semi_natural_scaled = center_scale(percent_semi_natural),
         percent_semi_natural_20km_scaled = center_scale(percent_semi_natural_20km)
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
                    log_park_size_scaled, log_IIC_scaled,
                    percent_semi_natural_scaled, percent_semi_natural_20km_scaled)

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_richness <- simmed_diversity[[1]]

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

a1 <- ggplot(city_data, aes(log_park_size_scaled , mean)) +
  geom_smooth(method = lm) +
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4) +
  ylab("Mean Species Richness") +
  xlab("Median log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

a2 <- ggplot(city_data, aes(log_IIC_scaled , mean)) +
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

a3 <- ggplot(city_data, aes(percent_semi_natural_scaled , mean)) +
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
library(bayesplot)
library(projpred)

#---------------------------------------------------------
m3.1 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + log_IIC_scaled + percent_semi_natural_scaled, 
                           data = city_data_richness)

summary(m3.1)
plot(m3.1)
pp_check(m3.1) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1, pars = "log_park_size_scaled", prob = 0.9)

mcmc_areas(as.matrix(m3.1),prob_outer = .95)
mcmc_pairs(as.matrix(m3.1),pars = c("log_park_size_scaled","log_IIC_scaled","percent_semi_natural_scaled"))

#-------------------------------------------------------------------------------
# Instead we could try to determine which variables add 
# most to the predictive power of the model

m3.1.0 <- rstanarm::stan_glm(mean ~ 1, 
                           data = city_data_richness)

(loo <- loo(m3.1))
(loo0 <- loo(m3.1.0))
loo_compare(loo0, loo)

#fitg_cv <- cv_varsel(m3.1, method='forward', cv_method='LOO')
#plot(fitg_cv, stats = c('elpd', 'rmse'))
vs <- varsel(m3.1, method = "L1", nterms_max = 3, nclusters_pred = 10,
             seed = 5555)
a4 <- plot(vs, stats = c('elpd', 'rmse'))

# And we get a LOO based recommendation for the model size to choose
(nsel <- suggest_size(vs, alpha=0.1))
(vsel <- ranking(vs)[[1]][1:nsel])

# Next we form the projected posterior for the chosen model.
projg <- project(vs, nv = nsel, ns = 4000)
projdraws <- as.matrix(projg)
round(colMeans(projdraws),1)
round(posterior_interval(projdraws),1)
a5 <- mcmc_areas(projdraws, pars=c("(Intercept)",vsel))

cowplot::plot_grid(a4, a5, ncol = 1, rel_heights = c(2,1))

#---------------------------------------------------------
m3.1.1 <- rstanarm::stan_glm(mean ~ log_park_size_scaled, 
                             data = city_data_richness)

plot(m3.1.1)
pp_check(m3.1.1) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1.1, pars = "log_park_size_scaled", prob = 0.9)
plot(m3.1.1, "areas", pars = "log_park_size_scaled", prob = 0.9)

#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# I used the mean response estimates to conduct the variable selection
# But we actually have uncertainty in those diversity estimates
# Let's propagate the uncretainty into the submodel chosen by varsel()

# get the full matrix of estimates richness values (not just the quantile summaries)
mean_richness <- simmed_diversity[[1]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 20
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(mean_richness)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_park_size_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

for(i in 1:n_models){
  city_data$response <- mean_richness[,i]
  
  fit <- rstanarm::stan_glm(response ~ log_park_size_scaled, 
                            data = city_data)
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_park_size_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
}

posterior_draws <- as.data.frame(cbind(intercept, log_park_size_scaled, sigma))
# plot the sample densities
mcmc_areas(posterior_draws, 
           pars = c("intercept", "log_park_size_scaled")) +
  labs(title = 
         "Posterior densities of samples from\nmodels fit 100 simulated communities") +
  theme_classic()


# plot on a predictive scale
base <- ggplot(city_data, aes(x = log_park_size_scaled, y = mean)) +
  ylab("Mean Species Richness") +
  xlab("Median log(Park Size)") +
  scale_color_manual(values=my_palette) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

n_draws <- 100 # draw n lines from the post-posterior
base + geom_abline(
  aes(intercept = intercept, slope = log_park_size_scaled), 
  data = sample_n(posterior_draws, n_draws), 
  color = "grey", 
  alpha = 0.5
) + 
  geom_abline(intercept = median(posterior_draws[,1]), slope = median(posterior_draws[,2]), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

