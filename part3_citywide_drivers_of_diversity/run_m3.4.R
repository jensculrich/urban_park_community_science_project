library(tidyverse)
library(rstanarm)
library(bayesplot)
library(projpred)

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
IIC_connectivity_data <- read.csv("./data/city_wide_data/landscape_connectivity_IIC.csv") %>%
  rename("city" = "city_names")
IIC_connectivity_data <- IIC_connectivity_data[1:22,] # last row got duplicated in processing
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
         percent_tree_scaled = center_scale(percent_tree),
         percent_semi_natural_scaled = center_scale(percent_semi_natural),
         log_total_area_scaled = center_scale(log(total_area_sqm))
         #,
         #percent_semi_natural_20km_scaled = center_scale(percent_semi_natural_20km)
  ) 

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity2.RDS")

# elements of simmed_diversity list are as follows:
mean_richness <- simmed_diversity[[10]] 

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm=TRUE)
mean_richness_quantiles_df <- as.data.frame(t(mean_richness_quantiles))
colnames(mean_richness_quantiles_df) <- c("lower90", "lower50",
                                          "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, mean_richness_quantiles_df)

#-------------------------------------------------------------------------------
# model relationship between city-wide predictors and mean richness

#-------------------------------------------------------------------------------
# first just model relationship between city-wide predictors and MEAN ESTIMATE of mean richness

## IMPORTANTLY this does not account for the uncertainty in our simulated richness predictions


#---------------------------------------------------------
m3.4 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + 
                             log_IIC_scaled + 
                             percent_semi_natural_scaled + 
                             log_total_area_scaled +
                             percent_tree, 
                           data = city_data)

summary(m3.4)
plot(m3.4)
pp_check(m3.4) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.4, pars = "log_park_size_scaled", prob = 0.9)

mcmc_areas(as.matrix(m3.4),prob_outer = .95)
mcmc_pairs(as.matrix(m3.4),pars = c("log_park_size_scaled","isolation_scaled","percent_semi_natural_scaled"))

#-------------------------------------------------------------------------------
# Instead we could try to determine which variables add 
# most to the predictive power of the model

fitg_cv <-  varsel(m3.4, method = "L1", nterms_max = 2, nclusters_pred = 10,
                   seed = 5555)
# model size suggested by the program
plot(fitg_cv, stats = c('elpd', 'rmse'))
# And we get a LOO based recommendation for the model size to choose
(nsel <- suggest_size(fitg_cv, alpha=0.5))
(vsel <- ranking(fitg_cv)[[1]][1:nsel])

#---------------------------------------------------------
m3.4.1 <- rstanarm::stan_glm(mean ~ 
                               log_IIC_scaled + 
                               percent_tree_scaled + 
                               log_total_area_scaled, 
                             data = city_data)

plot(m3.4.1)
pp_check(m3.4.1) # ?bayesplot::ppc_hist

# 80% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.4.1, pars = "percent_tree_scaled", prob = 0.9)
plot(m3.4.1, "areas", pars = "percent_tree_scaled", prob = 0.9)

#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# I used the mean response estimates to conduct the variable selection
# But we actually have uncertainty in those diversity estimates
# Let's propagate the uncretainty into the submodel chosen by varsel()

# get the full matrix of estimates richness values (not just the quantile summaries)
mean_richness <- simmed_diversity[[10]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 20
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(mean_richness)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_total_area_scaled <- vector(length = n_subsamples*n_models)
percent_tree_scaled <- vector(length = n_subsamples*n_models)
log_IIC_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

for(i in 1:n_models){
  city_data$response <- mean_richness[,i]
  
  fit <- rstanarm::stan_glm(response ~ log_IIC_scaled + 
                              percent_tree_scaled + 
                              log_total_area_scaled, 
                            data = city_data)
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_IIC_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  percent_tree_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  log_total_area_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,4]
  sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
}

posterior_draws <- as.data.frame(cbind(intercept, log_IIC_scaled, percent_tree_scaled, log_total_area_scaled, sigma))
# plot the sample densities
mcmc_areas <- mcmc_areas(posterior_draws, 
                         pars = c("intercept", #"log_IIC_scaled", 
                                  "percent_tree_scaled", "log_total_area_scaled")) +
  #labs(title = "Posterior densities of samples from\nmodels fit to 100 simulated communities") +
  theme_classic() +
  scale_x_continuous(name = "Posterior Model Estimate") +
  scale_y_discrete(labels = c("Intercept", 
                              #"log(IIC - Connectivity)", 
                              "% Tree Cover", "Total Area of City")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(angle = 45))

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
  "San Fransisco",
  "St. Louis",
  "Tampa"
)

# plot on a predictive scale
base <- ggplot(city_data, aes(x = log_total_area_scaled, y = mean)) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("Total Area of City") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

n_draws <- 100 # draw n lines from the post-posterior
base <- base + geom_abline(
  aes(intercept = intercept, slope = log_total_area_scaled), 
  data = sample_n(posterior_draws, n_draws), 
  color = "grey", 
  alpha = 0.5
) + 
  geom_abline(intercept = median(posterior_draws[,1]), slope = median(posterior_draws[,4]), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

# plot on a predictive scale
base2 <- ggplot(city_data, aes(x = log_IIC_scaled, y = mean)) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("(IIC - Connectivity)") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))

n_draws <- 100 # draw n lines from the post-posterior
base2 <- base2 + geom_abline(
  aes(intercept = intercept, slope = log_IIC_scaled), 
  data = sample_n(posterior_draws, n_draws), 
  color = "grey", 
  alpha = 0.5
) + 
  geom_abline(intercept = median(posterior_draws[,1]), slope = median(posterior_draws[,2]), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

# plot on a predictive scale
base3 <- ggplot(city_data, aes(x = percent_tree_scaled, y = mean)) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("% Tree Cover") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))

n_draws <- 100 # draw n lines from the post-posterior
base3 <- base3 + geom_abline(
  aes(intercept = intercept, slope = percent_tree_scaled), 
  data = sample_n(posterior_draws, n_draws), 
  color = "grey", 
  alpha = 0.5
) + 
  geom_abline(intercept = median(posterior_draws[,1]), slope = median(posterior_draws[,3]), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

m3.4_plot <- cowplot::plot_grid(mcmc_areas, base, base3, ncol = 3, rel_widths = c(1, 1))

saveRDS(m3.4_plot, "./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.4_plot.rds")

# get legend
base4 <- ggplot(city_data, aes(x = log_total_area_scaled, y = mean)) +
  ylab("Total Number of Species\nOccurring in City Parks") +
  xlab("Total Area of City") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14)) + 
  geom_abline(intercept = median(posterior_draws[,1]), slope = median(posterior_draws[,4]), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

legend <- cowplot::get_legend(base4)
plot(legend)

saveRDS(legend, "./part3_citywide_drivers_of_diversity/figures/m3_plots/m3_legend.rds")
