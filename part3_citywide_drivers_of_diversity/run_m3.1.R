library(tidyverse)
library(rstanarm)
library(bayesplot)

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

#my_palette <- viridis::viridis(n=n_cities+2, option = "turbo")
#my_palette <- my_palette[3:(n_cities+2)] # remove the really dark colours

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
  select(city, percent_grass_shrub, percent_tree, percent_semi_natural)

# join the data
city_data <- park_size_data %>%
  left_join(., connectivity_data, by ="city") %>%
  left_join(., IIC_connectivity_data, by ="city") %>%
  left_join(., landcover_data, by = "city") 

city_data <- city_data[order(city_data$city), ]

city_data <- city_data %>%
  cbind(., city_factor = seq(1:n_cities)) %>% 
  mutate(log_IIC = log(IIC), 
         log_park_size_scaled = center_scale(median_log_park_size),
         percent_tree_scaled = center_scale(percent_tree),
         percent_grassshrub_scaled = center_scale(percent_grass_shrub),
         log_IIC_scaled = center_scale(log_IIC),
         isolation_scaled = center_scale(mean_isolation),
         percent_semi_natural_scaled = center_scale(percent_semi_natural)
  ) 

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity2.RDS")

# elements of simmed_diversity list are as follows:
#list(mean_richness, median_richness, median_richness_prop,
#mean_prop_disturbance_avoidant, mean_prop_edge_avoidant, mean_prop_disturbance_or_edge_avoidant,
#beta_diversity, beta_repl, beta_richdif,
#gamma_diversity, gamma_diversity_prop)
mean_richness <- simmed_diversity[[2]] # 2 is median not the mean

#  calculate Means and CI's for the diversity metrics for each city
mean_richness_quantiles <- apply(mean_richness, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
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
m3.1 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + 
                             log_IIC_scaled + 
                             percent_tree_scaled + 
                             percent_grassshrub_scaled, 
                           data = city_data)

summary(m3.1)
plot(m3.1)
pp_check(m3.1) # ?bayesplot::ppc_hist

# 90% interval of estimated reciprocal_dispersion parameter
posterior_interval(m3.1, pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"),
                                  prob = 0.9)
posterior_interval(m3.1, pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"),
                  prob = 0.5)

mcmc_areas(as.matrix(m3.1),prob_outer = .95)
mcmc_pairs(as.matrix(m3.1),pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"))

#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# There was some uncertainty in those diversity estimates
# Let's propagate the uncertainty into the posthoc regression model

# get the full matrix of estimates richness values (not just the quantile summaries)
mean_richness <- simmed_diversity[[2]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 20
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(mean_richness)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_park_size_scaled <- vector(length = n_subsamples*n_models)
percent_tree_scaled <- vector(length = n_subsamples*n_models)
percent_grassshrub_scaled <- vector(length = n_subsamples*n_models)
log_IIC_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

for(i in 1:n_models){
  city_data$response <- mean_richness[,i]
  
  fit <- rstanarm::stan_glm(response ~ log_park_size_scaled + 
                              percent_tree_scaled +
                              percent_grassshrub_scaled +
                              log_IIC_scaled, 
                            data = city_data)
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_park_size_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  percent_tree_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  percent_grassshrub_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,4]
  log_IIC_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
  sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,6]
}

posterior_draws <- as.data.frame(cbind(
  intercept, log_park_size_scaled, percent_tree_scaled, 
  percent_grassshrub_scaled, log_IIC_scaled, sigma))
# plot the sample densities
mcmc_areas <- mcmc_areas(posterior_draws, 
           pars = c("log_park_size_scaled", 
                    "percent_tree_scaled", "percent_grassshrub_scaled", "log_IIC_scaled")) +
  #labs(title = "Posterior Densities of Retained Predictors") +
  theme_classic() +
  scale_x_continuous(name = "Posterior Model Estimate") +
  scale_y_discrete(labels = c("Median log(Park Size)", "% Tree Cover", "% Grass/Shrub Cover", "log(IIC)")) +
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

##------------------------------------------------------------------------------
# now plot the first trend

sd_size <- sd(city_data$median_log_park_size)
mean_size <- mean(city_data$median_log_park_size)
pred = seq(-2, 2, length.out=100)
x <- pred * sd_size + mean_size

pred_data <- as.data.frame(cbind(x, pred))

n_draws <- 100 # draw n lines from the post-posterior
predictions <- matrix(nrow = nrow(pred_data), ncol = n_draws)
sampled_posterior <- sample_n(posterior_draws, n_draws)

# for each value of pred data
for(i in 1:nrow(predictions)){
  # draw a potential relationship, and predict the outcome given the pred value
  for(j in 1:n_draws){
    
    predictions[i,j] <- sampled_posterior[j, 1] + sampled_posterior[j,2] * pred_data[i,2]
    
  }
}

y <- rowMeans(predictions)

new_dat <- as.data.frame(cbind(pred_data, y, predictions))

# plot on a predictive scale
base <- ggplot(new_dat, aes(x = x, y = y)) +
  ylab("Median Park Species Richness") +
  xlab(expression(paste("Median log(Park Size (m"^2, "))"))) +
  theme_classic() + 
  ylim(5, 50) +
  xlim(c(min(x), max(x))) +
  theme(legend.position = c(0.025, 0.975), # x=1 (right), y=0 (bottom)
        legend.justification = c(0, 1), # Justify the bottom-right corner of the legend box to these coordinates
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

line_layers <- list()
for (i in 4:ncol(new_dat)) {
  temp_data <- data.frame(x = new_dat$x, y = new_dat[,i])
  # Create a geom_line layer and add it to the list
  line_layers[[i-3]] <- geom_line(data = temp_data, colour="grey", alpha = 0.5)
}

base <- base + line_layers + 
  geom_line(size = 2, colour = "black")

base <- base + 
  geom_errorbar(data = city_data, aes(x=median_log_park_size, y=mean, ymin=lower50, ymax=upper50, colour=percent_tree), size=2) +
  geom_errorbar(data = city_data, aes(x=median_log_park_size, y=mean, ymin = lower90, ymax=upper90, colour=percent_tree), size=1) +
  geom_point(data = city_data, aes(x=median_log_park_size, y=mean, colour=percent_tree), size = 4) +
  scale_colour_viridis_c(name="City-Wide Tree Cover", labels = scales::label_percent()) 

##------------------------------------------------------------------------------
# now plot the second trend

sd_tree <- sd(city_data$percent_tree)
mean_tree <- mean(city_data$percent_tree)
pred = seq(-1.5, 2, length.out=100) # tree cover can't be lower than a certain amount
x <- pred * sd_tree + mean_tree

pred_data <- as.data.frame(cbind(x, pred))

n_draws <- 100 # draw n lines from the post-posterior
predictions <- matrix(nrow = nrow(pred_data), ncol = n_draws)
sampled_posterior <- sample_n(posterior_draws, n_draws)

# for each value of pred data
for(i in 1:nrow(predictions)){
  # draw a potential relationship, and predict the outcome given the pred value
  for(j in 1:n_draws){
    
    predictions[i,j] <- sampled_posterior[j, 1] + sampled_posterior[j,3] * pred_data[i,2]
    
  }
}

y <- rowMeans(predictions)

new_dat <- as.data.frame(cbind(pred_data, y, predictions))

# plot on a predictive scale
base2 <- ggplot(new_dat, aes(x = x, y = y)) +
  ylab("Median Park Species Richness") +
  xlab("City-Wide Tree Cover") +
  ylim(5, 50) +
  scale_x_continuous(labels = scales::label_percent(), limits = c(min(x), max(x))) + 
  theme_classic() + 
  theme(legend.position = c(0.025, 0.975), # x=1 (right), y=0 (bottom)
        legend.justification = c(0, 1), # Justify the bottom-right corner of the legend box to these coordinates
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

line_layers <- list()
for (i in 4:ncol(new_dat)) {
  temp_data <- data.frame(x = new_dat$x, y = new_dat[,i])
  # Create a geom_line layer and add it to the list
  line_layers[[i-3]] <- geom_line(data = temp_data, colour="grey", alpha = 0.5)
}

base2 <- base2 + line_layers + 
  geom_line(size = 2, colour = "black")

base2 <- base2 + 
  geom_errorbar(data = city_data, aes(x=percent_tree, y=mean, ymin=lower50, ymax=upper50, colour=median_log_park_size), size=2) +
  geom_errorbar(data = city_data, aes(x=percent_tree, y=mean, ymin = lower90, ymax=upper90, colour=median_log_park_size), size=1) +
  geom_point(data = city_data, aes(x=percent_tree, y=mean, colour=median_log_park_size), size = 4) +
  scale_colour_viridis_c(name=expression(paste("Median log(Park Size (m"^2, "))"))) 


##------------------------------------------------------------------------------
# combine the panels

m3.1_plot <- cowplot::plot_grid(mcmc_areas, base, base2, ncol = 3,
                                labels = c("a)", "b)", "c)"), 
                                label_size = 16)
m3.1_plot

saveRDS(m3.1_plot, "./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.1_plot.rds")
