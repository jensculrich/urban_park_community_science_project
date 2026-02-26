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
  select(city, percent_grass_shrub, percent_tree, percent_semi_natural)
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
         percent_tree_scaled = center_scale(percent_tree),
         isolation_scaled = center_scale(mean_isolation),
         percent_semi_natural_scaled = center_scale(percent_semi_natural)
         #,
         #percent_semi_natural_20km_scaled = center_scale(percent_semi_natural_20km)
  ) 

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

mean_prop_disturbance_avoidant <- simmed_diversity[[3]] 

#  calculate Means and CI's for the diversity metrics for each city
mean_prop_disturbance_avoidant_quantiles <- apply(mean_prop_disturbance_avoidant, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
mean_prop_disturbance_avoidant_quantiles_df <- as.data.frame(t(mean_prop_disturbance_avoidant_quantiles))
colnames(mean_prop_disturbance_avoidant_quantiles_df) <- c("lower90", "lower50",
                                                           "mean", "upper50", "upper90")

#-------------------------------------------------------------------------------
# and then add to the city data df

city_data <- cbind(city_data, mean_prop_disturbance_avoidant_quantiles_df)

city_data <- city_data %>%
  # beta reg doesn't link boundary cases at absolute 0
  mutate(mean_adj = mean + 0.0001)

#---------------------------------------------------------
# full model

m3.2 <- rstanarm::stan_betareg(mean_adj ~ log_park_size_scaled + log_IIC_scaled + percent_semi_natural_scaled, 
                           data = city_data, link = "logit", link.phi = "log" )

summary(m3.2)
plot(m3.2)
pp_check(m3.2) # ?bayesplot::ppc_hist

mcmc_areas(as.matrix(m3.2),prob_outer = .95)
mcmc_pairs(as.matrix(m3.2),pars = c("log_park_size_scaled","log_IIC_scaled","percent_semi_natural_scaled"))

##---------------------------------------------------------------
## try to compare loo scores 

# can't do it on a beta regression unforunately so keep the outcome a linear term for now
m3.2_linear <- rstanarm::stan_glm(mean_adj ~ log_park_size_scaled + log_IIC_scaled + percent_semi_natural_scaled, 
                               data = city_data)

fitg_cv <-  varsel(m3.2_linear, method = "L1", nterms_max = 2, nclusters_pred = 10,
                   seed = 5555)
# model size suggested by the program
plot(fitg_cv, stats = c('elpd', 'rmse'))
# And we get a LOO based recommendation for the model size to choose
(nsel <- suggest_size(fitg_cv, alpha=0.5))
(vsel <- ranking(fitg_cv)[[1]][1:nsel])

#---------------------------------------------------------
# alternatively fit the reduced model suggested by vsel separately

m3.2.4 <- rstanarm::stan_betareg(mean_adj ~ log_IIC_scaled + log_park_size_scaled, 
                                 data = city_data, link = "logit", link.phi = "log")

summary(m3.2.4)
plot(m3.2.4)
pp_check(m3.2.4) # ?bayesplot::ppc_hist

#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# I used the mean response estimates to conduct the variable selection
# But we actually have uncertainty in those diversity estimates
# Let's propagate the uncretainty into the submodel chosen by varsel()

# get the full matrix of estimates of diversity values (not just the quantile summaries)
mean_prop_disturbance_avoidant <- simmed_diversity[[3]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 20
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(mean_prop_disturbance_avoidant)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_park_size_scaled <- vector(length = n_subsamples*n_models)
log_IIC_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

for(i in 1:n_models){
  city_data$response <- mean_prop_disturbance_avoidant[,i]
  
  fit <- rstanarm::stan_betareg(mean_adj ~ log_IIC_scaled + log_park_size_scaled, 
                                data = city_data, link = "logit", link.phi = "log" )
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_IIC_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  log_park_size_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,4]
}

posterior_draws <- as.data.frame(cbind(intercept, log_park_size_scaled, log_IIC_scaled, sigma))
# plot the sample densities
mcmc_areas <- mcmc_areas(posterior_draws, 
                         pars = c("intercept", "log_park_size_scaled", "log_IIC_scaled")) +
  #labs(title = "Posterior Densities of Retained Predictors") +
  theme_classic()  +
  scale_x_continuous(name = "Posterior Model Estimate") +
  scale_y_discrete(labels = c("Intercept", "Median log(Park Size)", "log(IIC - Connectivity)")) +
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

ilogit <- function(x){exp(x)/(1+exp(x))}

# plot on a predictive scale
base <- ggplot(city_data, aes(x = log_park_size_scaled, y = mean)) +
  ylab("(%) Disturbance Avoidant Species") +
  xlab("Median log(Park Size)") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


x_plot <- seq(-2, 2, length.out=22)
y_plot <- ilogit(median(posterior_draws[,1]) + median(posterior_draws[,2]) * x_plot)

#plot_data <- data.frame(x_plot, y_plot)

n_draws <- 100
rows <- sample_n(posterior_draws, n_draws)
for(i in 1:n_draws){
  
  y <- ilogit(rows[i,1] + rows[i,2] * x_plot)
  plot_data <- data.frame(x_plot, y)
  base <- base + 
    geom_line(data=plot_data,aes(x_plot, y), 
              color = "grey", alpha = 0.2)
}


base <- base + 
  geom_line(aes(x_plot, y_plot), 
              color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

# plot on a predictive scale
base2 <- ggplot(city_data, aes(x = log_IIC_scaled, y = mean)) +
  ylab("(%) Disturbance Avoidant Species") +
  xlab("log(IIC - Connectivity)") +
  scale_color_manual(values=my_palette, labels = city_names_labels, name="City") + 
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


x_plot2 <- seq(-3.5, 2, length.out=22)
y_plot2 <- ilogit(median(posterior_draws[,1]) + median(posterior_draws[,3]) * x_plot)

#plot_data <- data.frame(x_plot, y_plot)

n_draws <- 100
rows <- sample_n(posterior_draws, n_draws)
for(i in 1:n_draws){
  
  y <- ilogit(rows[i,1] + rows[i,3] * x_plot2)
  plot_data <- data.frame(x_plot2, y)
  base2 <- base2 + 
    geom_line(data=plot_data,aes(x_plot2, y), 
              color = "grey", alpha = 0.2)
}


base2 <- base2 + 
  geom_line(aes(x_plot2, y_plot2), 
            color = "black", size = 2) + 
  geom_errorbar(aes(ymin = lower50, ymax=upper50, colour=city), size=2) +
  geom_errorbar(aes(ymin = lower90, ymax=upper90, colour=city), size=1) +
  geom_point(aes(colour=city), size = 4)

m3.2_plot <- cowplot::plot_grid(mcmc_areas, base, base2, ncol = 3)

saveRDS(m3.2_plot, "./part3_citywide_drivers_of_diversity/figures/m3_plots/m3.2_plot.rds")
