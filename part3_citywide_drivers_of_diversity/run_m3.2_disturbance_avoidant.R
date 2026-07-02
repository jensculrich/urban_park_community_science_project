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

#-------------------------------------------------------------------------------
# get the city covariate data

city_data <- read.csv("./data/city_wide_data/derived_city_wide_data.csv")

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

#---------------------------------------------------------
# full model

m3.2 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + 
                                 log_IIC_scaled + 
                                 percent_tree_scaled + 
                                 percent_grassshrub_scaled, 
                           data = city_data, 
                           family = Gamma(link="log"))

summary(m3.2)
posterior_interval(m3.2, pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"),
                   prob = 0.5)
# percent_grassshrub_scaled, log_IIC_scaled,  has 50% BCI overlapping with zero
# i.e., has little clear effect on the diversity responses
# in addition removing them does not introduce a confound according to a causal diagram of the system
# i.e., we do not expect a causal pathway from other predictors through these ones into diversity
# and so it makes sense to remove them from the causal model.
plot(m3.2)
pp_check(m3.2) #+ xlim(c(0, 0.01)) # ?bayesplot::ppc_hist
#mcmc_areas(as.matrix(m3.2),prob_outer = .95)
#mcmc_pairs(as.matrix(m3.2),pars = c("log_park_size_scaled","log_IIC_scaled","percent_semi_natural_scaled"))

#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# There was some uncertainty in those diversity estimates
# Let's propagate the uncertainty into the posthoc regression model

# get the full matrix of estimates of diversity values (not just the quantile summaries)
mean_disturbance_or_edge_avoidant <- simmed_diversity[[14]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 40
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(mean_disturbance_or_edge_avoidant)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_park_size_scaled <- vector(length = n_subsamples*n_models)
percent_tree_scaled <- vector(length = n_subsamples*n_models)
percent_grassshrub_scaled <- vector(length = n_subsamples*n_models)
log_IIC_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

set.seed(1)

for(i in 1:n_models){
  city_data$response <- mean_disturbance_or_edge_avoidant[,i]
  #city_data$response_adj <- city_data$response + 0.0001 # beta reg doesn't link boundary cases at absolute 0
  
  
  fit <- rstanarm::stan_glm(response ~ log_park_size_scaled + 
                              percent_tree_scaled, 
                            data = city_data, 
                            family = Gamma(link="log"))
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_park_size_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  #log_IIC_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  percent_tree_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  #percent_grassshrub_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
  #sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
}

posterior_draws <- as.data.frame(cbind(
  intercept, log_park_size_scaled, percent_tree_scaled, 
  percent_grassshrub_scaled, log_IIC_scaled, sigma))

quantile(posterior_draws$intercept, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$log_park_size_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$percent_tree_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$percent_grassshrub_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$log_IIC_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))


# plot the sample densities
mcmc_areas <- mcmc_areas(posterior_draws, 
                         pars = c("intercept",
                                  "log_park_size_scaled",
                                  "percent_tree_scaled"
                                  #"percent_grassshrub_scaled",
                                  #"log_IIC_scaled"
                                  )) +
  #labs(title = "Posterior Densities of Retained Predictors") +
  theme_classic() +
  scale_x_continuous(name = "Posterior Model Estimate (log-scaled)") +
  scale_y_discrete(labels = c("Intercept",
                              "Median log(Park Size)", 
                              "% Tree Cover"
                              #"% Grass/Shrub Cover", 
                              #"log(IIC)"
                              )) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(angle = 45))

mcmc_areas
