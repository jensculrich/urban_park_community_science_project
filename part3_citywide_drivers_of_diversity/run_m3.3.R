library(tidyverse)
library(rstanarm)
library(bayesplot)

center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

ilogit <- function(x){exp(x)/(1+exp(x))}

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
         log_IIC_scaled = center_scale(log_IIC),
         percent_tree_scaled = center_scale(percent_tree),
         percent_grassshrub_scaled = center_scale(percent_grass_shrub),
         isolation_scaled = center_scale(mean_isolation),
         percent_semi_natural_scaled = center_scale(percent_semi_natural)
  ) 

#-------------------------------------------------------------------------------
# get the city-wide diversity predictions

simmed_diversity <- readRDS("./part3_citywide_drivers_of_diversity/simmed_diversity.RDS")

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
# model relationship between city-wide predictors and mean richness

#-------------------------------------------------------------------------------
# first just model relationship between city-wide predictors and MEAN ESTIMATE of mean richness

## IMPORTANTLY this does not account for the uncertainty in our simulated richness predictions


#---------------------------------------------------------
m3.3 <- rstanarm::stan_betareg(mean ~ log_park_size_scaled + 
                             log_IIC_scaled + 
                             percent_tree_scaled + 
                             percent_grassshrub_scaled, 
                           data = city_data, link = "logit", link.phi = "log")

summary(m3.3)
plot(m3.3)
posterior_interval(m3.3, pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"),
                   prob = 0.5)
pp_check(m3.3) # fit is really not good!

m3.3 <- rstanarm::stan_glm(mean ~ log_park_size_scaled + 
                                 log_IIC_scaled + 
                                 percent_tree_scaled + 
                                 percent_grassshrub_scaled, 
                               data = city_data)

summary(m3.3)
plot(m3.3)
posterior_interval(m3.3, pars = c("log_park_size_scaled","log_IIC_scaled","percent_tree_scaled", "percent_grassshrub_scaled"),
                   prob = 0.5)
pp_check(m3.3) # fit is much better!
# percent_grassshrub_scaled has 50% BCI overlapping with zero
# i.e., has little clear effect on the diversity responses
# in addition removing them does not introduce a confound according to a causal diagram of the system
# i.e., we do not expect a causal pathway from other predictors through these ones into diversity
# and so it makes sense to remove them from the causal model.


#---------------------------------------------------------
# Finally, incorporate uncertainty - 
# There was some uncertainty in those diversity estimates
# Let's propagate the uncertainty into the posthoc regression model

# get the full matrix of estimates richness values (not just the quantile summaries)
beta_diversity <- simmed_diversity[[7]] 

# how many samples to take from each model fit to each draw of the response data?
n_subsamples <- 40
# how many models to fit? I'll fit one for every simulated set of diversity responses
n_models <- ncol(beta_diversity)
# empty vectors of param values to fill
intercept <- vector(length = n_subsamples*n_models)
log_park_size_scaled <- vector(length = n_subsamples*n_models)
percent_tree_scaled <- vector(length = n_subsamples*n_models)
percent_grassshrub_scaled <- vector(length = n_subsamples*n_models)
log_IIC_scaled <- vector(length = n_subsamples*n_models)
sigma <- vector(length = n_subsamples*n_models)

set.seed(1)
for(i in 1:n_models){
  city_data$response <- beta_diversity[,i]
  
  fit <- rstanarm::stan_glm(response ~ log_park_size_scaled + 
                                  log_IIC_scaled + 
                                  percent_tree_scaled, #+ 
                                  #percent_grassshrub_scaled
                                data = city_data)
  
  draws <- as.data.frame(fit)
  
  sample_rows <- sample(1:nrow(draws), n_subsamples, replace = FALSE)
  
  intercept[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,1]
  log_park_size_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,2]
  log_IIC_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,3]
  percent_tree_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,4]
  #percent_grassshrub_scaled[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
  sigma[(n_subsamples*(i-1)+1):((n_subsamples*(i-1))+n_subsamples)] <- draws[sample_rows,5]
}

posterior_draws <- as.data.frame(cbind(
  intercept, log_park_size_scaled, percent_tree_scaled, 
  percent_grassshrub_scaled, log_IIC_scaled, sigma))

quantile(posterior_draws$intercept, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$log_park_size_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$percent_tree_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$percent_grassshrub_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
quantile(posterior_draws$log_IIC_scaled, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))

##-----------------------------------------------------------------------------
# make a classic counterfactual interval band plot
# now plot the first trend

sd_size <- sd(city_data$median_log_park_size)
mean_size <- mean(city_data$median_log_park_size)
pred = seq(-2, 2, length.out=100)
x <- pred * sd_size + mean_size

pred_data <- as.data.frame(cbind(x, pred))

n_draws <- 2000 # draw n lines from the post-posterior
predictions <- matrix(nrow = nrow(pred_data), ncol = n_draws)
sampled_posterior <- sample_n(posterior_draws, n_draws)

# for each value of pred data
for(i in 1:nrow(predictions)){
  # draw a potential relationship, and predict the outcome given the pred value
  for(j in 1:n_draws){
    
    predictions[i,j] <- #ilogit(
      sampled_posterior[j, 1] + sampled_posterior[j,2] * pred_data[i,2]
    #)
  }
}

#  calculate Means and CI's for the diversity metrics for each city
prediction_quantiles <- apply(predictions, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
prediction_quantiles_df <- as.data.frame(t(prediction_quantiles))
colnames(prediction_quantiles_df) <- c("lower90", "lower50",
                                       "mean", "upper50", "upper90")

prediction_quantiles_df <- as.data.frame(cbind(pred_data, prediction_quantiles_df))

ylim_lower <- 0.15
temp <- as.data.frame(cbind(park_size_data$median_log_park_size, ylim_lower))
# plot on a predictive scale
p <- ggplot(prediction_quantiles_df, aes(x = x, y = mean)) +
  ylab("Community Dissimilarity\n(Jaccard Index)") +
  theme_classic() + 
  geom_ribbon(aes(ymin=lower90, ymax=upper90), alpha=0.25) +
  geom_ribbon(aes(ymin=lower50, ymax=upper50), alpha = 0.5) +
  geom_line() +
  ylim(ylim_lower, 0.4) +
  geom_point(data=temp, aes(V1, ylim_lower), shape = "|", size = 10, colour="#A25050") +
  scale_x_continuous(name=expression(paste("Median log(Park Size (m"^2, "))")), 
                     limits=c(min(x), max(x))) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
p

##------------------------------------------------------------------------------
# now plot the second trend

sd_IIC <- sd(city_data$log_IIC)
mean_IIC <- mean(city_data$log_IIC)
pred = seq(-3.25, 2, length.out=100)
x <- pred * sd_IIC + mean_IIC

pred_data <- as.data.frame(cbind(x, pred))

n_draws <- 2000 # draw n lines from the post-posterior
predictions <- matrix(nrow = nrow(pred_data), ncol = n_draws)
sampled_posterior <- sample_n(posterior_draws, n_draws)

# for each value of pred data
for(i in 1:nrow(predictions)){
  # draw a potential relationship, and predict the outcome given the pred value
  for(j in 1:n_draws){
    
    predictions[i,j] <- #ilogit(
      sampled_posterior[j, 1] + sampled_posterior[j,5] * pred_data[i,2]#)
    
  }
}

#  calculate Means and CI's for the diversity metrics for each city
prediction_quantiles <- apply(predictions, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
prediction_quantiles_df <- as.data.frame(t(prediction_quantiles))
colnames(prediction_quantiles_df) <- c("lower90", "lower50",
                                       "mean", "upper50", "upper90")

prediction_quantiles_df <- as.data.frame(cbind(pred_data, prediction_quantiles_df))

ylim_lower <- 0.15
temp <- as.data.frame(cbind(log(IIC_connectivity_data$IIC), ylim_lower))
# plot on a predictive scale
p2 <- ggplot(prediction_quantiles_df, aes(x = x, y = mean)) +
  ylab("Community Dissimilarity\n(Jaccard Index)") +
  xlab("log(IIC - Connectivity)") +
  ylim(ylim_lower, 45) +
  scale_x_continuous(limits = c(min(x), max(x))) + 
  theme_classic() + 
  geom_ribbon(aes(ymin=lower90, ymax=upper90), alpha=0.25) +
  geom_ribbon(aes(ymin=lower50, ymax=upper50), alpha = 0.5) +
  geom_line() +
  ylim(ylim_lower, 0.4) +
  geom_point(data=temp, aes(V1, ylim_lower), shape = "|", size = 10, colour="#A25050") +
  theme(legend.position = c(0.025, 0.975), # x=1 (right), y=0 (bottom)
        legend.justification = c(0, 1), # Justify the bottom-right corner of the legend box to these coordinates
        legend.text = element_text(size=14),
        legend.title = element_text(size=16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
p2

figure5.2 <- cowplot::plot_grid(p, p2, ncol = 2,
                                labels = c("c)", "d)"), 
                                label_size = 16)
figure5.2

saveRDS(figure5.2, "./part3_citywide_drivers_of_diversity/figures/m3_plots/figure5.2.rds")


# plot the sample densities
mcmc_areas <- mcmc_areas(posterior_draws, 
                         pars = c("log_park_size_scaled", 
                                  "percent_tree_scaled", 
                                  #"percent_grassshrub_scaled", 
                                  "log_IIC_scaled")) +
  #labs(title = "Posterior Densities of Retained Predictors") +
  theme_classic() +
  scale_x_continuous(name = "Posterior Model Estimate (logit-scaled)") +
  scale_y_discrete(labels = c("Median log(Park Size)", "% Tree Cover", 
                              #"% Grass/Shrub Cover", 
                              "log(IIC)")) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 14))
mcmc_areas