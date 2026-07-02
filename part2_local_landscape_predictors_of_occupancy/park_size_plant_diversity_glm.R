library(tidyverse)
library(rstanarm)
library(bayesplot)

# get prepared urban park site data
my_data <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))
site_data <- my_data$site_data %>%
  select(city, log_total_green_space_area, log_total_green_space_area_scaled,
         log_n_plant_genera, plant_genera_density_scaled) %>%
  mutate(city = as.factor(city)) %>%
  mutate(n_plant_genera = as.integer(exp(log_n_plant_genera))) %>%
  rename("park_size" = "log_total_green_space_area_scaled")

# raw data
plot(site_data$n_plant_genera ~ site_data$park_size)
# log scale plant genera
plot(site_data$log_n_plant_genera ~ site_data$park_size)

hist(site_data$n_plant_genera)

# using a poisson regression for the relationship between park size and the number of flowering plant genera
options(mc.cores=4)
m_plant_div <- rstanarm::stan_glmer(n_plant_genera ~ park_size + 
                             (1|city), 
                           data = site_data, 
                           family = poisson(link="log"))

summary(m_plant_div)
# combine draws into a df
posterior_draws <- as.data.frame(m_plant_div)
quantile(posterior_draws$park_size, prob = c(0.05, 0.25, 0.5, 0.75, 0.95))
plot(m_plant_div, pars = c("park_size"))
pp_check(m_plant_div)

##-----------------------------------------------------------------------------
# make a classic counterfactual interval band plot
# now plot the first trend

sd_size <- sd(site_data$log_total_green_space_area)
mean_size <- mean(site_data$log_total_green_space_area)

pred = seq(-2.5, 2.5, length.out=100)
x <- pred * sd_size + mean_size

pred_data <- as.data.frame(pred)

n_draws <- 1000 # draw n lines from the post-posterior
predictions <- matrix(nrow = nrow(pred_data), ncol = n_draws)
sampled_posterior <- sample_n(posterior_draws, n_draws)

# for each value of pred data
for(i in 1:nrow(predictions)){
  # draw a potential relationship, and predict the outcome given the pred value
  for(j in 1:n_draws){
    
    predictions[i,j] <- exp(sampled_posterior[j, 1] + sampled_posterior[j,2] * pred_data$pred[i])
    
  }
}

#  calculate Means and CI's for the diversity metrics for each city
prediction_quantiles <- apply(predictions, MARGIN = 1, FUN = quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
prediction_quantiles_df <- as.data.frame(t(prediction_quantiles))
colnames(prediction_quantiles_df) <- c("lower90", "lower50",
                                       "mean", "upper50", "upper90")

prediction_quantiles_df <- as.data.frame(cbind(pred_data, prediction_quantiles_df, x))

ylim_lower <- 0
#temp <- as.data.frame(cbind(park_size_data$median_log_park_size, ylim_lower))
# plot on a predictive scale
p <- ggplot(prediction_quantiles_df, aes(x = x, y = mean)) +
  ylab("Number of Plant Genera Detected") +
  theme_classic() + 
  geom_point(data=site_data, aes(log_total_green_space_area, n_plant_genera), size = 1, alpha = 0.5, colour="#A25050") +
  geom_ribbon(aes(ymin=lower90, ymax=upper90), alpha=0.25) +
  geom_ribbon(aes(ymin=lower50, ymax=upper50), alpha = 0.5) +
  geom_line() +
  ylim(ylim_lower, 500) +
  scale_x_continuous(name=expression(paste("log(Park Size (m"^2, "))")), 
                     limits=c(min(x), max(x))) +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))
p

