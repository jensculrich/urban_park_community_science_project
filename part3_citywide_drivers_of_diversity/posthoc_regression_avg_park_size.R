# We want to know whether city-wide characteristics predict variation in importance of
# park-specific characteristics (landscape isolation and landscape veg cover)
# and also the overall occurrence rate

library(tidyverse)
library(cmdstanr)

# predictor center scaling function
center_scale <- function(x) {
  (x - mean(x)) / sd(x)
}

##------------------------------------------------------------------------------
# get model fit

## get param estimates from m2.1
stan_out_m2.1 <- readRDS(
  "./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1_jan2.rds")

# summarise all variables with default and additional summary measures
estimates <- stan_out_m2.1$draws(
  variables = c(
    "psi_0", 
    "sigma_psi_species",
    "sigma_psi_city",
    "mu_psi_wingspan",
    "sigma_psi_wingspan",
    "mu_psi_park_size",
    "sigma_psi_park_size",
    "mu_psi_tree_cover",
    "sigma_psi_tree_cover",
    "mu_psi_plant_diversity",
    "sigma_psi_plant_diversity",
    "mu_psi_landscape_isolation",
    "sigma_psi_landscape_isolation",
    "mu_psi_landscape_grassherb",
    "sigma_psi_landscape_grassherb",
    "mu_psi_landscape_woody",
    "sigma_psi_landscape_woody",
    
    "p0", 
    "sigma_p_species",
    "sigma_p_city",
    "p_city_detections",
    "p_wingspan",
    "p_feature_diversity",
    "p_ease_of_id",
    "delta0",
    "delta_regional_cluster",
    "sigma_p_species_date",
    "epsilon0",
    "epsilon_regional_cluster",
    "sigma_p_species_date_sq",
    
    # city effects
    "psi_city",
    "psi_wingspan",
    "psi_park_size",
    "psi_tree_cover",
    "psi_plant_diversity",
    "psi_landscape_isolation",
    "psi_landscape_grassherb",
    "psi_landscape_woody",
    "p_city"),
  
  format = "draws_matrix"
)

# clear space
rm(stan_out_m2.1)
gc()

##------------------------------------------------------------------------------
# read and modify the city wide data

# get the city covariate data
city_data <- read.csv("./data/city_wide_data/all_cities_average_park_size_classified_parks_only.csv")

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

city_data <- city_data %>%
  filter(city %in% city_names) %>%
  cbind(., city_factor = seq(1:length(n_cities)))

city_data <- city_data %>%
  mutate(log_avg_park_size = log(average_park_size_sqm),
         log_park_size_scaled = center_scale(log_avg_park_size))

##------------------------------------------------------------------------------
# post hoc regression psi_0

# get indices for random effects distributions for any particular city
first_psi <- which( colnames(estimates)=="psi_city[1]" )

# for now just take the means
mean_psi0 <- vector(length = n_cities)
  
for(city_number in 1:n_cities){
  mean_psi0[city_number] <- mean(estimates[,city_number + (first_psi - 1)])
}  
  
df <- cbind(city_data, mean_psi0) 

ggplot(df, aes(log_avg_park_size , mean_psi0)) +
  geom_point(aes(colour=city), size = 3) +
  geom_smooth(method = lm)

##------------------------------------------------------------------------------
# post hoc regression psi_park_szie

# get indices for random effects distributions for any particular city
first_psi_size <- which( colnames(estimates)=="psi_park_size[1]" )

# for now just take the means
mean_psi_size <- vector(length = n_cities)

for(city_number in 1:n_cities){
  mean_psi_size[city_number] <- mean(estimates[,city_number + (first_psi_size - 1)])
}  

df <- cbind(city_data, mean_psi_size) 

ggplot(df, aes(log_avg_park_size , mean_psi_size)) +
  geom_point(aes(colour=city), size = 3)  +
  geom_smooth(method = lm)

##------------------------------------------------------------------------------
# post hoc regression psi_park_isolation

# get indices for random effects distributions for any particular city
first_psi_isolation <- which( colnames(estimates)=="psi_landscape_isolation[1]" )

# for now just take the means
mean_psi_isolation <- vector(length = n_cities)

for(city_number in 1:n_cities){
  mean_psi_isolation[city_number] <- mean(estimates[,city_number + (first_psi_isolation - 1)])
}  

df <- cbind(city_data, mean_psi_isolation) 

ggplot(df, aes(log_avg_park_size , mean_psi_isolation)) +
  geom_point(aes(colour=city), size = 3)  +
  geom_smooth(method = lm)

##------------------------------------------------------------------------------
# post hoc regression psi_park_isolation

# get indices for random effects distributions for any particular city
first_psi_grassherb <- which( colnames(estimates)=="psi_landscape_grassherb[1]" )

# for now just take the means
mean_psi_grassherb <- vector(length = n_cities)

for(city_number in 1:n_cities){
  mean_psi_grassherb[city_number] <- mean(estimates[,city_number + (first_psi_grassherb - 1)])
}  

df <- cbind(city_data, mean_psi_grassherb) 

ggplot(df, aes(log_avg_park_size , mean_psi_grassherb)) +
  geom_point(aes(colour=city), size = 3) +
  geom_smooth(method = lm)
