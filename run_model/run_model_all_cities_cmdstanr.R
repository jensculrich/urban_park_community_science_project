
# select a region
regions <- c(
  "northeast",
  "southeast",
  "texas",
  "california",
  "all"
)

region <- regions[5]

# list of city names

# all
if(region == regions[5]){
  city_names <- c(
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
  )
}


my_data <- readRDS( paste0("./run_model/prepped_data/prepped_data_", region, ".rds"))

# prepare to fit occupncy model
# library(rstan)
library(cmdstanr)

# data to feed to the model
# detection data
V <- my_data$V # detections (1==detected)
V_NA <- my_data$V_NA # NAs (0==no known survey effort made)

# covariate data
species_info <- my_data$species_info
n_species <- my_data$n_species # number of species
species <- as.integer(as.factor(species_info$species))
site_data <- my_data$site_data
n_sites <- my_data$n_sites # number of sites
sites <- site_data$multicity_site_id
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
years <- seq(1, n_years_minus1)
n_surveys <- my_data$n_surveys
surveys <- sequence(n_surveys)
surveys <- (surveys - mean(surveys)) / sd(surveys)

migratory <- species_info$migratory
migratory[is.na(migratory)] <- 0

## predictors
# species
feature_diversity <- species_info$featureDiversity_scaled
ease_of_id <- species_info$research_grade_proportion_scaled
wingspan <- species_info$aveWingspan_scaled
# site
park_size <- site_data$log_total_green_space_area_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
isolation <- site_data$log_isolation_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
total_detections_by_city <- site_data$total_detections
contributor_detections_by_city <- site_data$max_detections_by_contributor_scaled
city <- as.integer(as.factor(unique(site_data$city)))
n_cities <- length(unique(city))

temp <- site_data %>%
  select(city, total_detections, mean_detections_by_contributor_scaled, max_detections_by_contributor_scaled) %>%
  group_by(city) %>%
  slice(1)

## ranges
#ranges <- my_data$ranges 

# extra stuff
R <- my_data$R
species_integer_vector <- my_data$species_integer_vector
city_integer_vector <- my_data$city_integer_vector
multicity_site_id_vector <- my_data$multicity_site_id_vector
city_id_vector <- my_data$city_id_vector
site_survey_year_vector <- my_data$site_survey_year_vector
prev_index_vector <- my_data$prev_index_vector
confirmed_occurrence <- my_data$confirmed_occurrence
species_cluster_id_vector <- my_data$species_cluster_integer_vector
n_species_clusters <- length(unique(species_cluster_id_vector))
regional_cluster_id_vector <- my_data$region_cluster_integer_vector
n_regional_clusters <- length(unique(regional_cluster_id_vector))


stan_data <- list(R = R, n_surveys = n_surveys, surveys = surveys,
                  V = V, V_NA = V_NA, site_survey_year_vector = site_survey_year_vector,
                  n_species = n_species, species = species, species_integer_vector = species_integer_vector,
                  n_sites = n_sites, sites = sites, multicity_site_id_vector = multicity_site_id_vector,
                  n_cities = n_cities, city = city, city_id_vector = city_id_vector,
                  feature_diversity = feature_diversity, ease_of_id = ease_of_id, wingspan = wingspan,
                  park_size = park_size, isolation = isolation, total_detections_by_city = total_detections_by_city,
                  confirmed_occurrence = confirmed_occurrence, prev_index_vector = prev_index_vector, 
                  species_cluster_id_vector = species_cluster_id_vector, n_species_clusters = n_species_clusters,
                  regional_cluster_id_vector = regional_cluster_id_vector, n_regional_clusters = n_regional_clusters         
) 

## Parameters monitored 
params <- c(
  
  "psi1_0", 
  "sigma_psi1_species",
  "sigma_psi1_city",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation",
  
  "gamma0", 
  "sigma_gamma_species",
  "sigma_gamma_city",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation",
  
  "phi0", 
  "sigma_phi_species",
  "sigma_phi_city",
  "mu_phi_wingspan",
  "sigma_phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_isolation",
  "sigma_phi_isolation",
  
  "p0", 
  "sigma_p_species",
  "sigma_p_city",
  "p_total_detections",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  #"mu_p_species_date",
  "delta0",
  "delta_regional_cluster",
  "sigma_p_species_date",
  #"mu_p_species_date_sq",
  "epsilon0",
  "epsilon_regional_cluster",
  "sigma_p_species_date_sq",
  
  # city effects
  "psi1_city",
  "psi1_wingspan",
  "psi1_park_size",
  "psi1_isolation",
  "gamma_city",
  "gamma_wingspan",
  "gamma_park_size",
  "gamma_isolation",
  "phi_city",
  "phi_wingspan",
  "phi_park_size",
  "phi_isolation",
  "p_city",
  
  # species effects and PPC
  #"W_species_rep",
  "psi1_species",
  "gamma_species", "phi_species",
  "p_species"
)

# MCMC settings
n_iterations <- 150
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- parallel::detectCores()
delta = 0.97

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0.5, 1),
       sigma_psi1_city = runif(1, 0.5, 1),
       sigma_psi1_wingspan = runif(1, 0.5, 1),
       sigma_psi1_park_size  = runif(1, 0.5, 1),
       sigma_psi1_isolation  = runif(1, 0.5, 1),
       psi1_wingspan = runif(1, -1, 1),
       mu_psi1_park_size = runif(1, 0, 1),
       gamma0 = runif(1, -3, -2),
       sigma_gamma_species = runif(1, 0.5, 1),
       sigma_gamma_city = runif(1, 0.5, 1),
       sigma_gamma_wingspan = runif(1, 0.5, 1),
       sigma_gamma_park_size = runif(1, 0.5, 1),
       sigma_gamma_isolation = runif(1, 0.5, 1),
       gamma_wingspan = runif(1, 1, 2),
       mu_gamma_park_size = runif(1, 0.5, 1),
       phi0 = runif(1, 2, 3),
       sigma_phi_species= runif(1, 0.5, 1),
       sigma_phi_city= runif(1, 0.5, 1),
       sigma_phi_wingspan= runif(1, 0.5, 1),
       sigma_phi_park_size= runif(1, 0.5, 1),
       sigma_phi_isolation= runif(1, 0.5, 1),
       phi_wingspan = runif(1, -1, 0),
       mu_phi_park_size = runif(1, 0.5, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0.5, 1),
       sigma_p_city = runif(1, 0.5, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       sigma_p_species_date = runif(1,0.5, 1),
       sigma_p_species_date_sq = runif(1, 0.5, 1)
  )
)


## --------------------------------------------------
### Run model

# choose a model
#stan_model <- "./models/dynamic_occupancy_model.stan"
stan_model <- cmdstan_model("./models/dynamic_occupancy_model_all_cities.stan")

## Call Stan from R
stan_out <- stan_model$sample(
  data = stan_data, 
  refresh = 50,
  init = inits, 
  chains = n_chains, 
  parallel_chains = n_cores,
  iter_sampling = n_iterations, 
  iter_warmup = n_burnin, 
  thin = n_thin,
  seed = 1,
  adapt_delta=delta)

# save the object
stan_out$save_object(file = "stan_out.rds")
# or read the object back into the environment
#stan_out <- readRDS("./model_outputs/stan_out_dec7.rds")


# summarise all variables with default and additional summary measures
temp <- as.data.frame(stan_out$summary(
  variables = c(
    "psi1_0", 
    "sigma_psi1_species",
    "sigma_psi1_city",
    "mu_psi1_wingspan",
    "sigma_psi1_wingspan",
    "mu_psi1_park_size",
    "sigma_psi1_park_size",
    "mu_psi1_isolation",
    "sigma_psi1_isolation",
    
    "gamma0", 
    "sigma_gamma_species",
    "sigma_gamma_city",
    "mu_gamma_wingspan",
    "sigma_gamma_wingspan",
    "mu_gamma_park_size",
    "sigma_gamma_park_size",
    "mu_gamma_isolation",
    "sigma_gamma_isolation",
    
    "phi0", 
    "sigma_phi_species",
    "sigma_phi_city",
    "mu_phi_wingspan",
    "sigma_phi_wingspan",
    "mu_phi_park_size",
    "sigma_phi_park_size",
    "mu_phi_isolation",
    "sigma_phi_isolation",
    
    "p0", 
    "sigma_p_species",
    "sigma_p_city",
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
    "psi1_city",
    "psi1_wingspan",
    "psi1_park_size",
    "psi1_isolation",
    "gamma_city",
    "gamma_wingspan",
    "gamma_park_size",
    "gamma_isolation",
    "phi_city",
    "phi_wingspan",
    "phi_park_size",
    "phi_isolation",
    "p_city"),
  
  posterior::default_summary_measures(),
  extra_quantiles = ~posterior::quantile2(., probs = c(0.25, .75))
))

stan_out$summary(variables = c(
  "psi1_0", 
  "sigma_psi1_species",
  "sigma_psi1_city",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation",
  
  "gamma0", 
  "sigma_gamma_species",
  "sigma_gamma_city",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation",
  
  "phi0", 
  "sigma_phi_species",
  "sigma_phi_city",
  "mu_phi_wingspan",
  "sigma_phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_isolation",
  "sigma_phi_isolation",
  
  "p0", 
  "sigma_p_species",
  "sigma_p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "delta0",
  "delta_regional_cluster",
  "sigma_p_species_date",
  "epsilon0",
  "epsilon_regional_cluster",
  "sigma_p_species_date_sq"
))

library(bayesplot)
mcmc_trace(stan_out$draws(), pars = c(
  "psi1_0", 
  "sigma_psi1_species",
  "sigma_psi1_city",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation"
))

mcmc_trace(stan_out$draws(), pars = c(
  "gamma0", 
  "sigma_gamma_species",
  "sigma_gamma_city",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation"
))

mcmc_trace(stan_out$draws(), pars = c(
  "p_city_detections"
))
