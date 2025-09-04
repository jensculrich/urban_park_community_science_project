
# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data_multicity.R")

# list of city names

# northeast
city_names <- c(
  "Boston", 
  "NYC", 
  "Philadelphia"
  #"Washington",
  
)

# southeast
city_names <- c(
  #"Atlanta"
  "Charlotte",
  "Dallas",
  "Houston",
  "Raleigh"
)

# southwest
city_names <- c(
  "LA",
  #"Pheonix",
  #"Riverside",
  "SD",
  "SF"
)

# midwest
city_names <- c(
  #"Chicago",
  #"Des Moines"
  #"Detroit", 
  "Minneapolis",
  #"Omaha",
  #"St. Louis"
)

# now choose a city (enter the number of the city)
#city <- city_names[2]

min_species_detections <- 2 # binary park/year/species detections
min_species_for_community_sampling_event = 1 
family_sampling = TRUE # Should enter either TRUE or FALSE 
# family_sampling:
# if false infer sampling event for all butterflies if any butterflies detected
# if true only infer sampling event for butterflies in same family as any butterflies detected

my_data <- prep_data(city_names,
                     min_species_detections,
                     min_species_for_community_sampling_event,
                     family_sampling
)

saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_southeast.rds"))
my_data <- readRDS( paste0("./run_model/prepped_data/prepped_data_southeast.rds"))

# prepare to fit occupncy model
library(rstan)

# data to feed to the model
# detection data
V <- my_data$V_detections # detections (1==detected)
V_NA <- my_data$V_NA # NAs (0==no known survey effort made)

# covariate data
species_info <- my_data$species_info
n_species <- my_data$n_species # number of species
species <- as.integer(as.factor(species_info$species))
site_data <- my_data$site_data
n_sites <- my_data$n_sites # number of sites
sites <- site_data$new_id
n_years <- my_data$n_years # number of surveys
n_years_minus1 <- n_years - 1
years <- seq(1, n_years_minus1)
n_surveys <- my_data$n_surveys
surveys <- sequence(n_surveys)
surveys <- (surveys - mean(surveys)) / sd(surveys)

## predictors
# species
feature_diversity <- species_info$featureDiversity_scaled
ease_of_id <- species_info$research_grade_proportion_scaled
wingspan <- species_info$aveWingspan_scaled
# site
park_size <- site_data$log_total_green_space_area_scaled
connectivity <- site_data$connectivity_scaled
city <- as.integer(as.factor(site_data$city))
n_cities <- length(city_names)

stan_data <- c("V", "V_NA", "species", "sites", "years", "surveys", 
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys",
               "feature_diversity", "ease_of_id", "wingspan",
               "park_size", "connectivity", "city", "n_cities"
) 

## Parameters monitored 
params <- c(
  
  "psi1_0", 
  "sigma_psi1_species",
  "psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_connectivity",
  "sigma_psi1_connectivity",
  "psi1_park_size",
  "psi1_connectivity",

  "gamma0", 
  "sigma_gamma_species",
  "gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_connectivity",
  "sigma_gamma_connectivity",
  "gamma_park_size",
  "gamma_connectivity",

  "phi0", 
  "sigma_phi_species",
  "phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_connectivity",
  "sigma_phi_connectivity",
  "phi_park_size",
  "phi_connectivity",

  "p0", 
  "sigma_p_species",
  #"p_city",
  #"sigma_p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "mu_p_species_date",
  "sigma_p_species_date",
  "mu_p_species_date_sq",
  "sigma_p_species_date_sq",
  
  "W_species_rep",
  "psi1_species",
  "gamma_species", "phi_species",
  "p_species"
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.97

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 1, 2),
       psi1_wingspan = runif(1, -1, 1),
       #psi1_park_size = runif(1, -1, 1),
       gamma0 = runif(1, -3, -2),
       gamma_wingspan = runif(1, 1, 2),
       #gamma_park_size = runif(1, 0, 1),
       #gamma_connectivity = runif(1, -2, -1),
       phi0 = runif(1, 2, 3),
       phi_wingspan = runif(1, -1, 0),
       #phi_park_size = runif(1, 0, 1),
       #phi_connectivity = runif(1, 1, 2),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 1, 2),
       #sigma_p_site = runif(1, 0, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       mu_p_species_date_sq = runif(1, -1, 0),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)


## --------------------------------------------------
### Run model

stan_model <- "./models/dynamic_occupancy_model_simple_multicity.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                 data = stan_data, 
                 init = inits, 
                 pars = params,
                 chains = n_chains, iter = n_iterations, 
                 warmup = n_burnin, thin = n_thin,
                 seed = 1,
                 control=list(adapt_delta=delta),
                 open_progress = FALSE,
                 cores = n_cores)

saveRDS(stan_out, paste0("./model_outputs/stan_out_southeast_2km_connectivity_family_50buffers_simple.rds"))

#stan_out <- readRDS( paste0("./model_outputs/stan_out_NYC_2km_connectivity_family_100buffers.rds"))

print(stan_out, digits = 3, 
      pars = c(
        "psi1_0", 
        "sigma_psi1_species",
        "psi1_wingspan",
        "mu_psi1_park_size",
        "sigma_psi1_park_size",
        "mu_psi1_connectivity",
        "sigma_psi1_connectivity",
        "psi1_park_size",
        "psi1_connectivity",
        
        "gamma0", 
        "sigma_gamma_species",
        "gamma_wingspan",
        "mu_gamma_park_size",
        "sigma_gamma_park_size",
        "mu_gamma_connectivity",
        "sigma_gamma_connectivity",
        "gamma_park_size",
        "gamma_connectivity",
        
        "phi0", 
        "sigma_phi_species",
        "phi_wingspan",
        "mu_phi_park_size",
        "sigma_phi_park_size",
        "mu_phi_connectivity",
        "sigma_phi_connectivity",
        "phi_park_size",
        "phi_connectivity",
        
        "p0", 
        "sigma_p_species",
        "p_wingspan",
        "p_feature_diversity",
        "p_ease_of_id",
        "mu_p_species_date",
        "sigma_p_species_date",
        "mu_p_species_date_sq",
        "sigma_p_species_date_sq",
      ))

print(stan_out, digits = 3, 
      pars = c("p_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "psi1_0", 
  "sigma_psi1_species",
  "psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_connectivity",
  "sigma_psi1_connectivity",
  "psi1_park_size",
  "psi1_connectivity"
  
))

traceplot(stan_out, pars = c(
 
  "gamma0", 
  "sigma_gamma_species",
  "gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_connectivity",
  "sigma_gamma_connectivity",
  "gamma_park_size",
  "gamma_connectivity",
  
  "phi0", 
  "sigma_phi_species",
  "phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_connectivity",
  "sigma_phi_connectivity",
  "phi_park_size",
  "phi_connectivity"
))

traceplot(stan_out, pars = c(
  
  "p0", 
  "sigma_p_species",
  "p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "mu_p_species_date",
  "sigma_p_species_date",
  "mu_p_species_date_sq",
  "sigma_p_species_date_sq"
  
))

