
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
  )
}


# or choose one city
# city_names <- "Philadelphia"

min_species_detections <- 1 # binary park/year/species detections
min_site_years_w_detection <- 1 # remove parks never surveyed across repeat years
min_species_for_community_sampling_event <- 1 # if 1 species detected, any other species in same fam could have been 
family_sampling <- TRUE # Should enter either TRUE or FALSE 
remove_outlier_parks <- TRUE # remove very small parks
write_city_data_csv <- TRUE
# family_sampling:
# if false infer sampling event for all butterflies if any butterflies detected
# if true only infer sampling event for butterflies in same family as any butterflies detected

# for community sampling events inferred by [taxonomic] family, source this file:
source("./part2_local_landscape_predictors_of_occupancy/run_model/prep_data.R")

my_data <- prep_data(city_names,
                     min_species_detections,
                     min_site_years_w_detection,
                     min_species_for_community_sampling_event,
                     family_sampling,
                     remove_outlier_parks,
                     write_city_data_csv
)


saveRDS(my_data, paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_", region, ".rds"))
my_data <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data_", region, ".rds"))


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

## predictors
# species
feature_diversity <- species_info$featureDiversity_scaled
ease_of_id <- species_info$research_grade_proportion_scaled
wingspan <- species_info$aveWingspan_scaled
# site
park_size <- site_data$log_total_green_space_area_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
tree_cover <- site_data$tree_cover_scaled_across_all_cities
plant_diversity <- site_data$plant_genera_density_scaled_across_all_cities
landscape_isolation <- site_data$log_isolation_scaled_across_all_cities # scaled_2 is scaled to only parks being modeled
landscape_grassherb <- site_data$proportion_landscape_grassherb_scaled_across_all_cities
landscape_woody <- site_data$proportion_landscape_woody_scaled_across_all_cities
city <- as.integer(as.factor(unique(site_data$city)))
n_cities <- length(unique(city))

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
species_city_cluster <- my_data$species_city_cluster
n_species_city_clusters <- length(unique(species_city_cluster))

# plot
library(tidyverse)
ggplot(site_data, aes(
  x = log_total_green_space_area, y = log(isolation), colour = city)) +
  geom_point()

ggplot(site_data, aes(
  x = log_total_green_space_area, y = total_detections_by_city, colour = city)) +
  geom_point()

# prepare to fit occupncy model
library(rstan)

stan_data <- c("R", "n_surveys", "surveys", 
               "V", "V_NA", "site_survey_year_vector",
               "n_species", "species", "species_integer_vector",
               "n_sites", "sites", "multicity_site_id_vector",
               "n_cities","city", "city_id_vector",
               "feature_diversity", "ease_of_id", "wingspan",
               "park_size", "isolation", 
               "confirmed_occurrence", "prev_index_vector", 
               "species_cluster_id_vector", "n_species_clusters",
               "regional_cluster_id_vector", "n_regional_clusters",
               "species_city_cluster", "n_species_city_clusters"
               
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
  
  "p0", 
  "sigma_p_species",
  "sigma_p_city",
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
       sigma_psi1_city = runif(1, 0, 1),
       sigma_psi1_wingspan = runif(1, 0, 1),
       sigma_psi1_park_size  = runif(1, 0, 1),
       sigma_psi1_isolation  = runif(1, 0, 1),
       psi1_wingspan = runif(1, -1, 1),
       mu_psi1_park_size = runif(1, 0, 1),
       gamma0 = runif(1, -3, -2),
       sigma_gamma_species = runif(1, 0, 1),
       sigma_gamma_city = runif(1, 0, 1),
       sigma_gamma_wingspan = runif(1, 0, 1),
       sigma_gamma_park_size = runif(1, 0, 1),
       sigma_gamma_isolation = runif(1, 0, 1),
       gamma_wingspan = runif(1, 1, 2),
       mu_gamma_park_size = runif(1, 0, 1),
       phi0 = runif(1, 2, 3),
       sigma_phi_species= runif(1, 0, 1),
       sigma_phi_city= runif(1, 0, 1),
       sigma_phi_wingspan= runif(1, 0, 1),
       sigma_phi_park_size= runif(1, 0, 1),
       sigma_phi_isolation= runif(1, 0, 1),
       phi_wingspan = runif(1, -1, 0),
       mu_phi_park_size = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 1, 2),
       sigma_p_city = runif(1, 0, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       #mu_p_species_date = runif(1, -1, 1),
       sigma_p_species_date = runif(1, 0, 1),
       #mu_p_species_date_sq = runif(1, -1, 0),
       sigma_p_species_date_sq = runif(1, 0, 1)
  )
)


## --------------------------------------------------
### Run model

# choose a model
#stan_model <- "./models/dynamic_occupancy_model.stan"
stan_model <- "./models/dynamic_occupancy_model_all_cities.stan"

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


saveRDS(stan_out, paste0("./model_outputs/stan_out_", region, "3.rds"))

# read old data
#stan_out <- readRDS( paste0("./model_outputs/stan_out_", region, "_dec5.rds"))

# print outputs
print(stan_out, digits = 3, 
      pars = c(
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
        "mu_p_species_date",
        "sigma_p_species_date",
        "mu_p_species_date_sq",
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
        "p_city"
      ))

print(stan_out, digits = 3, 
      pars = c("psi1_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "psi1_0", 
  "psi1_city",
  "sigma_psi1_city",
  "sigma_psi1_species",
  "mu_psi1_wingspan",
  "sigma_psi1_wingspan",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation",
  "psi1_park_size",
  "psi1_isolation"
  
))

traceplot(stan_out, pars = c(
  "gamma0", 
  "gamma_city",
  "sigma_gamma_city",
  "sigma_gamma_species",
  "mu_gamma_wingspan",
  "sigma_gamma_wingspan",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation",
  "gamma_park_size",
  "gamma_isolation"
))  

traceplot(stan_out, pars = c(
  "phi0", 
  "phi_city",
  "sigma_phi_city",
  "sigma_phi_species",
  "mu_phi_wingspan",
  "sigma_phi_wingspan",
  "mu_phi_park_size",
  "sigma_phi_park_size",
  "mu_phi_isolation",
  "sigma_phi_isolation",
  "phi_park_size",
  "phi_isolation"
))

traceplot(stan_out, pars = c(
  "p0", 
  "sigma_p_species",
  "p_city",
  "sigma_p_city",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  #"mu_p_species_date",
  "sigma_p_species_date",
  #"mu_p_species_date_sq",
  "sigma_p_species_date_sq"
))

