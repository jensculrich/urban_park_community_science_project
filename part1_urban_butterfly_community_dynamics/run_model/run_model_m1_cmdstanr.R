# I used cmdstanr to fit the models through the computing cluster.
# This code tells the model how to start running.
library(cmdstanr)

## -----------------------------------------------------------------------------
### Prepare some data to feed to the dynamic occupancy model

# city names for cities that we will be looking at:
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

# There are some choices that can be made here that may the analysis
# These flexible commands allow the analysis to be rerun under different conditions:
min_species_detections <- 1 # only include species detected n or more times (binary park/year/species detections)
min_site_years_w_detection <- 2 # 2 == remove parks never surveyed across more than one years (i.e., site must be surveyed in 2 or more years)
min_species_for_community_sampling_event <- 1 # if n species detected, any other species in same taxonomic group could have been 
family_sampling <- TRUE # TRUE == non-detection inferred only if species from the same family were detecteted (as opposed to any butterfly species)
remove_outlier_parks <- TRUE # remove very small parks from the model fitting procedure
write_city_data_csv <- FALSE # save some info about the cities and the detections/parks in them

# source this file which contains the prep data function:
source("./part1_urban_butterfly_community_dynamics/run_model/prep_data.R")
# and get some data (takes a minute or two to process)
my_data <- prep_data(city_names,
                     min_species_detections,
                     min_site_years_w_detection,
                     min_species_for_community_sampling_event,
                     family_sampling,
                     remove_outlier_parks,
                     write_city_data_csv
)

# save the data
saveRDS(my_data, paste0("./part1_urban_butterfly_community_dynamics/run_model/prepped_data/prepped_data_m1.rds"))
# alternatively load in some data that you've already prepared
my_data <- readRDS( paste0("./part1_urban_butterfly_community_dynamics/run_model/prepped_data/prepped_data_", region, ".rds"))

## -----------------------------------------------------------------------------
### unpack the data and bundle into a list for stan

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

## covariate data
# species covariates
feature_diversity <- species_info$featureDiversity_scaled
ease_of_id <- species_info$research_grade_proportion_scaled
wingspan <- species_info$aveWingspan_scaled
migratory <- species_info$migratory

# site or city covariates
park_size <- site_data$log_total_green_space_area_scaled_across_all_cities 
connectivity <- site_data$log_isolation_scaled_across_all_cities * (-1) 
total_detections_by_city <- site_data$total_detections
city <- as.integer(as.factor(unique(site_data$city)))
n_cities <- length(unique(city))

## extra stuff the model needs to be able to run properly
R <- my_data$R # total number of unique (city by species by site) detection/non-detection events
species_integer_vector <- my_data$species_integer_vector # vector length(R) indicating which species sampled
city_integer_vector <- my_data$city_integer_vector # vector length(R) indicating which species sampled
multicity_site_integer_vector <- my_data$multicity_site_integer_vector # vector length(R) indicating which species sampled
regional_cluster_id_vector <- my_data$region_cluster_integer_vector # vector length(R) indicating which regional cluster sampled
  # regional clusters == e.g., Northeastern cities, Southeastern cities, southern california cities, etc.
n_regional_clusters <- length(unique(regional_cluster_id_vector)) # n regional clusters
city_id_vector <- my_data$city_id_vector # vector of unique city IDs
species_cluster_id_vector <- my_data$species_cluster_integer_vector # vector of unique speciesXcluster IDs 
  # eg species cluster id is the same for species A in Boston and NYC, but different in Atlanta (which is in a different cluster)
n_species_clusters <- length(unique(species_cluster_id_vector)) # n species clusters
species_city_id_vector <- my_data$species_city_cluster # vector of unique speciesXcity IDs
  # eg species cluster is different for species A whether in Boston or NYC or Atlanta 
n_species_city_clusters <- length(unique(species_city_id_vector)) # n speciesXcity IDs
site_survey_year_vector <- my_data$site_survey_year_vector # is it the first time a site is sampled or some later year?
prev_index_vector <- my_data$prev_index_vector # indicates the row in the df containing the most recent previous sampling 
  # event for a species being sampled at some site in year > 1 (so we can estimate colonization/persistence from prev state)
confirmed_occurrence <- my_data$confirmed_occurrence # 1 == species was detected > 1 in some month of the year at a site
  # i.e., we have confirmed that it is present. Forks the likelihood into a simpler logistic rather than marginalized outcome.

# bundle the data into a list
stan_data <- list(R = R, n_surveys = n_surveys, surveys = surveys,
                  V = V, V_NA = V_NA, site_survey_year_vector = site_survey_year_vector,
                  n_species = n_species, species = species, species_integer_vector = species_integer_vector,
                  n_sites = n_sites, sites = sites, multicity_site_integer_vector = multicity_site_integer_vector,
                  n_cities = n_cities, city = city, city_id_vector = city_id_vector,
                  feature_diversity = feature_diversity, ease_of_id = ease_of_id, wingspan = wingspan, migratory = migratory,
                  park_size = park_size, isolation = isolation, total_detections_by_city = total_detections_by_city,
                  confirmed_occurrence = confirmed_occurrence, prev_index_vector = prev_index_vector, 
                  species_cluster_id_vector = species_cluster_id_vector, n_species_clusters = n_species_clusters,
                  regional_cluster_id_vector = regional_cluster_id_vector, n_regional_clusters = n_regional_clusters,
                  species_city_id_vector = species_city_id_vector,   n_species_city_clusters = n_species_city_clusters      
) 

## -----------------------------------------------------------------------------
### Set stan settings

# MCMC settings
n_iterations <- 1000 # n post warmup samples
n_thin <- 2 # save every other sample (to reduce the save size of the output)
n_warmup <- 500 # n warmup samples
n_chains <- 4 # run 4 parallel chains
n_cores <- parallel::detectCores()
delta = 0.97 # increase the adapt delta somewhat above the default (slows down the warmup step, but helps traverse more difficult surface)

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
init_generate <- function(chain_id)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0.5, 1),
       sigma_psi1_city = runif(1, 0.5, 1),
       sigma_psi1_park_size  = runif(1, 0.5, 1),
       sigma_psi1_isolation  = runif(1, 0.5, 1),
       psi1_wingspan = runif(1, -1, 1),
       mu_psi1_park_size = runif(1, 0, 1),
       gamma0 = runif(1, -3, -2),
       sigma_gamma_species = runif(1, 0.5, 1),
       sigma_gamma_city = runif(1, 0.5, 1),
       sigma_gamma_park_size = runif(1, 0.5, 1),
       sigma_gamma_isolation = runif(1, 0.5, 1),
       gamma_wingspan = runif(1, 1, 2),
       mu_gamma_park_size = runif(1, 0.5, 1),
       phi0 = runif(1, 2, 3),
       sigma_phi_species= runif(1, 0.5, 1),
       sigma_phi_city= runif(1, 0.5, 1),
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


## -----------------------------------------------------------------------------
### Run model

# choose a model
stan_model <- cmdstan_model(
  "./part1_urban_butterfly_community_dynamics/models/dynamic_occupancy_model_all_cities.stan")

## Call Stan from R
stan_out <- stan_model$sample(
  data = stan_data, 
  refresh = 50,
  init = inits, 
  chains = n_chains, 
  parallel_chains = n_cores,
  iter_sampling = n_iterations, 
  iter_warmup = n_warmup, 
  thin = n_thin,
  seed = 1,
  adapt_delta=delta)

# save the object
stan_out$save_object(file = "./part1_urban_butterfly_community_dynamics/model_outputs/stan_out.rds")
# or read the object back into the environment
#stan_out <- readRDS("./model_outputs/stan_out_dec7.rds")


## -----------------------------------------------------------------------------
### Get some model diagnostics

## get param estimates from the region
stan_out <- readRDS(
  "./part1_urban_butterfly_community_dynamics/model_outputs/stan_out_apr2.rds")
gc()

stan_out$diagnostic_summary()

rhats <- stan_out$summary(c("psi1_0", 
                   "sigma_psi1_species",
                   "sigma_psi1_city",
                   "mu_psi1_park_size",
                   "sigma_psi1_park_size",
                   "mu_psi1_isolation",
                   "sigma_psi1_isolation",
                   
                   "gamma0", 
                   "sigma_gamma_species",
                   "sigma_gamma_city",
                   "mu_gamma_park_size",
                   "sigma_gamma_park_size",
                   "mu_gamma_isolation",
                   "sigma_gamma_isolation",
                   
                   "phi0", 
                   "sigma_phi_species",
                   "sigma_phi_city",
                   "mu_phi_park_size",
                   "sigma_phi_park_size",
                   "mu_phi_isolation",
                   "sigma_phi_isolation",
                   
                   "p0", 
                   "sigma_p_species",
                   "sigma_p_city",
                   "p_city_detections",
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
                   "p_species"), "rhat")

mcmc_rhat_hist(rhats$rhat) +
  ggplot2::ggtitle("Rhats for all m1 parameters")
  
library(bayesplot)
mcmc_trace(stan_out$draws(), pars = c(
  "psi1_0", 
  "sigma_psi1_species",
  "sigma_psi1_city",
  "mu_psi1_park_size",
  "sigma_psi1_park_size",
  "mu_psi1_isolation",
  "sigma_psi1_isolation"
))

mcmc_trace(stan_out$draws(), pars = c(
  "gamma0", 
  "sigma_gamma_species",
  "sigma_gamma_city",
  "mu_gamma_park_size",
  "sigma_gamma_park_size",
  "mu_gamma_isolation",
  "sigma_gamma_isolation"
))

mcmc_trace(stan_out$draws(), pars = c(
  "p_city_detections"
))
