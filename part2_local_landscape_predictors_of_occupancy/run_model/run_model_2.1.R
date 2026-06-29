# run model 2.1

# list of city names
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

# save prepped data
#saveRDS(my_data, paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))
# OR just read in previously prepared data and fit the model
my_data <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/run_model/prepped_data/prepped_data.rds"))


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
park_size <- site_data$log_total_green_space_area_scaled 
tree_cover <- site_data$tree_cover_scaled
plant_diversity <- site_data$plant_genera_density_scaled
landscape_connectivity <- site_data$log_connectivity_scaled
landscape_grassherb <- site_data$proportion_landscape_grassherb_scaled
landscape_woody <- site_data$proportion_landscape_woody_scaled
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


# prepare to fit occupancy model
library(cmdstanr)

stan_data <- list(R = R, n_surveys = n_surveys, surveys = surveys,
                  V = V, V_NA = V_NA, site_survey_year_vector = site_survey_year_vector,
                  n_species = n_species, species = species, species_integer_vector = species_integer_vector,
                  n_sites = n_sites, sites = sites, multicity_site_id_vector = multicity_site_id_vector,
                  n_cities = n_cities, city = city, city_id_vector = city_id_vector,
                  feature_diversity = feature_diversity, ease_of_id = ease_of_id, wingspan = wingspan, migratory = migratory,
                  tree_cover = tree_cover, plant_diversity = plant_diversity,
                  landscape_grassherb = landscape_grassherb, landscape_woody = landscape_woody,
                  park_size = park_size, landscape_connectivity = landscape_connectivity, total_detections_by_city = total_detections_by_city,
                  confirmed_occurrence = confirmed_occurrence, prev_index_vector = prev_index_vector, 
                  species_cluster_id_vector = species_cluster_id_vector, n_species_clusters = n_species_clusters,
                  regional_cluster_id_vector = regional_cluster_id_vector, n_regional_clusters = n_regional_clusters,
                  species_city_id_vector = species_city_id_vector,   n_species_city_clusters = n_species_city_clusters 
) 

# MCMC settings
n_iterations <- 1000
n_thin <- 2
n_burnin <- 500
n_chains <- 4
n_cores <- parallel::detectCores()
delta = 0.97

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
init_generate <- function(chain_id)
  
  list(psi_0 = runif(1, -1, 1),
       sigma_psi_species = runif(1, 0.5, 1),
       sigma_psi_city = runif(1, 0.5, 1),
       sigma_psi_wingspan = runif(1, 0.5, 1),
       sigma_psi_park_size  = runif(1, 0.5, 1),
       sigma_psi_tree_cover  = runif(1, 0.5, 1),
       sigma_psi_landscape_connectivity  = runif(1, 0.5, 1),
       sigma_psi_landscape_grassherb  = runif(1, 0.5, 1),
       sigma_psi_landscape_woody  = runif(1, 0.5, 1),
       mu_psi_wingspan = runif(1, -1, 1),
       mu_psi_park_size = runif(1, 0, 1),
       mu_psi_tree_cover= runif(1, 0, 1),
       mu_psi_landscape_connectivity = runif(1, 0, 1),
       mu_psi_landscape_grassherb = runif(1, 0, 1),
       mu_psi_landscape_woody = runif(1, 0, 1),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0.5, 1),
       sigma_p_city = runif(1, 0.5, 1),
       p_wingspan = runif(1, -1, 1),
       p_feature_diversity = runif(1, -1, 1),
       p_ease_of_id = runif(1, -1, 1),
       sigma_p_species_date = runif(1,0.5, 1),
       sigma_p_species_date_sq = runif(1, 0.5, 1)
       
  )


## --------------------------------------------------
### Run model

# choose a model
#stan_model <- "./models/dynamic_occupancy_model.stan"
stan_model <- cmdstan_model("./models/occupancy_model_m2.1.stan")

## Call Stan from R
stan_out <- stan_model$sample(
  data = stan_data, 
  refresh = 50,
  init = init_generate, 
  chains = n_chains, 
  parallel_chains = n_cores,
  iter_sampling = n_iterations, 
  iter_warmup = n_burnin, 
  thin = n_thin,
  seed = 1,
  adapt_delta=delta)

# save the object
stan_out$save_object(file = "stan_out_m2.1.rds")

# read old data
#stan_out <- readRDS( paste0("./part2_local_landscape_predictors_of_occupancy/model_outputs/stan_out_m2.1.rds"))

stan_out$diagnostic_summary()

rhats <- stan_out$summary(c("psi_0", 
                            "sigma_psi_species",
                            "sigma_psi_city",
                            "psi_wingspan",
                            "psi_migratory",
                            "mu_psi_park_size",
                            "sigma_psi_park_size",
                            "mu_psi_tree_cover",
                            "sigma_psi_tree_cover",
                            "mu_psi_plant_diversity",
                            "sigma_psi_plant_diversity",
                            "mu_psi_landscape_connectivity",
                            "sigma_psi_landscape_connectivity",
                            "mu_psi_landscape_grassherb",
                            "sigma_psi_landscape_grassherb",
                            "mu_psi_landscape_woody",
                            "sigma_psi_landscape_woody",
                            
                            "p0", 
                            "sigma_p_species",
                            "sigma_p_city",
                            "p_city_detections",
                            "p_wingspan",
                            "p_migratory",
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
                            "psi_landscape_connectivity",
                            "psi_landscape_grassherb",
                            "psi_landscape_woody",
                            "p_city"), 
                          "rhat")

library(bayesplot)

mcmc_rhat_hist(rhats$rhat) +
  ggplot2::ggtitle("Rhats for all m2.1 parameters")

mcmc_trace(stan_out$draws(), pars = c(
  "psi_0", 
  "sigma_psi_species",
  "sigma_psi_city",
  "mu_psi_park_size",
  "sigma_psi_park_size",
  "mu_psi_tree_cover",
  "sigma_psi_tree_cover",
  "mu_psi_plant_diversity",
  "sigma_psi_plant_diversity",
  "mu_psi_landscape_connectivity",
  "sigma_psi_landscape_connectivity",
  "mu_psi_landscape_grassherb",
  "sigma_psi_landscape_grassherb",
  "mu_psi_landscape_woody",
  "sigma_psi_landscape_woody"
))