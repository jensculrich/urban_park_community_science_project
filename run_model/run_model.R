# for community sampling events inferred by [taxonomic] order, source this file:
# source("./run_model/prep_data_iNat_w_community_sampling_events.R")

# for community sampling events inferred by [taxonomic] family, source this file:
source("./run_model/prep_data.R")

min_species_detections <- 2 # binary park/year/species detections
min_species_for_community_sampling_event = 1 
family_sampling = FALSE # Should enter either TRUE or FALSE 
# family_sampling:
# if false infer sampling event for all butterflies if any butterflies detected
# if true only infer sampling event for butterflies in same family as any butterflies detected

my_data <- prep_data(min_species_detections,
                     min_species_for_community_sampling_event,
                     family_sampling
                     )

#saveRDS(my_data, paste0("./run_model/prepped_data/prepped_data_LA_family.RDS"))

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
park_size <- site_data$total_green_space_area
connectivity <- site_data$avg_dist_1000m_scaled
tree_cover <- site_data$tree_cover_scaled

stan_data <- c("V", "V_NA", "species", "sites", "years", "surveys", 
               "n_species", "n_sites", "n_years", "n_years_minus1", "n_surveys",
               "feature_diversity", "ease_of_id", "wingspan",
               "park_size", "connectivity", "tree_cover"
) 

## Parameters monitored 
params <- c("psi1_0", 
           "sigma_psi1_species",
           "psi1_wingspan",
           "psi1_park_size",
           "psi1_connectivity",
           "psi1_tree_cover",
           
           "gamma0", 
           "sigma_gamma_species",
           "gamma_wingspan",
           "gamma_park_size",
           "gamma_connectivity",
           "gamma_tree_cover",
           
           "phi0", 
           "sigma_phi_species",
           "phi_wingspan",
           "phi_park_size",
           "phi_connectivity",
           "phi_tree_cover",
           
           "p0", 
           "sigma_p_species",
           "p_wingspan",
           "p_feature_diversity",
           "p_ease_of_id",
           "mu_p_species_date",
           "sigma_p_species_date",
           "mu_p_species_date_sq",
           "sigma_p_species_date_sq",
           
           "W_species_rep",
           "psi1_species", "gamma_species", "phi_species", "p_species"
)

# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- n_chains
delta = 0.95

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  
  list(psi1_0 = runif(1, -1, 1),
       sigma_psi1_species = runif(1, 0, 1),
       psi1_wingspan = runif(1, -1, 1),
       psi1_park_size = runif(1, -1, 1),
       gamma0 = runif(1, -4, -3),
       phi0 = runif(1, 2, 3),
       p0 = runif(1, -1, 1),
       sigma_p_species = runif(1, 0, 1),
       sigma_p_site = runif(1, 0, 1),
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

stan_model <- "./models/dynamic_occupancy_model_stricter_priors.stan"

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

saveRDS(stan_out, "./model_outputs/stan_out_2km_connectivity_order.rds")
#stan_out <- readRDS("./model_outputs/stan_out.rds")

print(stan_out, digits = 3, 
      pars = c("psi1_0", 
               "sigma_psi1_species",
               "psi1_wingspan",
               "psi1_park_size",
               "psi1_connectivity",
               "psi1_tree_cover",
               
               "gamma0", 
               "sigma_gamma_species",
               "gamma_wingspan",
               "gamma_park_size",
               "gamma_connectivity",
               "gamma_tree_cover",
               
               "phi0", 
               "sigma_phi_species",
               "phi_wingspan",
               "phi_park_size",
               "phi_connectivity",
               "phi_tree_cover",
               
               "p0", 
               "sigma_p_species",
               "p_wingspan",
               "p_feature_diversity",
               "p_ease_of_id",
               "mu_p_species_date",
               "sigma_p_species_date",
               "mu_p_species_date_sq",
               "sigma_p_species_date_sq"
      ))

print(stan_out, digits = 3, 
      pars = c("p_species"
      ))

# traceplots
traceplot(stan_out, pars = c(
  "psi1_0", 
  "sigma_psi1_species",
  "psi1_wingspan",
  "psi1_park_size",
  
  "gamma0", 
  "sigma_gamma_species",
  "gamma_wingspan",
  "gamma_park_size",
  
  "phi0", 
  "sigma_phi_species",
  "phi_wingspan",
  "phi_park_size"
))

traceplot(stan_out, pars = c(
  "p0", 
  "sigma_p_species",
  "p_wingspan",
  "p_feature_diversity",
  "p_ease_of_id",
  "mu_p_species_date",
  "sigma_p_species_date",
  "mu_p_species_date_sq"
))

